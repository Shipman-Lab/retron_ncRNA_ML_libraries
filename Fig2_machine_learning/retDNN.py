import os
import numpy as np
import random

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F
import torchvision
from torchvision import transforms
import scipy.stats as stats

import ray
from ray import tune
from ray.tune import CLIReporter
from ray.tune.schedulers import ASHAScheduler

def _get_conv1d_out_length(l_in, padding, dilation, kernel, stride):
    return int((l_in + 2*padding - dilation*(kernel-1) - 1)/stride)+1

def _get_1h_seq(seq):
    nt_code = {
    'N': [0,0,0,0],
    'A': [1,0,0,0],
    'C': [0,1,0,0],
    'G': [0,0,1,0],
    'T': [0,0,0,1]}
    
    seq_1h = []
    for nt in seq:
        seq_1h.append(nt_code[nt])
    return np.array(seq_1h)

def _rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def pearsonr_pt(x, y):
    mean_x = torch.mean(x)
    mean_y = torch.mean(y)
    xm = x.sub(mean_x)
    ym = y.sub(mean_y)
    r_num = xm.dot(ym)
    r_den = torch.norm(xm, 2) * torch.norm(ym, 2)
    r_val = r_num / r_den
    return r_val

class retDataset(Dataset):
    def __init__(self, seq_list, score_dict, seq_dict, max_len=200, left_pad_max=8):
        """
        retDataset is the data seq_idrator for the retDNN model
        Arguments
        ---------
        seq_list:    a list of sequence names that are queried
        score_dict:  a dict of sequence names as keys and scores as values
        seq_dict:    a dict of sequence names as keys and sequences themselves as values
        """
        self.seq_list = seq_list
        self.score_dict = score_dict
        self.seq_dict = seq_dict
        self.max_len = max_len
        self.left_pad_max = left_pad_max
    
    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        
        seq_id = self.seq_list[idx]        
        score = self.score_dict[seq_id]
        seq = self.seq_dict[seq_id]
        
        ##Get sequence, pad it to max_len, and augment data by moving a random len
        lp = random.randint(0,self.left_pad_max)
        seq = _get_1h_seq('N'*(lp) + str(seq) + 'N'*(self.max_len-len(seq)-lp))
        #seq = _get_1h_seq(str(seq) + 'N'*(self.max_len-len(seq)))

        seq = np.array(seq).swapaxes(0,1)
        
        return (seq_id, torch.from_numpy(seq.astype(np.float32)),
                torch.from_numpy(np.array([score]).astype(np.float32)))
    
    def __len__(self):
        return len(self.seq_list)

class retDNNModel(nn.Module):
    def __init__(self, n_channel=4, max_len=200, 
                 conv1kc=64, conv1ks=10, conv1st=1, conv1pd=10, pool1ks=10, pool1st=5 , pdrop1=0.2, #conv_block_1 parameters
                 conv2kc=32, conv2ks=5 , conv2st=1, conv2pd=5 , pool2ks=5 , pool2st=5 , pdrop2=0.2, #conv_block_2 parameters
                 convdc = 3, convdkc=4 , convdks=3 , #dilation block parameters
                 fchidden = 10, pdropfc=0.5, final=1, #fully connected parameters
                 
                 opt="Adam", loss="mse", lr=1e-3, momentum=0.9, weight_decay=1e-3
                ):
        super(retDNNModel, self).__init__()
        
        self.convdc = convdc
        self.opt = opt
        self.loss = loss
        self.momentum = momentum
        self.weight_decay = weight_decay
        self.lr = lr
        
        ## CNN + dilated CNN
        self.conv_block_1 = nn.Sequential(
            nn.Conv1d(n_channel, conv1kc, kernel_size=conv1ks, stride=conv1st, padding=conv1pd),
            nn.LeakyReLU(),
            nn.BatchNorm1d(conv1kc),
            nn.MaxPool1d(kernel_size=pool1ks, stride=pool1st),
            nn.Dropout(p=pdrop1),
        )
        conv_block_1_out_len = _get_conv1d_out_length(max_len, conv1pd, 1, conv1ks, conv1st) #l_in, padding, dilation, kernel, stride
        mpool_block_1_out_len = _get_conv1d_out_length(conv_block_1_out_len, 0, 1, pool1ks, pool1st)
        
        self.conv_block_2 = nn.Sequential(
            nn.Conv1d(conv1kc, conv2kc, kernel_size=conv2ks, stride=conv2st, padding=conv2pd),
            nn.LeakyReLU(),
            nn.BatchNorm1d(conv2kc),
            nn.MaxPool1d(kernel_size=pool2ks, stride=pool2st),
            nn.Dropout(p=pdrop2),
        )
        conv_block_2_out_len = _get_conv1d_out_length(mpool_block_1_out_len, conv2pd, 1, conv2ks, conv2st)
        mpool_block_2_out_len = _get_conv1d_out_length(conv_block_2_out_len, 0, 1, pool2ks, pool2st)
        
        self.dilations = nn.ModuleList()
        for i in range(convdc):
            self.dilations.append(nn.Sequential(
                nn.Conv1d(conv2kc, conv2kc, kernel_size=convdks, padding=2**i, dilation=2**i),
                nn.LeakyReLU(),
                nn.BatchNorm1d(conv2kc))
            )
        dilation_blocks_lens = []
        for i in range(convdc):
            dilation_blocks_lens.append( _get_conv1d_out_length(mpool_block_2_out_len, 2**i, 2**i, convdks, 1)  * conv2kc) #(l_in, padding, dilation, kernel, stride):
        
        total_length = mpool_block_2_out_len * conv2kc #+ sum(dilation_blocks_lens)
        #print(total_length, mpool_block_2_out_len, conv2kc, dilation_blocks_lens, convdks)
        
        self.fc = nn.Sequential(
            nn.Dropout(p=pdropfc),
            nn.Linear(total_length, fchidden),
            nn.ReLU(),
            nn.Dropout(p=pdropfc),
            nn.Linear(fchidden, final)
        )
        
    def forward(self, seq):
        seq = self.conv_block_1(seq)
        seq = self.conv_block_2(seq)
        current = seq
        for i in range(self.convdc):
            residual = self.dilations[i](current)
            current = current.add(residual)
        #res = current.transpose(1, 2)
        res = torch.flatten(current, 1)
        res = self.fc(res)
        return res
        
    def compile(self, device='cpu'):
        self.to(device)
        if self.opt=="Adam":
            self.optimizer = optim.Adam(self.parameters(), lr=self.lr)
        if self.opt=="SGD":
            self.optimizer = optim.SGD(self.parameters(), lr=self.lr, momentum = self.momentum)
        if self.opt=="Adagrad":
            self.optimizer = optim.Adagrad(self.parameters(), lr=self.lr, weight_decay = self.weight_decay)
        if self.loss=="mse":
            self.loss_fn = F.mse_loss
            
    def train_model(self, train_loader, val_loader=None, epochs=100, clip=10, device='cpu', modelfile='models/best_ret_checkpoint.pt', logfile = None, tuneray=False, verbose=True):
        train_losses = []
        train_Rs = []
        valid_losses = []
        valid_Rs = []
        best_model=-1
        for epoch in range(epochs):
            training_loss = 0.0
            train_R = 0.0
            
            for batch_idx, batch in enumerate(train_loader):
                self.train()
                self.optimizer.zero_grad()
                seqid,seq_X,y = batch
                seq_X = seq_X.to(device)
                y = y.to(device)
                out = self(seq_X)
                loss = self.loss_fn(out,y)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.parameters(), clip)
                #torch.nn.utils.clip_grad_value_(self.parameters(), clip)
                self.optimizer.step()
                training_loss = loss.data.item() * seq_X.size(0) #divide by two since every sequence is paired F/R
                train_R = pearsonr_pt(out[:,0],y[:,0]).to('cpu').detach().numpy()
                
                train_losses.append(training_loss)
                train_Rs.append(train_R)
                
                if verbose and batch_idx%10==0:
                    print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}\tR: {:.6f}'.format(
                            epoch, batch_idx * train_loader.batch_size, len(train_loader.dataset), 100. * batch_idx * float(train_loader.batch_size) / len(train_loader.dataset),
                            loss.item(), train_R))
                    
            if val_loader:
                target_list = []
                pred_list = []
                valid_loss = 0.0
                valid_R = 0.0
                self.eval()
                for batch_idx, batch in enumerate(val_loader):
                    seqid,seq_X,y = batch
                    seq_X = seq_X.to(device)
                    y = y.to(device)
                    out = self(seq_X)
                    loss = self.loss_fn(out[:,0],y[:,0])
                    valid_loss += loss.data.item() * seq_X.size(0)
                    pred_list.append(out.to('cpu').detach().numpy())
                    target_list.append(y.to('cpu').detach().numpy())
                targets = np.concatenate(target_list)
                preds = np.concatenate(pred_list)
                valid_R = stats.pearsonr(targets[:,0], preds[:,0])[0]
                valid_loss /= len(val_loader.dataset)

                valid_losses.append(valid_loss)
                valid_Rs.append(valid_R)
                
                if tuneray:
                    with tune.checkpoint_dir(step=epoch) as checkpoint_dir:
                        path = os.path.join(checkpoint_dir, "checkpoint")
                        torch.save(
                            (self.state_dict(), self.optimizer.state_dict()), path)

                    tune.report(loss=valid_loss, corel=valid_R)

                if verbose:
                    print('Validation: Loss: {:.6f}\tR: {:.6f}'.format(valid_loss, valid_R))
                if logfile:
                    logfile.write('Validation: Loss: {:.6f}\tR: {:.6f}\n'.format(valid_loss, valid_R))

                if (valid_R>best_model):
                    best_model = valid_R
                    if modelfile:
                        print('Best model updated.')
                        self.save_model(modelfile)
                
                
        return {'train_loss':train_losses, 'train_Rs':train_Rs, 'valid_loss':valid_losses, 'valid_Rs':valid_Rs}
    
    def save_model(self, filename):
        torch.save(self.state_dict(), filename)
    
    def load_model(self, filename):
        self.load_state_dict(torch.load(filename))


