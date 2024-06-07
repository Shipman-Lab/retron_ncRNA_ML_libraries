from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import xlsxwriter
import numpy as np
import pandas as pd
import random

library_length = 250

first_fasta = 12000 # starting fasta number for site 3
# skpp primers for site 3
current_skpp_fow = "GGTCGAGCCGGAACT"
current_skpp_rev = "TCTGGGTGCGCATCC"
nts = ['G','A','T','C']
random.seed(35)

# import current synthesis parts
with open('REp_site3_parts_synth.fasta') as fasta_file:
    fasta_nums = []
    sequences = []
    chopped_lengths = []
    for fasta_index, seq_record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        fasta_nums.append(first_fasta + fasta_index)
        seq = str(seq_record.seq)
        # chop off random nucleotides at end of shorter library members
        # chop off forward primer
        fow_primer_pos = seq.find(current_skpp_fow)
        if fow_primer_pos == -1:
            import pdb
            pdb.set_trace()
        chopped_seq = seq[fow_primer_pos + len(current_skpp_fow):]
        # chop off rev primer
        rev_primer_pos = chopped_seq.find(current_skpp_rev)
        if rev_primer_pos == -1:
            import pdb
            pdb.set_trace()
        chopped_seq = chopped_seq[0:rev_primer_pos]
        sequences.append(chopped_seq)
        chopped_lengths.append(len(chopped_seq))

# find all unique sequence lengths & set 5% bins
unique_seq_lengths = np.unique(chopped_lengths)
current_min = unique_seq_lengths.max()
current_max = unique_seq_lengths.max()
length_bins = [current_max]

while not current_min < unique_seq_lengths.min():
    current_min = np.ceil(current_max - 0.05*current_max)
    length_bins.append(current_min)
    current_max = current_min

length_bins.reverse()

# make a primer set to length bin dictionary:
skpp_dict = {}
record_dict_forward = list(SeqIO.parse("../skpp15-forward.fasta", "fasta"))
record_dict_reverse = list(SeqIO.parse("../skpp15-reverse.fasta", "fasta"))
# for site 3, start again on index 1
start_skpp_index = 1
for bin_index, lower_bin_length in enumerate(length_bins[:]):
    fow_primer = str(record_dict_forward[start_skpp_index + bin_index].seq)
    rev_primer = str(record_dict_reverse[start_skpp_index + bin_index].seq.reverse_complement())
    skpp_dict[lower_bin_length] = [fow_primer, rev_primer]

# now loop through bins & make full sequences
synth_df = pd.DataFrame({"fasta_num": fasta_nums,
                         "seq": sequences,
                         "seq_length_before_fill": chopped_lengths,
                         "bin": np.NaN,
                         "skpp_primer_fow": np.NaN})
synth_df = synth_df.sort_values(by="seq_length_before_fill", axis=0, ascending=True)
overlap_df = pd.DataFrame(columns=["fasta_num", "seq", "seq_length_before_fill", "bin","skpp_primer_fow"])

fasta_num = np.max(fasta_nums) + 1
for bin_index, lower_bin_length in enumerate(length_bins[:-1]):
    bin_seq_indices = (synth_df["seq_length_before_fill"] > length_bins[bin_index]) &\
                      (synth_df["seq_length_before_fill"] <= length_bins[bin_index+1])
    binned_seqs = synth_df.loc[bin_seq_indices, :]
    # 10% overlap with the next sublibrary
    num_overlap = int(np.floor(sum(bin_seq_indices)/10))
    bin_stop_row = synth_df.index.get_loc(synth_df[bin_seq_indices].iloc[-1].name)
    new_bin_stop_row = bin_stop_row + num_overlap
    # import pdb
    # pdb.set_trace()
    overlap_seq = synth_df.iloc[bin_stop_row+1:new_bin_stop_row+1, :].copy()

    synth_df.loc[bin_seq_indices, "bin"] = bin_index
    synth_df.loc[bin_seq_indices, "skpp_primer_fow"] = skpp_dict[lower_bin_length][0]

    overlap_seq.loc[:, "bin"] = bin_index
    overlap_seq.loc[:, "skpp_primer_fow"] = skpp_dict[lower_bin_length][0]
    overlap_seq.loc[:, "fasta_num"] = np.NaN

    for seq_index in binned_seqs.index:
        seq = synth_df.loc[seq_index, "seq"]
        # add on primers at the beginning & end
        primer_seq = skpp_dict[lower_bin_length][0] + seq + skpp_dict[lower_bin_length][1]
        # add random nucleotides until the sequence is 250 long
        random_end = ''.join(random.choice(nts) for i in range(library_length - len(primer_seq)))
        final_seq = primer_seq + random_end 
        synth_df.loc[seq_index, "seq"] = final_seq
    for seq_index in overlap_seq.index:
        seq = overlap_seq.loc[seq_index, "seq"]
        # add on primers at the beginning & end
        primer_seq = skpp_dict[lower_bin_length][0] + seq + skpp_dict[lower_bin_length][1]
        # add random nucleotides until the sequence is 250 long
        random_end = ''.join(random.choice(nts) for i in range(library_length - len(primer_seq)))
        final_seq = primer_seq + random_end 
        overlap_seq.loc[seq_index, "seq"] = final_seq
        overlap_seq.loc[seq_index, "fasta_num"] = fasta_num
        fasta_num += 1
    # append overlap seq onto synth df & reset indices
    # import pdb
    # pdb.set_trace()
    overlap_df = pd.concat([overlap_df, overlap_seq])

synth_df = pd.concat([synth_df, overlap_df], ignore_index=True)

synthesis_list = []
for seq_index in synth_df.index:
    synthesis_list.append(SeqRecord(Seq(synth_df.loc[seq_index, "seq"]),
                          id=str(synth_df.loc[seq_index, "fasta_num"]), description=''))

# write fasta site 2
SeqIO.write(synthesis_list,"REp_site3_parts_synth_sublibraries_OVERLAP.fasta", "fasta")

synth_df.to_excel("Rep_site3_parts_synth_sublibraries_OVERLAP.xlsx")







