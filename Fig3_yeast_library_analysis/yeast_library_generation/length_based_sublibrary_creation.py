from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import xlsxwriter
import numpy as np
import pandas as pd
import random

library_length = 250
# first_fasta = 2931 #starting fasta number for site 1
first_fasta = 7356 # starting fasta number for site 2
current_skpp_fow = "CGATCGCCCTTGGTG"
current_skpp_rev = "ACGCCGGCTAAACC"
nts = ['G','A','T','C']
random.seed(120)

# import current synthesis parts
# with open('REp_site1_parts_synth.fasta') as fasta_file:
with open('REp_site2_parts_synth.fasta') as fasta_file:
    fasta_nums = []
    sequences = []
    chopped_lengths = []
    for fasta_index, seq_record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        fasta_nums.append(first_fasta + fasta_index)
        seq = str(seq_record.seq)
        # chop off random nucleotides at end of shorter library members
        # chop off forward primer
        fow_primer_pos = seq.find(current_skpp_fow)
        chopped_seq = seq[fow_primer_pos + len(current_skpp_fow):]
        # chop off rev primer
        rev_primer_pos = chopped_seq.find(current_skpp_rev)
        chopped_seq = chopped_seq[0:rev_primer_pos]
        sequences.append(chopped_seq)
        chopped_lengths.append(len(chopped_seq))

# find all unique sequence lengths & set 15% bins
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
# for site 1, start on skpp 2 (index 1)
# start_skpp_index = 1
# for site 2, start on skpp 12 (index 11)
start_skpp_index = 11
for bin_index, lower_bin_length in enumerate(length_bins[:]):
    fow_primer = str(record_dict_forward[start_skpp_index + bin_index].seq)
    rev_primer = str(record_dict_reverse[start_skpp_index + bin_index].seq.reverse_complement())
    skpp_dict[lower_bin_length] = [fow_primer, rev_primer]

# now loop through bins & make full sequences
synth_df = pd.DataFrame({"fasta_num": fasta_nums, "seq": sequences, "seq_length_before_fill": chopped_lengths})

for bin_index, lower_bin_length in enumerate(length_bins[:-1]):
    bin_seq_indices = (synth_df["seq_length_before_fill"] > length_bins[bin_index]) &\
                      (synth_df["seq_length_before_fill"] <= length_bins[bin_index+1])
    binned_seqs = synth_df.loc[bin_seq_indices, :]
    for seq_index in binned_seqs.index:
        seq = synth_df.loc[seq_index, "seq"]
        # add on primers at the beginning & end
        primer_seq = skpp_dict[lower_bin_length][0] + seq + skpp_dict[lower_bin_length][1]
        # add random nucleotides until the sequence is 250 long
        random_end = ''.join(random.choice(nts) for i in range(library_length - len(primer_seq)))
        final_seq = primer_seq + random_end 
        synth_df.loc[seq_index, "seq"] = final_seq

# repeat one of the libraries but barcode centered :)
# (I'll do the second shortest library for the most dramatic effect)
bin_seq_indices = (synth_df["seq_length_before_fill"] > length_bins[1]) &\
                  (synth_df["seq_length_before_fill"] <= length_bins[2])
binned_seqs = synth_df.loc[bin_seq_indices, :]
addl_df = pd.DataFrame({"fasta_num": np.NaN * np.empty((sum(bin_seq_indices))),
                        "seq": np.NaN * np.empty((sum(bin_seq_indices))),
                        "seq_length_before_fill": np.NaN * np.empty((sum(bin_seq_indices)))})
final_fasta = np.max(synth_df["fasta_num"])

for fasta_addition, seq_index in enumerate(binned_seqs.index):
    fasta_num = int(final_fasta) + fasta_addition + 1
    seq = sequences[seq_index]
    # add on primers at the beginning & end
    primer_seq = skpp_dict[np.max(length_bins)][0] + seq + skpp_dict[np.max(length_bins)][1]
    # add random nucleotides until the sequence is 250 long
    random_nts = ''.join(random.choice(nts) for i in range(library_length - len(primer_seq)))
    # divide random_nts in half
    len_diff = library_length - len(primer_seq)
    front_random = random_nts[:int(np.ceil(len_diff/2))]
    back_random =random_nts[int(np.ceil(len_diff/2)):]
    final_seq = front_random + primer_seq + back_random
    addl_df.loc[fasta_addition, "fasta_num"] = fasta_num
    addl_df.loc[fasta_addition, "seq"] = final_seq

synth_df = pd.concat([synth_df, addl_df], ignore_index=True)

synthesis_list = []
for seq_index in synth_df.index:
    synthesis_list.append(SeqRecord(Seq(synth_df.loc[seq_index, "seq"]),
                          id=str(synth_df.loc[seq_index, "fasta_num"]), description=''))

# #write fasta site 1
# SeqIO.write(synthesis_list,"REp_site1_parts_synth_sublibraries.fasta", "fasta")

# synth_df.to_excel("Rep_site1_parts_synth_sublibraries.xlsx")

# write fasta site 2
SeqIO.write(synthesis_list,"REp_site2_parts_synth_sublibraries.fasta", "fasta")

synth_df.to_excel("Rep_site2_parts_synth_sublibraries.xlsx")







