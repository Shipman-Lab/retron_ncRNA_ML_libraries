#get access to old libraries with pickle
import sys,os
import pickle
from Bio import SeqIO
import xlsxwriter
import fuzzysearch
import re
from collections import Counter
from Bio.Seq import Seq


"""Line Arguments"""
run_number = sys.argv[1] #MiSeq run number for Data_Path
condition = sys.argv[2] #file name (minus .fastq) 
running_location = sys.argv[3] #either local or otherwise
retron = sys.argv[4] #86
library_number = sys.argv[5] #library, note: not currently set up for pooling
extension_nuc = sys.argv[6] #extended nucleotide


"""Globals"""
wd = os.getcwd()
if running_location == 'local':
    user_profile = os.environ ['USERPROFILE']
    Data_Path = #to be set by user
fastq_reads = '%s/%s.fastq' % (Data_Path,condition)

#pull in counts using fasta sequences, then map them onto the experiments (to avoid duplicates within a library) using these:
with open("%s/part_list_dict_synth.p" % wd, 'rb') as fp:
    ordered_part_list_dict = pickle.load(fp)
with open("%s/part_dict_synth.p" % wd, 'rb') as fp:
    ordered_part_dict = pickle.load(fp)

count = 0
count_dict = {}

total_seqs = 0
match_count = 0
no_first_T = []
no_post_flank = []
no_lib_match = []

"""Run"""
#Create Results folder
newpath = ((r'%s/%s_Results') % (Data_Path,condition))
if not os.path.exists(newpath): os.makedirs(newpath)
#pull in synthesis_parts
part_counts_trimmed = {} #trimmed seq is key: fasta number, sequence_name, matches
parts_trimmed_list = []
for part_name in ordered_part_list_dict[library_number]:
	temp_part = ordered_part_dict[library_number][part_name][3]
	if temp_part != 'variable_sequence_repeat':
		part_trimmed = ordered_part_dict[library_number][part_name][1]
		part_counts_trimmed[part_trimmed] = (ordered_part_dict[library_number][part_name][0], [])
		parts_trimmed_list.append(part_trimmed)
for seq_record in SeqIO.parse(fastq_reads, "fastq"):
    total_seqs += 1
    if str(seq_record.seq)[0] == 'T':      
        if extension_nuc == 'A':
            seq_trimmed = re.split('AAAAAAAA', str(seq_record.seq)[1:])
        elif extension_nuc == 'C':
             seq_trimmed = re.split('CCCCC', str(seq_record.seq)[1:])
        elif extension_nuc == 'G':
             seq_trimmed = re.split('GGGGG', str(seq_record.seq)[1:])
        if len(seq_trimmed)>1: 
            # test_seqs2.append(seq_trimmed[0])
            rev_seq = Seq(seq_trimmed[0]).reverse_complement()
            try:
                part_counts_trimmed[rev_seq][1].append(seq_record)
                match_count += 1
            except KeyError:
                no_lib_match.append(str(rev_seq)) 
        else: no_post_flank.append(seq_record)    
    else: no_first_T.append(seq_record)
c = Counter(no_lib_match)

"""Output"""
#write excel
workbook = xlsxwriter.Workbook('%s/%s_Results/msd.xlsx' % (Data_Path,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'fasta #')
worksheet.write(0,1,'variable sequence')
worksheet.write(0,2,'total counts')
#add data
row = 1
col = 0
for part in parts_trimmed_list:
    worksheet.write(row,col,part_counts_trimmed[part][0], bold)
    worksheet.write(row,col+1,part)
    worksheet.write(row,col+2,len(part_counts_trimmed[part][1]))
    row += 1
for seq in c.most_common():
    if seq[1] >= 10:
        worksheet.write(row,col+1,seq[0])
        worksheet.write(row,col+2,seq[1])
        row += 1
workbook.close()

print '%s total sequences' % total_seqs
print '%s matches found' % match_count
print '%s do not contain an initial T' % len(no_first_T)
print '%s do not have an appropriate nucleotide extension' % len(no_post_flank)
print '%s were trimmed, but have no match in the target set' % len(no_lib_match)
