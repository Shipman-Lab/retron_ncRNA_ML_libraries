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


"""Globals"""
wd = os.getcwd()
if running_location == 'local':
    user_profile = os.environ ['USERPROFILE']
    Data_Path = #to be set by user
    Blast_Data_Path = #to be set by user
fastq_reads = '%s/%s.fastq' % (Data_Path,condition)

#pull in counts using fasta sequences, then map them onto the experiments (to avoid duplicates within a library) using these:
with open("%s/part_list_dict_synth.p" % wd, 'rb') as fp:
    ordered_part_list_dict = pickle.load(fp)
with open("%s/part_dict_synth.p" % wd, 'rb') as fp:
    ordered_part_dict = pickle.load(fp)

trim86part_L = 'GTCGCCAG'
trim86part_R = 'CTGGCGAC'
if retron == '86':
    trim_parts = (trim86part_L,trim86part_R)


BC_dict = {} #loop BC is key: [0]fasta number, [1]total matches (seqs)
BC_list = []

total_seqs = 0
no_flanking_seqs = []
no_BC_match = []

"""Run"""
#Create Results folder
newpath = ((r'%s/%s_Results') % (Data_Path,condition))
if not os.path.exists(newpath): os.makedirs(newpath)
#pull in synthesis_parts
for part_name in ordered_part_list_dict[library_number]:
	temp_part = ordered_part_dict[library_number][part_name][1] #msd
	if temp_part != 'variable_sequence_repeat':
		loop_trimmed = re.search('%s(.*)%s' % (trim_parts[0],trim_parts[1]), temp_part)
		BC_dict[loop_trimmed.group(1)] = (ordered_part_dict[library_number][part_name][0], []) #adds fasta
		BC_list.append(loop_trimmed.group(1))
#pull in sequences
for seq_record in SeqIO.parse(fastq_reads, "fastq"):
    total_seqs += 1
    temp_record = seq_record.reverse_complement()
    loop_trimmed = re.search('%s(.*)%s' % (trim_parts[0],trim_parts[1]), str(temp_record.seq))
    if loop_trimmed:
        try:
            BC_dict[loop_trimmed.group(1)][1].append(str(temp_record))
        except KeyError:
            no_BC_match.append(loop_trimmed.group(1))
    else: no_flanking_seqs.append(temp_record)

print '%s total sequences' % total_seqs
print '%s do not contain matching loop flanking sequences' % len(no_flanking_seqs)
print '%s do not contain a matching barcode' % len(no_BC_match)

unmatched_BC_counter = Counter(no_BC_match)

"""Output"""
#write excel
workbook = xlsxwriter.Workbook('%s/%s_Results/Loop_Library_PCRd.xlsx' % (Data_Path,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'fasta #')
worksheet.write(0,1,'BC')
worksheet.write(0,2,'total counts')
#add data
row = 1
col = 0
for BC in BC_list:
    worksheet.write(row,col,BC_dict[BC][0], bold)
    worksheet.write(row,col+1,BC)
    worksheet.write(row,col+2,len(BC_dict[BC][1]))
    row += 1

for seq in unmatched_BC_counter.most_common():
    if seq[1] >= 10:
        worksheet.write(row,col+1,seq[0])
        worksheet.write(row,col+2,seq[1])
        row += 1
workbook.close()
