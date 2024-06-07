#get access to old libraries with pickle
import sys,os
import pickle
from Bio import SeqIO
import xlsxwriter
import fuzzysearch
import re
from collections import Counter

"""Line Arguments"""
run_number = sys.argv[1] #MiSeq run number for Data_Path
condition = sys.argv[2] #file name (minus .fastq) 
running_location = sys.argv[3] #either local or otherwise
retron = sys.argv[4] #86
library_number = sys.argv[5] #library, note: not currently set up for pooling
L_lig_site = sys.argv[6] #e.g. ec861L
R_lig_site = sys.argv[7]


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

ec86_Fow_Primer = 'TTATGCTAGGTGATGCAGCGGATTTCATGAAAG'
ec86_Rev_Primer_revcomp = 'GGTGCGCATACGGAATCTTATCATAGTTAAATG'
lig_site_dict = {'ec86_1L': 'GCATTGAA',
				 'ec86_1R': 'GTAAGGGT',
				 'ec86_2R': 'ACTTTCAT',
				 'ec86_3L': 'GTCGCCAG',
                 'ec86_3R': 'CTGGCGAC',
                 'ec86_2L': 'CACCCTTA'}

BsaI_L = 'GGTCTCA'
BsaI_R = 'AGAGACC'

seq_count = 0

match_count = 0
no_lig_site_seqs = []
no_synth_match_seqs = []
lig_sites = [lig_site_dict[L_lig_site], lig_site_dict[R_lig_site]]

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
		bsa_trimmed = re.search('%s(.*)%s' % (BsaI_L,BsaI_R), temp_part)
		part_counts_trimmed[bsa_trimmed.group(1)[4:-4]] = (ordered_part_dict[library_number][part_name][0], part_name, [])
		parts_trimmed_list.append(bsa_trimmed.group(1)[4:-4])
for seq_record in SeqIO.parse(fastq_reads, "fastq"):
    seq_count += 1
    temp_record = seq_record.reverse_complement()
    lig_trimmed = re.search('%s(.*)%s' % (lig_sites[0],lig_sites[1]), str(temp_record.seq))
    if lig_trimmed:
        try:
            part_counts_trimmed[lig_trimmed.group(1)][2].append(str(seq_record.seq)[:5])
            match_count += 1
        except KeyError:
            no_synth_match_seqs.append(lig_trimmed.group(1))
    else: no_lig_site_seqs.append(temp_record)
c = Counter(no_synth_match_seqs)

"""Output"""
#write excel
workbook = xlsxwriter.Workbook('%s/%s_Results/Plasmid_Library.xlsx' % (Data_Path,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'fasta #')
worksheet.write(0,1,'variable sequence')
worksheet.write(0,2,'total counts')
worksheet.write(0,3,'part_name')
#add data
row = 1
col = 0
for part in parts_trimmed_list:
    worksheet.write(row,col,part_counts_trimmed[part][0], bold)
    worksheet.write(row,col+1,part)
    worksheet.write(row,col+2,len(part_counts_trimmed[part][2]))
    worksheet.write(row,col+3,part_counts_trimmed[part][1])
    row += 1
for seq in c.most_common():
    if seq[1] >= 10:
        worksheet.write(row,col+1,seq[0])
        worksheet.write(row,col+2,seq[1])
        row += 1
workbook.close()

print 'Total sequences analyzed: %s' % seq_count
print 'Matched Sequences: %s' %match_count
print 'No ligation sites: %s' %len(no_lig_site_seqs)
print 'Ligation sites, but no match in synth library: %s' %len(no_synth_match_seqs)
