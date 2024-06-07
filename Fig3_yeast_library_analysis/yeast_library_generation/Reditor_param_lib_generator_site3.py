"""Import Modules"""
import sys,os
import openpyxl
from openpyxl import load_workbook 
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
from Bio import SeqUtils
from Bio.Data import CodonTable
import itertools
from collections import Counter
import xlsxwriter
import random

"""Arguments"""
library_length = 250
first_fasta = 12000 #starting fasta number
skpp = 'skpp15-4' #skpp set to use

"""Globals"""
user_profile = os.environ ['USERPROFILE']
wd = os.getcwd()
nts = ['G','A','T','C']
PAM_positions = [102, 110, 118, 126, 134]
#site1
site = 'AGAGCCCCCCTGTAAGTCAGCACTCCTCAGAGGATCCCACCCGCCTTCTGCTTTCCTGTTGGGACCCTGATTCACACACCTCGTGACTGTCTACGCCACGTAGGACACATGGTTATCTGGGGTAAAGGCTGTACGGCATCCCTCCATTGTTCACTCTGGGGGTTCATCACCTCAACCATTATGTCTCAGCTCCAAACAAAATCTCCTCTGCTTTCATTTAAGCTCCTT'
#site gRNAs as dictionary where the key is the position of the PAM GG
site_gRNAs= {102:['TGACTGTCTACGCCACGT','PAM_-13'], 
			 110:['TACGCCACGTAGGACACA','PAM_-5'], 
			 118:['GTAGGACACATGGTTATC','PAM_-3'], 
			 126:['CATGGTTATCTGGGGTAA','PAM_11'], 
			 134:['TCTGGGGTAAAGGCTGTA','PAM_19']}
#different chassis, split into two around the donor, except for 22-25, which are split into three (predonor, extra left stem, and postdonor)
chassis= {1:['TGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCCGTTTTTTCTGACGTAAGGGTGCGCA'], #1_ec86wt msr/msd (RNA) bare_true_wt
		  2:['TGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA'], #2_ec86wt msr/msd (RNA) bare_CRISPEY_wt
		  3:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #3_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27
		  4:['AAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTT'], #4_ec86wt msr/msd (RNA) bare_CRISPEY_Pext23
		  5:['TTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAA'], #5_ec86wt msr/msd (RNA) bare_CRISPEY_Pext19
		  6:['CGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACG'], #6_ec86wt msr/msd (RNA) bare_CRISPEY_Pext16
		  7:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','GGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #7_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_d136
		  8:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #8_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_d139
		  9:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAACCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #9_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_d141
		  10:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCTGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #10_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_144T_old11
		  11:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACCCGTATCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #11_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_147A_old12
		  12:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','TATCATGTAATACTCTACGGCGTAAGGGTGCGCATACGGAATCTTATCA'], #12_compN0_328
		  13:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','ATATCGGCGTAATTTCGGTTCGTAAGGGTGCGCATACGGAATCTTATCA'], #13_compN0_64294
		  14:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','ATCTATCGTGATATTCAGATCGTAAGGGTGCGCATACGGAATCTTATCA'], #14_compN0_185887
		  15:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','ATATGGACACGGTGAGCACTCGTAAGGGTGCGCATACGGAATCTTATCA'], #15_compN0_232445
		  16:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','ATCATTTACGGTACGATCAGCGTAAGGGTGCGCATACGGAATCTTATCA'], #16_compN0_101299
		  17:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAATCCAACATTCACTGCGTAAGGGTGCGCATACGGAATCTTATCA'], #17_compN1_446844
		  18:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','GGATCACTCACTGTATACCGCGTAAGGGTGCGCATACGGAATCTTATCA'], #18_compN1_40329
		  19:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACTGTAAATAAACCGCGTAAGGGTGCGCATACGGAATCTTATCA'], #19_compN1_461043
		  20:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','GCGGGCGGTGTAACGTTCCACGTAAGGGTGCGCATACGGAATCTTATCA'], #20_compN1_25575
		  21:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','AGGAAACTGTAGATAACCCGCGTAAGGGTGCGCATACGGAATCTTATCA'], #21_compN1_470993
		  22:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','TG','CCAGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #22_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_stem8_old13
		  23:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','TGTT','AACCAGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #23_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27 stem10_old14
		  24:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','TGTTGG','CCAACCAGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA'], #24_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_stem12_old15
		  25:['TGATAAGATTCCGTATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACTGTCTGTTTTCCT','TGTTGGAA','AGCCAACCAGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCA']} #25_ec86wt msr/msd (RNA) bare_CRISPEY_Pext27_stem14_old16

#import the amplification primer sets
primer_dict = {}
record_dict_forward = SeqIO.index("%s/skpp15-forward.fasta" % wd, "fasta")
record_dict_reverse = SeqIO.index("%s/skpp15-reverse.fasta" % wd, "fasta")
for key in record_dict_forward:
	primer_dict[key[:-2]] = (record_dict_forward[key],record_dict_reverse[key[:-2]+'-R'])

record_dict_bcs = SeqIO.index("%s/10base_BCs_hamming5_GC_40_70_4501.fasta" % wd, "fasta")
barcode_list = []
for key in record_dict_bcs:
	barcode_list.append(str(record_dict_bcs[key].seq))
barcodes_used = []

part_dict = {}
part_fasta_list = []

"""Defs"""
def grab_a_barcode():
	bc = random.choice(barcode_list)
	barcode_list.remove(bc)
	barcodes_used.append(bc)	
	return bc
def make_part(barcode, donor, gRNA, chassis_number):
	primer_f = str(primer_dict[skpp][0].seq)
	primer_r = str(primer_dict[skpp][1].seq.reverse_complement())
	if len(chassis[chassis_number]) == 2:
		synth_part = primer_f+'CACCTGCATCATCCT'+donor+chassis[chassis_number][1]+gRNA+'GTTTATCAGCAGGTG'+primer_r
		ncRNA = chassis[chassis_number][0]+donor+chassis[chassis_number][1]+gRNA+'GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTT'
	elif len(chassis[chassis_number]) == 3:
		synth_part = primer_f+'CACCTGCATCATCCT'+chassis[chassis_number][1]+donor+chassis[chassis_number][2]+gRNA+'GTTTATCAGCAGGTG'+primer_r
		ncRNA = chassis[chassis_number][0]+chassis[chassis_number][1]+donor+chassis[chassis_number][2]+gRNA+'GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTT'
	if len(synth_part)>library_length:
		print "length error, gRNA: %s, chassis_number: %s" % (gRNA, chassis_number)	
	rando_end = ''.join(random.choice(nts) for i in range(library_length-len(synth_part)))
	synthesis_whole = synth_part+rando_end
	return ncRNA, synthesis_whole

#################
###### Run ######
#################
fasta = first_fasta
#generate gRNA/donor combos
for PAM_pos in PAM_positions:
	recoded_site = site[:PAM_pos]+'CT'+site[PAM_pos+2:]
	PAM_pos_int = int(site_gRNAs[PAM_pos][1].split('_')[1])
	
	#make 112 length donors
	for i in [-28,-14,0,19,38]: #center of the donors
		DONOR_name = 'DONORf_length_112_center_%s' % i
		for ch in range(1,22): #just 21 for the longest donor, use 1,26 for all the others
			part_name = '%s_%s_chassis_%s' % (site_gRNAs[PAM_pos][1],DONOR_name,ch)
			bc = grab_a_barcode()
			barcoded_recoded_site = recoded_site[:114]+bc+recoded_site[114:]
			donor = barcoded_recoded_site[114-5-56+i:114+5+56+i] #the 56 is half 112, change for other length donors
			parts = make_part(bc, donor, site_gRNAs[PAM_pos][0], ch)
			part_dict[fasta] = [part_name,bc,PAM_pos_int,112,'F',i,ch,site_gRNAs[PAM_pos][0],donor,parts[0],parts[1]]
			#part_dict[fasta] = part_name, barcode, PAM position, donor length, donor direction, donor center, chassis number, gRNA, donor, ncRNA (with gRNA and tracr), synthesis part
			fasta += 1

	#make 94 length donors
	for i in [-20,-10,0,14,28]: #center of the donors
		#make forward donors
		DONOR_name = 'DONORf_length_94_center_%s' % i
		for ch in range(1,26): #just 21 for the longest donor, use 1,26 for all the others
			part_name = '%s_%s_chassis_%s' % (site_gRNAs[PAM_pos][1],DONOR_name,ch)
			bc = grab_a_barcode()
			barcoded_recoded_site = recoded_site[:114]+bc+recoded_site[114:]
			donor = barcoded_recoded_site[114-5-47+i:114+5+47+i] #the 56 is half 112, change for other length donors
			parts = make_part(bc, donor, site_gRNAs[PAM_pos][0], ch)
			part_dict[fasta] = [part_name,bc,PAM_pos_int,94,'F',i,ch,site_gRNAs[PAM_pos][0],donor,parts[0],parts[1]]
			#part_dict[fasta] = part_name, barcode, PAM position, donor length, donor direction, donor center, chassis number, gRNA, donor, ncRNA (with gRNA and tracr), synthesis part
			fasta += 1
		#make reverse donors
		DONOR_name = 'DONORr_length_94_center_%s' % i
		for ch in range(1,26): #just 21 for the longest donor, use 1,26 for all the others
			part_name = '%s_%s_chassis_%s' % (site_gRNAs[PAM_pos][1],DONOR_name,ch)
			bc = grab_a_barcode()
			barcoded_recoded_site = recoded_site[:114]+bc+recoded_site[114:]
			donor = str(Seq(barcoded_recoded_site[114-5-47+i:114+5+47+i]).reverse_complement()) #the 56 is half 112, change for other length donors
			parts = make_part(bc, donor, site_gRNAs[PAM_pos][0], ch)
			part_dict[fasta] = [part_name,bc,PAM_pos_int,94,'R',i,ch,site_gRNAs[PAM_pos][0],donor,parts[0],parts[1]]
			#part_dict[fasta] = part_name, barcode, PAM position, donor length, donor direction, donor center, chassis number, gRNA, donor, ncRNA (with gRNA and tracr), synthesis part
			fasta += 1		

	#make 78 length donors
	for i in [-12,-6,0,10,20]: #center of the donors
		#make forward donors
		DONOR_name = 'DONORf_length_78_center_%s' % i
		for ch in range(1,26): #just 21 for the longest donor, use 1,26 for all the others
			part_name = '%s_%s_chassis_%s' % (site_gRNAs[PAM_pos][1],DONOR_name,ch)
			bc = grab_a_barcode()
			barcoded_recoded_site = recoded_site[:114]+bc+recoded_site[114:]
			donor = barcoded_recoded_site[114-5-39+i:114+5+39+i] #the 56 is half 112, change for other length donors
			parts = make_part(bc, donor, site_gRNAs[PAM_pos][0], ch)
			part_dict[fasta] = [part_name,bc,PAM_pos_int,78,'F',i,ch,site_gRNAs[PAM_pos][0],donor,parts[0],parts[1]]
			#part_dict[fasta] = part_name, barcode, PAM position, donor length, donor direction, donor center, chassis number, gRNA, donor, ncRNA (with gRNA and tracr), synthesis part
			fasta += 1

	#make 64 length donors
	for i in [-5,-2,0,7,14]: #center of the donors
		#make forward donors
		DONOR_name = 'DONORf_length_64_center_%s' % i
		for ch in range(1,26): #just 21 for the longest donor, use 1,26 for all the others
			part_name = '%s_%s_chassis_%s' % (site_gRNAs[PAM_pos][1],DONOR_name,ch)
			bc = grab_a_barcode()
			barcoded_recoded_site = recoded_site[:114]+bc+recoded_site[114:]
			donor = barcoded_recoded_site[114-5-32+i:114+5+32+i] #the 56 is half 112, change for other length donors
			parts = make_part(bc, donor, site_gRNAs[PAM_pos][0], ch)
			part_dict[fasta] = [part_name,bc,PAM_pos_int,64,'F',i,ch,site_gRNAs[PAM_pos][0],donor,parts[0],parts[1]]
			#part_dict[fasta] = part_name, barcode, PAM position, donor length, donor direction, donor center, chassis number, gRNA, donor, ncRNA (with gRNA and tracr), synthesis part
			fasta += 1
		#make reverse donors
		DONOR_name = 'DONORr_length_64_center_%s' % i
		for ch in range(1,26): #just 21 for the longest donor, use 1,26 for all the others
			part_name = '%s_%s_chassis_%s' % (site_gRNAs[PAM_pos][1],DONOR_name,ch)
			bc = grab_a_barcode()
			barcoded_recoded_site = recoded_site[:114]+bc+recoded_site[114:]
			donor = str(Seq(barcoded_recoded_site[114-5-32+i:114+5+32+i]).reverse_complement()) #the 56 is half 112, change for other length donors
			parts = make_part(bc, donor, site_gRNAs[PAM_pos][0], ch)
			part_dict[fasta] = [part_name,bc,PAM_pos_int,64,'R',i,ch,site_gRNAs[PAM_pos][0],donor,parts[0],parts[1]]
			#part_dict[fasta] = part_name, barcode, PAM position, donor length, donor direction, donor center, chassis number, gRNA, donor, ncRNA (with gRNA and tracr), synthesis part
			fasta += 1	

	#make 54 length donors
	for i in [-4,-2,0,5,10]: #center of the donors
		#make forward donors
		DONOR_name = 'DONORf_length_54_center_%s' % i
		for ch in range(1,26): #just 21 for the longest donor, use 1,26 for all the others
			part_name = '%s_%s_chassis_%s' % (site_gRNAs[PAM_pos][1],DONOR_name,ch)
			bc = grab_a_barcode()
			barcoded_recoded_site = recoded_site[:114]+bc+recoded_site[114:]
			donor = barcoded_recoded_site[114-5-27+i:114+5+27+i] #the 56 is half 112, change for other length donors
			parts = make_part(bc, donor, site_gRNAs[PAM_pos][0], ch)
			part_dict[fasta] = [part_name,bc,PAM_pos_int,54,'F',i,ch,site_gRNAs[PAM_pos][0],donor,parts[0],parts[1]]
			#part_dict[fasta] = part_name, barcode, PAM position, donor length, donor direction, donor center, chassis number, gRNA, donor, ncRNA (with gRNA and tracr), synthesis part
			fasta += 1

synthesis_list = []
for fasta_part in range(first_fasta,fasta):
	synthesis_list.append(SeqRecord(Seq(part_dict[fasta_part][10]),id=str(fasta_part),description=''))

"""Output"""
#save part dictionary to use for analysis
with open("%s/part_dictionary_site3.p" % wd, 'wb') as fp:
	pickle.dump(part_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)

#write fasta
SeqIO.write(synthesis_list,"%s/REp_site3_parts_synth.fasta" % wd, "fasta")

#write excel
workbook = xlsxwriter.Workbook('REditor_param_site3_summary.xlsx')
bold = workbook.add_format({'bold': True})
monospace = workbook.add_format({'font_name': 'Lucida Sans Typewriter'})
worksheet1 = workbook.add_worksheet('parts')
#Add titles (row, col: zero referenced)
worksheet1.write(0,0,'fasta', bold)
worksheet1.write(0,1,'part_name', bold)
worksheet1.write(0,2,'barcode', bold)
worksheet1.write(0,3,'PAM_position', bold)
worksheet1.write(0,4,'donor_length', bold)
worksheet1.write(0,5,'donor_direction', bold)
worksheet1.write(0,6,'donor_center', bold)
worksheet1.write(0,7,'chassis_number', bold)
worksheet1.write(0,8,'gRNA', bold)
worksheet1.write(0,9,'donor', bold)
worksheet1.write(0,10,'ncRNA (with gRNA and tracr)', bold)
worksheet1.write(0,11,'sythesis_part', bold)
#Add data
row = 1
col = 0
for fasta_part in range(first_fasta,fasta):
	worksheet1.write(row,col,fasta_part)
	worksheet1.write(row,col+1,part_dict[fasta_part][0])
	worksheet1.write(row,col+2,part_dict[fasta_part][1],monospace)
	worksheet1.write(row,col+3,part_dict[fasta_part][2])
	worksheet1.write(row,col+4,part_dict[fasta_part][3])
	worksheet1.write(row,col+5,part_dict[fasta_part][4])
	worksheet1.write(row,col+6,part_dict[fasta_part][5])
	worksheet1.write(row,col+7,part_dict[fasta_part][6])
	worksheet1.write(row,col+8,part_dict[fasta_part][7],monospace)
	worksheet1.write(row,col+9,part_dict[fasta_part][8],monospace)
	worksheet1.write(row,col+10,part_dict[fasta_part][9],monospace)
	worksheet1.write(row,col+11,part_dict[fasta_part][10],monospace)
	row += 1
workbook.close()

