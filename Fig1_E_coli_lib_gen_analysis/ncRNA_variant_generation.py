"""Import Modules"""
import sys,os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pyperclip
import xlsxwriter
import random
import itertools
import difflib
import pickle

"""Arguments"""

"""Globals"""
user_profile = os.environ ['USERPROFILE']
wd = os.getcwd()
primer_dict = {}
nts = ['G','A','T','C']
ligation_sites_dict = {'86_1L':'TGAAT', #adding base annotated as msd, but not actually in msd
					 '86_2L':'CTTA',
					 '86_1R':'GTAA',
					 '86_2R':'ACTT',
					 '86_3R':'CTGG',
					 '83_2L':'TTTA',
					 '83_3R':'TTGC',
					 '86_3L': 'CCAG'}
BsaI_L = 'GGTCTCA'
BsaI_R = 'AGAGACC'
ec86_wt_msd = 'CTGAGTTACTGTCTGTTTTCCTTGTTGGAACGGAGAGCATCGCCTGATGCTCTCCGAGCCAACCAGGAAACCCGTTTTTTCTGAC' #actually 85 bases, consistent with sequencing, in RNA direction
ec86_wt_msr = 'ATGCGCACCCTTAGCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACA' #note: wt msr is actually b3v2 version (diff last base)
ec86_wt_msr_noprimer = 'GCGAGAGGTTTATCATTAAGGTCAACCTCTGGATGTTGTTTCGGCATCCTGCATTGAATCTGAGTTACA' #note: wt msr is actually b3v2 version (diff last base)
ec86_b3v2_frag_left = 'AGCTGTTTGTCGCCAG'
ec86_b3v2_frag_right = 'CTGGCGACAACCCGTTTTTTCTGAC'
fasta_number = 1
part_dict = {}
part_list_dict = {}
library_list = []
synthesis_list = []
synthesis_lists_sub = {}
expt_dict = {}
record_dict_bcs = SeqIO.index("%s/9base_BCs.fasta" % wd, "fasta")
barcode_list = []
for key in record_dict_bcs:
	barcode_list.append(str(record_dict_bcs[key].seq))
rule_name = ['sw','sw_comp','pp','pp_comp'] #strong/weak, purine/pyrimidine, comp means compensating mutations on the other side of the stem to keep structure
rule_dicts = [{'G':'C', 'A':'T', 'T':'A', 'C':'G'},{'G':'C', 'A':'T', 'T':'A', 'C':'G'},{'G':'A', 'A':'G', 'T':'C', 'C':'T'},{'G':'A', 'A':'G', 'T':'C', 'C':'T'}] #strong/weak, purine/pyrimidine


"""Defs"""
def make_sides(retron,primer,lig_sites):
	"""input is the amplication primer name and four base overlap for ligation sites
		output is a list of left and right sides"""
	left = str(primer_dict[primer][0].seq)+BsaI_L+ligation_sites_dict[retron+'_'+lig_sites[0]]
	right = ligation_sites_dict[retron+'_'+lig_sites[1]]+BsaI_R+str(primer_dict[primer][1].seq.reverse_complement())
	return (left,right)

def assemble_parts(sides,variable_part,library,fasta_number):
	"""input is the variable inside and outside parts of the sequence
		output is the assembled part with random nucleotides added 3' up to 150 bases"""
	synthesis_part = sides[0]+variable_part+sides[1]
	if len(synthesis_part)>170:
		print "length error, library %s" % library
	rando_end = ''.join(random.choice(nts) for i in range(170-len(synthesis_part)))
	synthesis_whole = synthesis_part+rando_end
	synthesis_lists_sub[library].append(SeqRecord(Seq(synthesis_whole),id=str(fasta_number)))
	synthesis_list.append(SeqRecord(Seq(synthesis_whole),id=str(fasta_number),description=''))
	return synthesis_whole

def add_library_to_dicts(library):
	"""input is the library name
		adds it to the various dictionaries"""	
	part_list_dict[library] = []
	part_dict[library] = {}
	synthesis_lists_sub[library] = []

def add_wt_msd(library,part_prefix,sides,retron,fasta_number):
	part_list_dict[library].append(part_prefix+'_wt')
	synth_part = assemble_parts(sides,retron,library,fasta_number)
	part_dict[library][part_prefix+'_wt'] = (fasta_number,retron,'unaltered',synth_part)

def add_wt_msr(library,part_prefix,sides,msr_full,msr_noprimer,frag_left,frag_right,fasta_number,retron):
	part_list_dict[library].append(part_prefix+'_wt')
	bc = random.choice(barcode_list)
	barcode_list.remove(bc)
	variable_part = msr_noprimer+frag_left+bc
	if retron == '86':
		msd = variable_part[59:]+frag_right
		msr = msr_full[:13]+variable_part[:69]			
	synth_part = assemble_parts(sides,variable_part,library,fasta_number)
	part_dict[library][part_prefix+'_wt'] = (fasta_number,msd,msr,synth_part)
	return variable_part

"""Run"""
#import the amplification primer sets
record_dict_forward = SeqIO.index("%s/skpp15-forward.fasta" % wd, "fasta")
record_dict_reverse = SeqIO.index("%s/skpp15-reverse.fasta" % wd, "fasta")
for key in record_dict_forward:
	primer_dict[key[:-2]] = (record_dict_forward[key],record_dict_reverse[key[:-2]+'-R'])

#################
### Library 1 ###
#################
library = 'L1'
library_list.append(library) #library name
retron = '86'
lig_sites = ('1L','1R')
primer = 'skpp15-25'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L1'
expt = '86--msd:change every base of wt msd to every other base (one at a time) + wild-type'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
var_fasta[ec86_wt_msd] = add_wt_msd(library,part_prefix,sides,ec86_wt_msd,fasta_number) #all libraries should have wt as a reference
fasta_number+=1
#variable part
for c,base in enumerate(ec86_wt_msd):
	avail_nts = ['G','A','T','C']
	avail_nts.remove(base)
	for new_base in avail_nts:
		part_name = part_prefix+'_'+str(c+1)+new_base
		variable_part = ec86_wt_msd[:c]+new_base+ec86_wt_msd[c+1:]
		msd = variable_part #in this case variable part is the msd
		msr = 'unaltered' #no changes to msr
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)	

#################
### Library 3 ###
#################
library = 'L3'
library_list.append(library) #library name
retron = '86'
lig_sites = ('1L','1R')
primer = 'skpp15-27'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L3'
expt = '86--msd:delete bases one at a time, two at a time, three at a time, four at a time, and five at a time (careful to avoid ambiguities) + wild-type'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
var_fasta[ec86_wt_msd] = add_wt_msd(library,part_prefix,sides,ec86_wt_msd,fasta_number) #all libraries should have wt as a reference
fasta_number+=1
#variable part
for length in range(1,6):
	for i in range(0,len(ec86_wt_msd)-length+1):
		part_name = part_prefix+'_'+str(i+1)+'_d'+str(length)
		variable_part = ec86_wt_msd[:i]+ec86_wt_msd[i+length:]
		msd = variable_part #in this case variable part is the msd
		msr = 'unaltered' #no changes to msr
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1		
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)	

#################
### Library 4 ###
#################
library = 'L4'
library_list.append(library) #library name
retron = '86'
lig_sites = ('1L','1R')
primer = 'skpp15-28'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L4'
expt = '86--msd:insert bases, one at a time, adding each nucleotide (careful to avoid duplicates), three at a time and five at a time with with three different sets each + wild-type'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
var_fasta[ec86_wt_msd] = add_wt_msd(library,part_prefix,sides,ec86_wt_msd,fasta_number) #all libraries should have wt as a reference
fasta_number+=1
#variable part
for addition in ['G','A','T','C','CAG','ACA','CTG','TTGAG','TCGGT','ATGTC']:
	for i in range(0,len(ec86_wt_msd)+1):
		part_name = part_prefix+'_'+str(i+1)+'_'+addition
		variable_part = ec86_wt_msd[:i]+addition+ec86_wt_msd[i:]
		msd = variable_part #in this case variable part is the msd
		msr = 'unaltered' #no changes to msr
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1		
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)	

#################
### Library 5 ###
#################
library = 'L5'
library_list.append(library) #library name
retron = '86'
lig_sites = ('1L','1R')
primer = 'skpp15-29'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L5'
expt = '86--msd:break the hairpin and make compensatory changes in a sliding window of 5 bases + wild-type (keep strong/weak (sw: A-T --> T-A) or purine/pyrimidine (pp: A-T --> G-C))'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
var_fasta[ec86_wt_msd] = add_wt_msd(library,part_prefix,sides,ec86_wt_msd,fasta_number) #all libraries should have wt as a reference
fasta_number+=1
#variable part
for rule, swap in zip(rule_name, rule_dicts):
	for i in range(0,25-4): #stem is 25 bases
		part_name = part_prefix+'_'+str(i+1)+'_'+rule
		hairpin_start = 16 #actualy the nucleotide before the first hairpin nuc
		hairpin_end = 70 #actually the nucleotide before the last hairpin nuc
		hairpin_length = 25
		variable_list = [ec86_wt_msd[0:hairpin_start+i]] #up to changes in stem
		for nuc in range(0,5): #range of 5 nucleotides
			variable_list.append(swap[ec86_wt_msd[hairpin_start+i+nuc]]) #swapped bases on rising stem
		variable_list.append(ec86_wt_msd[hairpin_start+i+5:hairpin_end-i-5]) #from first swap to second swap
		if rule[3:] == 'comp':
			for nuc in reversed(range(1,6)): #range of 5 nucleotides
				variable_list.append(swap[ec86_wt_msd[hairpin_end-i-nuc]]) #swapped bases on falling stem
		else:
			for nuc in reversed(range(1,6)): #range of 5 nucleotides
				variable_list.append(ec86_wt_msd[hairpin_end-i-nuc]) #keep original bases on falling stem (break hairpin)
		variable_list.append(ec86_wt_msd[hairpin_end-i:])
		variable_part = ''.join(variable_list)
		# variable_part = ec86_wt_msd[0:16+i]+swap[ec86_wt_msd[16+i]]+swap[ec86_wt_msd[17+i]]+swap[ec86_wt_msd[18+i]]+swap[ec86_wt_msd[19+i]]+swap[ec86_wt_msd[20+i]]+ec86_wt_msd[21+i:65-i]+swap[ec86_wt_msd[65-i]]+swap[ec86_wt_msd[66-i]]+swap[ec86_wt_msd[67-i]]+swap[ec86_wt_msd[68-i]]+swap[ec86_wt_msd[69-i]]+ec86_wt_msd[70-i:]
		msd = variable_part #in this case variable part is the msd
		msr = 'unaltered' #no changes to msr
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1		
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)

#################
### Library 6 ###
#################
library = 'L6'
library_list.append(library) #library name
retron = '86'
lig_sites = ('2L','3R')
primer = 'skpp15-30'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L6'
expt = '86--msr: change every base to every other base + wildtype'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
wt_variable_seq = add_wt_msr(library,part_prefix,sides,ec86_wt_msr,ec86_wt_msr_noprimer,ec86_b3v2_frag_left,ec86_b3v2_frag_right,fasta_number,retron)
var_fasta[wt_variable_seq] = fasta_number
fasta_number+=1
#variable part
for c,base in enumerate(ec86_wt_msr_noprimer[0:58]):
	avail_nts = ['G','A','T','C']
	avail_nts.remove(base)	
	for new_base in avail_nts:
		part_name = part_prefix+'_'+str(c+1)+new_base
		bc = random.choice(barcode_list)
		barcode_list.remove(bc)
		variable_part = ec86_wt_msr_noprimer[:c]+new_base+ec86_wt_msr_noprimer[c+1:]+ec86_b3v2_frag_left+bc
		msd = variable_part[59:]+ec86_b3v2_frag_right
		msr = ec86_wt_msr[:13]+variable_part[:69]
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)	

#################
### Library 7 ###
#################
library = 'L7'
library_list.append(library) #library name
retron = '86'
lig_sites = ('2L','3R')
primer = 'skpp15-31'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L7'
expt = '86--msr: delete bases one at a time, two at a time, three at a time, four at a time, and five at a time (careful to avoid ambiguities) + wild-type'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
wt_variable_seq = add_wt_msr(library,part_prefix,sides,ec86_wt_msr,ec86_wt_msr_noprimer,ec86_b3v2_frag_left,ec86_b3v2_frag_right,fasta_number,retron)
var_fasta[wt_variable_seq] = fasta_number
fasta_number+=1
#variable part
for length in range(1,6):
	for i in range(0,59-length): #does not delete duplexed msr/msd bases because it could interfere with PCRing barcodes
		part_name = part_prefix+'_'+str(i+1)+'_d'+str(length)
		bc = random.choice(barcode_list)
		barcode_list.remove(bc)
		variable_part = ec86_wt_msr_noprimer[:i]+ec86_wt_msr_noprimer[i+length:]+ec86_b3v2_frag_left+bc
		msd = variable_part[59-length:]+ec86_b3v2_frag_right
		msr = ec86_wt_msr[:13]+variable_part[:69-length]
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)	

#################
### Library 8 ###
#################
library = 'L8'
library_list.append(library) #library name
retron = '86'
lig_sites = ('2L','3R')
primer = 'skpp15-32'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L8'
expt = '86--msr: insert bases, one at a time, adding each nucleotide (careful to avoid duplicates), three at a time and five at a time with with three different sets each + wild-type'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
wt_variable_seq = add_wt_msr(library,part_prefix,sides,ec86_wt_msr,ec86_wt_msr_noprimer,ec86_b3v2_frag_left,ec86_b3v2_frag_right,fasta_number,retron)
var_fasta[wt_variable_seq] = fasta_number
fasta_number+=1
#variable part
for addition in ['G','A','T','C','CAG','ACA','CTG']:
	for i in range(0,59-length): #does not delete duplexed msr/msd bases because it could interfere with PCRing barcodes
		part_name = part_prefix+'_'+str(i+1)+'_'+addition
		bc = random.choice(barcode_list)
		barcode_list.remove(bc)
		variable_part = ec86_wt_msr_noprimer[:i]+addition+ec86_wt_msr_noprimer[i:]+ec86_b3v2_frag_left+bc
		msd = variable_part[59+len(addition):]+ec86_b3v2_frag_right
		msr = ec86_wt_msr[:13]+variable_part[:69+len(addition)]
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)

#################
### Library 9 ###
#################
library = 'L9'
library_list.append(library) #library name
retron = '86'
lig_sites = ('2L','3R')
primer = 'skpp15-33'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L9'
expt = '86--msr: all permutations to the three loop bases in the critical hairpin (wild-type sequencing being one of them), lengthen and shorten elongate the critical hairpin'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
wt_variable_seq = add_wt_msr(library,part_prefix,sides,ec86_wt_msr,ec86_wt_msr_noprimer,ec86_b3v2_frag_left,ec86_b3v2_frag_right,fasta_number,retron)
var_fasta[wt_variable_seq] = fasta_number
fasta_number+=1
#variable part
for perm in itertools.product('GATC',repeat=3): #this is the loop changing part
	perm_str = ''.join(perm)
	part_name = part_prefix+'_loop_'+perm_str
	bc = random.choice(barcode_list)
	barcode_list.remove(bc)
	variable_part = ec86_wt_msr_noprimer[:38]+perm_str+ec86_wt_msr_noprimer[41:]+ec86_b3v2_frag_left+bc
	msd = variable_part[59:]+ec86_b3v2_frag_right
	msr = ec86_wt_msr[:13]+variable_part[:69]
	part_list_dict[library].append(part_name)
	try: 
		int_fasta = var_fasta[variable_part] #indicates duplicate part
		part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
	except KeyError:
		var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
		synth_part = assemble_parts(sides,variable_part,library,fasta_number)
		part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
		fasta_number+=1	
for i in range(-4,0): #make shorter versions
	part_name = part_prefix+'_stem_d'+str(i)
	bc = random.choice(barcode_list)
	barcode_list.remove(bc)
	variable_part = ec86_wt_msr_noprimer[:38+i]+'TTT'+ec86_wt_msr_noprimer[41-i:]+ec86_b3v2_frag_left+bc
	msd = variable_part[59+i+i:]+ec86_b3v2_frag_right
	msr = ec86_wt_msr[:13]+variable_part[:69+i+i]
	part_list_dict[library].append(part_name)
	try: 
		int_fasta = var_fasta[variable_part] #indicates duplicate part
		part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
	except KeyError:
		var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
		synth_part = assemble_parts(sides,variable_part,library,fasta_number)
		part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
		fasta_number+=1
for stem_nuc in [['G','C'],['GC','GC']]: #make longer versions, add GC to stem
	part_name = part_prefix+'_stem_+'+stem_nuc[0]
	bc = random.choice(barcode_list)
	barcode_list.remove(bc)
	variable_part = ec86_wt_msr_noprimer[:38]+stem_nuc[0]+'TTT'+stem_nuc[1]+ec86_wt_msr_noprimer[41:]+ec86_b3v2_frag_left+bc
	msd = variable_part[59+len(stem_nuc[0])+len(stem_nuc[1]):]+ec86_b3v2_frag_right
	msr = ec86_wt_msr[:13]+variable_part[:69+len(stem_nuc[0])+len(stem_nuc[1])]
	part_list_dict[library].append(part_name)
	try: 
		int_fasta = var_fasta[variable_part] #indicates duplicate part
		part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
	except KeyError:
		var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
		synth_part = assemble_parts(sides,variable_part,library,fasta_number)
		part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
		fasta_number+=1	
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)
	
#################
### Library 10 ##
#################
library = 'L10'
library_list.append(library) #library name
retron = '86'
lig_sites = ('2L','3R')
primer = 'skpp15-34'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L10'
expt = '86--msr: break the critical hairpin and make compensatory changes in a sliding window of 4 bases + wild-type'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
wt_variable_seq = add_wt_msr(library,part_prefix,sides,ec86_wt_msr,ec86_wt_msr_noprimer,ec86_b3v2_frag_left,ec86_b3v2_frag_right,fasta_number,retron)
var_fasta[wt_variable_seq] = fasta_number
fasta_number+=1
#variable part
for rule, swap in zip(rule_name, rule_dicts):
	for i in range(0,8-3): #stem is 8 bases
		part_name = part_prefix+'_'+str(i+1)+'_'+rule
		bc = random.choice(barcode_list)
		barcode_list.remove(bc)
		hairpin_start = 30 #actualy the nucleotide before the first hairpin nuc
		hairpin_end = 49 #actually the nucleotide before the last hairpin nuc?
		hairpin_length = 8
		variable_list = [ec86_wt_msr_noprimer[0:hairpin_start+i]]
		for nuc in range(0,4): #range of 4 nucleotides
			variable_list.append(swap[ec86_wt_msr_noprimer[hairpin_start+i+nuc]]) #swapped bases on rising stem
		variable_list.append(ec86_wt_msr_noprimer[hairpin_start+i+4:hairpin_end-i-4]) #from first swap to second swap
		if rule[3:] == 'comp':
			for nuc in reversed(range(1,5)): #range of 4 nucleotides
				variable_list.append(swap[ec86_wt_msr_noprimer[hairpin_end-i-nuc]]) #swapped bases on falling stem
		else:
			for nuc in reversed(range(1,5)): #range of 4 nucleotides
				variable_list.append(ec86_wt_msr_noprimer[hairpin_end-i-nuc]) #keep original bases on falling stem (break hairpin)
		variable_list.append(ec86_wt_msr_noprimer[hairpin_end-i:])
		variable_frag = ''.join(variable_list)
		variable_part = variable_frag+ec86_b3v2_frag_left+bc	
		msd = variable_part[59:]+ec86_b3v2_frag_right
		msr = ec86_wt_msr[:13]+variable_part[:69]
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)				

#################
### Library 11 ##
#################
library = 'L11'
library_list.append(library) #library name
retron = '86'
lig_sites = ('2L','3R')
primer = 'skpp15-35'
sequence_numbers = [fasta_number] #first, last 
part_prefix = '86_r2_L11'
expt = '86--msr: break the other hairpin and make compensatory changes in a sliding window of 4 bases + wild-type (region 2)'
sides = make_sides(retron,primer,lig_sites)
add_library_to_dicts(library)
var_fasta = {} #this is a dictionary of all the variable parts in this library with their fasta numbers, to avoid duplicate parts within a sublibrary
wt_variable_seq = add_wt_msr(library,part_prefix,sides,ec86_wt_msr,ec86_wt_msr_noprimer,ec86_b3v2_frag_left,ec86_b3v2_frag_right,fasta_number,retron)
var_fasta[wt_variable_seq] = fasta_number
fasta_number+=1
#variable part
for rule, swap in zip(rule_name, rule_dicts):
	for i in range(0,11-3): #stem is 8 bases
		part_name = part_prefix+'_'+str(i+1)+'_'+rule
		bc = random.choice(barcode_list)
		barcode_list.remove(bc)
		hairpin_start = 3 #actualy the nucleotide before the first hairpin nuc
		hairpin_end = 30 #actually the nucleotide before the last hairpin nuc?
		hairpin_length = 11
		variable_list = [ec86_wt_msr_noprimer[0:hairpin_start+i]]
		for nuc in range(0,4): #range of 4 nucleotides
			variable_list.append(swap[ec86_wt_msr_noprimer[hairpin_start+i+nuc]]) #swapped bases on rising stem
		variable_list.append(ec86_wt_msr_noprimer[hairpin_start+i+4:hairpin_end-i-4]) #from first swap to second swap
		if rule[3:] == 'comp':
			for nuc in reversed(range(1,5)): #range of 4 nucleotides
				variable_list.append(swap[ec86_wt_msr_noprimer[hairpin_end-i-nuc]]) #swapped bases on falling stem
		else:
			for nuc in reversed(range(1,5)): #range of 4 nucleotides
				variable_list.append(ec86_wt_msr_noprimer[hairpin_end-i-nuc]) #keep original bases on falling stem (break hairpin)
		variable_list.append(ec86_wt_msr_noprimer[hairpin_end-i:])
		variable_frag = ''.join(variable_list)
		variable_part = variable_frag+ec86_b3v2_frag_left+bc	
		msd = variable_part[59:]+ec86_b3v2_frag_right
		msr = ec86_wt_msr[:13]+variable_part[:69]
		part_list_dict[library].append(part_name)
		try: 
			int_fasta = var_fasta[variable_part] #indicates duplicate part
			part_dict[library][part_name] = (var_fasta[variable_part],msd,msr,'variable_sequence_repeat')
		except KeyError:
			var_fasta[variable_part] = fasta_number #KeyError if no duplicate yet
			synth_part = assemble_parts(sides,variable_part,library,fasta_number)
			part_dict[library][part_name] = (fasta_number,msd,msr,synth_part)
			fasta_number+=1
sequence_numbers.append(fasta_number-1) #first, last 
expt_dict[library] = (retron,lig_sites,primer,sequence_numbers,expt)			

"""Output"""
#write fasta
SeqIO.write(synthesis_list,"%s/ncRNA_variant_generation.fasta" % wd, "fasta") #no duplicates within libraries
#export dictionaries
with open("%s/part_list_dict_synth.p" % wd, 'wb') as fp:
	pickle.dump(part_list_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
with open("%s/part_dict_synth.p" % wd, 'wb') as fp:
	pickle.dump(part_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
#write excel with sublibraries and descriptions
workbook = xlsxwriter.Workbook('ncRNA_variant_generation.xlsx')
bold = workbook.add_format({'bold': True})
monospace = workbook.add_format({'font_name': 'Lucida Sans Typewriter'})
for key in expt_dict: #key is library
	worksheet = workbook.add_worksheet(expt_dict[key][2])
	#Add titles (row, col: zero referenced)
	worksheet.write(0,0,'Retron:')
	worksheet.write(0,1,expt_dict[key][0])
	worksheet.write(1,0,'ligation sites:')
	worksheet.write(1,1,str(expt_dict[key][1]))
	worksheet.write(2,0,'primer:')
	worksheet.write(2,1,expt_dict[key][2])
	worksheet.write(3,0,'sequence_numbers:')
	worksheet.write(3,1,str(expt_dict[key][3]))
	worksheet.write(4,0,'Experiment:')
	worksheet.write(4,1,expt_dict[key][4],bold)
	worksheet.write(6,0,'fasta id:', bold)
	worksheet.write(6,1,'part name:',bold)
	worksheet.write(6,2,'msd:',bold)
	worksheet.write(6,3,'msr:',bold)
	worksheet.write(6,4,'synthesized part:',bold)
	row = 7
	col = 0
	for part in part_list_dict[key]:
		worksheet.write(row,col,part_dict[key][part][0]) #fasta
		worksheet.write(row,col+1,part) #part name
		worksheet.write(row,col+2,part_dict[key][part][1],monospace) #msd
		worksheet.write(row,col+3,part_dict[key][part][2],monospace) #msr
		worksheet.write(row,col+4,part_dict[key][part][3],monospace) #synth part
		row +=1
workbook.close()
#write excel with each synthesized part in a row (for Agilent)
workbook2 = xlsxwriter.Workbook('ncRNA_variant_synth.xlsx')
worksheet = workbook2.add_worksheet()
row = 0
col = 0
for part in synthesis_list:
	worksheet.write(row,col,str(part.seq))
	row +=1
workbook2.close()	