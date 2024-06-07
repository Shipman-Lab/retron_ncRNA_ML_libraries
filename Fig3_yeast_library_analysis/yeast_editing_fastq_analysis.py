# this is a script to process all fastqs in a subfolder named "data"
# it will save two outputs per fastq
#
# one output is fastq_name_map_dict.p
#   this pickle contains a dictionary with four numbers:
#       total number of reads processed
#       total number of reads mapped (based on the fuzzy search parameters)
#       total number of reads containing a barcode
#       total reads with code errors (currently only containing 2 barcodes)
#
# the other output is a csv of a dataframe
#   the index is the barcode
#   the value is the count of that barcode identified in the fastq

import pandas as pd
import numpy as np
import os
import time
import subprocess
import sys
from Bio import SeqIO
from match_read_and_barcode import match_read_and_barcode, match_read_and_barcode_plasmid, process_barcode_direction
from input_output_helper_functions import find_filepaths, unzip_if_necessary

# how to run the script
# python3 -m yeast_editing_fastq_analysis run_path file_key_path miseq_name
# to access the lab Hive (on a Mac): /Volumes/Shipman-Lab/BaseSpace/...

# get input arguments
args = sys.argv
if len(args) != 4:
    raise ValueError("Not enough arguments! The function expects"
                     "'python yeast_editing_fastq_analysis run_path file_key miseq_name")
run_path = args[1]
file_key_path = args[2]
miseq_name = args[3]

# check for uncommitted changes & get a version number
status = subprocess.check_output(["git", "status"])
if "Changes not staged for commit" in str(status, 'utf-8').strip():
    raise ValueError("Uncommitted changes - please commit before running")
git_short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
git_short_hash = str(git_short_hash, "utf-8").strip()

if not os.path.isdir(run_path):
    raise ValueError("Run path is not a directory - make sure lab Hive is mounted")

# get barcodes & direction for each site
site1_part_df = pd.read_excel("REditor_param_site1_summary.xlsx")
site2_part_df = pd.read_excel("REditor_param_site2_summary.xlsx")
site3_part_df = pd.read_excel("REditor_param_site3_summary.xlsx")
barcodes_site1 = site1_part_df["barcode"]
barcode_dir_site1 = site1_part_df["donor_direction"]
barcodes_site2 = site2_part_df["barcode"]
barcode_dir_site2 = site2_part_df["donor_direction"]
barcodes_site3 = site3_part_df["barcode"]
barcode_dir_site3 = site3_part_df["donor_direction"]
all_barcodes = np.unique(np.concatenate([barcodes_site1, barcodes_site2, barcodes_site3]))

# make final output dataframes for each site & load in file key
file_key = pd.read_excel(file_key_path)
multiindex = pd.MultiIndex.from_tuples(tuple(zip(file_key["site"],
                                                 file_key["rep"],
                                                 file_key["type"],
                                                 file_key["timepoint"])),
                                       names=("site", "rep", "type", "timepoint"))
full_run_qual_df = pd.DataFrame(index=multiindex, columns=["total_num_reads",
                                                           "total_mapped_reads",
                                                           "total_barcode_reads",
                                                           "total_reads_w_code_errors"]).sort_index()
full_barcode_df = pd.DataFrame(index=multiindex, columns=all_barcodes).sort_index()
filepath_df = pd.DataFrame(index=multiindex, columns=["full_run_path"]).sort_index()

# put filepath into the qual df:
for run_id in file_key["file_id"]:
    full_file_paths = find_filepaths(run_path, run_id)
    unzipped_file_path = unzip_if_necessary(full_file_paths, run_id)
    file_key_row = file_key.loc[file_key["file_id"] == run_id, :]
    filepath_df.loc[(file_key_row["site"],
                          file_key_row["rep"],
                          file_key_row["type"],
                          file_key_row["timepoint"]),
                         "full_run_path"] = unzipped_file_path

for count, index in enumerate(full_run_qual_df.index):
    print("-" * 60)
    print("working on %s (%i/%i)" % (index, count, len(full_run_qual_df.index)))
    start_time = time.time()
    site = index[0]
    sample_type = index[2]
    if site == 1:
        barcodes = barcodes_site1
        directions = barcode_dir_site1
        L_fuzzy_string = "CGAATCGCGGCGTCTTGGTA"
        R_fuzzy_string = "GCACGGCAGACTGGCTCACT"
    elif site == 2:
        barcodes = barcodes_site2
        directions = barcode_dir_site2
        L_fuzzy_string = "AAACTAACGGTCAGCAGGAC"
        R_fuzzy_string = "TTGCGGTCGACAGGATGAGC"
    elif site == 3:
        barcodes = barcodes_site3
        directions = barcode_dir_site3
        L_fuzzy_string = "GCCACGTAGGACACATGGTT"
        R_fuzzy_string = "ATCTGGGGTAAAGGCTGTAC"
    else:
        raise ValueError("Site is not 1, 2, or 3")
    # the plasmids have forward & reverse donors, so I need a slightly different
    # pipeline for them
    barcodes_with_dir = process_barcode_direction(barcodes, directions)
    if sample_type == "plasmid":
        with open(filepath_df.loc[index, "full_run_path"].values[0]) as handle:
            dir_barcode_donor_df = pd.DataFrame(index=barcodes, data={"dir_barcode": barcodes_with_dir,
                                                                      "dir": np.array(directions)})
            output_dict = {"total_num_reads": 0,
                           "total_mapped_reads": 0,
                           "total_barcode_reads": 0,
                           "total_reads_w_code_errors": 0}
            full_barcode_df.loc[index, barcodes] = 0
            for record in SeqIO.parse(handle, "fastq"):
                output_dict["total_num_reads"] += 1
                seq = record.seq
                barcode_dict = match_read_and_barcode_plasmid(seq, dir_barcode_donor_df, L_fuzzy_string,
                                                              R_fuzzy_string)
                barcode_match = barcode_dict["barcode_match"]
                if barcode_dict["does_it_map"] & (barcode_match == ""):
                    output_dict["total_mapped_reads"] += 0
                if barcode_dict["does_it_map"]:
                    output_dict["total_mapped_reads"] += 1
                if barcode_match != "":
                    full_barcode_df.loc[index, barcode_match] += 1
                    output_dict["total_barcode_reads"] += 1
                if barcode_dict["code_error"] != "None":
                    output_dict["total_reads_w_code_errors"] += 1
                for column in output_dict.keys():
                    full_run_qual_df.loc[index, column] = output_dict[column]
    else:
        # import pdb 
        # pdb.set_trace()
        with open(filepath_df.loc[index, "full_run_path"].values[0]) as handle:
            output_dict = {"total_num_reads": 0,
                           "total_mapped_reads": 0,
                           "total_barcode_reads": 0,
                           "total_reads_w_code_errors": 0}
            full_barcode_df.loc[index, barcodes] = 0
            for record in SeqIO.parse(handle, "fastq"):
                output_dict["total_num_reads"] += 1
                seq = record.seq
                barcode_dict = match_read_and_barcode(seq, barcodes, L_fuzzy_string, R_fuzzy_string)
                barcode_match = barcode_dict["barcode_match"]
                if barcode_dict["does_it_map"]:
                    output_dict["total_mapped_reads"] += 1
                if barcode_match != "":
                    full_barcode_df.loc[index, barcode_match] += 1
                    output_dict["total_barcode_reads"] += 1
                if barcode_dict["code_error"] != "None":
                    output_dict["total_reads_w_code_errors"] += 1
                for column in output_dict.keys():
                    full_run_qual_df.loc[index, column] = output_dict[column]
    full_run_qual_df.to_excel("%s_run_qual_df_v%s.xlsx" % (miseq_name, git_short_hash))
    full_barcode_df.to_excel("%s_barcode_df_v%s.xlsx" % (miseq_name, git_short_hash))
    print(
        "processing took %s seconds" % (time.time() - start_time))
