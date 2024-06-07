"""
BEFORE BEGINNING: mount the Shipman-Lab hive drive on your Mac (what is currently supported by this script)

TO RUN:
    change everything between the two hashes to match your current run
        (note: you must include your run information in the file_key.xlsx file tracked by GitHub)
    commit this file after saving & changes - the code will not run with uncommitted changes
    navigate to this folder in terminal
    python3 -m edit_site_analysis_script
"""

##Import
from Bio import SeqIO
import pandas as pd
from collections import Counter
import time
import subprocess
import gzip
import shutil
import time
import os as os
import numpy as np
from edit_site_analysis_functions import extract_and_match, find_file

## CHANGE EVERYTHING (IF NEEDED) BETWEEN THE TWO HASHES
# run path must currently point to the fastq generation folder: known bug
# if the run path is on your Desktop
# run_path = os.path.expanduser("~/Volumes/Shipman-Lab/BaseSpace/msDMP_01-705234547")
# if the run path is on the hive
run_path = r"C:\BaseSpace"
run_name = "msAGK_25"
# this check is to get rid of the silent error where the script can't
# find your files because the Hive isn't mounted
# needs to be changed to the correct path if you're not on a Mac
if not os.path.isdir(run_path):
    raise ValueError("Hive not mounted")
## END CHANGE REGION

# check git hash & that there are no uncommitted changes
status = subprocess.check_output(["git", "status"])
if "Changes not staged for commit" in str(status, 'utf-8').strip():
    raise ValueError("Uncommitted changes - please commit before running")
git_short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
git_short_hash = str(git_short_hash, "utf-8").strip()

# load in file key
# if you don't care about some of these, just leave them blank in the file key
file_key = pd.read_excel("phage_AGK_25_filekey.xlsx")\
           [["run_id", "info", "plasmid", "wt_nt", "edit_nt",
             "L_inside", "R_inside", "L_outside", "R_outside"]]
outcome_df = file_key.copy()
# assumes you have 9 or fewer replicates per type!!!
outcome_df[["wt", "edited", "unmatched_region", "unmatched_edit_nt"]] = np.NaN
run_ids = outcome_df["run_id"].dropna()

for fastq_name in run_ids:
    filepath = find_file(fastq_name, run_name, run_path)
    print(filepath)
    print("working on %s" %fastq_name)
    start_time = time.time()
    outcomes_dict = {'wt':0, 'edited':0, 'unmatched_region':0, 'unmatched_edit_nt':0}
    all_reads_str = []
    read_counter = []
    L_outside = outcome_df.loc[outcome_df["run_id"] == fastq_name, "L_outside"].values[0]
    R_outside = outcome_df.loc[outcome_df["run_id"] == fastq_name, "R_outside"].values[0]
    L_inside = outcome_df.loc[outcome_df["run_id"] == fastq_name, "L_inside"].values[0]
    R_inside = outcome_df.loc[outcome_df["run_id"] == fastq_name, "R_inside"].values[0]
    wt_nt = outcome_df.loc[outcome_df["run_id"] == fastq_name, "wt_nt"].values[0]
    edited_nt = outcome_df.loc[outcome_df["run_id"] == fastq_name, "edit_nt"].values[0]
    # need to go into the folder & unzip the file
    with gzip.open(filepath, "rt") as handle:
        for seq_record in SeqIO.parse(handle, "fastq"): 
            all_reads_str.append(str(seq_record.seq))
        read_counter = Counter(all_reads_str)
        for read in read_counter:
            read_result = extract_and_match(read, L_outside, R_outside, L_inside,
                                            R_inside, wt_nt, edited_nt)
            if read_result is None:
                import pdb
                pdb.set_trace()
                read_result = extract_and_match(read, L_outside, R_outside, L_inside,
                                                R_inside, wt_nt, edited_nt)
            outcomes_dict[read_result] += read_counter[read]
    # put into output df
    index = outcome_df.index[outcome_df["run_id"] == fastq_name]
    if len(index) != 1:
        raise ValueError("There is more than one row in the outcome df with the same MiSeq run id")
    index = index[0]
    outcome_df.loc[index, ["wt", "edited", "unmatched_region", "unmatched_edit_nt"]] = outcomes_dict
    print("---  processing took %s seconds ---" % (time.time() - start_time))
    outcome_df.to_excel("%s_summary_df_vers"%(run_name) + str(git_short_hash) + ".xlsx")