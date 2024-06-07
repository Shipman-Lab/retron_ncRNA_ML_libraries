import fuzzysearch
import os
from pathlib import Path

def extract_and_match(read, L_outside, R_outside, L_inside, R_inside, wt_nt, edited_nt):
    left_outside = fuzzysearch.find_near_matches(L_outside, read, max_l_dist=4)
    right_outside = fuzzysearch.find_near_matches(R_outside, read, max_l_dist=4)
    if len(left_outside) == 1 and len(right_outside) == 1:
        # allow two mismatches for the PAM edit
        left_inside = fuzzysearch.find_near_matches(L_inside, read, max_l_dist=2)
        right_inside = fuzzysearch.find_near_matches(R_inside, read, max_l_dist=2)
        # check all combinations of matches, because sometimes the inside search seq
        # happens in the bad end of the read and then there are more than two matches >:(
        if (len(left_inside) >= 1) & (len(right_inside) >= 1):
            for left_match in left_inside:
                for right_match in right_inside:
                    var_nt = read[left_match.end:right_match.start]
                    if len(var_nt) == 1:
                        if var_nt == wt_nt:
                            return 'wt'
                        elif var_nt == edited_nt:
                            return 'edited'
                        else:
                            return 'unmatched_edit_nt'
            return 'unmatched_edit_nt'
        else:
            return 'unmatched_region'
    else:
        return 'unmatched_region'

def find_file(file_name, run, basespace_folder):
    miseq_folder_names = os.listdir(basespace_folder)
    miseq_folder_dict = {}
    for folder in miseq_folder_names:
        miseq_folder_dict[folder.split('-')[0]] = folder
    file_folder = os.path.join(basespace_folder, miseq_folder_dict[run])
    print(file_folder)
    for file in Path(file_folder).rglob(file_name + "*.fastq.gz"):
        return file