import fuzzysearch
import numpy as np
from Bio.Seq import Seq

def match_read_and_barcode(read, barcode_list, L_fuzzy_string, R_fuzzy_string):
    """
    This is a wrapper that should process a single read for GENOME SITES ONLY
    This function takes in a single read, list of barcodes, and left/right search
    strings to return which (if any) barcode matches, as well as a general map
    and any errors encountered in the overall processing. 

    Inputs:
    read:                 sequence as a string
    L_fuzzy_string:       left search string
    R_fuzzy_string:       right search string
    barcode_list:

    barcode_match:        string for matching barcode or empty string for no match
    does_it_map:          bool for if the read maps to the right region
    code_error:           string output for if there were any errors
                          processing the read (like double barcodes, etc.)
    """
    [mapped, L_fuzzy_end_index, R_fuzzy_start_index] = does_it_map(read, L_fuzzy_string, R_fuzzy_string)
    if mapped == False:
        return {"barcode_match": "",
        "does_it_map": mapped,
        "code_error": "None"}

    chopped_read = chop_read(read, L_fuzzy_end_index, R_fuzzy_start_index)
    [matched_barcode, code_error] = find_barcode(chopped_read, barcode_list)

    return {"barcode_match": matched_barcode,
            "does_it_map": mapped,
            "code_error": code_error}

def match_read_and_barcode_plasmid(read, dir_barcode_donor_df, L_fuzzy_string_F, R_fuzzy_string_F):
    """
    This is a wrapper that should process a single read for PLASMID READS
    This function takes in a single read and a DataFrame with barcodes as the index
    and full donor sequences as the column under "donor"

    Inputs:
    read:                 sequence as a string
    barcode_donor_df:     DataFrame with all barcodes (original), barcodes in donor direction ("dir_barcode"),
                          and full donor sequence ("donor")

    barcode_match:        string for matching barcode or empty string for no match (original)
    does_it_map:          bool for if the read maps to the right region
    code_error:           string output for if there were any errors
                          processing the read (like double barcodes, etc.)
    """
    # import pdb 
    # pdb.set_trace()
    
    map_dict = does_it_map_fow_rev(read, L_fuzzy_string_F, R_fuzzy_string_F)

    if map_dict["mapped"] == False:
        return {"barcode_match": "",
                "does_it_map": False,
                "code_error": "None"}

    dir_barcode_only = dir_barcode_donor_df.loc[dir_barcode_donor_df["dir"] == map_dict["dir"], :]
    found_barcodes = [x for x in dir_barcode_only["dir_barcode"] if x in map_dict["chopped_read"]]

    if len(found_barcodes) == 0:
        return {"barcode_match": "",
                "does_it_map": True,
                "code_error": "None"}

    # there's an issue with found_barcode_indices, where it always takes the first barcode (not the two correct ones)
    if len(found_barcodes) == 2:
        return {"barcode_match": "",
                "does_it_map": True,
                "code_error": "more_than_one_barcode"}
    else:
        matched_barcode = dir_barcode_donor_df.loc[dir_barcode_donor_df["dir_barcode"] == found_barcodes[0], :].index[0]

    return {"barcode_match": matched_barcode,
            "does_it_map": True,
            "code_error": "None"}

def does_it_map_fow_rev(read, L_fuzzy_string_F, R_fuzzy_string_F):
    """
    This function takes in a read and figures out if it maps in the forward
    or reverse direction
    """
    L_fuzzy_string_R = str(Seq(L_fuzzy_string_F).reverse_complement())
    R_fuzzy_string_R = str(Seq(R_fuzzy_string_F).reverse_complement())

    # first check which fuzzy strings map --> forward or reverse
    [mapped_F, L_fuzzy_end_index_F, R_fuzzy_start_index_F] = does_it_map(read, L_fuzzy_string_F, R_fuzzy_string_F)
    [mapped_R, L_fuzzy_end_index_R, R_fuzzy_start_index_R] = does_it_map(read, R_fuzzy_string_R, L_fuzzy_string_R)

    if mapped_F == True:
        mapped = True
        direction = "F"
        chopped_read = chop_read(read, L_fuzzy_end_index_F, R_fuzzy_start_index_F)
    elif mapped_R == True:
        mapped = True
        direction = "R"
        chopped_read = chop_read(read, L_fuzzy_end_index_R, R_fuzzy_start_index_R)
    elif (mapped_R == True) & (mapped_F == True):
        print("both forward & reverse fuzzy searches mapped to read")
        mapped = False
        direction = ""
        chopped_read = ""
    else:
        mapped = False
        direction = ""
        chopped_read = ""

    return {"mapped": mapped, "dir": direction, "chopped_read": chopped_read}


def does_it_map(read, L_fuzzy_string, R_fuzzy_string):
    """
    This function takes a single read & searches for sequences
    near, but not on the barcoded region using a quite fuzzy search

    Inputs:
    read:                 sequence as a string
    L_fuzzy_string:       left search string
    R_fuzzy_string:       right search string

    Outputs:
    does_it_map:          bool for if the read maps to the right region
    L_fuzzy_end_index:    index at which the left fuzzy string match ends
                          (for cutting up the read later) 
    R_fuzzy_start_index:  index at which the right fuzzy string match starts
                          (for cutting up the read later)  
    """
    L_match = fuzzysearch.find_near_matches(L_fuzzy_string, read, max_l_dist=2)
    R_match = fuzzysearch.find_near_matches(R_fuzzy_string, read, max_l_dist=2)

    if (len(L_match) != 1) | (len(R_match) != 1):
        does_it_map = False
        L_fuzzy_end_index = np.NaN
        R_fuzzy_start_index = np.NaN
    else:
        does_it_map = True
        L_fuzzy_end_index = L_match[0].end - 1 # weird artifact of fuzzysearch puts the match 1 up
        R_fuzzy_start_index = R_match[0].start

    return [does_it_map, L_fuzzy_end_index, R_fuzzy_start_index]



def chop_read(read, L_fuzzy_end_index, R_fuzzy_start_index):
    """
    Takes a read and chops it based on the L & R fuzzy matches
    """
    chopped_read = read[L_fuzzy_end_index+1:R_fuzzy_start_index]
    return chopped_read



def find_barcode(chopped_read, barcode_list):
    """
    Finds matching barcode of barcode list with no fuzziness

    code_error returns an error if there are two or more barcodes
    in a read
    """
    found_barcodes = [x for x in barcode_list if x in chopped_read]

    if len(found_barcodes) == 1:
        matched_barcode = found_barcodes[0]
        code_error = "None"
    elif len(found_barcodes) == 0:
        matched_barcode = ""
        code_error = "None"
    else:
        matched_barcode = ""
        code_error = "more_than_one_barcode"

    return [matched_barcode, code_error]

def process_barcode_direction(barcode_list, barcode_dir):
    """
    Reverse complements barcodes that are in the reverse direction
    (for plasmid coverage analysis)
    """
    barcode_correct_dir_list = np.array(barcode_list)
    barcode_list = np.array(barcode_list)
    barcode_dir = np.array(barcode_dir)
    for barcode_index, barcode in enumerate(barcode_list):
        if barcode_dir[barcode_index] == "R":
            barcode_correct_dir_list[barcode_index] = str(Seq(barcode_list[barcode_index]).reverse_complement())

    return barcode_correct_dir_list