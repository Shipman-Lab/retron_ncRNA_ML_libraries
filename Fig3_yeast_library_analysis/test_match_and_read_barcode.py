from match_read_and_barcode import *
import numpy as np
import pandas as pd
from Bio.Seq import Seq

def test_does_it_map():
    perfect_site_read = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                         "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTGGTAGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                         "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                         "ATTCGAGAC"
    read_w_changed_PAMs = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                           "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTCTTAGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                           "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                           "ATTCGAGAC"
    edited_read = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                   "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTCTTAGAAGTCCAAGGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                   "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                   "ATTCGAGAC"
    plasmid_read = "TGAATCTGAGTTACTGTCTGTTTTCCTCACTCTCCATCCTGCGTACGCCATCCCCTTCGATGCGAGTGCGA"\
                    "ATCGCGGCGTCTTGGTAGAAGTCCAAGGCACGGCAGACTCTCTCACTGGACGAAGAGTGCATCCGAGAACCA"\
                    "GAAACAAGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCAATTAACGAAATTGCCCC"\
                    "AGTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAA"

    L_fuzzy_site = "CGAATCGCGGCGTCTTGGTA"
    R_fuzzy_site = "GCACGGCAGACTGGCTCACT"

    # for the plasmid, we need to choose fuzzy search parameters in the constant part of the ncRNA
    # barcode is 88 nt from the start of the read
    L_fuzzy_plasmid = "GAGTTACTGTCTGTTTTCCT" # this is on the primer, but there's not any room unfort
    R_fuzzy_plasmid = "AGGAAACCCGTTTCTTCTGA"

    barcode_used = "GAAGTCCAAG"

    # perfect read tests
    perfect_read_result = does_it_map(perfect_site_read, L_fuzzy_site, R_fuzzy_site)
    perfect_read_plasmid_fuzzy = does_it_map(perfect_site_read, L_fuzzy_plasmid, R_fuzzy_plasmid)
    assert perfect_read_result == [True, 113, 114]
    assert perfect_read_plasmid_fuzzy == [False, np.NaN, np.NaN]

    # changed PAM
    PAM_read_result = does_it_map(read_w_changed_PAMs, L_fuzzy_site, R_fuzzy_site)
    PAM_read_plasmid_fuzzy = does_it_map(read_w_changed_PAMs, L_fuzzy_plasmid, R_fuzzy_plasmid)
    assert PAM_read_result == [True, 113, 114]
    assert PAM_read_plasmid_fuzzy == [False, np.NaN, np.NaN]

    # edited read
    edited_read_result = does_it_map(edited_read, L_fuzzy_site, R_fuzzy_site)
    edited_read_plasmid_fuzzy = does_it_map(edited_read, L_fuzzy_plasmid, R_fuzzy_plasmid)
    assert edited_read_result == [True, 113, 124]
    assert edited_read_plasmid_fuzzy == [False, np.NaN, np.NaN]

    # plasmid read
    plasmid_read_result = does_it_map(plasmid_read, L_fuzzy_site, R_fuzzy_site)
    plasmid_read_plasmid_fuzzy = does_it_map(plasmid_read, L_fuzzy_plasmid, R_fuzzy_plasmid)
    assert plasmid_read_plasmid_fuzzy == [True, 26, 149]


def test_chop_read():
    perfect_site_read = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                        "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTGGTAGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                        "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                        "ATTCGAGAC"
    edited_read = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                   "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTCTTAGAAGTCCAAGGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                   "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                   "ATTCGAGAC"
    plasmid_read = "TGAATCTGAGTTACTGTCTGTTTTCCTCACTCTCCATCCTGCGTACGCCATCCCCTTCGATGCGAGTGCGA"\
                    "ATCGCGGCGTCTTGGTAGAAGTCCAAGGCACGGCAGACTCTCTCACTGGACGAAGAGTGCATCCGAGAACCA"\
                    "GAAACAAGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCAATTAACGAAATTGCCCC"\
                    "AGTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAA"


    assert chop_read(perfect_site_read, 113, 114) == ""
    assert chop_read(edited_read, 113, 124) == "GAAGTCCAAG"
    assert chop_read(plasmid_read, 26, 149) == "CACTCTCCATCCTGCGTACGCCATCCCCTTCGATGCGAGTGCGA"\
                                               "ATCGCGGCGTCTTGGTAGAAGTCCAAGGCACGGCAGACTCTCTCA"\
                                               "CTGGACGAAGAGTGCATCCGAGAACCAGAAACA"

def test_find_barcode():
    perfect_site_read = ""
    edited_read = "GAAGTCCAAG"
    plasmid_read = "CACTCTCCATCCTGCGTACGCCATCCCCTTCGATGCGAGTGCGA"\
                   "ATCGCGGCGTCTTGGTAGAAGTCCAAGGCACGGCAGACTCTCTCA"\
                   "CTGGACGAAGAGTGCATCCGAGAACCAGAAACA"
    plasmid_double_barcode = "CACTCTCTACCAAAGCCCATCCTGCGTACGCCATCCCCTTCGATGCGAGTGCGA"\
                             "ATCGCGGCGTCTTGGTAGAAGTCCAAGGCACGGCAGACTCTCTCA"\
                             "CTGGACGAAGAGTGCATCCGAGAACCAGAAACA"
    barcode_list = ["GAAGTCCAAG", "CTACCAAAGC"]

    assert find_barcode(perfect_site_read, barcode_list) == ["", "None"]
    assert find_barcode(edited_read, barcode_list) == ["GAAGTCCAAG", "None"]
    assert find_barcode(plasmid_read, barcode_list) == ["GAAGTCCAAG", "None"]
    assert find_barcode(plasmid_double_barcode, barcode_list) == ["", "more_than_one_barcode"]

def test_match_read_and_barcode():
    perfect_site_read = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                         "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTGGTAGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                         "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                         "ATTCGAGAC"
    read_w_changed_PAMs = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                           "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTCTTAGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                           "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                           "ATTCGAGAC"
    edited_read = "AGCTAGACTTGTTTACTTTGTATTTATTAGTATTTGCCTCACATGAGAACACACACTCTCCATCCTGCGTAC"\
                   "GCCATCCCCTTCGATGCGAGTGCGAATCGCGGCGTCTTCTTAGAAGTCCAAGGCACGGCAGACTGGCTCACTGGACGAAGAGT"\
                   "GCATCCGAGAACCAGAAACAAGCAACGTTCCAGTGAGTTGTTCCACACATCCTTTAAAGTTTAGCGTATAGAAA"\
                   "ATTCGAGAC"
    plasmid_read = "TGAATCTGAGTTACTGTCTGTTTTCCTCACTCTCCATCCTGCGTACGCCATCCCCTTCGATGCGAGTGCGA"\
                    "ATCGCGGCGTCTTGGTAGAAGTCCAAGGCACGGCAGACTCTCTCACTGGACGAAGAGTGCATCCGAGAACCA"\
                    "GAAACAAGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCATACGGAATCTTATCAATTAACGAAATTGCCCC"\
                    "AGTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAA"

    L_fuzzy_site = "CGAATCGCGGCGTCTTGGTA"
    R_fuzzy_site = "GCACGGCAGACTGGCTCACT"
    L_fuzzy_plasmid = "GAGTTACTGTCTGTTTTCCT" # this is on the primer, but there's not any room unfort
    R_fuzzy_plasmid = "AGGAAACCCGTTTCTTCTGA"
    barcode_list = ["GAAGTCCAAG", "CTACCAAAGC"]

    assert match_read_and_barcode(perfect_site_read, barcode_list, L_fuzzy_site, R_fuzzy_site) == {"barcode_match": "",
                                                                                                   "does_it_map": True,
                                                                                                   "code_error": "None"}
    assert match_read_and_barcode(read_w_changed_PAMs, barcode_list, L_fuzzy_site, R_fuzzy_site) == {"barcode_match": "",
                                                                                                     "does_it_map": True,
                                                                                                     "code_error": "None"}
    assert match_read_and_barcode(edited_read, barcode_list, L_fuzzy_site, R_fuzzy_site) == {"barcode_match": "GAAGTCCAAG",
                                                                                             "does_it_map": True,
                                                                                             "code_error": "None"}
    assert match_read_and_barcode(plasmid_read, barcode_list, L_fuzzy_plasmid, R_fuzzy_plasmid) == {"barcode_match": "GAAGTCCAAG",
                                                                                                    "does_it_map": True,
                                                                                                    "code_error": "None"}

def test_match_and_read_barcode_plasmid():
    read_fow = "TGAATCTGAGTTACTGTCTGTTTTCCTTGCGAATCGCCTCGTCTTGGTACACGTAACTGGCACGGCAGACTGGCTCACTAGGAACCGTTTCT"\
                "TCTGACGTAAGGGTGCGCATACGGAATCTTATCACGCGGCGTCTTGGTAGCGAAATTGCGCACGGCAGACTCTCTCACTATCATTTACGGTA"\
                "CGATCAGCGTAAGGGTGCGCATACGGAATCTTATCACGATGCGAGTGCGAATCGGTTTCAGAG"
    read_rev =  "CTCTGAAACCGATTCGCACTCGCATCGTGATAAGATTCCGTATGCGCACCCTTACGCTGATCGTACCGTAAATGATAGTGAGAGAGTCTGCC"\
                "GTGCGCAATTTCGCTACCAAGACGCCGCGTGATAAGATTCCGTATGCGCACCCTTACGTCAGAAGAAACGGTTCCTAGTGAGCCAGTCTGCCG"\
                "TGCCACGTAACTGTACCAAGACGAGGCGATTCGCAAGGAAAACAGACAGTAACTCAGATTCA"
    L_fuzzy_site_fow = "TGCGAATCGCCTCGTCTTGGTA"
    R_fuzzy_site_fow = "GCACGGCAGACTGGCTCACTAG"

    barcode_df = pd.DataFrame(index=["CACGTAACTG", "CAGTTACGTG"],
                              data={"dir_barcode": ["CACGTAACTG", "CACGTAACTG"], "dir": ["F", "R"]})

    assert match_read_and_barcode_plasmid(read_fow, barcode_df, L_fuzzy_site_fow, R_fuzzy_site_fow) == {"barcode_match": "CACGTAACTG",
                                                                    "does_it_map": True,
                                                                    "code_error": "None"}
    # import pdb
    # pdb.set_trace()
    assert match_read_and_barcode_plasmid(read_rev, barcode_df, L_fuzzy_site_fow, R_fuzzy_site_fow) == {"barcode_match": "CAGTTACGTG",
                                                                    "does_it_map": True,
                                                                    "code_error": "None"}

def test_does_it_map_fow_rev():
    read_fow = "TGAATCTGAGTTACTGTCTGTTTTCCTTGCGAATCGCCTCGTCTTGGTACACGTAACTGGCACGGCAGACTGGCTCACTAGGAACCGTTTCT"\
                "TCTGACGTAAGGGTGCGCATACGGAATCTTATCACGCGGCGTCTTGGTAGCGAAATTGCGCACGGCAGACTCTCTCACTATCATTTACGGTA"\
                "CGATCAGCGTAAGGGTGCGCATACGGAATCTTATCACGATGCGAGTGCGAATCGGTTTCAGAG"
    read_rev =  str(Seq(read_fow).reverse_complement())

    L_fuzzy_site_fow = "TGCGAATCGCCTCGTCTTGGTA"
    R_fuzzy_site_fow = "GCACGGCAGACTGGCTCACTAG"

    assert does_it_map_fow_rev(read_fow, L_fuzzy_site_fow, R_fuzzy_site_fow) == {"mapped": True,
                            "dir": "F", "chopped_read": "CACGTAACTG"}
    assert does_it_map_fow_rev(read_rev, L_fuzzy_site_fow, R_fuzzy_site_fow) == {"mapped": True,
                            "dir": "R", "chopped_read": "CAGTTACGTG"}
    assert does_it_map_fow_rev("ATGCA", L_fuzzy_site_fow, R_fuzzy_site_fow) == {"mapped": False,
                            "dir": "", "chopped_read": ""}


