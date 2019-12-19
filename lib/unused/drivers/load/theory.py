""" library of reader functions for the theory file
"""

import autoparse.find as apf
from ptt import read_inp_str
from ptt import end_section_wname2
from ptt import build_keyword_dct
from keywords import THY_REQUIRED_KEYWORDS, THY_SUPPORTED_KEYWORDS


THEORY_INP = 'inp/theory.dat'


def read_theory_file():
    """ Read the theory input file into a string
    """
    return read_inp_str(THEORY_INP)


# FUNCTION TO READ IN A STRING FOR A SPECIFIC LEVEL #

def read_theory_sections(thy_inp_str):
    """ species input
    """
    # Obtain the species string
    thy_sections = apf.all_captures(
        end_section_wname2('level'), thy_inp_str)
    # Make sure some section has been defined
    assert thy_sections is not None

    # Build dictionary of theory methods
    thy_methods = {}
    for section in thy_sections:
        name = section[0]
        keyword_dct = build_thy_keyword_dct(section[1])
        thy_methods[name] = keyword_dct

    return thy_methods


def build_thy_keyword_dct(thy_str):
    """ Build a dictionary for all the theory keywords
    """
    keyword_dct = build_keyword_dct(thy_str)
    assert keyword_dct
    assert check_thy_dct(keyword_dct)

    return keyword_dct


def check_thy_dct(dct):
    """ Make sure the theory dictionary keywords are all correct
    """
    keys = dct.keys()
    req_keys_def = all(key in keys for key in THY_REQUIRED_KEYWORDS)
    sup_keys_def = all(key in THY_SUPPORTED_KEYWORDS for key in keys)
    return bool(req_keys_def and sup_keys_def)
