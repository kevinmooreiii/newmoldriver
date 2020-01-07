""" library of reader functions for the theory file
"""

import autoparse.find as apf
from lib.load.ptt import read_inp_str
from lib.load.ptt import end_section_wname2
from lib.load.ptt import build_keyword_dct
from lib.load.keywords import THY_REQUIRED_KEYWORDS, THY_SUPPORTED_KEYWORDS


THEORY_INP = 'inp/theory.dat'


def build_thy_dct(job_path):
    """ species input
    """
    # Obtain the species string
    thy_str = read_inp_str(job_path, THY_INP)
    thy_sections = apf.all_captures(
        end_section_wname2('level'), thy_str)
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
