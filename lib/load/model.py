""" Library of reader functions for the model file
"""

import autoparse.find as apf
from lib.load import ptt
from lib.load.keywords import MODEL_SUPPORTED_DCT

MODEL_INP = 'inp/models.dat'


# FUNCTION TO READ IN A STRING FOR A SPECIFIC LEVEL #

def read_models_sections(job_path):
    """ species input
    """
    mod_str = ptt.read_inp_str(job_path, MODEL_INP)
    # Obtain the species string
    model_sections = apf.all_captures(
        ptt.end_section_wname2('model'), mod_str)
    # Make sure some section has been defined
    assert model_sections is not None

    # Build dictionary of models methods
    model_methods = {}
    for section in model_sections:
        name = section[0]
        keyword_dct = build_model_keyword_dct(section[1])
        model_methods[name] = keyword_dct

    return model_methods


def build_model_keyword_dct(model_str):
    """ Build a dictionary for all the models keywords
    """
    # Grab the various sections required for each model
    pf_str = apf.first_capture(ptt.paren_section('pf'), model_str)
    es_str = apf.first_capture(ptt.paren_section('es'), model_str)
    etrans_str = apf.first_capture(ptt.paren_section('etransfer'), model_str)
    options_str = apf.first_capture(ptt.paren_section('options'), model_str)
    assert pf_str is not None
    assert es_str is not None
    assert options_str is not None

    # Get the dictionary for each section and check them
    pf_dct = ptt.build_keyword_dct(pf_str)
    es_dct = ptt.build_keyword_dct(es_str)
    etransfer_dct = ptt.build_keyword_dct(etrans_str)
    options_dct = ptt.build_keyword_dct(options_str)
    # assert check_model_dct(keyword_dct)

    # Combine dcts into single model dct
    model_dct = {}
    model_dct['pf'] = pf_dct
    model_dct['es'] = es_dct
    model_dct['etransfer'] = etransfer_dct
    model_dct['options'] = options_dct

    return model_dct


def check_model_dct(dct):
    """ Make sure the models dictionary keywords are all correct
    """
    proper_def = True
    for key, val in dct.items():
        if key in MODEL_SUPPORTED_DCT:
            if val not in MODEL_SUPPORTED_DCT[key]:
                proper_def = False
        else:
            proper_def = False

    return proper_def
