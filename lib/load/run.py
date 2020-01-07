""" library of reader functions for the run.dat file
"""

import autoparse.find as apf
from lib.load.ptt import read_inp_str
from lib.load.ptt import parse_idx_inp
from lib.load.ptt import paren_section
from lib.load.ptt import end_section
from lib.load.ptt import keyword_pattern
from lib.load.ptt import end_section_wname2
from lib.load.ptt import build_keyword_dct
from lib.load.ptt import remove_empty_lines
from lib.load.keywords import RUN_INP_REQUIRED_KEYWORDS

RUN_INP = 'inp/run.dat'



# PARSE THE INPUT SECTION OF THE FILE #
def build_run_inp_dct(job_path):
    """ Build a dictionary for all the theory keywords
    """
    run_str = read_inp_str(job_path, RUN_INP)
    keyword_dct = build_keyword_dct(inp_block(run_str))
    assert keyword_dct
    assert check_run_keyword_dct(keyword_dct)

    return keyword_dct


def inp_block(inp_str):
    """ Read the string that has the global model information
    """
    return remove_empty_lines(apf.first_capture(end_section('input'), inp_str))


def check_run_keyword_dct(dct):
    """ Make sure the theory dictionary keywords are all correct
    """
    keys = dct.keys()
    req_keys_def = all(key in keys for key in RUN_INP_REQUIRED_KEYWORDS)
    return req_keys_def


# PARSE THE OBJ SECTION OF THE FILE #
def objects_lst(job_path):
    """ Get the sections for the run block
    """
    run_str = read_inp_str(job_path, RUN_INP)
    obj_str = object_block(run_str)
    # Read one of a set of objects to run calcs on (only one supported)
    pes_block_str = apf.first_capture(paren_section('pes'), obj_str)
    spc_block_str = apf.first_capture(paren_section('spc'), obj_str)
    # ts_block_str = apf.first_capture(paren_section('ts'), section_str)
    # wells_block_str = apf.first_capture(paren_section('wells'), section_str)
    if pes_block_str is not None:
        run_lst = get_pes_idxs(remove_empty_lines(pes_block_str))
    elif spc_block_str is not None:
        run_lst = get_spc_idxs(remove_empty_lines(spc_block_str))
    # elif ts_block_str is not None:
    #     obj_str = ts_block_str
    # elif ts_block_str is not None:
    #     obj_str = ts_block_str
    # elif wells_block_str is not None:
    #     obj_str = wells_block_str
    else:
        raise ValueError
    return run_lst


def object_block(inp_str):
    """ Read the string that has the global model information
    """
    return remove_empty_lines(apf.first_capture(end_section('obj'), inp_str))


# put a length criterion which can be done
def get_pes_idxs(pes_str):
    """ Determine the indices corresponding to what PES and channels to run
    """
    run_pes = []
    for line in pes_str.splitlines():
        [pes_idxs, chn_idxs, proc] = line.strip().split(';')
        pes_lst = parse_idx_inp(pes_idxs)
        chn_lst = parse_idx_inp(chn_idxs)
        proc = proc.strip()
        for pes in pes_lst:
            for chn in chn_lst:
                run_pes.append([pes, chn, proc])

    return run_pes


def get_spc_idxs(pes_str):
    """ Determine the indices corresponding to what species to run
    """
    run_spc = []
    for line in pes_str.splitlines():
        [spc_idxs, proc] = line.strip().split(';')
        spc_lst = parse_idx_inp(spc_idxs)
        proc = proc.strip()
        for spc in spc_lst:
            run_spc.append([spc, proc])

    return run_spc


# PARSE THE PROC SECTION OF THE FILE #
def read_proc_sections(run_inp_str):
    """ species input
    """
    # Obtain the species string
    proc_sections = apf.all_captures(
        end_section_wname2('proc'), run_inp_str)
    # Make sure some section has been defined
    assert proc_sections is not None

    # Build dictionary of procedures
    procs = {}
    for section in proc_sections:
        name = section[0]
        keyword_dct = build_proc_keyword_dct(section[1])
        procs[name] = keyword_dct

    return procs


def build_proc_keyword_dct(proc_str):
    """ Build a dictionary for all the theory keywords
    """
    # keyword_dct = build_keyword_dct(thy_str)
    # assert keyword_dct
    # assert check_thy_dct(keyword_dct)

    jobs_str = apf.first_capture(paren_section('jobs'), proc_str)
    es_tsks_str = apf.first_capture(paren_section('es_tsks'), proc_str)
    options_str = apf.first_capture(paren_section('options'), proc_str)
    model_str = apf.first_capture(keyword_pattern('model'), proc_str)

    # Maybe do the checking in other functions cuz theres a lot..
    # Check for requirements on jobs; change str as needed
    assert jobs_str is not None
    jobs = [line.strip() for line in jobs_str.splitlines()
            if line.strip() != '']

    # Check for requirements on es_tsks; change str as needed
    if 'es' in jobs:
        assert es_tsks_str is not None
    # Really need a way to get the es_tsk list

    # Check for requirements on models; change str as needed
    if 'rates' in jobs or 'thermo' in jobs or 'fit' in jobs:
        assert model_str is not None

    # Check for requirements on options; change str as needed
    options_dct = build_keyword_dct(options_str)

    # Add the parts to the proc dictionary
    proc_dct = {}
    proc_dct['jobs'] = jobs
    proc_dct['model'] = model_str
    proc_dct['es_tsks'] = es_tsks_str
    proc_dct['options'] = options_dct

    return proc_dct
