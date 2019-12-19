""" Functions for setting the Electronic Structure tasks for
    an outlined procedure
"""

from keywords import ES_TSK_SUPPORTED_LST


def set_es_tsks(es_tsk_str, model_dct, saddle=False):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    # Set the task list using either the given models or supplied list
    if 'models' in es_tsk_str:
        tsk_lst = es_tsks_from_models(model_dct, saddle=saddle)
    else:
        tsk_lst = es_tsks_from_lst(es_tsk_str)

    return tsk_lst


def es_tsks_from_lst(es_tsks_str):
    """ Take the es tsk list string from input and set the tasks
        Right now, we presume the tasks given in the file are correct
    """
    # Split the string into different strings of keywords
    tsk_lst = [line.strip() for line in es_tsks_str.splitlines()
               if line.strip() != '']

    # Ensure that all the tasks are in the supported tasks
    assert check_es_tsks_supported(tsk_lst)

    # Hard to get the full tsk_lst for both function
    # tsk inp_geo out_info overwrite
    # Take the overwrite opt keyword plus el method from model

    return tsk_lst


def es_tsks_from_models(model_dct, saddle=False):  # barrierless=False):
    """ Set a list of tasks using the model dictionary for
        setting up the electronic structure calculations
    """

    # Make the assumption that we will pass saddle to these functions?

    # Initialize task lst based on saddle or not
    if saddle:
        tsk_lst = ['find_ts']
    else:
        tsk_lst = ['find_geom']

    # Append conformer things
    tsk_lst.append('conf_samp')
    tsk_lst.append('conf_hess')

    # Add tasks based on model choice
    if model_dct['sym'] == 'sampling':
        tsk_lst.append('sym_samp')
    if model_dct['tors'] == '1dhr':
        tsk_lst.append('hr_scan')
    if model_dct['vib'] == 'vpt2':
        tsk_lst.append('conf_vpt2')

    # Add an energy calculation task
    tsk_lst.append('conf_energy')

    # Add in the tasks for running ircs
    if model_dct['ts_sadpt'] == 'vtst':
        tsk_lst.append('run_irc')

    # Ensure that all the tasks are in the supported tasks
    assert check_es_tsks_supported(tsk_lst)

    return tsk_lst


def check_es_tsks_supported(es_tsk_lst):
    """ Check to see if the list of es tasks are supported by the code
    """
    return all(tsk in ES_TSK_SUPPORTED_LST for tsk in es_tsk_lst)
