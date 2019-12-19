""" drivers for thermochemistry evaluations
"""

import lib.drivers.thermo as lthmdriver
import lib.filesystem.build as lfs


REF_CALLS = {"basic": "get_basic",
             "cbh0": "get_cbhzed",
             "cbh1": "get_cbhone",
             "cbh2": "get_cbhtwo"}


def run(tsk_info_lst, spc_dct, spc_queue, ref, run_prefix, save_prefix,
        ene_coeff=[1.], options=[True, True, True, False]):
    """ main driver for thermo run
    """

    # Set Boolean values for various options
    runes = options[0]     # Run electronic structure theory
    runmess = options[1]   # Run mess or just make the mess input file
    runthermo = options[2] if runmess else False  # Get NASA polynomial
    # everypf = options[3]   # Make PF for all comps, only w/ final info

    # Fix any issues in tsk_list
    tsk_info_lst = lthmdriver.fix(tsk_info_lst)

    # Add reference molecules
    full_queue = lthmdriver.add_spc_refs_to_lst(spc_queue, spc_dct, ref)

    # Prepare filesystem
    lfs.build_filesystem(run_prefix, save_prefix)

    # Run all necessary ES calculations needed for partition functions
    if runes:
        lthmdriver.run_es_for_thermo(
            full_queue, tsk_info_lst, spc_dct, run_prefix, save_prefix)

    # Determine the ES Theory needed for partition functions
    pf_thy = lthmdriver.set_pf_theory_levels()
    spc_model, pf_levels, ref_levels, glob_pf_str = pf_thy

    # Calculate the partition functions using MESS
    if runmess:

        # Collate PF spc input; initializing enes for each species
        spc_dct = lthmdriver.get_pf_inp_and_ene(
            full_queue, spc_dct,
            spc_model, pf_levels,
            save_prefix)

        # Run MESS
        harm_thy_info = pf_levels[1]
        lthmdriver.run_messpf(
            spc_queue, spc_dct, harm_thy_info, glob_pf_str)

    # Calculate heats-of-formation and NASA polynomials
    if runthermo:

        # pf levels
        lthmdriver.set_pf_levels(
            tsk_info_lst, spc_dct, full_queue, save_prefix, ene_coeff)

        # Need to get zpe for reference molecules too
        lthmdriver.get_hf0k(spc_dct, spc_queue, ref)

        # Write the CHEMKIN strings for NASA plynomial
        lthmdriver.get_nasa_polynomial(
            spc_dct, spc_model, spc_queue,
            pf_levels, ref_levels)
