""" driver for rate constant evaluations
"""

import lib.drivers.ktp as lktpdriver
import lib.filesystem.build as lfs


def run(tsk_info_lst, es_dct, spc_dct, rct_names_lst, prd_names_lst,
        run_prefix, save_prefix, ene_coeff=[1.],
        vdw_params=[False, False, True],
        options=[True, True, True, False],
        etrans=[200.0, 0.85, 15.0, 57.0, 200.0, 3.74, 5.5, 28.0],
        pst_params=[1.0, 6],
        rad_rad_ts='vtst'):
    """ main driver for generation of full set of rate constants on a single PES
    """

    # Prepare filesystem
    lfs.build_filesystem(run_prefix, save_prefix)

    # Set Boolean values for various options
    runes = options[0]     # Run electronic structure theory
    runmess = options[1]   # Run mess or just make the mess input file
    runrates = options[2] if runmess else False  # Get Rate Constants

    # Build lists of tasks for electronic structure
    spc_tsk_lst, ts_tsk_lst = lktpdriver.build_spc_list()

    # Run all necessary ES calculations needed for partition functions
    if runes:
        lktpdriver.init_es_run()

    # Form the reaction list (RUNS ES FOR TSs?)
    lktpdriver.ts_tsks()

    # Set the models for MESS
    lktpdriver.determine_mess_models()
    lktpdriver.determine_model_thy_levels()

    # Write MESS Reaction Files
    lktpdriver.write_mess_reaction_files()

    if runmess:
        # Run the Rates with MESS
        lktpdriver.run_mess_rates()

    if runrates:
        # Read the rate constants
        lfit.read_rates()

        # Fit the rate constants
        lfit.fit_rates()
