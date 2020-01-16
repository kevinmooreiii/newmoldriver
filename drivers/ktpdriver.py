""" driver for rate constant evaluations
"""

import autofile.fs
from drivers import mech as lmech
from routines.pf.calc import zpe as calczpe
from routines.pf.fit import fit_rates
from routines.pf import rates as messrates
from lib.filesystem import build as fbuild
from lib.runner import rates as raterunner


def run(pes_formula,
        spc_dct,
        thy_dct,
        tsk_info_lst,
        rct_names_lst,
        prd_names_lst,
        model_dct,
        run_inp_dct,
        run_rates=True,
        run_fits=True):
    """ main driver for generation of full set of rate constants on a single PES
    """
    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']
    etrans = model_dct['etransfer']
    ene_coeff = model_dct['options']['ene_coeff']
    temps = model_dct['options']['temps']
    pressures = model_dct['options']['pressures']
    pst_params = model_dct['options']['pst_params']
    multi_info = model_dct['options']['multi_info']
    assess_pdep = model_dct['options']['assess_pdep']

    # Prepare filesystem
    fbuild.prefix_filesystem(run_prefix, save_prefix)
    spc_save_fs = autofile.fs.species(save_prefix)

    # Get the levels in lists from the user
    pf_levels = lmech.set_es_model_info(model_dct['es'], thy_dct)
    pf_model = lmech.set_pf_model_info(model_dct['pf'])

    # Run the rates
    if run_rates:

        # Write the strings for the MESS input file
        rct_ichs = spc_dct['ts_0']['rxn_ichs'][0]
        header_str, energy_trans_str = messrates.rate_headers(
            rct_ichs, temps, pressures, **etrans)

        # Write the MESS strings for all the PES channels
        print('Starting mess file preparation for {}:'.format(pes_formula)
        idx_dct = {}
        well_str, bim_str, ts_str = lmech.write_channel_mess_strs(
            spc_dct, rxn_lst, pes_formula,
            pf_model, pf_levels, multi_info, pst_params,
            spc_save_fs, save_prefix, idx_dct)

        # Run mess to produce rate output
        mess_path = raterunner.run_rates(
            header_str, energy_trans_str, well_str, bim_str, ts_str,
            spc_dct['ts_0'], pf_levels[0],
            spc_dct['ts_0']['rxn_fs'][3])

    # Fit rate output to modified Arrhenius forms, print in ChemKin format
    if run_fits:
        fit_rates(spc_dct, pes_formula, idx_dct,
                  pf_levels, pf_model,
                  ene_str, mess_path, assess_pdep)
