""" driver for rate constant evaluations
"""

import autofile.fs
from drivers import mech as lmech
from routines.pf.calc import zpe as calczpe
from routines.pf.fit import fit_rates
from routines.pf import rates as messrates
from lib.filesystem import build as fbuild
from lib.runner import rates as raterunner


def run(spc_dct,
        thy_dct,
        tsk_info_lst,
        rct_names_lst,
        prd_names_lst,
        model_dct,
        run_inp_dct,
        rate_jobs=(True, True)):
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

    # Determine whether to just run MESS or run Rates+Fits
    runrates = rate_jobs[0]
    runfits = rate_jobs[1]

    # Get the reactant and product speceies in a list to be run
    spc_queue = lmech.form_spc_queue(rct_names_lst, prd_names_lst)
    full_queue = spc_queue + [spc for spc in spc_dct if 'ts_' in spc]

    # Form the reaction list
    rxn_lst = lmech.format_run_rxn_lst(rct_names_lst, prd_names_lst)

    # Format the task info list
    _, ts_tsk_lst = lmech.format_tsk_lst(tsk_info_lst)

    # Get the levels in lists from the user
    pf_levels = lmech.set_es_model_info(model_dct['es'], thy_dct)
    pf_model = lmech.set_pf_model_info(model_dct['pf'])

    # Run the rates
    if runrates:

        # Collect ground energies and zero-point energies
        for spc in full_queue:
            spc_info = (spc_dct[spc]['ich'],
                        spc_dct[spc]['chg'],
                        spc_dct[spc]['mul'])
            print('SPC CHECK')
            print(spc)
            print(spc_dct)
            if 'ts_' in spc:
                spc_dct[spc] = lmech.set_sadpt_info(
                    ts_tsk_lst, spc_dct, spc, thy_dct,
                    run_prefix, save_prefix,
                    kickoff=(0.1, False))
                spc_save_path = spc_dct[spc]['rxn_fs'][3]
                saddle = True
                save_path = spc_save_path
                # Set the pes_formula using the original zma
                pes_formula = lmech.set_pes_formula(spc_dct)
                print('Starting mess file prep for {}:'.format(pes_formula))
            else:
                spc_save_fs.leaf.create(spc_info)
                spc_save_path = spc_save_fs.leaf.path(spc_info)
                saddle = False
                save_path = save_prefix

            # Cacluate the ZPVE and set in the dict
            spc_dct[spc]['zpe'], _ = calczpe.get_zpe(
                spc, spc_dct[spc], spc_save_path, pf_levels, ts_model)

            # Set ene and string
            spc_dct[spc]['ene'] = lmech.get_high_energy(
                ts_tsk_lst, thy_dct, spc_info, save_path, saddle, ene_coeff)
            ene_str = lmech.get_ckin_ene_lvl_str(
                ts_tsk_lst, thy_dct, ene_coeff)

        # Write the strings for the MESS input file
        rct_ichs = spc_dct['ts_0']['rxn_ichs'][0]
        header_str, energy_trans_str = messrates.rate_headers(
            rct_ichs, temps, pressures, *etrans)

        # Write the MESS strings for all the PES channels
        mess_strs = ['', '', '']
        idx_dct = {}
        well_str, bim_str, ts_str = lmech.write_channel_mess_strs(
            spc_dct, rxn_lst, pes_formula,
            ts_model, pf_levels, multi_info, pst_params,
            spc_save_fs, save_prefix,
            idx_dct, mess_strs)

        # Run mess to produce rate output
        mess_path = raterunner.run_rates(
            header_str, energy_trans_str, well_str, bim_str, ts_str,
            spc_dct['ts_0'], pf_levels[0],
            spc_dct['ts_0']['rxn_fs'][3])

    if runfits:
        # Fit rate output to modified Arrhenius forms, print in ChemKin format
        fit_rates(spc_dct, pes_formula, idx_dct,
                  pf_levels, ts_model,
                  ene_str, mess_path, assess_pdep)
