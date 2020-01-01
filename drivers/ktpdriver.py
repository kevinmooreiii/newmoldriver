""" driver for rate constant evaluations
"""

import automol.inchi
import automol.geom
import autofile.fs

# Calling the new libs
from lib.submission import substr
from lib.filesystem import build as fbuild
from lib.calc import zpe as calczpe
from lib.runner import rates as raterunner
from routines.pf.fit import fit_rates
from routines.pf import rates as messrates
import mech as lmech


TEMPS = [500., 550., 600., 650., 700., 750., 800.,
         850., 900., 950., 1000., 1050., 1100., 1150.,
         1200., 1250., 1300., 1350., 1400., 1450., 1500.,
         1550., 1600., 1650., 1700., 1750., 1800., 1850.,
         1900., 1950., 2000.]
PRESS = [0.03, 0.1, 0.3, 1., 3., 10., 30., 100.]

KICKOFF_SIZE = 0.1
KICKOFF_BACKWARD = False


def run(spc_dct, tsk_info_lst, rct_names_lst, prd_names_lst,
        run_prefix, save_prefix, ene_coeff=[1.],
        options=[True, True, True, False],
        etrans=[200.0, 0.85, 15.0, 57.0, 200.0, 3.74, 5.5, 28.0],
        pst_params=[1.0, 6]):
    """ main driver for generation of full set of rate constants on a single PES
    """

    # Prepare filesystem
    fbuild.prefix_filesystem(run_prefix, save_prefix)
    spc_save_fs = autofile.fs.species(save_prefix)

    # Determine whether to just run MESS or run Rates+Fits
    runmess = options[2]
    runrates = options[3] if runmess else False
    runfits = True

    # Get the reactant and product speceies in a list to be run
    spc_queue = lmech.form_spc_queue(rct_names_lst, prd_names_lst)
    full_queue = spc_queue + [spc for spc in spc_dct if 'ts_' in spc]

    # Form the reaction list
    rxn_lst = lmech.format_run_rxn_lst(rct_names_lst, prd_names_lst)

    # Format the task info list
    _, ts_tsk_lst = lmech.format_tsk_lst(tsk_info_lst)

    # Get the levels in lists from the user
    pf_levels, ref_levels, ts_model = lmech.set_model_info(ts_tsk_lst)

    # Run the rates
    if runrates:

        # Collect ground energies and zero-point energies
        for spc in full_queue:
            spc_info = (spc_dct[spc]['ich'],
                        spc_dct[spc]['chg'],
                        spc_dct[spc]['mul'])
            if 'ts_' in spc:
                spc_dct[spc] = lmech.set_sadpt_info(
                    ts_tsk_lst, spc_dct, spc, run_prefix, save_prefix)
                spc_save_path = spc_dct[spc]['rxn_fs'][3]
                saddle = True
                save_path = spc_save_path
            else:
                spc_save_fs.leaf.create(spc_info)
                spc_save_path = spc_save_fs.leaf.path(spc_info)
                saddle = False
                save_path = save_prefix

            # Set the pes_formula using the original zma
            pes_formula = lmech.set_pes_formula(spc_dct)
            print('Starting mess file prep for {}:'.format(pes_formula))

            # Cacluate the ZPVE and set in the dict
            zpe, _ = calczpe.get_zpe(
                spc, spc_dct[spc], spc_save_path, pf_levels, ts_model)
            spc_dct[spc]['zpe'] = zpe

            # Set ene and string
            spc_dct[spc]['ene'] = lmech.get_high_energy(
                ts_tsk_lst, spc_info, save_path, saddle, ene_coeff)
            ene_str = lmech.get_ckin_ene_lvl_str(
                ts_tsk_lst, ene_coeff)

        # Collect formula and header string for the PES
        tsname_0 = 'ts_0'
        rct_ichs = spc_dct[tsname_0]['rxn_ichs'][0]
        multi_info = ['molpro2015', 'caspt2', 'cc-pVDZ', 'RR']
        # multi_info = ['molpro2015', 'caspt2', 'cc-pVTZ', 'RR']

        # Write the strings for the MESS input file
        header_str, energy_trans_str = messrates.rate_headers(
            rct_ichs, TEMPS, PRESS, *etrans)

        mess_strs = ['', '', '']
        idx_dct = {}
        first_ground_ene = 0.
        species = messrates.make_all_species_data(
            rxn_lst, spc_dct, save_prefix, ts_model, pf_levels,
            substr.PROJROT)
        for idx, rxn in enumerate(rxn_lst):
            tsname = 'ts_{:g}'.format(idx)
            tsform = automol.geom.formula(
                automol.zmatrix.geometry(spc_dct[tsname]['original_zma']))
            if tsform != pes_formula:
                print('Reaction ist contains reactions on different potential',
                      'energy surfaces: {} and {}'.format(tsform, pes_formula))
                print('Will proceed to construct only {}'.format(pes_formula))
                continue
            mess_strs, first_ground_ene = messrates.make_channel_pfs(
                tsname, rxn, species, spc_dct, idx_dct, mess_strs,
                first_ground_ene, spc_save_fs, ts_model, pf_levels,
                multi_info, substr.PROJROT,
                pst_params=pst_params)
        well_str, bim_str, ts_str = mess_strs
        ts_str += '\nEnd\n'
        print(well_str)
        print(bim_str)
        print(ts_str)

        # Run mess to produce rate output
        mess_path = raterunner.run_rates(
            header_str, energy_trans_str, well_str, bim_str, ts_str,
            spc_dct[tsname_0], ref_levels[4],
            spc_dct[tsname_0]['rxn_fs'][3])

    if runfits:
        # Fit rate output to modified Arrhenius forms, print in ChemKin format
        fit_rates(spc_dct, pes_formula, idx_dct,
                  pf_levels, ref_levels, ts_model,
                  ene_str, mess_path)
