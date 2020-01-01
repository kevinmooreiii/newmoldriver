""" driver for rate constant evaluations
"""

import automol.inchi
import automol.geom
import esdriver
import autofile.fs

# Calling the new libs
from lib.submission import substr
from lib.filesystem import build as fbuild
from lib.filesystem import inf as finf
from lib.filesystem import path as fpath
from lib.filesystem import read as fread
from lib.calc import zpe as calczpe
from lib.runner import rates as raterunner
from lib.reaction.ts import ts_class
from routines.pf.fit import fit_rates
from routines.pf import rates as messrates
from routines.pf import get_high_level_energy
import mech as lmech


TEMPS = [500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.]
PRESS = [0.03, 0.1, 0.3, 1., 3., 10., 30., 100.]

KICKOFF_SIZE = 0.1
KICKOFF_BACKWARD = False


def run(spc_dct, tsk_info_lst, rct_names_lst, prd_names_lst,
        run_prefix, save_prefix, ene_coeff=[1.],
        vdw_params=[False, False, True],
        options=[True, True, True, False],
        etrans=[200.0, 0.85, 15.0, 57.0, 200.0, 3.74, 5.5, 28.0],
        pst_params=[1.0, 6],
        rad_rad_ts='vtst',
        mc_nsamp=[True, 10, 1, 3, 100]):
    """ main driver for generation of full set of rate constants on a single PES
    """

    # Prepare prefix filesystem
    fbuild.prefix_filesystem(run_prefix, save_prefix)

    # Determine options
#    runes = options[0]  # run electronic structure theory (True/False)
    runmess = options[2]  # run mess (True) / only make mess input file (False)
    runrates = options[3]
    if not runmess:
        runrates = False

    # First run ESDriver for the species on the PES so that
    # exothermicity info is available to sort reactions by exothermicity
    spc_queue = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[rxn])
        rxn_spc.extend(list(prd_names_lst[rxn]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)

    spc_tsk_lst = []
    ts_tsk_lst = []
    ts_tsk = False
    for tsk in tsk_info_lst:
        if 'find_ts' in tsk[0]:
            ts_tsk = True
        if ts_tsk:
            ts_tsk_lst.append(tsk)
        else:
            spc_tsk_lst.append(tsk)

    # Form the reaction list
    rxn_lst = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_lst.append(
            {'species': [], 'reacs': list(rct_names_lst[rxn]), 'prods':
             list(prd_names_lst[rxn])})

    # Run the rates
    print('RATES TEST')
    print(runrates)
    print(spc_dct)
    if runrates:
        # for spc in spc_dct:
        #     print('PES FORMULA ASSESS TEST')
        #     print(spc_dct[spc])
        #     print(spc)
        #     if 'original_zma' in spc_dct[spc]:
        #         pes_formula = automol.geom.formula(
        #             automol.zmatrix.geometry(spc_dct[spc]['original_zma']))
        #         print('Starting mess file preparation for {}:'.format(
        #             pes_formula))
        #         break

        # Figure out the model and theory levels for the MESS files
        geo_lvl = ''
        harm_lvl = ''
        anharm_lvl = ''
        tors_lvl = ''
        sym_lvl = ''
        geo_lvl_ref = ''
        harm_lvl_ref = ''
        anharm_lvl_ref = ''
        tors_lvl_ref = ''
        sym_lvl_ref = ''

        ts_model = ['RIGID', 'HARM', '']
        for tsk in ts_tsk_lst:
            if 'samp' in tsk[0] or 'find' in tsk[0]:
                geo_lvl = tsk[1]
                geom = True
                if 'find' in tsk[0]:
                    geo_lvl_ref = geo_lvl
            if 'grad' in tsk[0] or 'hess' in tsk[0]:
                harm_lvl = tsk[1]
                harm_lvl_ref = tsk[2]
                if 'grad' in tsk[0]:
                    grad = True
                if 'hess' in tsk[0]:
                    hess = True
                if not geom:
                    ene_lvl = tsk[1]
                    geo_lvl = tsk[1]
            if 'hr' in tsk[0] or 'tau' in tsk[0]:
                # print('found')
                tors_lvl = tsk[1]
                tors_lvl_ref = tsk[2]
                if 'md' in tsk[0]:
                    ts_model[0] = 'MDHR'
                if 'tau' in tsk[0]:
                    ts_model[0] = 'TAU'
                else:
                    ts_model[0] = '1DHR'
            if 'anharm' in tsk[0] or 'vpt2' in tsk[0]:
                anharm_lvl = tsk[1]
                anharm_lvl_ref = tsk[2]
                ts_model[1] = 'ANHARM'
                if not hess:
                    geo_lvl = tsk[1]
            if 'sym' in tsk[0]:
                sym_lvl = tsk[1]
                sym_lvl_ref = tsk[2]
                if 'samp' in tsk[0]:
                    ts_model[2] = 'SAMPLING'
                if '1DHR' in tsk[0]:
                    ts_model[2] = '1DHR'

        geo_thy_info_ref = finf.get_thy_info(geo_lvl_ref)
        harm_thy_info = finf.get_thy_info(harm_lvl)
        tors_thy_info = None
        anharm_thy_info = None
        sym_thy_info = None
        harm_ref_thy_info = None
        tors_ref_thy_info = None
        anharm_ref_thy_info = None
        sym_ref_thy_info = None
        if tors_lvl:
            tors_thy_info = finf.get_thy_info(tors_lvl)
        if anharm_lvl:
            anharm_thy_info = finf.get_thy_info(anharm_lvl)
        if sym_lvl:
            sym_thy_info = finf.get_thy_info(sym_lvl)
        if harm_lvl_ref:
            harm_ref_thy_info = finf.get_thy_info(harm_lvl_ref)
        if tors_lvl_ref:
            tors_ref_thy_info = finf.get_thy_info(tors_lvl_ref)
        if anharm_lvl_ref:
            anharm_ref_thy_info = finf.get_thy_info(anharm_lvl_ref)
        if sym_lvl_ref:
            sym_ref_thy_info = finf.get_thy_info(sym_lvl_ref)
        pf_levels = [
            harm_thy_info, tors_thy_info,
            anharm_thy_info, sym_thy_info]
        ref_levels = [
            harm_ref_thy_info, tors_ref_thy_info,
            anharm_ref_thy_info, sym_ref_thy_info]

        # Collect ground energies and zero-point energies
        spc_save_fs = autofile.fs.species(save_prefix)
        ts_queue = []
        for spc in spc_dct:   # Have toget them for the TS too
            if 'ts_' in spc:
                ts_queue.append(spc)
            # if spc in ts_found:
            #     ts_queue.append(spc)
        print('getting ready for zpe:')
        for spc in spc_queue + ts_queue:
            print('SPC TEST')
            print(spc_dct[spc])
            spc_info = (
                spc_dct[spc]['ich'], spc_dct[spc]['chg'], spc_dct[spc]['mul'])
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
            print('SPC_CHK')
            print(spc_dct[spc])
            print(spc_save_path)
            # Set the pes_formula using the original zma
            for spc_2 in spc_dct:
                # print('PES FORMULA ASSESS TEST')
                # print(spc_dct[spc])
                # print(spc)
                if 'original_zma' in spc_dct[spc_2]:
                    pes_formula = automol.geom.formula(
                        automol.zmatrix.geometry(spc_dct[spc_2]['original_zma']))
                    print('Starting mess file preparation for {}:'.format(
                        pes_formula))
                    break
            print('ZPE TEST')
            print(spc_dct[spc])
            zpe, _ = calczpe.get_zpe(
                spc, spc_dct[spc], spc_save_path, pf_levels, ts_model)
            spc_dct[spc]['zpe'] = zpe
            ene_strl = []
            ene_lvl = ''
            ene_lvl_ref = ''
            ene_idx = 0
            spc_dct[spc]['ene'] = 0.
            ene_str = '! energy level:'
            # print('looking at ts tasks')
            for tsk in ts_tsk_lst:
                if 'ene' in tsk[0]:
                    if ene_idx > len(ene_coeff)-1:
                        print('Warning - an insufficient energy coefficient list was provided')
                        break
                    ene_lvl = tsk[1]
                    ene_lvl_ref = tsk[2]
                    ene_ref_thy_info = finf.get_thy_info(ene_lvl_ref)
                    ene_thy_info = finf.get_thy_info(ene_lvl)
                    ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
                        ene_coeff[ene_idx], ene_thy_info[3], ene_thy_info[1], ene_thy_info[2],
                        ene_ref_thy_info[3], ene_ref_thy_info[1], ene_ref_thy_info[2]))
                    ene = get_high_level_energy(
                        spc_info=spc_info,
                        thy_low_level=ene_ref_thy_info,
                        thy_high_level=ene_thy_info,
                        save_prefix=save_path,
                        saddle=saddle)
                    # print('ene test:', ene_idx, ene_coeff[ene_idx], ene)
                    spc_dct[spc]['ene'] += ene*ene_coeff[ene_idx]
                    ene_idx += 1
        ene_str += '!               '.join(ene_strl)

        # Collect formula and header string for the PES
        tsname_0 = 'ts_0'
        rct_ichs = spc_dct[tsname_0]['rxn_ichs'][0]
        header_str, energy_trans_str = messrates.rate_headers(
            rct_ichs, TEMPS, PRESS, *etrans)
        multi_info = ['molpro2015', 'caspt2', 'cc-pVDZ', 'RR']
        # multi_info = ['molpro2015', 'caspt2', 'cc-pVTZ', 'RR']

        mess_strs = ['', '', '']
        idx_dct = {}
        first_ground_ene = 0.
        species = messrates.make_all_species_data(
            rxn_lst, spc_dct, save_prefix, ts_model, pf_levels,
            substr.PROJROT)
        for idx, rxn in enumerate(rxn_lst):
            tsname = 'ts_{:g}'.format(idx)
            # if spc_dct[ts]['rad_rad']:
            tsform = automol.geom.formula(
                automol.zmatrix.geometry(spc_dct[tsname]['original_zma']))
            if tsform != pes_formula:
                print('Reaction list contains reactions on different potential energy',
                      'surfaces: {} and {}'.format(tsform, pes_formula))
                print('Will proceed to construct only {}'.format(pes_formula))
                continue
            mess_strs, first_ground_ene = messrates.make_channel_pfs(
                tsname, rxn, species, spc_dct, idx_dct, mess_strs,
                first_ground_ene, spc_save_fs, ts_model, pf_levels,
                multi_info, substr.PROJROT,
                pst_params=pst_params)
            # print(idx_dct)
        well_str, bim_str, ts_str = mess_strs
        ts_str += '\nEnd\n'
        print(well_str)
        print(bim_str)
        print(ts_str)

        # run mess to produce rate output
        mess_path = raterunner.run_rates(
            header_str, energy_trans_str, well_str, bim_str, ts_str,
            spc_dct[tsname_0], geo_thy_info_ref,
            spc_dct[tsname_0]['rxn_fs'][3])

        # Fit rate output to modified Arrhenius forms, print in ChemKin format
        fit_rates(spc_dct, pes_formula, idx_dct,
                  pf_levels, ref_levels, ts_model,
                  ene_str, mess_path)
