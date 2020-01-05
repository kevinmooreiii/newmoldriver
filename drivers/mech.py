""" estoktp driver functions with no
better home
"""

import os
import sys
import numpy
import automol.geom
from drivers import esdriver
from drivers import ktpdriver

from routines.pf import get_high_level_energy
from routines.pf import rates as messrates
from lib.phydat import phycon
from lib.submission import substr
from lib.filesystem import inf as finf
from lib.filesystem import path as fpath
from lib.filesystem import read as fread
from lib.load import geo
from lib.reaction import rxnid


def run_driver(pes_dct, pesnums_lst, channels, connchnls_lst,
               spc_dct, cla_dct,
               tsk_info_lst,
               run_prefix, save_prefix,
               ene_coeff=[1.],
               vdw_params=[False, False, True],
               options=[True, True, True, False],
               etrans=[200.0, 0.85, 15.0, 57.0, 200.0, 3.74, 5.5, 28.0],
               pst_params=[1.0, 6],
               rad_rad_ts='vtst',
               hind_inc=30.0,
               mc_nsamp=[True, 10, 1, 3, 100],
               multi_info=['molpro2015', 'caspt2', 'cc-pVDZ', 'RR'],
               temps=[500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.],
               pressures=[0.03, 0.1, 0.3, 1., 3., 10., 30., 100.],
               assess_pdep=[0.3, 3.0, [500., 1000.0]],
               kickoff=(0.1, False),
               driver='es_spc'):
    """ Run the ktp driver for the PESs
    """

    for pes_idx, pes in enumerate(pes_dct, start=1):
        # Only run the PESs that the user has specified
        if pes_idx in pesnums_lst:
            # Build the names list
            pes_rct_names_lst = pes_dct[pes]['rct_names_lst']
            pes_prd_names_lst = pes_dct[pes]['prd_names_lst']
            pes_rxn_name_lst = pes_dct[pes]['rxn_name_lst']
            if isinstance(channels, str):
                if channels == 'all':
                    print(len(pes_rxn_name_lst))
                    pes_chns = numpy.arange(len(pes_rxn_name_lst)+1)
                elif '-' in channels:
                    start, end = channels.split('-')
                    pes_chns = numpy.arange(int(start), int(end)+1)
                elif '[' in channels:
                    nums = channels.replace('[', '').replace(']', '')
                    pes_chns = [int(num) for num in nums.split(',')]

            # Print out the channel that is being run
            for cidx, cvals in enumerate(connchnls_lst[pes_idx-1].values()):
                print('PES{}_{}: {} surface'.format(
                    str(pes_idx), str(cidx+1), pes))
                for cval in cvals:
                    print(pes_rxn_name_lst[cval])

            # Run the KTPDriver over the channels
            for cidx, cvals in enumerate(connchnls_lst[pes_idx-1].values()):
                print('ktp on PES{}_{}: {} for following channels...'.format(
                    str(pes_idx), str(cidx+1), pes))
                run_pes = False
                rct_names_lst = []
                prd_names_lst = []
                rxn_name_lst = []
                for chn_idx, _ in enumerate(pes_rxn_name_lst):
                    if chn_idx+1 in pes_chns and chn_idx in cvals:
                        run_pes = True
                        rct_names_lst.append(pes_rct_names_lst[chn_idx])
                        prd_names_lst.append(pes_prd_names_lst[chn_idx])
                        rxn_name_lst.append(pes_rxn_name_lst[chn_idx])
                        print('running channel {}: {} = {}'.format(
                            str(chn_idx+1),
                            ' + '.join(pes_rct_names_lst[chn_idx]),
                            ' + '.join(pes_prd_names_lst[chn_idx])))
                if run_pes:
                    rxn_lst = []
                    for rxn, _ in enumerate(rct_names_lst):
                        rxn_lst.append(
                            {'species': [],
                             'reacs': list(rct_names_lst[rxn]),
                             'prods': list(prd_names_lst[rxn])})
                    ts_idx = 0
                    for idx, rxn in enumerate(rxn_lst):
                        reacs = rxn['reacs']
                        prods = rxn['prods']
                        tsname = 'ts_{:g}'.format(ts_idx)
                        spc_dct[tsname] = {}
                        rname = rxn_name_lst[idx]
                        rname_eq = '='.join(rname.split('=')[::-1])
                        if rname in cla_dct:
                            spc_dct[tsname]['given_class'] = cla_dct[rname]
                        elif rname_eq in cla_dct:
                            spc_dct[tsname]['given_class'] = cla_dct[rname_eq]
                            reacs = rxn['prods']
                            prods = rxn['reacs']
                        else:
                            spc_dct[tsname]['given_class'] = None
                        if reacs and prods:
                            spc_dct[tsname]['reacs'] = reacs
                            spc_dct[tsname]['prods'] = prods
                        spc_dct[tsname]['ich'] = ''
                        ts_chg = 0
                        for rct in rct_names_lst[idx]:
                            print(spc_dct[rct])
                            ts_chg += spc_dct[rct]['chg']
                        spc_dct[tsname]['chg'] = ts_chg
                        mul_low, _, rad_rad = rxnid.ts_mul_from_reaction_muls(
                            rct_names_lst[idx], prd_names_lst[idx], spc_dct)
                        spc_dct[tsname]['mul'] = mul_low
                        spc_dct[tsname]['rad_rad'] = rad_rad
                        spc_dct[tsname]['hind_inc'] = (
                            hind_inc * phycon.DEG2RAD)
                        ts_idx += 1

                    # Run the appropriate driver
                    if driver == 'es_spc':
                        spc_esdriver(
                            rct_names_lst, prd_names_lst, tsk_info_lst,
                            spc_dct, run_prefix, save_prefix, vdw_params,
                            rad_rad_ts=rad_rad_ts,
                            mc_nsamp=mc_nsamp)
                    if driver == 'es_rxn':
                        rxn_esdriver(
                            rct_names_lst, prd_names_lst, tsk_info_lst,
                            spc_dct, run_prefix, save_prefix, vdw_params,
                            rad_rad_ts=rad_rad_ts,
                            mc_nsamp=mc_nsamp,
                            kickoff=kickoff)
                    elif driver == 'ktp':
                        ktpdriver.run(
                            spc_dct,
                            tsk_info_lst,
                            pes_rct_names_lst,
                            pes_prd_names_lst,
                            run_prefix,
                            save_prefix,
                            ene_coeff=ene_coeff,
                            options=options,
                            etrans=etrans,
                            pst_params=pst_params,
                            multi_info=multi_info,
                            temps=temps,
                            pressures=pressures,
                            assess_pdep=assess_pdep)


def get_user_input():
    """ Read the options from the command line
    """
    # Set mechanism and type to be read based on user input
    data_path = sys.argv[1]
    mechanism_name = sys.argv[2]
    mech_type = sys.argv[3]
    mech_path = os.path.join(data_path, 'data', mechanism_name)
    mech_file = 'mech.json'

    # Set further parameters for what reactions and PESs to be run
    if len(sys.argv) > 4:
        pesnums = sys.argv[4]
    if len(sys.argv) > 5:
        channels = sys.argv[5]
    print('pesnums and params.channels:', pesnums, channels)

    return data_path, mech_path, mech_type, mech_file, pesnums, channels


def build_geom_dct(data_path):
    """ Obtain dct containing geometries to use as input
    """
    geom_path = os.path.join(data_path, 'data', 'geoms')
    geom_dct = geo.geometry_dictionary(geom_path)
    return geom_dct


def etrans_lst(params):
    """ set the etrans list
    """
    return [params.EXP_FACTOR,
            params.EXP_POWER,
            params.EXP_CUTOFF,
            params.EPS1,
            params.EPS2,
            params.SIG1,
            params.SIG2,
            params.MASS1]


def spc_esdriver(rct_names_lst, prd_names_lst, tsk_info_lst,
                 spc_dct, run_prefix, save_prefix, vdw_params,
                 rad_rad_ts='pst',
                 mc_nsamp=[True, 10, 1, 3, 100]):
    """ Call the ESDriver Routines for species in the spc dct
    """

    # Get the reactant and product speceies in a list to be run
    spc_queue = form_spc_queue(rct_names_lst, prd_names_lst)
    spc_run_lst = format_run_spc_lst(spc_queue)

    # Format the task info list
    spc_tsk_lst, _ = format_tsk_lst(tsk_info_lst)

    # Execute ESDriver
    esdriver.run(
        spc_tsk_lst, spc_run_lst, spc_dct,
        run_prefix, save_prefix, vdw_params,
        rad_rad_ts=rad_rad_ts,
        mc_nsamp=mc_nsamp)


def rxn_esdriver(rct_names_lst, prd_names_lst, tsk_info_lst,
                 spc_dct, run_prefix, save_prefix, vdw_params,
                 rad_rad_ts='pst',
                 mc_nsamp=[True, 10, 1, 3, 100],
                 kickoff=(0.1, False)):
    """ Call the ESDriver Routines for reactions in the spc dct
    """
    # Format the task info list
    _, ts_tsk_lst = format_tsk_lst(tsk_info_lst)

    # Form the reaction list
    rxn_lst = format_run_rxn_lst(rct_names_lst, prd_names_lst)

    # Add addtional dictionary items for all the TSs
    # This presumes that es has been run previously for species list
    # to produce energies in save file system
    if ts_tsk_lst:
        print('\nBegin transition state prep')
        for spc in spc_dct:
            if 'ts_' in spc:
                spc_dct[spc] = set_sadpt_info(
                    ts_tsk_lst, spc_dct, spc,
                    run_prefix, save_prefix,
                    kickoff)
        print('End transition state prep\n')

        # Execute ESDriver
        esdriver.run(
            ts_tsk_lst, rxn_lst, spc_dct,
            run_prefix, save_prefix, vdw_params,
            rad_rad_ts=rad_rad_ts,
            mc_nsamp=mc_nsamp,
            kickoff=kickoff)


def form_spc_queue(rct_names_lst=(), prd_names_lst=()):
    """ form the species queue from tht elist
    """
    spc_queue = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[rxn])
        rxn_spc.extend(list(prd_names_lst[rxn]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)

    return spc_queue


def form_spc_queue_2(rxn_lst):
    """ second spc queue; redundant to one above
    """
    spc_queue = []
    for _, rxn in enumerate(rxn_lst):
        reacs = rxn['reacs']
        prods = rxn['prods']
        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)
    spc_queue = list(dict.fromkeys(spc_queue))
    return spc_queue


def format_tsk_lst(tsk_info_lst):
    """ Format the input task list appropriate for spc or rxns
    """
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

    return spc_tsk_lst, ts_tsk_lst


def format_run_spc_lst(spc_queue):
    """ get the list of species to run
    """
    return [{'species': spc_queue, 'reacs': [], 'prods': []}]


def format_run_rxn_lst(rct_names_lst, prd_names_lst):
    """ Get the lst of reactions to be run
    """
    run_lst = []
    for rxn, _ in enumerate(rct_names_lst):
        run_lst.append(
            {'species': [], 'reacs': list(rct_names_lst[rxn]), 'prods':
             list(prd_names_lst[rxn])})

    return run_lst


def set_sadpt_info(ts_tsk_lst, spc_dct, spc, run_prefix, save_prefix, kickoff):
    """ set the saddle point dct with info
    """
    es_ini_key = ts_tsk_lst[0][2]
    es_run_key = ts_tsk_lst[0][1]
    ini_thy_info = finf.get_es_info(es_ini_key)
    thy_info = finf.get_es_info(es_run_key)

    # Generate rxn data, reorder if necessary, and put in spc_dct for given ts
    rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul = finf.rxn_info(
        save_prefix, spc, spc_dct, thy_info, ini_thy_info)
    spc_dct[spc]['rxn_ichs'] = rxn_ichs
    spc_dct[spc]['rxn_chgs'] = rxn_chgs
    spc_dct[spc]['rxn_muls'] = rxn_muls
    spc_dct[spc]['low_mul'] = low_mul
    spc_dct[spc]['high_mul'] = high_mul

    # Generate rxn_fs from rxn_info stored in spc_dct
    [kickoff_size, kickoff_backward] = kickoff
    rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path = fpath.get_rxn_fs(
        run_prefix, save_prefix, spc_dct[spc])
    spc_dct[spc]['rxn_fs'] = [
        rxn_run_fs,
        rxn_save_fs,
        rxn_run_path,
        rxn_save_path]
    rct_zmas, prd_zmas, rct_cnf_save_fs, prd_cnf_save_fs = fread.get_zmas(
        spc_dct[spc]['reacs'], spc_dct[spc]['prods'], spc_dct,
        ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, substr.PROJROT)
    ret = rxnid.ts_class(
        rct_zmas, prd_zmas, spc_dct[spc]['rad_rad'],
        spc_dct[spc]['mul'], low_mul, high_mul,
        rct_cnf_save_fs, prd_cnf_save_fs, spc_dct[spc]['given_class'])
    ret1, ret2 = ret
    if ret1:
        [rxn_class, spc_zma,
         dist_name, brk_name, grid,
         frm_bnd_key, brk_bnd_key,
         tors_names, update_guess] = ret1
        spc_dct[spc]['class'] = rxn_class
        spc_dct[spc]['grid'] = grid
        spc_dct[spc]['tors_names'] = tors_names
        spc_dct[spc]['original_zma'] = spc_zma
        dist_info = [dist_name, 0., update_guess, brk_name]
        spc_dct[spc]['dist_info'] = dist_info
        spc_dct[spc]['frm_bnd_key'] = frm_bnd_key
        spc_dct[spc]['brk_bnd_key'] = brk_bnd_key
        # Adding in the rct and prd zmas for vrctst
        spc_dct[spc]['rct_zmas'] = rct_zmas
        spc_dct[spc]['prd_zmas'] = prd_zmas
        if ret2:
            spc_dct[spc]['bkp_data'] = ret2
        else:
            spc_dct[spc]['bkp_data'] = None
    else:
        spc_dct[spc]['class'] = None
        spc_dct[spc]['bkp_data'] = None

    return spc_dct[spc]


def set_model_info(ts_tsk_lst):
    """ Set the model info
    """
    # Initialize the levels
    geo_lvl = ''
    harm_lvl = ''
    anharm_lvl = ''
    tors_lvl = ''
    sym_lvl = ''
    # Initialize the reference levels
    geo_lvl_ref = ''
    harm_lvl_ref = ''
    anharm_lvl_ref = ''
    tors_lvl_ref = ''
    sym_lvl_ref = ''

    # Set model levels
    ts_model = ['RIGID', 'HARM', '']
    geom = False
    hess = False
    for tsk in ts_tsk_lst:
        if 'samp' in tsk[0] or 'find' in tsk[0]:
            geo_lvl = tsk[1]
            geom = True
            if 'find' in tsk[0]:
                geo_lvl_ref = geo_lvl
        if 'grad' in tsk[0] or 'hess' in tsk[0]:
            harm_lvl = tsk[1]
            harm_lvl_ref = tsk[2]
            if 'hess' in tsk[0]:
                hess = True
            if not geom:
                geo_lvl = tsk[1]
        if 'hr' in tsk[0] or 'tau' in tsk[0]:
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

    # Set the theory info objects
    harm_thy_info = finf.get_thy_info(harm_lvl)
    tors_thy_info = None
    anharm_thy_info = None
    sym_thy_info = None
    geo_ref_thy_info = finf.get_thy_info(geo_lvl_ref)
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

    # Combine levels into a list
    pf_levels = [
        harm_thy_info, tors_thy_info,
        anharm_thy_info, sym_thy_info]
    ref_levels = [
        harm_ref_thy_info, tors_ref_thy_info,
        anharm_ref_thy_info, sym_ref_thy_info,
        geo_ref_thy_info]

    return pf_levels, ref_levels, ts_model


def set_pes_formula(spc_dct):
    """ Set pes formula using zma
    """
    for spc_2 in spc_dct:
        print('TEST PES FORM')
        print(spc_2)
        print(spc_dct[spc_2])
        if 'original_zma' in spc_dct[spc_2]:
            pes_formula = automol.geom.formula(
                automol.zmatrix.geometry(spc_dct[spc_2]['original_zma']))
            print('Starting mess file preparation for {}:'.format(
                pes_formula))
            break
    return pes_formula


def get_high_energy(ts_tsk_lst, spc_info, save_path, saddle, ene_coeff):
    """
    """
    spc_ene = 0.0
    ene_idx = 0
    for tsk in ts_tsk_lst:
        if 'ene' in tsk[0]:
            if ene_idx > len(ene_coeff)-1:
                print('Warning - an insufficient ',
                      'energy coefficient list was provided')
                break
            ene_lvl = tsk[1]
            ene_lvl_ref = tsk[2]
            ene_ref_thy_info = finf.get_thy_info(ene_lvl_ref)
            ene_thy_info = finf.get_thy_info(ene_lvl)
            ene = get_high_level_energy(
                spc_info=spc_info,
                thy_low_level=ene_ref_thy_info,
                thy_high_level=ene_thy_info,
                save_prefix=save_path,
                saddle=saddle)
            spc_ene += ene*ene_coeff[ene_idx]
            ene_idx += 1

    return spc_ene


def get_ckin_ene_lvl_str(ts_tsk_lst, ene_coeff):
    """ Write the comment lines for the enrgy lvls for ckin
    """
    ene_strl = []
    ene_idx = 0
    ene_str = '! energy level:'
    for tsk in ts_tsk_lst:
        if 'ene' in tsk[0]:
            if ene_idx > len(ene_coeff)-1:
                print('Warning - an insufficient ',
                      'energy coefficient list was provided')
                break
            ene_lvl = tsk[1]
            ene_lvl_ref = tsk[2]
            ene_ref_thy_info = finf.get_thy_info(ene_lvl_ref)
            ene_thy_info = finf.get_thy_info(ene_lvl)
            ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
                ene_coeff[ene_idx],
                ene_thy_info[3],
                ene_thy_info[1],
                ene_thy_info[2],
                ene_ref_thy_info[3],
                ene_ref_thy_info[1],
                ene_ref_thy_info[2]))
            ene_idx += 1
    ene_str += '!               '.join(ene_strl)

    return ene_str


def write_channel_mess_strs(spc_dct, rxn_lst, pes_formula,
                            ts_model, pf_levels, multi_info, pst_params,
                            spc_save_fs, save_prefix,
                            idx_dct, mess_strs):
    """ Write all the MESS input file strings for the reaction channels
    """
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

    return well_str, bim_str, ts_str
