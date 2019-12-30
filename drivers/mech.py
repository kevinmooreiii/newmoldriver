""" estoktp driver functions with no
better home
"""

import os
import sys
import numpy
import ktpdriver

from routines import util
from lib.phydat import phycon


def run_ktp_driver(pes_dct, pesnums_lst, channels, connchnls_lst,
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
                   mc_nsamp=[True, 10, 1, 3, 100]):
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
                    nums = channels.replace('[', '').replace(']', '').split(',')
                    pes_chns = [int(num) for num in nums]

            # Print out the channel that is being run
            for cidx, cvals in enumerate(connchnls_lst[pes_idx-1].values()):
                print('PES{}_{}: {} surface'.format(
                    str(pes_idx), str(cidx+1), pes))
                for cval in cvals:
                    print(pes_rxn_name_lst[cval])

            # Run the KTPDriver over the channels
            for cidx, cvals in enumerate(connchnls_lst[pes_idx-1].values()):
                print('ktp on PES{}_{}: {} for the following channels...'.format(
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
                        if rxn_name_lst[idx] in cla_dct:
                            spc_dct[tsname]['given_class'] = cla_dct[rxn_name_lst[idx]]
                        elif '='.join(rxn_name_lst[idx].split('=')[::-1]) in cla_dct:
                            spc_dct[tsname]['given_class'] = cla_dct['='.join(rxn_name_lst[idx].split('=')[::-1])]
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
                        ts_mul_low, ts_mul_high, rad_rad = util.ts_mul_from_reaction_muls(
                            rct_names_lst[idx], prd_names_lst[idx], spc_dct)
                        spc_dct[tsname]['mul'] = ts_mul_low
                        spc_dct[tsname]['rad_rad'] = rad_rad
                        spc_dct[tsname]['hind_inc'] = (
                            hind_inc * phycon.DEG2RAD)
                        ts_idx += 1

                    # Run the ktpdriver
                    ktpdriver.run(
                        spc_dct,
                        tsk_info_lst,
                        pes_rct_names_lst,
                        pes_prd_names_lst,
                        run_prefix,
                        save_prefix,
                        ene_coeff=ene_coeff,
                        vdw_params=vdw_params,
                        options=options,
                        etrans=etrans,
                        pst_params=pst_params,
                        rad_rad_ts=rad_rad_ts,
                        mc_nsamp=mc_nsamp)


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
    geom_dct = util.geometry_dictionary(geom_path)
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


def run_spc_esdriver(rct_names_lst, prd_names_lst, tsk_info_lst,
                     spc_dct, run_prefix, save_prefix, vdw_params,
                     pst_params=[1.0, 6],
                     rad_rad_ts='pst'):
    """ Call the ESDriver Routines for species in the spc dct
    """
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

    runspecies = [{'species': spc_queue, 'reacs': [], 'prods': []}]
    esdriver.run(
        spc_tsk_lst, runspecies, spc_dct,
        run_prefix, save_prefix, vdw_params,
        pst_params=pst_params,
        rad_rad_ts=rad_rad_ts)


def run_rxn_esdriver(rct_names_lst, prd_names_lst, tsk_info_lst,
                     spc_dct, run_prefix, save_prefix, vdw_params,
                     pst_params=[1.0, 6],
                     rad_rad_ts='pst'):
    """ Call the ESDriver Routines for reactions in the spc dct
    """
    rxn_lst = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_lst.append(
            {'species': [], 'reacs': list(rct_names_lst[rxn]), 'prods':
             list(prd_names_lst[rxn])})

    ts_tsk_lst = []
    ts_tsk = False
    for tsk in tsk_info_lst:
        if 'find_ts' in tsk[0]:
            ts_tsk = True
        if ts_tsk:
            ts_tsk_lst.append(tsk)

    # Add addtional dictionary items for all the TSs
    # This presumes that es has been run previously for species list
    # to produce energies in save file system
    if ts_tsk_lst:
        print('\nBegin transition state prep')
        for ts in spc_dct:
            if 'ts_' in ts:
                # Exothermicity reordering requires electronic energy which requires theory level and
                # tsk_info
                es_ini_key = ts_tsk_lst[0][2]
                es_run_key = ts_tsk_lst[0][1]
                ini_thy_info = finf.get_es_info(es_ini_key)
                thy_info = finf.get_es_info(es_run_key)
                # generate rxn data, reorder if necessary, and put in spc_dct for given ts
                rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul = finf.rxn_info(
                    run_prefix, save_prefix, ts, spc_dct, thy_info, ini_thy_info)
                spc_dct[ts]['rxn_ichs'] = rxn_ichs
                spc_dct[ts]['rxn_chgs'] = rxn_chgs
                spc_dct[ts]['rxn_muls'] = rxn_muls
                spc_dct[ts]['low_mul'] = low_mul
                spc_dct[ts]['high_mul'] = high_mul
                # generate rxn_fs from rxn_info stored in spc_dct
                rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path = fpath.get_rxn_fs(
                    run_prefix, save_prefix, spc_dct[ts])
                spc_dct[ts]['rxn_fs'] = [rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path]

                rct_zmas, prd_zmas, rct_cnf_save_fs, prd_cnf_save_fs = fread.get_zmas(
                    spc_dct[ts]['reacs'], spc_dct[ts]['prods'], spc_dct,
                    ini_thy_info, save_prefix, run_prefix, KICKOFF_SIZE,
                    KICKOFF_BACKWARD, substr.PROJROT)
                ret = ts_class(
                    rct_zmas, prd_zmas, spc_dct[ts]['rad_rad'],
                    spc_dct[ts]['mul'], low_mul, high_mul,
                    rct_cnf_save_fs, prd_cnf_save_fs, spc_dct[ts]['given_class'])
                ret1, ret2 = ret
                if ret1:
                    rxn_class, ts_zma, dist_name, brk_name, grid, frm_bnd_key, brk_bnd_key, tors_names, update_guess = ret1
                    spc_dct[ts]['class'] = rxn_class
                    spc_dct[ts]['grid'] = grid
                    spc_dct[ts]['tors_names'] = tors_names
                    spc_dct[ts]['original_zma'] = ts_zma
                    dist_info = [dist_name, 0., update_guess, brk_name]
                    spc_dct[ts]['dist_info'] = dist_info
                    spc_dct[ts]['frm_bnd_key'] = frm_bnd_key
                    spc_dct[ts]['brk_bnd_key'] = brk_bnd_key
                    # Adding in the rct and prd zmas for vrctst
                    spc_dct[ts]['rct_zmas'] = rct_zmas
                    spc_dct[ts]['prd_zmas'] = prd_zmas
                    if ret2:
                        spc_dct[ts]['bkp_data'] = ret2
                    else:
                        spc_dct[ts]['bkp_data'] = None
                else:
                    spc_dct[ts]['class'] = None
                    spc_dct[ts]['bkp_data'] = None

        print('End transition state prep\n')

        # Run ESDriver
        ts_found = esdriver.run(
            ts_tsk_lst, rxn_lst, spc_dct,
            run_prefix, save_prefix, vdw_params)
