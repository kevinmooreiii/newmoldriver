""" electronic structure drivers
"""

from drivers import mech as lmech
import routines
from lib.filesystem import build as fbuild
from lib.filesystem import inf as finf
from lib.filesystem import check as fcheck
from lib.filesystem import path as fpath
from lib.submission import substr
from lib.submission import theolvls
from lib.reaction import ts as lts
from lib import msg


def run(tsk_info_lst, rxn_lst, spc_dct, run_prefix, save_prefix,
        vdw_params=(False, False, True),
        rad_rad_ts='vtst',
        mc_nsamp=(True, 10, 1, 3, 100),
        kickoff=(0.1, False)):
    """ driver for all electronic structure tasks
    """

    # Set the es_dct
    es_dct = theolvls.ES_DCT

    # Prepare species queue
    spc_queue = lmech.form_spc_queue_2(rxn_lst)

    # Prepare prefix filesystem
    fbuild.prefix_filesystem(run_prefix, save_prefix)

    # Loop over Tasks
    for tsk_info in tsk_info_lst:

        # Task and theory information
        [tsk, es_run_key, es_ini_key, overwrite] = tsk_info
        ini_thy_info = finf.get_es_info(es_ini_key)
        thy_info = finf.get_es_info(es_run_key)

        # If task is to find the transition state, find all TSs for rxn lst
        if tsk in ('find_ts', 'find_vdw'):
            for sadpt in spc_dct:
                if 'ts_' in sadpt:
                    msg.sadpt_tsk_msg(
                        tsk, sadpt, spc_dct, thy_info, ini_thy_info)
                    ts_info = (spc_dct[sadpt]['ich'],
                               spc_dct[sadpt]['chg'],
                               spc_dct[sadpt]['mul'])
                    if not spc_dct[sadpt]['class']:
                        print('skipping reaction because type =',
                              spc_dct[sadpt]['class'])
                        continue

                    # Find the transition state
                    [_, _,
                     rxn_run_path, rxn_save_path] = spc_dct[sadpt]['rxn_fs']
                    if 'ts' in tsk:
                        geo, _, final_dist = routines.es.find.find_ts(
                            spc_dct, spc_dct[sadpt], ts_info,
                            spc_dct[sadpt]['original_zma'],
                            spc_dct[sadpt]['class'],
                            spc_dct[sadpt]['dist_info'],
                            spc_dct[sadpt]['grid'],
                            spc_dct[sadpt]['bkp_data'], ini_thy_info, thy_info,
                            run_prefix, save_prefix,
                            rxn_run_path, rxn_save_path, overwrite,
                            rad_rad_ts=rad_rad_ts)

                        # Add an angle check which is added to spc dct for TS
                        angle = lts.check_angle(
                            spc_dct[sadpt]['original_zma'],
                            spc_dct[sadpt]['dist_info'],
                            spc_dct[sadpt]['class'])
                        spc_dct[sadpt]['dist_info'][1] = final_dist
                        spc_dct[sadpt]['dist_info'].append(angle)

                        # Add TS to species queue if TS is found
                        if not isinstance(geo, str):
                            print('Success, transition state',
                                  '{} added to species queue'.format(sadpt))
                            spc_queue.append(sadpt)
                    elif 'vdw' in tsk:
                        vdws = routines.es.wells.find_vdw(
                            sadpt, spc_dct, thy_info, ini_thy_info,
                            vdw_params,
                            es_dct[es_run_key]['mc_nsamp'], run_prefix,
                            save_prefix, 0.1, False,
                            substr.PROJROT, overwrite)
                        spc_queue.extend(vdws)
            continue

        # Loop over all species
        for spc in spc_queue:

            msg.tsk_msg(tsk, thy_info, ini_thy_info, spc)

            # Build the species filesys objs
            spc_fs = fpath.set_spc_fs(
                spc_dct, spc, run_prefix, save_prefix, rad_rad_ts)
            [spc_info,
             spc_run_fs, spc_save_fs,
             spc_run_path, spc_save_path] = spc_fs

            # Build the input and main run filesystem objects
            filesys, thy_level = fpath.set_fs(
                spc, spc_info, thy_info,
                spc_run_fs, spc_save_fs,
                spc_run_path, spc_save_path,
                setfs_chk=True, ini_fs=False)
            ini_filesys, ini_thy_level = fpath.set_fs(
                spc, spc_info, ini_thy_info,
                spc_run_fs, spc_save_fs,
                spc_run_path, spc_save_path,
                setfs_chk=bool(ini_thy_info[0] != 'input_geom'), ini_fs=True)

            # Run tasks
            if 'ts_' in spc:
                if any(string in tsk for string in ('samp', 'scan', 'geom')):
                    routines.es.geometry_generation(
                        tsk, spc_dct[spc], spc_info,
                        mc_nsamp,
                        ini_thy_level, thy_level,
                        ini_filesys, filesys, overwrite,
                        saddle=True, kickoff=kickoff)
                else:
                    routines.es.geometry_analysis(
                        tsk, thy_level, ini_filesys,
                        spc_info, spc_dct[spc], overwrite,
                        saddle=True, selection='min')
            else:
                if any(string in tsk for string in ('samp', 'scan', 'geom')):
                    if 'vdw_' not in spc:
                        routines.es.geometry_generation(
                            tsk, spc_dct[spc], spc_info,
                            mc_nsamp,
                            ini_thy_level, thy_level,
                            ini_filesys, filesys, overwrite,
                            saddle=False, kickoff=kickoff)
                    else:
                        routines.es.wells.fake_geo_gen(tsk, thy_level, filesys)
                else:
                    if 'conf' in tsk:
                        ini_cnf_save_fs = ini_filesys[5]
                        avail = fcheck.check_save(ini_cnf_save_fs, tsk, 'conf')
                    elif 'tau' in tsk:
                        ini_tau_save_fs = ini_filesys[7]
                        avail = fcheck.check_save(ini_tau_save_fs, tsk, 'tau')
                    elif 'scan' in tsk:
                        ini_scn_save_fs = ini_filesys[9]
                        avail = fcheck.check_save(ini_scn_save_fs, tsk, 'tau')
                    if avail:
                        routines.es.geometry_analysis(
                            tsk, thy_level, ini_filesys,
                            spc_info, spc_dct[spc], overwrite,
                            saddle=False, selection='min')
