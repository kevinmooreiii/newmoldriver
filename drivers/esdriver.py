""" electronic structure drivers
"""

import routines
from lib.filesystem import check as fcheck
from lib.filesystem import path as fpath
from lib import msg


def run(rxn_lst, spc_dct,
        es_tsk_lst,
        model_dct,
        run_options_dct, run_inp_dct):
    """ driver for all electronic structure tasks
    """
    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']
    # vdw_params = model_dct['options']['vdw_params']
    rad_rad_ts = model_dct['pf']['ts_barrierless']
    mc_nsamp = run_options_dct['mc_nsamp']
    kickoff = run_options_dct['kickoff']

    # Print the header message for the driver
    msg.program_header('es')

    # Species queue
    print('rxn_lst\n', rxn_lst)
    spc_queue = rxn_lst[0]['species']

    # Loop over Tasks
    for tsk_info in es_tsk_lst:

        # Task and theory information
        [tsk, thy_info, ini_thy_info, overwrite] = tsk_info

        # If task is to find the transition state, find all TSs for rxn lst
        if tsk in ('find_ts', 'find_vdw'):
            for sadpt in spc_dct:
                if 'ts_' in sadpt:
                    msg.sadpt_tsk_msg(
                        tsk, sadpt, spc_dct, thy_info, ini_thy_info)
                    if not spc_dct[sadpt]['class']:
                        print('skipping reaction because type =',
                              spc_dct[sadpt]['class'])
                        continue

                    # Find the transition state
                    print('es_dist_info\n', spc_dct[sadpt]['dist_info'])
                    if 'ts' in tsk:
                        geo, _, _ = routines.es.find.find_ts(
                            spc_dct, spc_dct[sadpt],
                            spc_dct[sadpt]['original_zma'],
                            ini_thy_info, thy_info,
                            run_prefix, save_prefix,
                            overwrite,
                            rad_rad_ts=rad_rad_ts)

                        # Add TS to species queue if TS is found
                        if not isinstance(geo, str):
                            print('Success, transition state',
                                  '{} added to species queue'.format(sadpt))
                            spc_queue.append(sadpt)
                    elif 'vdw' in tsk:
                        pass
                        # vdws = routines.es.wells.find_vdw(
                        #     sadpt, spc_dct, thy_info, ini_thy_info,
                        #     vdw_params,
                        #     thy_dct[es_run_key]['mc_nsamp'], run_prefix,
                        #     save_prefix, 0.1, False,
                        #     overwrite)
                        # spc_queue.extend(vdws)
            continue

        # Loop over all species
        for spc in spc_queue:

            print('spc in es queue', spc)
            msg.tsk_msg(tsk, thy_info, ini_thy_info, spc)

            # Build the input and main run filesystem objects
            filesys, thy_level = fpath.set_fs(
                spc_dct, spc, thy_info,
                run_prefix, save_prefix,
                setfs_chk=True, ini_fs=False)
            ini_filesys, ini_thy_level = fpath.set_fs(
                spc_dct, spc, ini_thy_info,
                run_prefix, save_prefix,
                setfs_chk=bool(ini_thy_info[0] != 'input_geom'),
                ini_fs=True)

            # Run tasks
            saddle = bool('ts_' in spc)
            vdw = bool('vdw' in spc)
            avail = True
            if any(string in tsk for string in ('samp', 'scan', 'geom')):
                if not vdw:
                    routines.es.geometry_generation(
                        tsk, spc_dct[spc], mc_nsamp,
                        ini_thy_level, thy_level, ini_filesys, filesys,
                        overwrite, saddle=False, kickoff=kickoff)
                else:
                    routines.es.wells.fake_geo_gen(tsk, thy_level, filesys)
            else:
                if 'conf' in tsk and not saddle:
                    ini_cnf_save_fs = ini_filesys[5]
                    avail = fcheck.check_save(ini_cnf_save_fs, tsk, 'conf')
                elif 'tau' in tsk and not saddle:
                    ini_tau_save_fs = ini_filesys[7]
                    avail = fcheck.check_save(ini_tau_save_fs, tsk, 'tau')
                elif 'scan' in tsk and not saddle:
                    ini_scn_save_fs = ini_filesys[9]
                    avail = fcheck.check_save(ini_scn_save_fs, tsk, 'tau')
                if avail:
                    routines.es.geometry_analysis(
                        tsk, thy_level, ini_filesys,
                        spc_dct[spc], overwrite,
                        saddle=saddle, selection='min')
