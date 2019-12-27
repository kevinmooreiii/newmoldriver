""" electronic structure drivers
"""

import automol.inchi
import automol.geom
import autofile.fs
import routines

# Calling the new libs
from lib.filesystem import build as fbuild
from lib.filesystem import inf as finf
from lib.reaction import wells as lwells
from lib.reaction import ts as lts
from lib import msg


def run(tsk_info_lst, es_dct, rxn_lst, spc_dct, run_prefix, save_prefix,
        vdw_params=[False, False, True],
        pst_params=[1.0, 6],
        rad_rad_ts='vtst'):
    """ driver for all electronic structure tasks
    """

    print("Tasks:\n", tsk_info_lst)
    # Prepare species queue
    spc_queue = []
    for _, rxn in enumerate(rxn_lst):
        reacs = rxn['reacs']
        prods = rxn['prods']

        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)

    spc_queue = list(dict.fromkeys(spc_queue))
    # removes duplicates

    # Prepare prefix filesystem
    fbuild.prefix_filesystem(run_prefix, save_prefix)

    # Loop over Tasks
    ts_found = []
    for tsk_info in tsk_info_lst:

        # Task information
        tsk = tsk_info[0]
        es_ini_key = tsk_info[2]
        es_run_key = tsk_info[1]
        overwrite = tsk_info[3]
        # Theory information
        ini_thy_info = finf.get_es_info(es_ini_key)
        thy_info = finf.get_es_info(es_run_key)

        # If task is to find the transition state, find all TSs for rxn lst
        if tsk in ('find_ts', 'find_vdw'):
            for ts in spc_dct:
                if 'ts_' in ts:
                    msg.sadpt_tsk_msg(
                        tsk, ts, spc_dct, thy_info, ini_thy_info)
                    ts_info = (spc_dct[ts]['ich'],
                               spc_dct[ts]['chg'],
                               spc_dct[ts]['mul'])
                    rxn_class = spc_dct[ts]['class']
                    if not rxn_class:
                        print('skipping reaction because type =', rxn_class)
                        continue
                    elif rad_rad_ts == 'pst':
                        print('skipping reaction because we are using PST')
                    ts_zma = spc_dct[ts]['original_zma']
                    dist_info = spc_dct[ts]['dist_info']
                    grid = spc_dct[ts]['grid']
                    bkp_data = spc_dct[ts]['bkp_data']
                    _, _, rxn_run_path, rxn_save_path = spc_dct[ts]['rxn_fs']
                    if 'ts' in tsk:
                        geo, _, final_dist = lts.find_ts(
                            spc_dct, spc_dct[ts], ts_info, ts_zma, rxn_class,
                            dist_info, grid, bkp_data, ini_thy_info, thy_info,
                            run_prefix, save_prefix, rxn_run_path,
                            rxn_save_path, overwrite,
                            pst_params=pst_params,
                            rad_rad_ts=rad_rad_ts)
                        spc_dct[ts]['dist_info'][1] = final_dist
                        angle = None
                        dist_name = dist_info[0]
                        if 'abstraction' in rxn_class or 'addition' in rxn_class:
                            brk_name = dist_info[3]
                            if dist_name and brk_name:
                                ts_bnd = automol.zmatrix.bond_idxs(
                                    ts_zma, dist_name)
                                brk_bnd = automol.zmatrix.bond_idxs(
                                    ts_zma, brk_name)
                                ang_atms = [0, 0, 0]
                                cent_atm = list(set(brk_bnd) & set(ts_bnd))
                                if cent_atm:
                                    ang_atms[1] = cent_atm[0]
                                    for idx in brk_bnd:
                                        if idx != ang_atms[1]:
                                            ang_atms[0] = idx
                                    for idx in ts_bnd:
                                        if idx != ang_atms[1]:
                                            ang_atms[2] = idx

                                geom = automol.zmatrix.geometry(ts_zma)
                                angle = automol.geom.central_angle(
                                    geom, *ang_atms)
                        spc_dct[ts]['dist_info'].append(angle)
                        if not isinstance(geo, str):
                            print('Success, transition state',
                                  '{} added to species queue'.format(ts))
                            spc_queue.append(ts)
                            ts_found.append(ts)
                    elif 'vdw' in tsk:
                        pass
                        # vdws = lwells.find_vdw(
                        #     ts, spc_dct, thy_info, ini_thy_info, ts_info,
                        #     vdw_params,
                        #     es_dct[es_run_key]['mc_nsamp'], run_prefix,
                        #     save_prefix, 0.1, False,
                        #     substr.PROJROT, overwrite)
                        # spc_queue.extend(vdws)
            continue

        # Loop over all species
        for spc in spc_queue:
            if 'ts_' in spc and rad_rad_ts != 'pst':
                msg.tsk_msg(tsk, thy_info, ini_thy_info, spc)
                spc_run_fs, spc_save_fs, spc_run_path, spc_save_path = spc_dct[spc]['rxn_fs']
                spc_info = finf.get_spc_info(spc_dct[spc])

            else:
                msg.tsk_msg(tsk, thy_info, ini_thy_info, spc)
                spc_info = finf.get_spc_info(spc_dct[spc])
                spc_run_fs = autofile.fs.species(run_prefix)
                spc_run_fs.leaf.create(spc_info)
                spc_run_path = spc_run_fs.leaf.path(spc_info)

                spc_save_fs = autofile.fs.species(save_prefix)
                spc_save_fs.leaf.create(spc_info)
                spc_save_path = spc_save_fs.leaf.path(spc_info)

            # add orb resti
            orb_restr = routines.util.orbital_restriction(
                spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_save_fs = autofile.fs.theory(spc_save_path)

            if 'ene' not in tsk and 'hess' not in tsk:
                if 'ts_' in spc:
                    thy_run_fs.leaf.create(thy_level[1:4])
                    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
                    thy_save_fs.leaf.create(thy_level[1:4])
                    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

                    thy_run_fs = autofile.fs.ts(thy_run_path)
                    thy_run_fs.trunk.create()
                    thy_run_path = thy_run_fs.trunk.path()

                    thy_save_fs = autofile.fs.ts(thy_save_path)
                    thy_save_fs.trunk.create()
                    thy_save_path = thy_save_fs.trunk.path()

                else:
                    thy_run_fs.leaf.create(thy_level[1:4])
                    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
                    thy_save_fs.leaf.create(thy_level[1:4])
                    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

                cnf_run_fs = autofile.fs.conformer(thy_run_path)
                cnf_save_fs = autofile.fs.conformer(thy_save_path)
                tau_run_fs = autofile.fs.tau(thy_run_path)
                tau_save_fs = autofile.fs.tau(thy_save_path)
                min_cnf_locs = routines.util.min_energy_conformer_locators(
                    cnf_save_fs)
                if min_cnf_locs:
                    min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
                    min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
                    scn_run_fs = autofile.fs.conformer(min_cnf_run_path)
                    scn_save_fs = autofile.fs.conformer(min_cnf_save_path)
                else:
                    scn_run_fs = None
                    scn_save_fs = None
            else:
                cnf_run_fs = None
                cnf_save_fs = None
                tau_run_fs = None
                tau_save_fs = None
                scn_run_fs = None
                scn_save_fs = None

            if ini_thy_info[0] != 'input_geom':
                orb_restr = routines.util.orbital_restriction(
                    spc_info, ini_thy_info)
                ini_thy_level = ini_thy_info[0:3]
                ini_thy_level.append(orb_restr)

                ini_thy_run_fs = autofile.fs.theory(spc_run_path)
                ini_thy_save_fs = autofile.fs.theory(spc_save_path)
                if 'ts_' in spc:
                    ini_thy_run_fs.leaf.create(ini_thy_level[1:4])
                    ini_thy_run_path = ini_thy_run_fs.leaf.path(
                        ini_thy_level[1:4])
                    ini_thy_save_fs.leaf.create(ini_thy_level[1:4])
                    ini_thy_save_path = ini_thy_save_fs.leaf.path(
                        ini_thy_level[1:4])

                    ini_thy_run_fs = autofile.fs.ts(ini_thy_run_path)
                    ini_thy_run_fs.trunk.create()
                    ini_thy_run_path = ini_thy_run_fs.trunk.path()

                    ini_thy_save_fs = autofile.fs.ts(ini_thy_save_path)
                    ini_thy_save_fs.trunk.create()
                    ini_thy_save_path = ini_thy_save_fs.trunk.path()

                else:
                    ini_thy_run_fs.leaf.create(ini_thy_level[1:4])
                    ini_thy_run_path = ini_thy_run_fs.leaf.path(
                        ini_thy_level[1:4])
                    ini_thy_save_fs.leaf.create(ini_thy_level[1:4])
                    ini_thy_save_path = ini_thy_save_fs.leaf.path(
                        ini_thy_level[1:4])

                ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
                ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)

                ini_tau_run_fs = autofile.fs.tau(ini_thy_run_path)
                ini_tau_save_fs = autofile.fs.tau(ini_thy_save_path)
                min_cnf_locs = routines.util.min_energy_conformer_locators(
                    ini_cnf_save_fs)
                if min_cnf_locs:
                    min_cnf_run_path = ini_cnf_run_fs.leaf.path(min_cnf_locs)
                    min_cnf_save_path = ini_cnf_save_fs.leaf.path(min_cnf_locs)
                    ini_scn_run_fs = autofile.fs.conformer(min_cnf_run_path)
                    ini_scn_save_fs = autofile.fs.conformer(min_cnf_save_path)
                else:
                    ini_scn_run_fs = None
                    ini_scn_save_fs = None

            else:
                ini_thy_run_fs = None
                ini_thy_run_path = None
                ini_thy_save_fs = None
                ini_thy_save_path = None
                ini_cnf_run_fs = None
                ini_cnf_save_fs = None
                ini_tau_run_fs = None
                ini_tau_save_fs = None
                ini_thy_level = ini_thy_info
                ini_scn_run_fs = None
                ini_scn_save_fs = None

            run_fs = autofile.fs.run(thy_run_path)
            run_fs.trunk.create()

            filesys = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                       cnf_run_fs, cnf_save_fs, tau_run_fs, tau_save_fs,
                       scn_run_fs, scn_save_fs, run_fs]

            ini_filesys = [ini_thy_run_fs, ini_thy_save_fs, ini_cnf_run_fs,
                           ini_cnf_save_fs, ini_tau_run_fs, ini_tau_save_fs,
                           ini_scn_run_fs, ini_scn_save_fs]

            # Run tasks
            if 'ts_' in spc:
                saddle = True
                if any(string in tsk for string in ('samp', 'scan', 'geom')):
                    routines.es.geometry_generation(
                        tsk, spc_dct[spc], spc_info,
                        es_dct[es_run_key]['mc_nsamp'],
                        ini_thy_level, thy_level,
                        ini_filesys, filesys, overwrite,
                        saddle=saddle)
                else:
                    selection = 'min'
                    routines.es.geometry_analysis(
                        tsk, thy_level, ini_filesys, selection,
                        spc_info, spc_dct[spc], overwrite)
            else:
                saddle = False
                if 'samp' in tsk or 'scan' in tsk or 'geom' in tsk:
                    if 'vdw_' not in spc:
                        routines.es.geometry_generation(
                            tsk, spc_dct[spc], spc_info,
                            es_dct[es_run_key]['mc_nsamp'],
                            ini_thy_level, thy_level,
                            ini_filesys, filesys, overwrite,
                            saddle=saddle)
                    else:
                        lwells.fake_geo_gen(
                            tsk, spc_dct[spc], es_dct[es_run_key], thy_level,
                            filesys, spc_info, overwrite)
                else:
                    selection = 'min'
                    if 'conf' in tsk:
                        min_cnf_locs = routines.util.min_energy_conformer_locators(
                            ini_cnf_save_fs)
                        if not min_cnf_locs:
                            msg.ini_info_noavail_msg(tsk)
                            continue
                        elif not ini_cnf_save_fs.leaf.file.geometry.exists(
                                 min_cnf_locs):
                            msg.ini_info_noavail_msg(tsk)
                            continue
                    elif 'tau' in tsk:
                        tau_locs = ini_tau_save_fs.leaf.existing()
                        if not tau_locs:
                            msg.ini_info_noavail_msg(tsk)
                            continue
                        elif not ini_tau_save_fs.leaf.file.geometry.exists(
                                 [tau_locs[0]]):
                            msg.ini_info_noavail_msg(tsk)
                            continue
                    elif 'scan' in tsk:
                        scn_locs = ini_scn_save_fs.leaf.existing()
                        if not scn_locs:
                            msg.ini_info_noavail_msg(tsk)
                            continue
                        elif not ini_scn_save_fs.leaf.file.geometry.exists(
                                 [scn_locs[0]]):
                            msg.ini_info_noavail_msg(tsk)
                            continue
                    routines.es.geometry_analysis(
                        tsk, thy_level, ini_filesys, selection,
                        spc_info, spc_dct[spc], overwrite, saddle=saddle)
    return ts_found
