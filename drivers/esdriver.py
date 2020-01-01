""" electronic structure drivers
"""

import automol.inchi
import automol.geom
import autofile.fs
import routines

# Calling the new libs
from lib.filesystem import build as fbuild
from lib.filesystem import inf as finf
from lib.filesystem import check as fcheck
from lib.submission import theolvls
from lib.reaction import wells as lwells
from lib.reaction import ts as lts
from lib import msg


def run(tsk_info_lst, rxn_lst, spc_dct, run_prefix, save_prefix,
        vdw_params=[False, False, True],
        pst_params=[1.0, 6],
        rad_rad_ts='vtst',
        mc_nsamp=[True, 10, 1, 3, 100]):
    """ driver for all electronic structure tasks
    """

    # Set the es_dct
    es_dct = theolvls.ES_DCT

    # Prepare species queue
    spc_queue = []
    for _, rxn in enumerate(rxn_lst):
        reacs = rxn['reacs']
        prods = rxn['prods']
        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)
    spc_queue = list(dict.fromkeys(spc_queue))

    # Prepare prefix filesystem
    fbuild.prefix_filesystem(run_prefix, save_prefix)

    # Loop over Tasks
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
            for sadpt in spc_dct:
                if 'ts_' in sadpt:
                    msg.sadpt_tsk_msg(
                        tsk, sadpt, spc_dct, thy_info, ini_thy_info)
                    ts_info = (spc_dct[sadpt]['ich'],
                               spc_dct[sadpt]['chg'],
                               spc_dct[sadpt]['mul'])
                    rxn_class = spc_dct[sadpt]['class']
                    if not rxn_class:
                        print('skipping reaction because type =', rxn_class)
                        continue
                    elif rad_rad_ts == 'pst':
                        print('using PST')

                    # Find the transition state
                    ts_zma = spc_dct[sadpt]['original_zma']
                    dist_info = spc_dct[sadpt]['dist_info']
                    grid = spc_dct[sadpt]['grid']
                    bkp_data = spc_dct[sadpt]['bkp_data']
                    _, _, rxn_run_path, rxn_save_path = spc_dct[sadpt]['rxn_fs']
                    if 'ts' in tsk:
                        geo, _, final_dist = lts.find_ts(
                            spc_dct, spc_dct[sadpt], ts_info, ts_zma, rxn_class,
                            dist_info, grid, bkp_data, ini_thy_info, thy_info,
                            run_prefix, save_prefix, rxn_run_path,
                            rxn_save_path, overwrite,
                            pst_params=pst_params,
                            rad_rad_ts=rad_rad_ts)

                        # Add an angle check which is added to spc dct for TS
                        angle = check_angle(ts_zma, dist_info, rxn_class)
                        spc_dct[sadpt]['dist_info'][1] = final_dist
                        spc_dct[sadpt]['dist_info'].append(angle)

                        # Add TS to species queue if TS is found
                        if not isinstance(geo, str):
                            print('Success, transition state',
                                  '{} added to species queue'.format(sadpt))
                            spc_queue.append(sadpt)
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

            msg.tsk_msg(tsk, thy_info, ini_thy_info, spc)
            # Build the species filesys objs
            spc_fs = set_spc_fs(spc_dct, spc, run_prefix, save_prefix, rad_rad_ts)
            spc_info, spc_run_fs, spc_save_fs, spc_run_path, spc_save_path = spc_fs

            # Build the main run filesystem objects
            setfs_chk = True
            # setfs_chk = bool('ene' not in tsk and 'hess' not in tsk)
            filesys, thy_level = set_fs(
                spc, spc_info, thy_info,
                spc_run_fs, spc_save_fs,
                spc_run_path, spc_save_path,
                setfs_chk, ini_fs=False)

            # Build the ini run filesystem objects
            setfs_chk = bool(ini_thy_info[0] != 'input_geom')
            ini_filesys, ini_thy_level = set_fs(
                spc, spc_info, ini_thy_info,
                spc_run_fs, spc_save_fs,
                spc_run_path, spc_save_path,
                setfs_chk, ini_fs=True)

            # Run tasks
            if 'ts_' in spc:
                if any(string in tsk for string in ('samp', 'scan', 'geom')):
                    routines.es.geometry_generation(
                        tsk, spc_dct[spc], spc_info,
                        mc_nsamp,
                        ini_thy_level, thy_level,
                        ini_filesys, filesys, overwrite,
                        saddle=True)
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
                            saddle=False)
                    else:
                        lwells.fake_geo_gen(
                            tsk, spc_dct[spc], es_dct[es_run_key], thy_level,
                            filesys, spc_info, overwrite)
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


def set_fs(spc, spc_info, thy_info,
           spc_run_fs, spc_save_fs,
           spc_run_path, spc_save_path,
           setfs_chk, ini_fs=False):
    """ set up filesystem """

    # Initialize the filesystems
    thy_run_fs = None
    thy_run_path = None
    thy_save_fs = None
    thy_save_path = None
    cnf_run_fs = None
    cnf_save_fs = None
    tau_run_fs = None
    tau_save_fs = None
    thy_level = thy_info
    scn_run_fs = None
    scn_save_fs = None

    # Build the theory filesys obj
    thy_level = thy_info
    orb_restr = routines.util.orbital_restriction(
        spc_info, thy_info)
    thy_level = thy_info[0:3]
    thy_level.append(orb_restr)
    thy_run_fs = autofile.fs.theory(spc_run_path)
    thy_save_fs = autofile.fs.theory(spc_save_path)

    if setfs_chk:
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

    filesys = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
               cnf_run_fs, cnf_save_fs, tau_run_fs, tau_save_fs,
               scn_run_fs, scn_save_fs]

    # Add run fs if needed
    if not ini_fs:
        filesys = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                   cnf_run_fs, cnf_save_fs, tau_run_fs, tau_save_fs,
                   scn_run_fs, scn_save_fs]
        run_fs = autofile.fs.run(thy_run_path)
        run_fs.trunk.create()
        filesys.append(run_fs)
    else:
        filesys = [thy_run_fs, thy_save_fs, cnf_run_fs,
                   cnf_save_fs, tau_run_fs, tau_save_fs,
                   scn_run_fs, scn_save_fs]

    return filesys, thy_level


def check_angle(ts_zma, dist_info, rxn_class):
    """ Check the angle to amend the dct
    """
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

    return angle


def set_spc_fs(spc_dct, spc, run_prefix, save_prefix, rad_rad_ts):
    """ set the species filesys objects
    """
    if 'ts_' in spc and rad_rad_ts != 'pst':
        spc_run_fs, spc_save_fs, spc_run_path, spc_save_path = spc_dct[spc]['rxn_fs']
        spc_info = finf.get_spc_info(spc_dct[spc])
    else:
        spc_info = finf.get_spc_info(spc_dct[spc])
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs.leaf.create(spc_info)
        spc_run_path = spc_run_fs.leaf.path(spc_info)

        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

    return spc_info, spc_run_fs, spc_save_fs, spc_run_path, spc_save_path
