""" electronic structure drivers
"""
import os
import automol.inchi
import automol.geom
import moldr
import autofile.fs
import scripts.es
from esdriver.load import load_params
from submission import substr

# driver libs
import lib.drivers.es as lesdriver
import lib.filesystem.build as lfs


KICKOFF_SIZE = 0.1
KICKOFF_BACKWARD = False

def run(tsk_info_lst, es_dct, rxn_lst, spc_dct, run_prefix, save_prefix,
        vdw_params=[False, False, True],
        pst_params=[1.0, 6],
        rad_rad_ts='vtst'):
    """ driver for all electronic structure tasks
    """

    # Print the user-specified task list to the log file
    lesdriver.print_tsk_lst(tsk_info_lst)

    # Build a queue of species to run calculations for
    lesdriver.build_species_queue(rxn_lst)

    # Prepare filesystem
    lfs.build_filesystem(run_prefix, save_prefix)

    # Loop over Tasks
    ts_found = []
    for tsk_info in tsk_info_lst:

        # Task and theory information
        [tsk, es_run_key, es_ini_key, overwrite] = tsk_info
        ini_thy_info = get_es_info(es_dct, es_ini_key)
        thy_info = get_es_info(es_dct, es_run_key)

        # If task is to find transition state, find all TSs for reaction list
        if tsk in ('find_ts', 'find_vdw'):
            for ts in spc_dct:
                find_ts(ts)
            continue

        # Loop over all species
        for spc in spc_queue:

            # Set the species run-save paths
            spc_run_fs, spc_save_fs, spc_run_path, spc_save_path = spc_fs(
                spc, rad_rad_ts)

            # Create the theory filesystem levels and objects
            thy_level = append_orb_restr_to_thy_lvl(thy_level, thy_info, spc_info)
            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_save_fs = autofile.fs.theory(spc_save_path)

            if 'ene' not in tsk and 'hess' not in tsk:
                thy_run_path, thy_save_path = thy_fs_paths(
                    spc, thy_run_fs, thy_save_fs, thy_level)

                cnf_run_fs = autofile.fs.conformer(thy_run_path)
                cnf_save_fs = autofile.fs.conformer(thy_save_path)
                tau_run_fs = autofile.fs.tau(thy_run_path)
                tau_save_fs = autofile.fs.tau(thy_save_path)
                min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
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

                # Create the theory filesystem levels and objects
                ini_thy_level = append_orb_restr_to_thy_lvl(
                    thy_level, ini_thy_info, spc_info)
                ini_thy_run_fs = autofile.fs.theory(spc_run_path)
                ini_thy_save_fs = autofile.fs.theory(spc_save_path)
                
                ini_thy_run_path, ini_thy_save_path = thy_fs_paths(
                    spc, ini_thy_run_fs, ini_thy_save_fs, ini_thy_level)

                ini_cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
                ini_cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)

                ini_tau_run_fs = autofile.fs.tau(ini_thy_run_path)
                ini_tau_save_fs = autofile.fs.tau(ini_thy_save_path)
                min_cnf_locs = moldr.util.min_energy_conformer_locators(ini_cnf_save_fs)
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

            fs = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                  cnf_run_fs, cnf_save_fs, tau_run_fs, tau_save_fs,
                  scn_run_fs, scn_save_fs, run_fs]

            ini_fs = [ini_thy_run_fs, ini_thy_save_fs, ini_cnf_run_fs,
                      ini_cnf_save_fs, ini_tau_run_fs, ini_tau_save_fs,
                      ini_scn_run_fs, ini_scn_save_fs]

            #Run tasks
            if 'ts_' in spc:
                if 'samp' in tsk or 'scan' in tsk or 'geom' in tsk:
                    geo = moldr.ts.reference_geometry(
                        spc_dct[spc], thy_level, ini_thy_level, fs, ini_fs,
                        spc_dct[spc]['dist_info'], overwrite)
                    if geo:
                        scripts.es.ts_geometry_generation(
                            tsk, spc_dct[spc], es_dct[es_run_key],
                            thy_level, fs, spc_info, overwrite)
                else:
                    selection = 'min'
                    scripts.es.ts_geometry_analysis(
                        tsk, thy_level, ini_fs, selection, spc_info, spc_dct[spc], overwrite)
            else:
                if 'samp' in tsk or 'scan' in tsk or 'geom' in tsk:
                    geo = moldr.geom.reference_geometry(
                        spc_dct[spc], thy_level, ini_thy_level, fs, ini_fs,
                        kickoff_size=KICKOFF_SIZE,
                        kickoff_backward=KICKOFF_BACKWARD,
                        projrot_script_str=substr.PROJROT,
                        overwrite=overwrite)
                    if geo:
                        if not 'vdw_' in spc:
                            scripts.es.geometry_generation(
                                tsk, spc_dct[spc], es_dct[es_run_key], thy_level,
                                fs, spc_info, overwrite)
                        else:
                            scripts.es.fake_geo_gen(
                                tsk, spc_dct[spc], es_dct[es_run_key], thy_level,
                                fs, spc_info, overwrite)
                else:
                    selection = 'min'
                    if 'conf' in tsk:
                        min_cnf_locs = moldr.util.min_energy_conformer_locators(ini_cnf_save_fs)
                        if not locs_found(min_cnf_locs, ini_cnf_save_fs, 'conf', tsk):
                            continue
                    elif 'tau' in tsk:
                        tau_locs = [ini_tau_save_fs.leaf.existing()[0]]
                        if not locs_found(tau_locs, ini_tau_save_fs, 'tau', tsk):
                            continue
                    elif 'scan' in tsk:
                        scn_locs = [ini_scn_save_fs.leaf.existing()[0]]
                        if not locs_found(tau_locs, ini_tau_save_fs, 'scan', tsk):
                            continue
                    scripts.es.geometry_analysis(tsk, thy_level, ini_fs,
                            selection, spc_info, overwrite)
    return ts_found


### Functions

def print_tsk_lst(tsk_info_lst):
    """ print tasks
    """
    print("Tasks:\n", tsk_info_lst)


def build_species_queue(rxn_lst):
    """ prepare species queue
    """
    spc_queue = []
    for _, rxn in enumerate(rxn_lst):
        reacs = rxn['reacs']
        prods = rxn['prods']

        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)

    # Removes duplicates?
    spc_queue = list(dict.fromkeys(spc_queue))

    return spc_queue


def find_ts(spc_dct, ts):
    """ find ts and associated vdw complexes
    """
    if 'ts_' in ts:
        print('Task {} \t for {} \t {}//{} \t {} = {}'.format(
            tsk, ts, '/'.join(thy_info), '/'.join(ini_thy_info),
            '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))
        ts_info = (spc_dct[ts]['ich'], spc_dct[ts]['chg'], spc_dct[ts]['mul'])
        rxn_class = spc_dct[ts]['class']
        if not rxn_class:
            print('skipping reaction because type =', rxn_class)
            continue
        elif rad_rad_ts == 'pst':
            print('skipping reaction because we are using PST')
            continue
        ts_zma = spc_dct[ts]['original_zma']
        dist_info = spc_dct[ts]['dist_info']
        grid = spc_dct[ts]['grid']
        bkp_data = spc_dct[ts]['bkp_data']
        _, _, rxn_run_path, rxn_save_path = spc_dct[ts]['rxn_fs']
        if 'ts' in tsk:
            geo, _, final_dist = scripts.es.find_ts(
                spc_dct, spc_dct[ts], ts_info, ts_zma, rxn_class,
                dist_info, grid, bkp_data, ini_thy_info, thy_info,
                run_prefix, save_prefix, rxn_run_path,
                rxn_save_path, overwrite,
                pst_params=pst_params,
                rad_rad_ts=rad_rad)

            spc_dct[ts]['dist_info'][1] = final_dist
            # Check the angle for addn and abst TS
            angle = check_abst_add_angle(dist_info, rxn_class, ts_zma)
            spc_dct[ts]['dist_info'].append(angle)
            if not isinstance(geo, str):
                print('Success, transition state {} added to species queue'.format(ts))
                spc_queue.append(ts)
                ts_found.append(ts)
        elif 'vdw' in tsk:
            vdws = scripts.es.find_vdw(
                ts, spc_dct, thy_info, ini_thy_info, ts_info, vdw_params,
                es_dct[es_run_key]['mc_nsamp'], run_prefix,
                save_prefix, KICKOFF_SIZE, KICKOFF_BACKWARD,
                substr.PROJROT, overwrite)
            spc_queue.extend(vdws)


def thy_fs_paths(spc, thy_run_fs, thy_save_fs, thy_level):
    """ build the theory run-save filesys paths
    """
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

    return thy_run_path, thy_save_path


def spc_fs(spc, spc_dct, thy_info, ini_thy_info, rad_rad_ts, run_prefix, save_prefix, tsk):
    """ build the species run-save filesystems
    """
    if 'ts_' in spc and rad_rad_ts != 'pst':
        print('\nTask {} \t {}//{} \t Species {}'.format(
            tsk, '/'.join(thy_info), '/'.join(ini_thy_info), spc))
        spc_run_fs, spc_save_fs, spc_run_path, spc_save_path = spc_dct[spc]['rxn_fs']
        spc_info = scripts.es.get_spc_info(spc_dct[spc])

    else:
        print('\nTask {} \t {}//{} \t Species {}: {}'.format(
            tsk, '/'.join(thy_info), '/'.join(ini_thy_info), spc,
            automol.inchi.smiles(spc_dct[spc]['ich'])))
        spc_info = scripts.es.get_spc_info(spc_dct[spc])
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs.leaf.create(spc_info)
        spc_run_path = spc_run_fs.leaf.path(spc_info)

        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

    return spc_run_fs, spc_save_fs, spc_run_path, spc_save_path


def append_orb_restr_to_thy_lvl(thy_level, thy_info, spc_info):
    """ append orb restricted to the theory level
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_level = thy_info[0:3]
    thy_level.append(orb_restr)

    return thy_level


def locs_found(locs, fs, string, tsk):
    """ Check if theory level has been set
    """
    found = True
    if not locs:
        print(
            'Initial level of theory for {} must be ',
            'run before {} '.format(string, tsk))
        found = False
    elif not fs.leaf.file.geometry.exists(locs):
        print(
            'Initial level of theory for {} must be ',
            'run before {} '.format(string, tsk))
        found = False
    return found


def check_abst_add_angle(dist_info, rxn_class, ts_zma):
    """ gets the angle for the TS for abstractions and additions
    """
    angle = None
    dist_name = dist_info[0]
    if 'abstraction' in rxn_class or 'addition' in rxn_class:
        brk_name = dist_info[3]
        if dist_name and brk_name:
            ts_bnd = automol.zmatrix.bond_idxs(ts_zma, dist_name)
            brk_bnd = automol.zmatrix.bond_idxs(ts_zma, brk_name)
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
            angle = automol.geom.central_angle(geom, *ang_atms)

    return angle


