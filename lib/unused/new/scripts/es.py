""" cycle over electronic structure calls
for species, TSs, and vdw species
"""
import os
import numpy
import automol
import elstruct
import autofile
import varecof_io
import moldr
from datalibs import phycon

# new libs
import lib.variational.vtst as lvtst
import lib.variational.vrctst as lvrctst
import lib.variational.setup as lvsetup
import lib.reaction.grid as lgrid


def run_energy(params, kwargs):
    """ energy for geometry in fiven fs directory
    """
    print('running task {}'.format('energy'))
    moldr.sp.run_energy(**params, **kwargs)


def run_grad(params, kwargs):
    """ gradient for geometry in given fs directory
    """
    print('running task {}'.format('gradient'))
    moldr.sp.run_gradient(**params, **kwargs)


def run_hess(params, kwargs):
    """ hessian for geometry in given fs directory
    """
    print('running task {}'.format('hessian'))
    moldr.sp.run_hessian(**params, **kwargs)


def run_vpt2(params, kwargs):
    """ second order vibrational perturbation theory for minimum energy conformer
    """
    print('running task {}'.format('vpt2'))
    print('anharm may not be working in moldr right now')
    moldr.sp.run_vpt2(**params, **kwargs)


def run_irc(params, kwargs):
    """ irc for geometry in given fs directory
    """
    print('running task {}'.format('irc'))
    moldr.ts.run_irc(**params, **kwargs)


def run_single_conformer(spc_info, thy_level, fs, overwrite, saddle=False, dist_info=[]):
    """ generate single optimized geometry for randomly sampled initial torsional angles
    """
    mc_nsamp = [False, 0, 0, 0, 0, 1]
    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2])
    thy_save_fs = fs[3]
    two_stage = False
    if saddle:
        two_stage = True
    moldr.conformer.conformer_sampling(
        spc_info=spc_info,
        thy_level=thy_level,
        thy_save_fs=thy_save_fs,
        cnf_run_fs=fs[4],
        cnf_save_fs=fs[5],
        script_str=sp_script_str,
        overwrite=overwrite,
        nsamp_par=mc_nsamp,
        saddle=saddle,
        dist_info=dist_info,
        two_stage=two_stage,
        **kwargs,
    )


def run_conf_samp(fs, params, opt_kwargs):
    """ generate a set of conformer geometries and energies via
    random sampling over torsional coordinates following by optimization
    """
    if params['nsamp_par'][0]:
        print('running task {} with abcd of {}'.format(
            'conformer sampling', ' '.join(
                [str(x) for x in params['nsamp_par'][1:]])))
    else:
        print('running task {} for {:g} points'.format(
            'conformer sampling', params['nsamp_par'][5]))
    params['thy_save_fs'] = fs[3]
    params['cnf_run_fs'] = fs[4]
    params['cnf_save_fs'] = fs[5]
    moldr.conformer.conformer_sampling(**params, **opt_kwargs)


def run_hr_scan(filesys, params, opt_kwargs):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('hr'))
    params['cnf_run_fs'] = filesys[4]
    params['cnf_save_fs'] = filesys[5]
    moldr.scan.hindered_rotor_scans(**params, **opt_kwargs)


def run_tau_sampling(filesys, params, opt_kwargs):
    """ energies, gradients, and hessians,
    for set of arbitrarily sampled torsional coordinates
    with all other coordinates optimized
    """

    params['tau_run_fs'] = filesys[6]
    params['tau_save_fs'] = filesys[7]
    params['thy_save_fs'] = filesys[3]
    moldr.tau.tau_sampling(**params, **opt_kwargs)

    del params['thy_save_fs']
    del params['nsamp_par']

    moldr.tau.run_tau_gradients(**params, **opt_kwargs)
    moldr.tau.run_tau_hessians(**params, **opt_kwargs)


def geometry_generation(tsk, spcdic, es_dct, thy_level, fs,
                        spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    print('Task in geometry_generation:', tsk)
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2])
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': opt_script_str,
              'overwrite': overwrite}
    choose_function = {'conf_samp': 'run_conf_samp',
                       'tau_samp': 'run_tau_samp',
                       'hr_scan': 'run_hr_scan'}

    if tsk in ['conf_samp', 'tau_samp']:
        params['nsamp_par'] = es_dct['mc_nsamp']
    elif tsk in ['hr_scan']:
        if 'hind_inc' in spcdic:
            params['scan_increment'] = spcdic['hind_inc']
        else:
            params['scan_increment'] = 30. * phycon.DEG2RAD

    if tsk in choose_function:
        eval(choose_function[tsk])(fs, params, opt_kwargs)


def ts_geometry_generation(tsk, spcdic, es_dct, thy_level, fs,
                           spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    # fs[3] = fs[11]
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2], saddle=True)
    print('tsk test in ts_geometry_generation:', tsk)
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': opt_script_str,
              'saddle' :  True,
              'tors_names': spcdic['tors_names'],
              'overwrite': overwrite}
    choose_function = {'conf_samp': 'run_conf_samp',
                       'tau_samp': 'run_tau_samp',
                       'hr_scan': 'run_hr_scan'}

    if tsk in ['conf_samp', 'tau_samp']:
        params['nsamp_par'] = es_dct['mc_nsamp']
        params['dist_info'] = spcdic['dist_info']
        params['two_stage'] = True
    if tsk in ['conf_samp']:
        params['rxn_class'] = spcdic['class']
    elif tsk in ['hr_scan']:
        if 'hind_inc' in spcdic:
            params['scan_increment'] = spcdic['hind_inc']
        else:
            params['scan_increment'] = 30. * phycon.DEG2RAD
        params['frm_bnd_key'] = spcdic['frm_bnd_key']
        params['brk_bnd_key'] = spcdic['brk_bnd_key']
        print('key test in ts_geom gen:', params['frm_bnd_key'], params['brk_bnd_key'])

    if tsk in choose_function:
        eval(choose_function[tsk])(fs, params, opt_kwargs)


def geometry_analysis(tsk, thy_level, ini_fs, selection, spc_info, overwrite):
    """ run the specified electronic structure task
    for a set of geometries
    """

    print('Task in geometry_analysis:', tsk)
    # specify the fs for the runs
    if 'conf' in tsk:
        run_dir = ini_fs[2]
        save_dir = ini_fs[3]
    elif 'tau' in tsk:
        run_dir = ini_fs[4]
        save_dir = ini_fs[5]
    elif 'hr' in tsk:
        run_dir = ini_fs[6]
        save_dir = ini_fs[7]
    else:
        return
    if isinstance(selection, str):
        if selection == 'all':
            locs_lst = save_dir.leaf.existing()
        elif selection == 'min':
            locs_lst = [moldr.util.min_energy_conformer_locators(save_dir)]
    elif isinstance(selection, int):
        locs_lst = moldr.util.geom_sort(save_dir)[0:selection]
    else:
        locs_lst = selection

    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2])
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': sp_script_str,
              'overwrite': overwrite}
    choose_function = {'conf_energy': 'run_energy',
                       'tau_energy': 'run_energy',
                       'hr_energy': 'run_energy',
                       'mep_energy': 'run_energy',
                       'conf_grad': 'run_grad',
                       'tau_grad': 'run_grad',
                       'hr_grad': 'run_grad',
                       'mep_grad': 'run_grad',
                       'conf_hess': 'run_hess',
                       'tau_hess': 'run_hess',
                       'hr_hess': 'run_hess',
                       'mep_hess': 'run_hess',
                       'conf_vpt2': 'run_vpt2',
                       'tau_vpt2': 'run_vpt2',
                       'conf_reopt': 'run_reopt',
                       'tau_reopt': 'run_reopt',
                       'hr_reopt': 'run_reopt',
                       'mep_reopt': 'run_reopt'}

    # cycle over the locations
    if tsk in choose_function:
        task_call = eval(choose_function[tsk])
        for locs in locs_lst:
            if locs:
                params['geo_run_fs'] = run_dir
                params['geo_save_fs'] = save_dir
                params['locs'] = locs
                task_call(params, kwargs)
            else:
                print('No initial geometry available for {} on {}'.format(
                    spc_info[0], '/'.join(thy_level[1:3])))


def ts_geometry_analysis(tsk, thy_level, ini_fs, selection, spc_info, spc_dic, overwrite):
    """ run the specified electronic structure task
    for a set of geometries
    """

    print('Task in ts geometry_analysis:', tsk)
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'overwrite': overwrite}
    # specify the fs for the runs
    if 'conf' in tsk:
        run_dir = ini_fs[2]
        save_dir = ini_fs[3]
    elif 'tau' in tsk:
        run_dir = ini_fs[4]
        save_dir = ini_fs[5]
    elif 'hr' in tsk:
        run_dir = ini_fs[6]
        save_dir = ini_fs[7]
        params['frm_bnd_key'] = spc_dic['frm_bnd_key']
        params['brk_bnd_key'] = spc_dic['brk_bnd_key']
        print('key test in ts_geom anal:', params['frm_bnd_key'], params['brk_bnd_key'])
    else:
        return
    if isinstance(selection, str):
        if selection == 'all':
            locs_lst = save_dir.leaf.existing()
        elif selection == 'min':
            locs_lst = [moldr.util.min_energy_conformer_locators(save_dir)]
    elif isinstance(selection, int):
        locs_lst = moldr.util.geom_sort(save_dir)[0:selection]
    else:
        locs_lst = selection
    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2], saddle=True)
    params['script_str'] = sp_script_str
    choose_function = {'conf_energy': 'run_energy',
                       'tau_energy': 'run_energy',
                       'hr_energy': 'run_energy',
                       'mep_energy': 'run_energy',
                       'conf_grad': 'run_grad',
                       'tau_grad': 'run_grad',
                       'hr_grad': 'run_grad',
                       'mep_grad': 'run_grad',
                       'conf_hess': 'run_hess',
                       'tau_hess': 'run_hess',
                       'hr_hess': 'run_hess',
                       'mep_hess': 'run_hess',
                       'conf_vpt2': 'run_vpt2',
                       'tau_vpt2': 'run_vpt2',
                       'conf_reopt': 'run_reopt',
                       'tau_reopt': 'run_reopt',
                       'hr_reopt': 'run_reopt',
                       'mep_reopt': 'run_reopt'}

    # cycle over the locations
    if tsk in choose_function:
        task_call = eval(choose_function[tsk])
        for locs in locs_lst:
            if locs:
                params['geo_run_fs'] = run_dir
                params['geo_save_fs'] = save_dir
                params['locs'] = locs
                task_call(params, kwargs)
            else:
                print('No initial geometry available for {} on {}'.format(
                    spc_info[0], '/'.join(thy_level[1:3])))


def rxn_info(run_prefix, save_prefix, ts, spc_dct, thy_info, ini_thy_info=None):
    """ prepare rxn info and reverse reactants and products if reaction is endothermic
    """
    rxn_ichs = [[], []]
    rxn_chgs = [[], []]
    rxn_muls = [[], []]
    print('\n TS for {}: {} = {}'.format(
        ts, '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))
    reacs = spc_dct[ts]['reacs']
    prods = spc_dct[ts]['prods']
    print('ts dct', spc_dct[ts])
    for spc in reacs:
        rxn_ichs[0].append(spc_dct[spc]['ich'])
        rxn_chgs[0].append(spc_dct[spc]['chg'])
        rxn_muls[0].append(spc_dct[spc]['mul'])
    for spc in prods:
        rxn_ichs[1].append(spc_dct[spc]['ich'])
        rxn_chgs[1].append(spc_dct[spc]['chg'])
        rxn_muls[1].append(spc_dct[spc]['mul'])
    # check direction of reaction
    try:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, thy_info)
    except:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ini_thy_info)
    print('reaction is {:.2f} endothermic'.format(rxn_exo*phycon.EH2KCAL))
    if rxn_exo > 0 and not spc_dct[ts]['given_class']:
        rxn_ichs = rxn_ichs[::-1]
        rxn_chgs = rxn_chgs[::-1]
        rxn_muls = rxn_muls[::-1]
        spc_dct[ts]['reacs'] = prods
        spc_dct[ts]['prods'] = reacs
        print('Reaction will proceed as {}: {} = {}'.format(
            ts, '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))
        # print('ts search will be performed in reverse direction')

    # set up the filesystem
    rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)
    low_mul = min(
        automol.mult.ts._low(rxn_muls[0]), automol.mult.ts._low(rxn_muls[1]))
    high_mul = max(
        automol.mult.ts._high(rxn_muls[0]), automol.mult.ts._high(rxn_muls[1]))

    return rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul


def get_rxn_fs(run_prefix, save_prefix, ts):
    """ get filesystems for a reaction
    """
    rxn_ichs = ts['rxn_ichs']
    rxn_chgs = ts['rxn_chgs']
    rxn_muls = ts['rxn_muls']
    ts_mul = ts['mul']

    rxn_run_fs = autofile.fs.reaction(run_prefix)
    rxn_run_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_run_path = rxn_run_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    rxn_ichs = tuple(map(tuple, rxn_ichs))
    rxn_chgs = tuple(map(tuple, rxn_chgs))
    rxn_muls = tuple(map(tuple, rxn_muls))
    rxn_save_fs = autofile.fs.reaction(save_prefix)
    rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_save_path = rxn_save_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    return rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path


def get_zmas(
        reacs, prods, spc_dct, ini_thy_info, save_prefix, run_prefix,
        kickoff_size, kickoff_backward, projrot_script_str):
    """get the zmats for reactants and products using the initial level of theory
    """
    if len(reacs) > 2:
        ich = spc_dct[reacs[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        reacs = reacs[:-1]
    elif len(prods) > 2:
        ich = spc_dct[prods[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        prods = prods[:-1]
    rct_geos, rct_cnf_save_fs_lst = get_geos(
        reacs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    prd_geos, prd_cnf_save_fs_lst = get_geos(
        prods, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
    prd_zmas = list(map(automol.geom.zmatrix, prd_geos))
    for geo in prd_geos:
        xyzs = automol.geom.coordinates(geo)
    if len(rct_zmas) > 2:
        rct_zmas.append(ichzma)
    if len(prd_zmas) > 2:
        prd_zmas.append(ichzma)
    return rct_zmas, prd_zmas, rct_cnf_save_fs_lst, prd_cnf_save_fs_lst


def get_geos(
        spcs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str):
    """get geos for reactants and products using the initial level of theory
    """
    spc_geos = []
    cnf_save_fs_lst = []
    for spc in spcs:
        spc_info = [spc_dct[spc]['ich'], spc_dct[spc]['chg'], spc_dct[spc]['mul']]
        orb_restr = moldr.util.orbital_restriction(spc_info, ini_thy_info)
        ini_thy_level = ini_thy_info[0:3]
        ini_thy_level.append(orb_restr)
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs.leaf.create(spc_info)
        spc_run_path = spc_run_fs.leaf.path(spc_info)
        ini_thy_save_fs = autofile.fs.theory(spc_save_path)
        ini_thy_save_path = ini_thy_save_fs.leaf.path(ini_thy_level[1:4])
        ini_thy_run_fs = autofile.fs.theory(spc_run_path)
        ini_thy_run_path = ini_thy_run_fs.leaf.path(ini_thy_level[1:4])
        cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        cnf_save_fs_lst.append(cnf_save_fs)
        cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
        min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
        if min_cnf_locs:
            geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        else:
            run_fs = autofile.fs.run(ini_thy_run_path)
            run_fs.trunk.create()
            tmp_ini_fs = [None, ini_thy_save_fs]
            tmp_fs = [spc_save_fs, spc_run_fs, ini_thy_save_fs, ini_thy_run_fs,
                      cnf_save_fs, cnf_run_fs, run_fs]
            geo = moldr.geom.reference_geometry(
                spc_dct[spc], ini_thy_level, ini_thy_level, tmp_fs, tmp_ini_fs,
                kickoff_size, kickoff_backward, projrot_script_str,
                overwrite=False)
        spc_geos.append(geo)
    return spc_geos, cnf_save_fs_lst


def ts_class(rct_zmas, prd_zmas, rad_rad, ts_mul, low_mul, high_mul, rct_cnf_save_fs_lst, prd_cnf_save_fs_lst, given_class):
    """ determine type of reaction and related ts info from the reactant and product z-matrices.
    Returns the type, the transition state z-matrix, the name of the coordinate to optimize,
    the grid of values for the initial grid search, the torsion names and symmetries, and
    whether or not to update the guess on successive steps.
    These parameters are set for both the initial and a backup evaluation for if the initial ts
    search fails.
    """

    # Force a trimolecular reaction to behave like a bimolecular.
    lrxn.rxnid.convert_termolec_to_bimolec(rct_zmas, prd_zmas)

    # Set the reaction type an get bond information
    brk_name, frm_bnd_key, brk_bnd_key = lrxn.rxnid.determine_reaction_type(
        rct_zmas, prd_zmas,
        rct_cnf_save_fs_lst, prd_cnf_save_fs_lst, rct_tors_names)

    # determine the grid for the preliminary grid search for all the different reaction types
    dist_coo, = automol.zmatrix.coordinates(ts_zma)[dist_name]
    syms = automol.zmatrix.symbols(ts_zma)
    bnd_len_key = tuple(sorted(map(syms.__getitem__, dist_coo)))

    grid, update_guess = lrxn.grid.build_grid(
        rtype, spin, ts_bnd_len, ts_zma, dist_name, npoints=None)

    if typ:
        ts_class_data = [
            typ, ts_zma, dist_name, brk_name, grid, frm_bnd_key, brk_bnd_key,
            tors_names, update_guess]
    else:
        ts_class_data = []
    if bkp_typ:
        bkp_ts_class_data = [
            bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess]
    else:
        bkp_ts_class_data = []

    return ts_class_data, bkp_ts_class_data


def find_ts(
        spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
        bkp_ts_class_data, ini_thy_info, thy_info, run_prefix, save_prefix,
        rxn_run_path, rxn_save_path, overwrite, attempt=1,
        pst_params=[1.0, 6],
        rad_rad_ts='vtst'):
    """ find the ts geometry
    """
    print('prepping ts scan for:', typ)

    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_info[0:2], saddle=True)

    orb_restr = moldr.util.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[0:3]
    ref_level.append(orb_restr)

    # Set up filesystem stuff
    thy_run_fs = autofile.fs.theory(rxn_run_path)
    thy_run_fs.leaf.create(ref_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(ref_level[1:4])

    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(ref_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(ref_level[1:4])

    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)

    ts_run_fs = autofile.fs.ts(thy_run_path)
    ts_run_fs.trunk.create()
    ts_run_path = ts_run_fs.trunk.path()
    run_fs = autofile.fs.run(ts_run_path)

    ts_save_fs = autofile.fs.ts(thy_save_path)
    ts_save_fs.trunk.create()
    ts_save_path = ts_save_fs.trunk.path()

    cnf_run_fs = autofile.fs.conformer(ts_run_path)
    cnf_save_fs = autofile.fs.conformer(ts_save_path)
    cnf_save_fs.trunk.create()

    # Unpack the dist_info list
    dist_name = dist_info[0]
    update_guess = dist_info[2]
    brk_name = dist_info[3]

    # Check if TS already is found, and determine if it fits original guess
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    
    # Check if there is a saddle point in the filesystem or find it
    if min_cnf_locs and not overwrite:
        final_dtst = lts.check_ts_in_fs()
    else:
        # Run for a barrierless reaction
        fs = [None, None, ts_run_fs, ts_save_fs,
              cnf_run_fs, cnf_save_fs, None, None,
              scn_run_fs, scn_save_fs, run_fs]
        print('running ts scan:')
        rad_rad = ('radical radical' in typ)
        low_spin = ('low spin' in typ)
        if rad_rad and low_spin and 'elimination' not in ts_dct['class']:
            # Run MEP Scan
            lvsetup.setup_multiref_mep()
            
            # Set a special active space for O2 otherwise, handle it below
            rcts = ts_dct['reacs']
            num_act_orb, num_act_elc = lvsetup.active_space(rcts)

            # Using rad_rad_ts model, run PST, VTST, VRC-TST
            if rad_rad_ts.lower() == 'vtst':
                lvtst.run_vtst_scan()
            elif rad_rad_ts.lower() == 'vrctst':
                lvrctst.calc_vrctst_rates()
            else:
                pass
        else:
            # Run scan and optimizations for saddle points
            lscan.run_scan()
            lts.opt_ts()

# additional functions that have moved entireley
# lwells.find_vdw()
# lwells.fake_conf()
# lwells.fake_geo_gen()
