"""
Functions to run electronic structure tasks
"""


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
        print('running task {} with abcd of {}'.format('conformer sampling', ' '.join(
            [str(x) for x in params['nsamp_par'][1:]])))
    else:
        print('running task {} for {:g} points'.format(
            'conformer sampling', params['nsamp_par'][5]))
    params['thy_save_fs'] = fs[3]
    params['cnf_run_fs'] = fs[4]
    params['cnf_save_fs'] = fs[5]
    moldr.conformer.conformer_sampling(**params, **opt_kwargs)


def run_hr_scan(fs, params, opt_kwargs):
    """ run a scan over the specified torsional coordinates
    """
    print('running task {}'.format('hr'))
    params['cnf_run_fs'] = fs[4]
    params['cnf_save_fs'] = fs[5]
    moldr.scan.hindered_rotor_scans(**params, **opt_kwargs)


def run_tau_sampling(fs, params, opt_kwargs):
    """ energies, gradients, and hessians,
    for set of arbitrarily sampled torsional coordinates
    with all other coordinates optimized
    """

    params['tau_run_fs'] = fs[6]
    params['tau_save_fs'] = fs[7]
    params['thy_save_fs'] = fs[3]
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


def ts_geometry_generation(tsk, spcdic, es_dct, thy_level, fs, spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
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
