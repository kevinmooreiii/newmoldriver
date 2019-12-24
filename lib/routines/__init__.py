""" routines modules
"""

from lib.routines.conformer import conformer_sampling
from lib.routines.conformer import run_conformers
from lib.routines.conformer import save_conformers
from lib.routines.conformer import int_sym_num_from_sampling
from lib.routines.conformer import symmetry_factor
from lib.routines.conformer import is_unique_stereo_dist_mat_energy
from lib.routines.conformer import are_torsions_same
from lib.routines.conformer import is_unique_tors_dist_mat_energy
from lib.routines.geom import reference_geometry
from lib.routines.geom import run_initial_geometry_opt
from lib.routines.geom import remove_imag
from lib.routines.geom import run_check_imaginary
from lib.routines.geom import run_kickoff_saddle
from lib.routines.geom import save_initial_geometry
from lib.routines.geom import projrot_frequencies
from lib.routines.scan import hindered_rotor_scans
from lib.routines.scan import run_scan
from lib.routines.scan import run_multiref_rscan
from lib.routines.scan import _run_1d_scan
from lib.routines.scan import _run_2d_scan
from lib.routines.scan import save_scan
from lib.routines.scan import infinite_separation_energy
from lib.routines.sp import run_energy
from lib.routines.sp import run_gradient
from lib.routines.sp import run_hessian
from lib.routines.sp import run_vpt2
from lib.routines.tau import tau_sampling
from lib.routines.tau import run_tau
from lib.routines.tau import save_tau
from lib.routines.tau import tau_pf_write
from lib.routines.ts import reference_geometry
from lib.routines.ts import cas_options_1
from lib.routines.ts import cas_options_2
from lib.routines.ts import multiref_wavefunction_guess
from lib.routines.ts import _run_irc
from lib.routines.util import run_qchem_par
from lib.routines.util import set_molpro_options_mat
from lib.routines.util import orbital_restriction
from lib.routines.util import geometry_dictionary
from lib.routines.util import min_energy_conformer_locators
from lib.routines.util import min_dist_conformer_zma
from lib.routines.util import min_dist_conformer_zma_geo
from lib.routines.util import locs_sort
from lib.routines.util import traj_sort
from lib.routines.util import nsamp_init
from lib.routines.util import reaction_energy
from lib.routines.util import reagent_energies
from lib.routines.util import ts_mul_from_reaction_muls
from lib.routines.util import run_script

__all__ = [
    'conformer_sampling',
    'run_conformers',
    'save_conformers',
    'int_sym_num_from_sampling',
    'symmetry_factor',
    'is_unique_stereo_dist_mat_energy',
    'are_torsions_same',
    'is_unique_tors_dist_mat_energy',
    'reference_geometry',
    'run_initial_geometry_opt',
    'remove_imag',
    'run_check_imaginary',
    'run_kickoff_saddle',
    'save_initial_geometry',
    'projrot_frequencies',
    'hindered_rotor_scans',
    'run_scan',
    'run_multiref_rscan',
    '_run_1d_scan',
    '_run_2d_scan',
    'save_scan',
    'infinite_separation_energy',
    'run_energy',
    'run_gradient',
    'run_hessian',
    'run_vpt2',
    'tau_sampling',
    'run_tau',
    'save_tau',
    'tau_pf_write',
    'reference_geometry',
    'cas_options_1',
    'cas_options_2',
    'multiref_wavefunction_guess',
    '_run_irc',
    'run_qchem_par',
    'set_molpro_options_mat',
    'orbital_restriction',
    'geometry_dictionary',
    'min_energy_conformer_locators',
    'min_dist_conformer_zma',
    'min_dist_conformer_zma_geo',
    'locs_sort',
    'traj_sort',
    'nsamp_init',
    'reaction_energy',
    'reagent_energies',
    'ts_mul_from_reaction_muls',
    'run_script'
]

import numpy
import automol
import elstruct
import autofile

# New libs
from lib.phydat import phycon
from lib.reaction import grid as rxngrid
from lib.reaction import rxnid
from lib.runner import driver
from lib import routines


# Dictionary of Electronic Structure Calculations
ES_TSKS = {
    # Conformers
    'conf_energy': 'run_energy',
    'conf_grad': 'run_grad',
    'conf_hess': 'run_hess',
    'conf_vpt2': 'run_vpt2',
    'conf_reopt': 'run_reopt',
    'conf_samp': 'run_conf_samp',
    # Tau Sampling
    'tau_energy': 'run_energy',
    'tau_grad': 'run_grad',
    'tau_hess': 'run_hess',
    'tau_vpt2': 'run_vpt2',
    'tau_reopt': 'run_reopt',
    'tau_samp': 'run_tau_samp',
    # Hindered Rotors
    'hr_energy': 'run_energy',
    'hr_grad': 'run_grad',
    'hr_hess': 'run_hess',
    'hr_reopt': 'run_reopt',
    'hr_scan': 'run_hr_scan',
    # MEP
    'mep_energy': 'run_energy',
    'mep_grad': 'run_grad',
    'mep_hess': 'run_hess',
    'mep_reopt': 'run_reopt'
}


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


def geometry_generation(tsk, spcdic, es_dct, thy_level, fs,
    spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    geo = moldr.geom.reference_geometry(
        spc_dct[spc], thy_level, ini_thy_level, fs, ini_fs,
        kickoff_size=KICKOFF_SIZE,
        kickoff_backward=KICKOFF_BACKWARD,
        projrot_script_str=substr.PROJROT,
        overwrite=overwrite)
    if geo:
        print('Task in geometry_generation:', tsk)
        _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(
            *thy_level[0:2])
        params = {'spc_info': spc_info,
                  'thy_level': thy_level,
                  'script_str': opt_script_str,
                  'overwrite': overwrite}

        if tsk in ['conf_samp', 'tau_samp']:
            params['nsamp_par'] = es_dct['mc_nsamp']
        elif tsk in ['hr_scan']:
            if 'hind_inc' in spcdic:
                params['scan_increment'] = spcdic['hind_inc']
            else:
                params['scan_increment'] = 30. * phycon.DEG2RAD

        if tsk in ES_TSKS:
            eval(ES_TSKS[tsk])(fs, params, opt_kwargs)


def ts_geometry_generation(tsk, spcdic, es_dct, thy_level, fs,
                           spc_info, overwrite):
    """ run an electronic structure task
    for generating a list of conformer or tau sampling geometries
    """
    geo = moldr.ts.reference_geometry(
        spc_dct[spc], thy_level, ini_thy_level, fs, ini_fs,
        spc_dct[spc]['dist_info'], overwrite)
    if geo:
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(
        *thy_level[0:2], saddle=True)
    print('tsk test in ts_geometry_generation:', tsk)
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': opt_script_str,
              'saddle':  True,
              'tors_names': spcdic['tors_names'],
              'overwrite': overwrite}

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
        print('key test in ts_geom gen:',
              params['frm_bnd_key'],
              params['brk_bnd_key'])

    if tsk in ES_TSKS:
        eval(ES_TSKS[tsk])(fs, params, opt_kwargs)


def geometry_analysis(tsk, thy_level, ini_fs, selection,
                      spc_info, overwrite):
    """ run the specified electronic structure task
    for a set of geometries
    """

    print('Task in geometry_analysis:', tsk)
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
    else:
        locs_lst = selection

    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(*thy_level[0:2])
    params = {'spc_info': spc_info,
              'thy_level': thy_level,
              'script_str': sp_script_str,
              'overwrite': overwrite}

    # cycle over the locations
    if tsk in ES_TSKS:
        task_call = eval(ES_TSKS[tsk])
        for locs in locs_lst:
            if locs:
                params['geo_run_fs'] = run_dir
                params['geo_save_fs'] = save_dir
                params['locs'] = locs
                task_call(params, kwargs)
            else:
                print('No initial geometry available for {} on {}'.format(
                    spc_info[0], '/'.join(thy_level[1:3])))


def ts_geometry_analysis(tsk, thy_level, ini_fs, selection,
                         spc_info, spc_dic, overwrite):
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
        print('key test in ts_geom anal:',
              params['frm_bnd_key'],
              params['brk_bnd_key'])
    else:
        return
    if isinstance(selection, str):
        if selection == 'all':
            locs_lst = save_dir.leaf.existing()
        elif selection == 'min':
            locs_lst = [moldr.util.min_energy_conformer_locators(save_dir)]
    else:
        locs_lst = selection
    sp_script_str, _, kwargs, _ = moldr.util.run_qchem_par(
        *thy_level[0:2], saddle=True)
    params['script_str'] = sp_script_str

    # cycle over the locations
    if tsk in ES_TSKS:
        task_call = eval(ES_TSKS[tsk])
        for locs in locs_lst:
            if locs:
                params['geo_run_fs'] = run_dir
                params['geo_save_fs'] = save_dir
                params['locs'] = locs
                task_call(params, kwargs)
            else:
                print('No initial geometry available for {} on {}'.format(
                    spc_info[0], '/'.join(thy_level[1:3])))



