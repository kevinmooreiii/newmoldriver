""" calculates certain quantities of interest using MESS+filesytem
"""

import automol
import elstruct
import autofile

# New Libs
from lib.phydat import phycon
from lib.runner import script
from lib.filesystem import minc as fsmin
from lib.filesystem import orb as fsorb
import routines.pf.messf.models as pfmodels
from routines.pf.messf.blocks import set_model_filesys


def get_high_level_energy(
        spc_info, thy_low_level, thy_high_level, save_prefix, saddle=False):
    """ get high level energy at low level optimized geometry
    """
    if saddle:
        spc_save_path = save_prefix
    else:
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_low_level)
    thy_low_level = thy_low_level[1:3]
    thy_low_level.append(orb_restr)

    ll_save_fs = autofile.fs.theory(spc_save_path)
    ll_save_path = ll_save_fs.leaf.path(thy_low_level)

    if saddle:
        ll_save_fs = autofile.fs.ts(ll_save_path)
        ll_save_fs.trunk.create()
        ll_save_path = ll_save_fs.trunk.path()

    cnf_save_fs = autofile.fs.conformer(ll_save_path)
    min_cnf_locs = fsmin.min_energy_conformer_locators(
        cnf_save_fs)
    if not min_cnf_locs:
        print('ERROR: No minimum conformer geometry for ',
              'this species {}'.format(spc_info[0]))
        return 0.0
    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)

    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_high_level)
    thy_high_level = thy_high_level[1:3]
    thy_high_level.append(orb_restr)

    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs.leaf.create(thy_high_level)

    min_ene = sp_save_fs.leaf.file.energy.read(thy_high_level)

    return min_ene


def get_zpe(spc, spc_dct, spc_save_path, pf_levels, spc_model):
    """ return the zpe for a given species according a specified set of
    partition function levels
    """
    spc_zpe, is_atom = get_zero_point_energy(
        spc, spc_dct, pf_levels, spc_model,
        elec_levels=[[0., 1]],
        sym_factor=1.0,
        save_prefix=spc_save_path)
    zpe_str = '{0:<8.2f}\n'.format(spc_zpe)
    if is_atom:
        zero_energy_str = 'End'
    else:
        zero_energy_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
        zero_energy_str += 'End'

    return spc_zpe, zero_energy_str


def get_zero_point_energy(
        spc, spc_dct_i, pf_levels, spc_model,
        elec_levels=((0., 1)), sym_factor=1.0,
        save_prefix='spc_save_path'):
    """ compute the ZPE including torsional and anharmonic corrections
    """

    # Set up input
    projrot_script_str = script.PROJROT
    spc_info = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])

    # Prepare the sets of file systems
    [_, _, harm_levels, _, _, tors_levels] = pf_levels
    print('spc_model')
    print(spc_model)
    tors_model, vib_model, _ = spc_model

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_prefix)

    # Set boolean to account for rad-rad reaction (not supported by vtst)
    rad_rad_ts = False
    saddle = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True
        saddle = True

    # Set the filesystem objects for various species models
    harmfs = set_model_filesys(
        thy_save_fs, spc_info, harm_levels, saddle=saddle)
    [harm_cnf_save_fs, _,
     harm_min_cnf_locs, _] = harmfs
    if tors_levels and not rad_rad_ts:
        torsfs = set_model_filesys(
            thy_save_fs, spc_info, tors_levels[0], saddle=saddle)
        [tors_cnf_save_fs, tors_cnf_save_path,
         tors_min_cnf_locs, tors_save_path] = torsfs

    # if saddle:
    if 'ts_' in spc:
        frm_bnd_key = spc_dct_i['frm_bnd_key']
        brk_bnd_key = spc_dct_i['brk_bnd_key']
    else:
        frm_bnd_key = []
        brk_bnd_key = []

    # Get reference harmonic
    harm_zpe = 0.0
    is_atom = False
    if not harm_min_cnf_locs:
        print('ERROR: No harmonic reference geometry for this species ',
              '{}'.format(spc_info[0]))
        return harm_zpe, is_atom
    harm_geo = harm_cnf_save_fs.leaf.file.geometry.read(harm_min_cnf_locs)
    if automol.geom.is_atom(harm_geo):
        zpe = 0.0
        is_atom = True
    else:
        hess = harm_cnf_save_fs.leaf.file.hessian.read(harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

        saddle = False
        mode_start = 6
        if 'ts_' in spc:
            mode_start = mode_start + 1
            saddle = True
        if automol.geom.is_linear(harm_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]

        harm_zpe = sum(freqs)*phycon.WAVEN2KCAL/2.

        # Determine the ZPVE based on the model
        print(vib_model)
        print(tors_model)
        print(rad_rad_ts)
        if (vib_model == 'harm' and tors_model == 'rigid') or rad_rad_ts:
            print('HARM_RIGID')
            zpe = harm_zpe
        elif vib_model == 'harm' and tors_model == '1dhr':
            print('HARM_1DHR')
            _, _, _, _, zpe, _ = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs, harm_cnf_save_fs,
                tors_min_cnf_locs, tors_cnf_save_fs,
                tors_save_path, tors_cnf_save_path,
                spc_dct_i, spc_info,
                frm_bnd_key, brk_bnd_key,
                sym_factor, elec_levels,
                projrot_script_str,
                saddle=saddle)
        elif vib_model == 'harm' and tors_model == 'mdhr':
            print('HARM and MDHR combination is not yet implemented')
        elif vib_model == 'harm' and tors_model == 'tau':
            print('HARM and TAU combination is not yet implemented')
        elif vib_model == 'vpt2' and tors_model == 'rigid':
            print('VPT2 and RIGID combination is not yet implemented')
        elif vib_model == 'vpt2' and tors_model == '1dhr':
            print('VPT2 and 1DHR combination is not yet implemented')
        elif vib_model == 'vpt2' and tors_model == 'tau':
            print('VPT2 and TAU combination is not yet implemented')

    return zpe, is_atom
