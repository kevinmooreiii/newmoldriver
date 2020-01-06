""" drivers
"""
import os
import numpy
from qcelemental import periodictable as ptab
import projrot_io
import automol
import elstruct
import autofile
import mess_io

# New Libs
from lib.phydat import phycon
from lib.runner.script import run_script
from lib.filesystem import orb as fsorb
from lib.filesystem import minc as fsmin
from routines.pf.messf import models as pfmodels


def species_block(
        spc, spc_dct_i, spc_info, spc_model, pf_levels, projrot_script_str,
        elec_levels=((0., 1)), sym_factor=1., save_prefix='spc_save_path'):
    """ prepare the species input for messpf
    """

    # Unpack the models and levels
    harm_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_prefix)

    # Set boolean to account for rad-rad reaction (not supported by vtst)
    rad_rad_ts = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True

    # Set the filesystem objects for various species models
    harmfs = set_model_filesys(
        thy_save_fs, spc_info, harm_level, saddle=('ts_' in spc))
    harm_cnf_save_fs, harm_cnf_save_path, harm_min_cnf_locs, harm_save_path = harmfs
    if sym_level:
        symfs = set_model_filesys(
            thy_save_fs, spc_info, sym_level, saddle=('ts_' in spc))
        sym_cnf_save_fs, sym_cnf_save_path, sym_min_cnf_locs, sym_save_path = symfs
    if tors_level and not rad_rad_ts:
        torsfs = set_model_filesys(
            thy_save_fs, spc_info, tors_level, saddle=('ts_' in spc))
        tors_cnf_save_fs, tors_cnf_save_path, tors_min_cnf_locs, tors_save_path = torsfs
    if vpt2_level:
        vpt2fs = set_model_filesys(
            thy_save_fs, spc_info, vpt2_level, saddle=('ts_' in spc))
        vpt2_cnf_save_fs, vpt2_cnf_save_path, vpt2_min_cnf_locs, vpt2_save_path = vpt2fs

    # Set additional info for a saddle point
    saddle = False
    dist_names = []
    tors_names = []
    if 'ts_' in spc:
        saddle = True
        tors_names = spc_dct_i['tors_names']
        mig = 'migration' in spc_dct_i['class']
        elm = 'elimination' in spc_dct_i['class']
        if mig or elm:
            dist_names.append(spc_dct_i['dist_info'][0])
            dist_names.append(spc_dct_i['dist_info'][3])

    # Set TS information
    frm_bnd_key, brk_bnd_key = get_bnd_keys(spc_dct_i, saddle)

    # Initialize electronic energy levels
    elec_levels = ini_elec_levels(spc_dct_i, spc_info)

    # Determine the species symmetry factor using the given model
    sym_factor = pfmodels.symmetry_factor(
        sym_model, spc_dct_i, spc_info, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs, tors_min_cnf_locs,
        sym_cnf_save_fs, sym_min_cnf_locs)

    # Build the species string and get the imaginary frequency
    if is_atom(harm_min_cnf_locs, harm_cnf_save_fs):
        # Get the mass needed for the MESS string; set imag to 0
        mass = atom_mass(harm_min_cnf_locs, harm_cnf_save_fs)
        imag = 0.0
        # Write the MESS string for an atom
        spc_str = mess_io.writer.atom(
            mass, elec_levels)
    else:
        if (vib_model == 'HARM' and tors_model == 'RIGID') or rad_rad_ts:
            geo, freqs, imag = pfmodels.vib_harm_tors_rigid(
                spc_info, harm_min_cnf_locs, harm_cnf_save_fs, saddle=saddle)
            hr_str = ""
            symf = sym_factor
        elif vib_model == 'HARM' and tors_model == '1DHR':
            geo, freqs, imag, hr_str, _, symf = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs, harm_cnf_save_fs,
                tors_min_cnf_locs, tors_cnf_save_fs,
                tors_save_path, tors_cnf_save_path,
                spc_dct_i, spc_info,
                frm_bnd_key, brk_bnd_key,
                sym_factor, elec_levels,
                projrot_script_str,
                saddle=saddle)
        elif vib_model == 'HARM' and tors_model == 'MDHR':
            print('HARM and MDHR combination is not yet implemented')
        elif vib_model == 'HARM' and tors_model == 'TAU':
            print('HARM and TAU combination is not yet implemented')
        elif vib_model == 'VPT2' and tors_model == 'RIGID':
            print('VPT2 and RIGID combination is not yet implemented')
        elif vib_model == 'VPT2' and tors_model == '1DHR':
            print('VPT2 and 1DHR combination is not yet implemented')
        elif vib_model == 'VPT2' and tors_model == 'TAU':
            print('VPT2 and TAU combination is not yet implemented')

        # Write the species string for the molecule
        core = mess_io.writer.core_rigidrotor(geo, symf)
        spc_str = mess_io.writer.molecule(
            core, freqs, elec_levels,
            hind_rot=hr_str)

    return spc_str, imag


def vtst_with_no_saddle_block(
        ts_dct, ts_label, reac_label, prod_label,
        spc_ene, rct_zpe, projrot_script_str,
        multi_info, elec_levels=((0., 1)), sym_factor=1.0):
    """ prepare the mess input string for a variational TS that does not have
    a saddle point. Do it by calling the species block for each grid point
    in the scan file system
    """

    ts_info = ['', ts_dct['chg'], ts_dct['mul']]
    orb_restr = fsorb.orbital_restriction(ts_info, multi_info)
    multi_level = multi_info[0:3]
    multi_level.append(orb_restr)

    rxn_save_path = ts_dct['rxn_fs'][3]
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(multi_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(multi_level[1:4])
    scn_save_fs = autofile.fs.scan(thy_save_path)

    # Read the scan save filesystem to get the molecular info
    sym_factor = 1.
    irc_pt_strs = []
    proj_rotors_str = ''
    pot = []

    elec_levels = [[0., ts_dct['mul']]]
    grid = ts_dct['grid']
    grid = numpy.append(grid[0], grid[1])
    dist_name = ts_dct['dist_info'][0]

    # Read the infinite separation energy
    inf_locs = [[dist_name], [1000.]]
    inf_sep_ene = scn_save_fs.leaf.file.energy.read(inf_locs)

    grid[::-1].sort()
    for idx, grid_val in enumerate(grid):
        # Set the filesystem locators for each grid point
        locs = [[dist_name], [grid_val]]

        # Get geometry, energy, and vibrational freqs
        if not scn_save_fs.leaf.file.geometry.exists(locs):
            continue
        else:
            geom = scn_save_fs.leaf.file.geometry.read(locs)
        if not scn_save_fs.leaf.file.energy.exists(locs):
            continue
        else:
            ene = scn_save_fs.leaf.file.energy.read(locs)
        if not scn_save_fs.leaf.file.hessian.exists(locs):
            continue
        else:
            hess = scn_save_fs.leaf.file.hessian.read(locs)
            scn_save_path = scn_save_fs.leaf.path(locs)
            freqs, _, _ = pfmodels.projrot_freqs_1(
                geom, hess, pot,
                proj_rotors_str, projrot_script_str,
                scn_save_path, saddle=True)

        # Calculate standard harmonic ZPVE
        zpe = sum(freqs)*phycon.WAVEN2KCAL/2.

        # Calcuate infinite separation ZPVE
        # Assumes the ZPVE = ZPVE(1st grid pt) as an approximation
        if idx == 0:
            rct_zpe = zpe

        # Calculate the reference energies
        erel = (ene - inf_sep_ene)*phycon.EH2KCAL
        erel_zpe_corr = erel + zpe - rct_zpe
        eref_abs = erel_zpe_corr + spc_ene

        # Iniialize the header of the string
        irc_pt_str = '!----------------------------------------------- \n'
        irc_pt_str += '! IRC Point {0}\n'.format(str(idx+1))

        # Write the MESS string for the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(
            geom, sym_factor, interp_emax='')
        irc_pt_str += mess_io.writer.species.molecule(
            core, freqs, elec_levels,
            hind_rot='', xmat=None, rovib_coups='', rot_dists='')

        # Append the zero energy for the molecule
        irc_pt_str += ('    ZeroEnergy[kcal/mol]      ',
                       '{0:<8.2f}\n'.format(eref_abs))
        if grid_val != grid[-1]:
            irc_pt_str += 'End \n'

        # Append string to list
        irc_pt_strs.append(irc_pt_str)

    # Write the MESS string for the entire variational section
    variational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, irc_pt_strs)

    return variational_str


# def vtst_saddle_block(scn_save_fs, geoms, frequencies, energies):
#     """ prepare the mess input string for a variational TS where there is a
#         saddle point on the MEP.
#         In this case, there is limited torsional information.
#     """
#
#     # Read scn save filesys to get enes, zpves, symnums
#     # Geometries, hessians, torsional potentials for each point on the MEP
#
#     # Determine the the number of points along the irc
#     nirc = 21
#
#     # Loop over all the points of the irc and build MESS strings
#     irc_pt_strings = []
#     for i in range(nirc):
#
#         # Iniialize the header of the string
#         irc_pt_string = '!-----------------------------------------------'
#         irc_pt_string += '! IRC Point {0}\n'.format(str(i+1))
#
#         # Write the molecule section for each irc point
#         core = mess_io.writer.mol_data.core_rigidrotor(
#             geom1, sym_factor, interp_emax='')
#         irc_pt_str += mess_io.writer.species.molecule(
#             core, freqs, elec_levels,
#             hind_rot='', xmat=None, rovib_coups='', rot_dists='')
#
#         # Append the zero point energy for the molecule
#         irc_pt_str += ('    ZeroEnergy[kcal/mol]      ',
#                        '{0:<8.2f}'.format(zero_energy))
#
#         # Append string to list
#         irc_pt_strings.append(irc_pt_string)
#
#     # Write the MESS string for the variational sections
#     variational_str = mess_io.writer.rxnchan.ts_variational(
#         ts_label, reac_label, prod_label, irc_pt_strings)
#
#     return variational_str


def pst_block(
        spc_dct_i, spc_dct_j, spc_model, pf_levels, projrot_script_str,
        spc_save_fs, elec_levels=((0., 1)), sym_factor=1.,
        pst_params=(1.0, 6)):
    """ prepare a Phase Space Theory species block
    """

    # Unpack the models and levels
    harm_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # prepare the four sets of file systems
    spc_info_i = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    spc_info_j = (spc_dct_j['ich'], spc_dct_j['chg'], spc_dct_j['mul'])
    spc_save_fs.leaf.create(spc_info_i)
    spc_save_fs.leaf.create(spc_info_j)
    save_path_i = spc_save_fs.leaf.path(spc_info_i)
    save_path_j = spc_save_fs.leaf.path(spc_info_j)

    # Set theory filesystem used throughout
    thy_save_fs_i = autofile.fs.theory(save_path_i)
    thy_save_fs_j = autofile.fs.theory(save_path_j)

    # Set the filesystem objects for the two species
    harmfs_i = set_model_filesys(
        thy_save_fs_i, spc_info_i, harm_level, saddle=False)
    harm_cnf_save_fs_i, harm_cnf_save_path_i, harm_min_cnf_locs_i, harm_save_path_i = harmfs_i
    harmfs_j = set_model_filesys(thy_save_fs_j, spc_info_j, harm_level, saddle=False)
    harm_cnf_save_fs_j, harm_cnf_save_path_j, harm_min_cnf_locs_j, harm_save_path_j = harmfs_j

    if sym_level:
        symfs_i = set_model_filesys(
            thy_save_fs_i, spc_info_i, sym_level, saddle=False)
        symfs_j = set_model_filesys(
            thy_save_fs_j, spc_info_j, sym_level, saddle=False)
        sym_cnf_save_fs_i, sym_cnf_save_path_i, sym_min_cnf_locs_i, sym_save_path_i = symfs_i
        sym_cnf_save_fs_j, sym_cnf_save_path_j, sym_min_cnf_locs_j, sym_save_path_j = symfs_j

    if tors_level:
        torsfs_i = set_model_filesys(thy_save_fs_i, spc_info_i, tors_level, saddle=False)
        torsfs_j = set_model_filesys(thy_save_fs_j, spc_info_j, tors_level, saddle=False)
        tors_cnf_save_fs_i, tors_cnf_save_path_i, tors_min_cnf_locs_i, tors_save_path_i = torsfs_i
        tors_cnf_save_fs_j, tors_cnf_save_path_j, tors_min_cnf_locs_j, tors_save_path_j = torsfs_j

    if vpt2_level:
        vpt2fs_i = set_model_filesys(thy_save_fs_i, spc_info_i, vpt2_level, saddle=False)
        vpt2fs_j = set_model_filesys(thy_save_fs_j, spc_info_j, vpt2_level, saddle=False)
        vpt2_cnf_save_fs_i, vpt2_cnf_save_path_i, vpt2_min_cnf_locs_i, vpt2_save_path_i = vpt2fs_i
        vpt2_cnf_save_fs_j, vpt2_cnf_save_path_j, vpt2_min_cnf_locs_j, vpt2_save_path_j = vpt2fs_j
    
    # Get the combined electronic energy levels
    elec_levels = combine_elec_levels(spc_dct_i, spc_dct_j)

    # Determine the species symmetry factor using the given model
    saddle = False
    dist_names = []
    tors_names = []
    frm_bnd_key = []
    brk_bnd_key = []
    sym_factor_i = pfmodels.symmetry_factor(
        sym_model, spc_dct_i, spc_info_i, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_i, tors_min_cnf_locs_i,
        sym_cnf_save_fs_i, sym_min_cnf_locs_i)
    sym_factor_j = pfmodels.symmetry_factor(
        sym_model, spc_dct_j, spc_info_j, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_j, tors_min_cnf_locs_j,
        sym_cnf_save_fs_j, sym_min_cnf_locs_j)

    # Get the stoichiometry
    stoich = pfmodels.get_stoich(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    spc_str = ''
    if vib_model == 'HARM' and tors_model == 'RIGID':
        geo_i, freqs_i, _ = pfmodels.vib_harm_tors_rigid(
            spc_info_i, harm_min_cnf_locs_i, harm_cnf_save_fs_i, saddle=False)
        geo_j, freqs_j, _ = pfmodels.vib_harm_tors_rigid(
            spc_info_j, harm_min_cnf_locs_j, harm_cnf_save_fs_j, saddle=False)
        freqs += freqs_i + freqs_j
        hind_rot_str = ""

    if vib_model == 'HARM' and tors_model == '1DHR':
        if is_atom(harm_min_cnf_locs_i, harm_cnf_save_fs_i):
            geo_i = harm_cnf_save_fs_i.leaf.file.geometry.read(
                harm_min_cnf_locs_i)
            freqs_i = []
            hr_str_i = ''
            symf_i = sym_factor_i
        else:
            geo_i, freqs_i, _, hr_str_i, _, symf_i = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_i, harm_cnf_save_fs_i,
                tors_min_cnf_locs_i, tors_cnf_save_fs_i,
                tors_save_path_i, tors_cnf_save_path_i,
                spc_dct_i, spc_info_i,
                frm_bnd_key, brk_bnd_key,
                sym_factor_i, elec_levels,
                projrot_script_str,
                hind_rot_geo=True,
                saddle=False)
        if is_atom(harm_min_cnf_locs_j, harm_cnf_save_fs_j):
            geo_j = harm_cnf_save_fs_j.leaf.file.geometry.read(
                harm_min_cnf_locs_j)
            freqs_j = []
            hr_str_j = ''
            symf_j = sym_factor_j
        else:
            geo_j, freqs_j, _, hr_str_j, _, symf_j = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_j, harm_cnf_save_fs_j,
                tors_min_cnf_locs_j, tors_cnf_save_fs_j,
                tors_save_path_j, tors_cnf_save_path_j,
                spc_dct_j, spc_info_j,
                frm_bnd_key, brk_bnd_key,
                sym_factor_j, elec_levels,
                projrot_script_str,
                hind_rot_geo=True,
                saddle=False)
        freqs = list(freqs_i) + list(freqs_j)
        hind_rot_str = hr_str_i + hr_str_j
        sym_factor = symf_i * symf_j

    # Write the MESS input strings
    core = mess_io.writer.core_phasespace(
        geo_i, geo_j, sym_factor, stoich,
        pot_prefactor=pst_params[0], pot_power_exp=pst_params[1])
    spc_str = mess_io.writer.molecule(
        core, freqs, elec_levels,
        hind_rot=hind_rot_str)

    return spc_str


def fake_species_block(
        spc_dct_i, spc_dct_j, spc_info_i, spc_info_j, spc_model,
        pf_levels, projrot_script_str,
        elec_levels=((0., 1)), sym_factor=1.,
        save_prefix_i='spc_save_path', save_prefix_j='spc_save_path'):
    """ prepare a fake species block corresponding to the
        van der Waals well between two fragments
    """
    harm_level, tors_level, _, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # prepare the four sets of file systems
    orb_restr = fsorb.orbital_restriction(
        spc_info_i, harm_level)
    har_levelp_i = harm_level[0:3]
    har_levelp_i.append(orb_restr)
    orb_restr = fsorb.orbital_restriction(
        spc_info_j, harm_level)
    har_levelp_j = harm_level[0:3]
    har_levelp_j.append(orb_restr)

    # Set theory filesystem used throughout
    thy_save_fs_i = autofile.fs.theory(save_prefix_i)
    thy_save_fs_j = autofile.fs.theory(save_prefix_j)

    # Set the filesystem objects for the two species
    harmfs_i = set_model_filesys(thy_save_fs_i, spc_info_i, harm_level, saddle=False)
    harm_cnf_save_fs_i, harm_cnf_save_path_i, harm_min_cnf_locs_i, harm_save_path_i = harmfs_i
    harmfs_j = set_model_filesys(thy_save_fs_j, spc_info_j, harm_level, saddle=False)
    harm_cnf_save_fs_j, harm_cnf_save_path_j, harm_min_cnf_locs_j, harm_save_path_j = harmfs_j

    if sym_level:
        symfs_i = set_model_filesys(thy_save_fs_i, spc_info_i, sym_level, saddle=False)
        symfs_j = set_model_filesys(thy_save_fs_j, spc_info_j, sym_level, saddle=False)
        sym_cnf_save_fs_i, sym_cnf_save_path_i, sym_min_cnf_locs_i, sym_save_path_i = symfs_i
        sym_cnf_save_fs_j, sym_cnf_save_path_j, sym_min_cnf_locs_j, sym_save_path_j = symfs_j

    if tors_level:
        torsfs_i = set_model_filesys(thy_save_fs_i, spc_info_i, tors_level, saddle=False)
        torsfs_j = set_model_filesys(thy_save_fs_j, spc_info_j, tors_level, saddle=False)
        tors_cnf_save_fs_i, tors_cnf_save_path_i, tors_min_cnf_locs_i, tors_save_path_i = torsfs_i
        tors_cnf_save_fs_j, tors_cnf_save_path_j, tors_min_cnf_locs_j, tors_save_path_j = torsfs_j

    spc_str = ''

    # Get the combined electronic energy levels
    elec_levels = combine_elec_levels(spc_dct_i, spc_dct_j)

    # Determine the species symmetry factor using the given model
    saddle = False
    dist_names = []
    tors_names = []
    frm_bnd_key = []
    brk_bnd_key = []
    sym_factor_i = pfmodels.symmetry_factor(
        sym_model, spc_dct_i, spc_info_i, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_i, tors_min_cnf_locs_i,
        sym_cnf_save_fs_i, sym_min_cnf_locs_i)
    sym_factor_j = pfmodels.symmetry_factor(
        sym_model, spc_dct_j, spc_info_j, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs_j, tors_min_cnf_locs_j,
        sym_cnf_save_fs_j, sym_min_cnf_locs_j)
    sym_factor = sym_factor_i * sym_factor_j

    # Get the freqs
    freqs = pfmodels.set_fake_freqs(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)
    geo = pfmodels.combine_geos_in_fake_well(
        harm_min_cnf_locs_i, harm_min_cnf_locs_j,
        harm_cnf_save_fs_i, harm_cnf_save_fs_j)

    if vib_model == 'HARM' and tors_model == 'RIGID':
        _, freqs_i, _ = pfmodels.vib_harm_tors_rigid(
            spc_info_i, harm_min_cnf_locs_i, harm_cnf_save_fs_i, saddle=saddle)
        _, freqs_j, _ = pfmodels.vib_harm_tors_rigid(
            spc_info_j, harm_min_cnf_locs_j, harm_cnf_save_fs_j, saddle=saddle)
        freqs += freqs_i + freqs_j
        hind_rot_str = ""

    if vib_model == 'HARM' and tors_model == '1DHR':
        if is_atom(harm_min_cnf_locs_i, harm_cnf_save_fs_i):
            freqs_i = []
            hr_str_i = ''
            symf_i = sym_factor_i
        else:
            _, freqs_i, _, hr_str_i, _, symf_i = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_i, harm_cnf_save_fs_i,
                tors_min_cnf_locs_i, tors_cnf_save_fs_i,
                tors_save_path_i, tors_cnf_save_path_i,
                spc_dct_i, spc_info_i,
                frm_bnd_key, brk_bnd_key,
                sym_factor_i, elec_levels,
                projrot_script_str,
                hind_rot_geo=True,
                saddle=False)
        if is_atom(harm_min_cnf_locs_j, harm_cnf_save_fs_j):
            freqs_j = []
            hr_str_j = ''
            symf_j = sym_factor_j
        else:
            _, freqs_j, _, hr_str_j, _, symf_j = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs_j, harm_cnf_save_fs_j,
                tors_min_cnf_locs_j, tors_cnf_save_fs_j,
                tors_save_path_j, tors_cnf_save_path_j,
                spc_dct_j, spc_info_j,
                frm_bnd_key, brk_bnd_key,
                sym_factor_j, elec_levels,
                projrot_script_str,
                hind_rot_geo=True,
                saddle=False)
        freqs += list(freqs_i) + list(freqs_j)
        hind_rot_str = hr_str_i + hr_str_j
        sym_factor = symf_i * symf_j

    core = mess_io.writer.core_rigidrotor(geo, sym_factor)
    spc_str = mess_io.writer.molecule(
        core, freqs, elec_levels,
        hind_rot=hind_rot_str)

    return spc_str


####################
# Helper functions #
####################

def set_model_filesys(thy_save_fs, spc_info, level, saddle=False):
    """ Gets filesystem objects for torsional calculations
    """
    # Set the level for the model
    levelp = level[0:3]
    levelp.append(fsorb.orbital_restriction(spc_info, level))

    # Get the save fileystem path
    save_path = thy_save_fs.leaf.path(levelp[1:4])
    if saddle:
        save_fs = autofile.fs.ts(save_path)
        save_fs.trunk.create()
        save_path = save_fs.trunk.path()

    # Get the fs object and the locs
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)

    # Get the save path for the conformers
    if min_cnf_locs:
        cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
    else:
        cnf_save_path = ''

    return cnf_save_fs, cnf_save_path, min_cnf_locs, save_path


def ini_elec_levels(spc_dct, spc_info):
    """ get initial elec levels
    """
    elec_levels = [[0., spc_info[2]]]
    if 'elec_levs' in spc_dct:
        elec_levels = spc_dct['elec_levs']

    return elec_levels


def combine_elec_levels(spc_dct_i, spc_dct_j):
    """ Put two elec levels together for two species
    """

    if 'elec_levs' in spc_dct_i:
        elec_levels_i = spc_dct_i['elec_levs']
    else:
        elec_levels_i = [[0., spc_dct_i['mul']]]
    if 'elec_levs' in spc_dct_j:
        elec_levels_j = spc_dct_j['elec_levs']
    else:
        elec_levels_j = [[0., spc_dct_j['mul']]]

    # Combine the energy levels
    init_elec_levels = []
    for _, elec_level_i in enumerate(elec_levels_i):
        for _, elec_level_j in enumerate(elec_levels_j):
            init_elec_levels.append(
                [elec_level_i[0]+elec_level_j[0],
                 elec_level_i[1]*elec_level_j[1]])

    # See if any levels repeat and thus need to be added together
    elec_levels = []
    for level in init_elec_levels:
        # Put level in in final list
        if level not in elec_levels:
            elec_levels.append(level)
        # Add the level to the one in the list
        else:
            idx = elec_levels.index(level)
            elec_levels[idx][1] += level[1]

    return elec_levels


def get_bnd_keys(spc_dct, saddle):
    """ get bond broken and formed keys for a transition state
    """
    if saddle:
        frm_bnd_key = spc_dct['frm_bnd_key']
        brk_bnd_key = spc_dct['brk_bnd_key']
    else:
        frm_bnd_key = []
        brk_bnd_key = []

    return frm_bnd_key, brk_bnd_key


def is_atom(har_min_cnf_locs, har_cnf_save_fs):
    """ Check if species is an atom
    """
    if har_min_cnf_locs is not None:
        har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
        print('This is an atom')
    return automol.geom.is_atom(har_geo)


def atom_mass(har_min_cnf_locs, har_cnf_save_fs):
    """ write the atom string
    """
    har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
    return ptab.to_mass(har_geo[0][0])
