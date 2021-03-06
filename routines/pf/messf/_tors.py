"""
  Functions handling hindered rotor model calculations
"""

import os
import numpy
from scipy.interpolate import interp1d
import automol
import mess_io
import projrot_io
import autofile

# New libs
from lib.phydat import phycon
from lib.runner import script


# MESS strings
def write_1dhr_tors_mess_strings(harm_geo, spc_info, spc_dct_i, ts_bnd, zma,
                                 tors_names, tors_grids, tors_sym_nums,
                                 tors_cnf_save_path, min_ene,
                                 saddle=False, hind_rot_geo=None):
    """ Gather the 1DHR torsional data and gather them into a MESS file
    """
    # Loop over the torsions
    hind_rot_str = ""
    proj_rotors_str = ""
    tors_info = zip(tors_names, tors_grids, tors_sym_nums)
    for tors_name, tors_grid, tors_sym in tors_info:

        # Read the hindered rotor potential
        pot = read_hr_pot(
            spc_info, tors_name, tors_grid,
            tors_cnf_save_path, min_ene)

        # Build potential lst from only successful calculations
        pot = hrpot_spline_fitter(pot)

        # Get the HR groups and axis for the rotor
        group, axis, atm_key = set_groups_ini(
            zma, tors_name, ts_bnd, saddle)
        if saddle:
            group, axis, pot = check_saddle_groups(
                zma, spc_dct_i, group, axis,
                pot, ts_bnd, tors_sym)
        group = list(numpy.add(group, 1))
        axis = list(numpy.add(axis, 1))
        if (atm_key+1) != axis[1]:
            axis.reverse()

        # Check for dummy transformations
        remdummy = check_dummy_trans(zma)

        # Write the MESS and ProjRot strings for the rotor
        hrgeo = harm_geo if hind_rot_geo else None
        hind_rot_str += mess_io.writer.rotor_hindered(
            group, axis, tors_sym, pot, remdummy=remdummy, geom=hrgeo)
        proj_rotors_str += projrot_io.writer.rotors(
            axis, group, remdummy=remdummy)

    return hind_rot_str, proj_rotors_str


def write_mdhr_tors_mess_strings(geom, spc_info, sym_num, spc_dct_i,
                                 ts_bnd, zma,
                                 tors_names, tors_grids, tors_sym_nums,
                                 tors_cnf_save_path, min_ene,
                                 saddle=False, hind_rot_geo=None):
    """ Gather the MDHR torsional data and gather them into a MESS file
    """
    # Loop over the torsions and get the int rot strings and potentials
    rotor_internal_str = ''
    hind_rot_potentials = []
    tors_info = zip(tors_names, tors_grids, tors_sym_nums)
    for tors_name, tors_grid, tors_sym in tors_info:

        # Read the hindered rotor potential (NEED NEW VERSION)
        pot = read_hr_pot(
            spc_info, tors_name, tors_grid,
            tors_cnf_save_path, min_ene)

        # Build potential lst from only successful calculations
        pot = hrpot_spline_fitter(pot)

        # Get the HR groups and axis for the rotor
        group, axis, atm_key = set_groups_ini(
            zma, tors_name, ts_bnd, saddle)
        if saddle:
            group, axis, pot = check_saddle_groups(
                zma, spc_dct_i, group, axis,
                pot, ts_bnd, sym_num)
        group = list(numpy.add(group, 1))
        axis = list(numpy.add(axis, 1))
        if (atm_key+1) != axis[1]:
            axis.reverse()

        # Add potential to master list
        hind_rot_potentials.append(pot)

        # Check for dummy transformations
        remdummy = check_dummy_trans(zma)

        # Write the MESS and ProjRot strings for the rotor
        rotor_internal_str += mess_io.writer.mol_data.rotor_internal(
            group, axis, tors_sym,
            rotor_id='', remdummy=remdummy,
            mass_exp_size=5, pot_exp_size=5,
            hmin=13, hmax=101,
            grid_size=100)

    # Write the MDHR potential file
    mdhr_str = write_mdhr_dat_file(hind_rot_potentials)

    return rotor_internal_str, mdhr_str


def write_mdhr_dat_file(potentials):
    """ Write a file containing the hindered rotor potentials
        Only writes the file for up to 4-dimensinal rotor
    """
    npts1 = len(potentials)
    npts2 = len(potentials[0]) if npts1 > 1 else 0
    npts3 = len(potentials[0][0]) if npts2 > 1 else 0
    npts4 = len(potentials[0][0][0]) if npts3 > 1 else 0

    # Write top line string with number of points in potential
    mdhr_str = '{0:>6d}'.format(npts1)
    if npts2 > 0:
        mdhr_str += '{0:>6d}'.format(npts2)
    if npts3 > 0:
        mdhr_str += '{0:>6d}'.format(npts3)
    if npts4 > 0:
        mdhr_str += '{0:>6d}'.format(npts4)

    # Add the nofreq line (need to know when to put the freqs)
    mdhr_str += '\n nofreq\n\n'

    # Write the strings with the potential values
    if npts2 == npts3 == npts4 == 0:
        for i in npts1:
            mdhr_str += (
                '{0:>6d}{1:>12.f}\n'.format(
                    i+1, potentials[i]))
    elif npts2 > 0 and npts3 == npts4 == 0:
        for i in npts1:
            for j in npts2:
                mdhr_str += (
                    '{0:>6d}{1:>6d}{2:>12.f}\n'.format(
                        i+1, j+1, potentials[i][j])
                )
    elif npts2 > 0 and npts3 > 0 and npts4 == 0:
        for i in npts1:
            for j in npts2:
                for k in npts3:
                    mdhr_str += (
                        '{0:>6d}{1:>6d}{2:>6d}{3:>12.f}\n'.format(
                            i+1, j+1, k+1, potentials[i][j][k])
                    )
    else:
        for i in npts1:
            for j in npts2:
                for k in npts3:
                    for m in npts4:
                        mdhr_str += (
                            '{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>12.f}\n'.format(
                                i+1, j+1, k+1, m+1, potentials[i][j][k][m])
                        )

    return mdhr_str


# Functions to handle setting up torsional defintion and potentials properly
def read_hr_pot(spc_info, tors_name, tors_grid, tors_cnf_save_path, min_ene):
    """ Get the potential for a hindered rotor
    """
    scn_save_fs = autofile.fs.scan(tors_cnf_save_path)
    locs_lst = []
    enes = []
    for grid_val in tors_grid:
        locs_lst.append([[tors_name], [grid_val]])
    for locs in locs_lst:
        if scn_save_fs.leaf.exists(locs):
            enes.append(scn_save_fs.leaf.file.energy.read(locs))
        else:
            enes.append(10.)
            print('ERROR: missing grid value for torsional potential ',
                  'of {}'.format(spc_info[0]))

    enes = numpy.subtract(enes, min_ene)
    pot = list(enes*phycon.EH2KCAL)

    return pot


def hrpot_spline_fitter(pot, thresh=-0.05):
    """ Get a physical hindered rotor potential via a series of spline fits
    """

    # Build a potential list from only successful calculations
    lpot = len(pot)+1
    idx_success = []
    pot_success = []
    pot.append(0.)
    for idx in range(lpot):
        if pot[idx] < 600.:
            idx_success.append(idx)
            pot_success.append(pot[idx])
    idx_success.append(lpot)
    pot_success.append(pot[0])
    pot_spl = interp1d(
        numpy.array(idx_success), numpy.array(pot_success), kind='cubic')
    for idx in range(lpot):
        pot[idx] = pot_spl(idx)

    # Do second spline fit of only positive values if any negative values found
    if any(val < thresh for val in pot):
        print('Found pot vals below',
              ' {0} kcal. Refit w/ positives'.format(thresh))
        print('Potential before spline:', pot)
        x_pos = numpy.array([i for i in range(lpot)
                             if pot[i] >= thresh])
        y_pos = numpy.array([pot[i] for i in range(lpot)
                             if pot[i] >= thresh])
        pos_pot_spl = interp1d(x_pos, y_pos, kind='cubic')
        pot_pos_fit = []
        for idx in range(lpot):
            pot_pos_fit.append(pos_pot_spl(idx))

        print('Potential after spline:', pot_pos_fit)
        # Perform second check to see if negative potentials have been fixed
        if any(val < thresh for val in pot_pos_fit):
            print('Found values below {0} kcal again.'.format(thresh),
                  ' Trying linear interp of positive vals')
            neg_idxs = [i for i in range(lpot) if pot_pos_fit[i] < thresh]
            clean_pot = []
            for i in range(lpot):
                if i in neg_idxs:
                    # Find the indices for positive vals around negative value
                    idx_0 = i - 1
                    while idx_0 in neg_idxs:
                        idx_0 = idx_0 - 1
                    for j in range(i, lpot):
                        if pot_pos_fit[j] >= thresh:
                            idx_1 = j
                            break
                    # Get a new value for this point on the potential by
                    # doing a linear interp of positives
                    interp_val = (
                        pot_pos_fit[idx_0] * (1.0-((i-idx_0)/(idx_1-idx_0))) +
                        pot_pos_fit[idx_1] * ((i-idx_0)/(idx_1-idx_0))
                    )
                    clean_pot.append(interp_val)
                else:
                    clean_pot.append(pot[i])
            final_potential = clean_pot.copy()

        else:
            final_potential = pot_pos_fit.copy()

    else:
        final_potential = pot.copy()

    print('Final potential in spline fitter:', final_potential)
    final_potential = final_potential[:-1]

    return final_potential


def set_groups_ini(zma, tors_name, ts_bnd, saddle):
    """ Set the initial set of groups
    """
    gra = automol.zmatrix.graph(zma, remove_stereo=True)
    coo_dct = automol.zmatrix.coordinates(zma, multi=False)
    axis = coo_dct[tors_name][1:3]
    atm_key = axis[1]
    if ts_bnd:
        for atm in axis:
            if atm in ts_bnd:
                atm_key = atm
                break
    group = list(
        automol.graph.branch_atom_keys(
            gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) - set(axis))
    if not group:
        for atm in axis:
            if atm != atm_key:
                atm_key = atm
        group = list(
            automol.graph.branch_atom_keys(
                gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) - set(axis))

    return group, axis, atm_key


def check_saddle_groups(zma, spc_dct_i, group, axis,
                        pot, ts_bnd, sym_num):
    """ Assess that hindered rotor groups and axes
    """
    n_atm = automol.zmatrix.count(zma)
    if 'addition' in spc_dct_i['class'] or 'abstraction' in spc_dct_i['class']:
        group2 = []
        ts_bnd1 = min(ts_bnd)
        ts_bnd2 = max(ts_bnd)
        for idx in range(ts_bnd2, n_atm):
            group2.append(idx)
        if ts_bnd1 in group:
            for atm in group2:
                if atm not in group:
                    group.append(atm)

    # Check to see if symmetry of XH3 rotor was missed
    if sym_num == 1:
        group2 = []
        for idx in range(n_atm):
            if idx not in group and idx not in axis:
                group2.append(idx)
        all_hyd = True
        symbols = automol.zmatrix.symbols(zma)
        hyd_count = 0
        for idx in group2:
            if symbols[idx] != 'H' and symbols[idx] != 'X':
                all_hyd = False
                break
            else:
                if symbols[idx] == 'H':
                    hyd_count += 1
        if all_hyd and hyd_count == 3:
            sym_num = 3
            lpot = int(len(pot)/3)
            potp = []
            potp[0:lpot] = pot[0:lpot]
            pot = potp

    return group, axis, pot


def check_dummy_trans(zma):
    """ check trans
    """
    atom_symbols = automol.zmatrix.symbols(zma)
    dummy_idx = []
    for atm_idx, atm in enumerate(atom_symbols):
        if atm == 'X':
            dummy_idx.append(atm_idx)
    remdummy = numpy.zeros(len(zma[0]))
    for dummy in dummy_idx:
        for idx, _ in enumerate(remdummy):
            if dummy < idx:
                remdummy[idx] += 1

    return remdummy


# Calculating certain quantities on the torsions
def calc_tors_freqs_zpe(tors_geo, sym_factor, elec_levels,
                        hind_rot_str, tors_save_path):
    """ Calculate the frequencies and ZPVES of the hindered rotors
        create a messpf input and run messpf to get tors_freqs and tors_zpes
    """
    dummy_freqs = [1000.]
    dummy_zpe = 0.0
    core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
    spc_str = mess_io.writer.molecule(
        core, dummy_freqs, elec_levels,
        hind_rot=hind_rot_str,
        )
    temp_step = 100.
    ntemps = 5
    zpe_str = '{0:<8.2f}\n'.format(dummy_zpe)
    zpe_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
    zpe_str += 'End\n'
    global_pf_str = mess_io.writer.global_pf(
        [], temp_step, ntemps, rel_temp_inc=0.001,
        atom_dist_min=0.6)
    spc_head_str = 'Species ' + ' Tmp'
    pf_inp_str = '\n'.join(
        [global_pf_str, spc_head_str,
         spc_str, zpe_str])

    bld_locs = ['PF', 0]
    bld_save_fs = autofile.fs.build(tors_save_path)
    bld_save_fs.leaf.create(bld_locs)
    pf_path = bld_save_fs.leaf.path(bld_locs)

    # run messpf
    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)
    pf_script_str = ("#!/usr/bin/env bash\n"
                     "export OMP_NUM_THREADS=10\n"
                     "messpf pf.inp pf.out >> stdout.log &> stderr.log")

    script.run_script(pf_script_str, pf_path)

    with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
        output_string = mess_file.read()

    # Read the freqs and zpes
    # tors_freqs = mess_io.reader.tors.freqs(output_string)
    tors_zpes = mess_io.reader.tors.zpves(output_string)

    # Calculate the torsional zpe
    tors_zpe = sum(tors_zpes) if tors_zpes else 0.0
    # tors_zpe_cor = 0.0
    # tors_zpe = 0.0
    # for (tors_freq, tors_1dhr_zpe) in zip(tors_freqs, tors_zpes):
    #     tors_zpe_cor += tors_1dhr_zpe - tors_freq*phycon.WAVEN2KCAL/2
    #     tors_zpe += tors_1dhr_zpe

    return tors_zpe


# Handle strucutral information about torsions
def get_tors_names(spc_dct_i, tors_cnf_save_fs, saddle=False):
    """ get the tors names
    """
    if saddle:
        tors_names = spc_dct_i['tors_names']
    else:
        if tors_cnf_save_fs.trunk.file.info.exists():
            inf_obj_s = tors_cnf_save_fs.trunk.file.info.read()
            tors_ranges = inf_obj_s.tors_ranges
            tors_ranges = autofile.info.dict_(tors_ranges)
            tors_names = list(tors_ranges.keys())
        else:
            print('No inf obj to identify torsional angles')

    return tors_names


def get_tors_grids(spc_dct_i, zma, tors_names, frm_bnd_key, brk_bnd_key):
    """ get tors parameters
    """
    # Prepare stuff
    if 'hind_inc' in spc_dct_i:
        scan_increment = spc_dct_i['hind_inc']
    else:
        scan_increment = 30. * phycon.DEG2RAD
    val_dct = automol.zmatrix.values(zma)

    # Set up torsional things
    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
        zma, tors_names, scan_increment,
        frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)
    tors_grids = [
        numpy.linspace(*linspace) + val_dct[name]
        for name, linspace in zip(tors_names, tors_linspaces)]

    return tors_grids


def get_tors_sym_nums(spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                      frm_bnd_key, brk_bnd_key, saddle=False):
    """ get tors parameters
    """
    zma = tors_cnf_save_fs.leaf.file.zmatrix.read(
        tors_min_cnf_locs)
    tors_names = get_tors_names(
        spc_dct_i, tors_cnf_save_fs, saddle=saddle)
    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
        zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key))

    return tors_sym_nums
