"""
  Build the species string based on the model
"""

import os
import numpy
from scipy.interpolate import interp1d
from qcelemental import periodictable as ptab
import elstruct
import automol
import mess_io
import projrot_io
import autofile

# New libs
from routines.es import conformer
from lib.phydat import phycon
from lib.runner import script


def vib_harm_tors_rigid(spc, spc_info, elec_levels, sym_factor,
                        harm_min_cnf_locs, harm_cnf_save_fs):
    """ Build the species string for a model of
        harmonic vibrational frequencies and rigid torsions
    """
    # Do the freqs obtain for two species for fake and pst
    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs.leaf.file.geometry.read(harm_min_cnf_locs)
        min_ene = harm_cnf_save_fs.leaf.file.energy.read(harm_min_cnf_locs)

        if automol.geom.is_atom(harm_geo):
            # Different for pst and fake
            print('This is an atom')
            mass = ptab.to_mass(harm_geo[0][0])
            spc_str = mess_io.writer.atom(
                mass, elec_levels)
        else:
            hess = harm_cnf_save_fs.leaf.file.hessian.read(harm_min_cnf_locs)
            freqs = elstruct.util.harmonic_frequencies(
                harm_geo, hess, project=False)
            mode_start = 6
            if 'ts_' in spc:
                mode_start = mode_start + 1
                imag_freq = freqs[0]
            if automol.geom.is_linear(harm_geo):
                mode_start = mode_start - 1
            freqs = freqs[mode_start:]

            hind_rot_str = ""

            core = mess_io.writer.core_rigidrotor(harm_geo, sym_factor)
            spc_str = mess_io.writer.molecule(
                core, freqs, elec_levels,
                hind_rot=hind_rot_str,
                )
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies ',
              'for species {}'.format(spc_info[0]))
        raise ValueError
        # spc_str = ''
        # imag_freq = 0.

    return spc_str, imag_freq


def vib_harm_tors_1dhr(harm_min_cnf_locs, harm_cnf_save_fs,
                       tors_min_cnf_locs, tors_cnf_save_fs,
                       tors_save_path, tors_cnf_save_path,
                       spc_dct_i, spc, spc_info,
                       frm_bnd_key, brk_bnd_key,
                       sym_factor, elec_levels,
                       projrot_script_str,
                       saddle=False):
    """ Build the species string for a model: Harm, 1DHR
    """
    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs.leaf.file.geometry.read(
            harm_min_cnf_locs)
        min_ene = harm_cnf_save_fs.leaf.file.energy.read(
            harm_min_cnf_locs)
        if automol.geom.is_atom(harm_geo):
            mass = ptab.to_mass(harm_geo[0][0])
            spc_str = mess_io.writer.atom(
                mass, elec_levels)
            imag_freq = 0.0
        else:
            hess = harm_cnf_save_fs.leaf.file.hessian.read(
                harm_min_cnf_locs)
            freqs = elstruct.util.harmonic_frequencies(
                harm_geo, hess, project=False)
            hind_rot_str = ""
            proj_rotors_str = ""

            if tors_min_cnf_locs is not None:

                # Get geometry for the torsional minimum
                zma = tors_cnf_save_fs.leaf.file.zmatrix.read(
                    tors_min_cnf_locs)
                tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(
                    tors_min_cnf_locs)

                # Set torsional stuff
                tors_names, ts_bnd = get_tors_names(
                    spc_dct_i, zma, tors_cnf_save_fs, saddle=saddle)
                tors_grids, tors_sym_nums = tors_params(
                    spc_dct_i, zma, tors_names, frm_bnd_key, brk_bnd_key)

                # Loop over the torsions
                pot = []
                tors_info = zip(tors_names, tors_grids, tors_sym_nums)
                for tors_name, tors_grid, sym_num in tors_info:

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
                            atm_key, ts_bnd, sym_num)
                    group = list(numpy.add(group, 1))
                    axis = list(numpy.add(axis, 1))
                    if (atm_key+1) != axis[1]:
                        axis.reverse()

                    # Check for dummy transformations
                    remdummy = check_dummy_trans(zma)

                    # Write the MESS and ProjRot strings for the rotor
                    hind_rot_str += mess_io.writer.rotor_hindered(
                        group, axis, sym_num, pot, remdummy=remdummy)
                    proj_rotors_str += projrot_io.writer.rotors(
                        axis, group, remdummy=remdummy)

                    # Divide total sym_factor by rotor sym number
                    sym_factor /= sym_num

                    # Increment index for the loop
                    # idx += 1

                # Calculate ZPVES of the hindered rotors
                tors_zpe = 0.0
                if saddle and tors_names is not None:
                    tors_zpe = calc_tors_freqs_zpe(
                        tors_geo, sym_factor, elec_levels,
                        hind_rot_str, tors_save_path)

                # Run one vers ProjRot to proj freqs for that version
                freqs1, imag_freq1, zpe_harm_no_tors = projrot_freqs_1(
                    tors_geo, hess, pot,
                    proj_rotors_str, projrot_script_str,
                    tors_save_path, saddle=False)

                # Now run the other version of ProjRot
                freqs2, imag_freq2, zpe_harm_no_tors_2, harm_zpe = projrot_freqs_2(
                        tors_save_path, pot, spc)

                # Determine freqs and imag_freqs
                freqs, imag_freq = determine_freqs(
                    freqs1, freqs2, imag_freq1, imag_freq2,
                    zpe_harm_no_tors, zpe_harm_no_tors_2,
                    harm_zpe, tors_zpe)

                # Write the MESS string
                core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str
                    )
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies',
              ' for species {}'.format(spc_info[0]))
        spc_str = ''
        imag_freq = 0.0

    return spc_str, imag_freq


def symmetry_factor(sym_model, spc_dct_i, spc_info, dist_names,
                    saddle, frm_bnd_key, brk_bnd_key, tors_names,
                    tors_cnf_save_fs, tors_min_cnf_locs,
                    sym_cnf_save_fs, sym_min_cnf_locs):
    """ Get the overall factor for a species
    """

    form_coords = []
    if 'sym' in spc_dct_i:
        sym_factor = spc_dct_i['sym']
        print('sym_factor from spc_dct_i:', sym_factor)
    else:
        if sym_model == 'SAMPLING':
            if not sym_min_cnf_locs:
                # Fix the return statement here
                print('ERROR: Reference geometry is missing for symmetry',
                      'for species {}'.format(spc_info[0]))
                return '', 0.
            sym_geo = sym_cnf_save_fs.leaf.file.geometry.read(sym_min_cnf_locs)
            sym_ene = sym_cnf_save_fs.leaf.file.energy.read(sym_min_cnf_locs)
            if dist_names:
                zma = tors_cnf_save_fs.leaf.file.zmatrix.read(
                    tors_min_cnf_locs)
                form_coords = list(
                    automol.zmatrix.bond_idxs(zma, dist_names[0]))
                form_coords.extend(list(dist_names[1]))
            sym_factor = conformer.symmetry_factor(
                sym_geo, sym_ene, sym_cnf_save_fs, saddle,
                frm_bnd_key, brk_bnd_key, form_coords, tors_names)
            print('sym_factor from conformer sampling:', sym_factor)
        elif sym_model == '1DHR':
            print('Warning: the 1DHR based symmetry number',
                  'has not yet been set up')
            sym_factor = 1
        else:
            print('Warning: no symmetry model requested,',
                  'setting symmetry factor to 1.0')
            sym_factor = 1

    return sym_factor


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

    script.util.run_script(pf_script_str, pf_path)

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
                        atm_key, ts_bnd, sym_num):
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


def get_tors_names(spc_dct_i, zma, tors_cnf_save_fs, saddle=False):
    """ get the tors names
    """
    if tors_cnf_save_fs.trunk.file.info.exists():
        inf_obj_s = tors_cnf_save_fs.trunk.file.info.read()
        tors_ranges = inf_obj_s.tors_ranges
        tors_ranges = autofile.info.dict_(tors_ranges)
        tors_names = list(tors_ranges.keys())
    else:
        print('No inf obj to identify torsional angles')

    # Get other info if it is a saddle point
    ts_bnd = None
    if saddle:
        dist_name = spc_dct_i['dist_info'][0]
        tors_names = spc_dct_i['tors_names']
        ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)

    return tors_names, ts_bnd


def tors_params(spc_dct_i, zma, tors_names, frm_bnd_key, brk_bnd_key):
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
    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
        zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key))

    return tors_grids, tors_sym_nums


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


def projrot_freqs_1(tors_geo, hess, pot,
                    proj_rotors_str, projrot_script_str,
                    save_path, saddle=False):
    """ Get frequencies from one version of ProjRot
    """
    coord_proj = 'cartesian'
    grad = ''
    # Write the string for the ProjRot input
    projrot_inp_str = projrot_io.writer.rpht_input(
        tors_geo, grad, hess, rotors_str=proj_rotors_str,
        coord_proj=coord_proj)

    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(save_path)
    bld_save_fs.leaf.create(bld_locs)
    path = bld_save_fs.leaf.path(bld_locs)
    print('Build Path for Partition Functions in species block')
    print(path)
    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_inp_str)

    script.run_script(projrot_script_str, path)

    freqs = []
    zpe_har_no_tors = 0.
    har_zpe = 0.
    if pot:
        rthrproj_freqs, _ = projrot_io.reader.rpht_output(
            path+'/hrproj_freq.dat')
        freqs = rthrproj_freqs
        zpe_har_no_tors = sum(freqs)*phycon.WAVEN2KCAL/2.
    rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
        path+'/RTproj_freq.dat')
    har_zpe = sum(rtproj_freqs)*phycon.WAVEN2KCAL/2.
    if not freqs:
        freqs = rtproj_freqs
    if saddle:
        if imag_freq:
            imag_freq = imag_freq[0]
        else:
            imag_freq = freqs[-1]
            freqs = freqs[:-1]

    return freqs, imag_freq, zpe_har_no_tors


def projrot_freqs_2(save_path, pot, spc):
    """ Get ProjRot frequencies via ProjRot 2
    """

    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(save_path)
    bld_save_fs.leaf.create(bld_locs)
    path = bld_save_fs.leaf.path(bld_locs)

    projrot_script_str2 = (
        "#!/usr/bin/env bash\n"
        "RPHt.exe >& /dev/null")
    script.run_script(projrot_script_str2, path)

    zpe_har_no_tors_2 = 0.0
    freqs_2 = []
    if pot:
        rthrproj_freqs_2, _ = projrot_io.reader.rpht_output(
            path+'/hrproj_freq.dat')
        freqs_2 = rthrproj_freqs_2
        zpe_har_no_tors_2 = sum(freqs_2)*phycon.WAVEN2KCAL/2.
    rtproj_freqs, imag_freq_2 = projrot_io.reader.rpht_output(
        path+'/RTproj_freq.dat')
    har_zpe = sum(rtproj_freqs)*phycon.WAVEN2KCAL/2.
    if not freqs_2:
        freqs_2 = rtproj_freqs
    if 'ts_' in spc:
        if imag_freq_2:
            imag_freq_2 = imag_freq_2[0]
        else:
            imag_freq_2 = freqs_2[-1]
            freqs_2 = freqs_2[:-1]

    return freqs_2, imag_freq_2, har_zpe, zpe_har_no_tors_2


def determine_freqs(freqs1, freqs2, imag_freq1, imag_freq2,
                    zpe_harm_no_tors, zpe_harm_no_tors_2,
                    harm_zpe, tors_zpe):
    """ get the freqs ftom two methods
    """
    harm_tors_zpe = harm_zpe - zpe_harm_no_tors
    harm_tors_zpe_2 = harm_zpe - zpe_harm_no_tors_2
    del_tors_zpe = harm_tors_zpe - tors_zpe
    del_tors_zpe_2 = harm_tors_zpe_2 - tors_zpe
    if del_tors_zpe <= del_tors_zpe_2:
        # zpe = zpe_harm_no_tors + tors_zpe
        freqs = freqs1
        imag_freq = imag_freq1
    else:
        # zpe = zpe_harm_no_tors_2 + tors_zpe
        freqs = freqs2
        imag_freq = imag_freq2

    return freqs, imag_freq


def set_fake_freqs(har_geo_i, har_geo_j,
                   har_min_cnf_locs_i, har_min_cnf_locs_j,
                   har_cnf_save_fs_i, har_cnf_save_fs_j):
    """ Set fake frequencies
    """
    freqs = [30, 50, 70, 100, 200]
    ntrans = 5
    is_atom_i = automol.geom.is_atom(har_geo_i)
    is_linear_i = automol.geom.is_linear(har_geo_i)
    is_atom_j = automol.geom.is_atom(har_geo_j)
    is_linear_j = automol.geom.is_linear(har_geo_i)
    if is_atom_i:
        ntrans = ntrans - 3
    if is_atom_j:
        ntrans = ntrans - 3
    if is_linear_i:
        ntrans = ntrans - 2
    if is_linear_j:
        ntrans = ntrans - 2
    if is_atom_i and is_atom_j:
        ntrans = 0
    freqs = freqs[0:ntrans]
    if not is_atom_i:
        hess_i = har_cnf_save_fs_i.leaf.file.hessian.read(
            har_min_cnf_locs_i)
        freqs_i = elstruct.util.harmonic_frequencies(
            har_geo_i, hess_i, project=False)
        mode_start = 6
        if automol.geom.is_linear(har_geo_i):
            mode_start = mode_start - 1
        freqs += freqs_i[mode_start:]
    if not is_atom_j:
        hess_j = har_cnf_save_fs_j.leaf.file.hessian.read(
            har_min_cnf_locs_j)
        freqs_j = elstruct.util.harmonic_frequencies(
            har_geo_j, hess_j, project=False)
        mode_start = 6
        if automol.geom.is_linear(har_geo_j):
            mode_start = mode_start - 1
        freqs += freqs_j[mode_start:]

    return freqs


def combine_geos_in_fake_well(har_geo_i, har_geo_j):
    """ put two geometries together in a fake well
    """
    max_z_i = max(atom[1][2] for atom in har_geo_i)
    min_z_j = min(atom[1][2] for atom in har_geo_j)
    har_geo = har_geo_i
    har_geo_j = automol.geom.translated(
        har_geo_j, [0., 0., max_z_i + 6. - min_z_j])
    har_geo += har_geo_j

    return har_geo
