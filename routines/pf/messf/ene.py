""" calculates certain quantities of interest using MESS+filesytem
"""
import os
import numpy
import projrot_io
import automol
import elstruct
import autofile
import mess_io

# New Libs
from lib.phydat import phycon
from lib.submission import substr
from lib.runner.script import run_script
from routines import util
from routines.pf.messf.models import hrpot_spline_fitter


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

    orb_restr = util.orbital_restriction(
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
    min_cnf_locs = util.min_energy_conformer_locators(
        cnf_save_fs)
    if not min_cnf_locs:
        print('ERROR: No minimum conformer geometry for ',
              'this species {}'.format(spc_info[0]))
        return 0.0
    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)

    orb_restr = util.orbital_restriction(
        spc_info, thy_high_level)
    thy_high_level = thy_high_level[1:3]
    thy_high_level.append(orb_restr)

    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs.leaf.create(thy_high_level)

    min_ene = sp_save_fs.leaf.file.energy.read(thy_high_level)

    return min_ene


def get_zero_point_energy(
        spc, spc_dct_i, pf_levels, spc_model, pf_script_str,
        elec_levels=((0., 1)), sym_factor=1.0,
        save_prefix='spc_save_path'):
    """ compute the ZPE including torsional and anharmonic corrections
    """

    projrot_script_str = substr.PROJROT

    spc_info = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    # prepare the sets of file systems
    har_level, tors_level, vpt2_level, _ = pf_levels
    tors_model, vib_model, _ = spc_model

    thy_save_fs = autofile.fs.theory(save_prefix)

    orb_restr = util.orbital_restriction(
        spc_info, har_level)
    har_levelp = har_level[0:3]
    har_levelp.append(orb_restr)

    har_save_path = thy_save_fs.leaf.path(har_levelp[1:4])
    print('har save path test')
    print(har_save_path)
    saddle = False
    if 'ts_' in spc:
        har_save_fs = autofile.fs.ts(har_save_path)
        har_save_fs.trunk.create()
        har_save_path = har_save_fs.trunk.path()
        saddle = True

    har_cnf_save_fs = autofile.fs.conformer(har_save_path)
    har_min_cnf_locs = util.min_energy_conformer_locators(har_cnf_save_fs)

    # Set boolean to account for rad-rad reaction (not supported by vtst)
    rad_rad_ts = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True

    if tors_level and not rad_rad_ts:
        orb_restr = util.orbital_restriction(
            spc_info, tors_level)
        tors_levelp = tors_level[0:3]
        tors_levelp.append(orb_restr)

        tors_save_path = thy_save_fs.leaf.path(tors_levelp[1:4])
        if 'ts_' in spc:
            tors_save_fs = autofile.fs.ts(tors_save_path)
            tors_save_fs.trunk.create()
            tors_save_path = tors_save_fs.trunk.path()

        tors_cnf_save_fs = autofile.fs.conformer(tors_save_path)
        tors_min_cnf_locs = util.min_energy_conformer_locators(
            tors_cnf_save_fs)
        if tors_min_cnf_locs:
            tors_cnf_save_path = tors_cnf_save_fs.leaf.path(tors_min_cnf_locs)

    if vpt2_level:
        orb_restr = util.orbital_restriction(
            spc_info, vpt2_level)
        vpt2_levelp = vpt2_level[0:3]
        vpt2_levelp.append(orb_restr)

        anh_save_path = thy_save_fs.leaf.path(vpt2_levelp[1:4])
        if 'ts_' in spc:
            anh_save_fs = autofile.fs.ts(anh_save_path)
            anh_save_fs.trunk.create()
            anh_save_path = anh_save_fs.trunk.path()

        anh_cnf_save_fs = autofile.fs.conformer(anh_save_path)
        # anh_min_cnf_locs = util.min_energy_conformer_locators(
        # anh_cnf_save_fs)

    if saddle:
        frm_bnd_key = spc_dct_i['frm_bnd_key']
        brk_bnd_key = spc_dct_i['brk_bnd_key']
    else:
        frm_bnd_key = []
        brk_bnd_key = []
    har_zpe = 0.0
    is_atom = False
    # get reference harmonic
    if not har_min_cnf_locs:
        print('ERROR: No harmonic reference geometry for this species ',
              '{}'.format(spc_info[0]))
        return har_zpe, is_atom
    har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
    if automol.geom.is_atom(har_geo):
        har_zpe = 0.0
        is_atom = True
    else:
        hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            har_geo, hess, project=False)

        mode_start = 6
        if 'ts_' in spc:
            mode_start = mode_start + 1
        if automol.geom.is_linear(har_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]

        har_zpe = sum(freqs)*phycon.WAVEN2KCAL/2.

    print('vib_model in zpe:', vib_model)
    print('tors_model in zpe:', tors_model)
    if (vib_model == 'HARM' and tors_model == 'RIGID') or rad_rad_ts:
        ret = har_zpe

    elif vib_model == 'HARM' and tors_model == '1DHR':

        zpe = har_zpe
        hind_rot_str = ""
        proj_rotors_str = ""
        tors_names = []
        if tors_min_cnf_locs is not None:
            if tors_cnf_save_fs.trunk.file.info.exists():
                inf_obj_s = har_cnf_save_fs.trunk.file.info.read()
                tors_ranges = inf_obj_s.tors_ranges
                tors_ranges = autofile.info.dict_(tors_ranges)
                tors_names = list(tors_ranges.keys())
            else:
                print('No inf obj to identify torsional angles')

            min_ene = tors_cnf_save_fs.leaf.file.energy.read(tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(tors_min_cnf_locs)
            zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)
            gra = automol.zmatrix.graph(zma, remove_stereo=True)
            tors_zpe_cor = 0.0
            coo_dct = automol.zmatrix.coordinates(zma, multi=False)
            # prepare axis, group, info
            scn_save_fs = autofile.fs.scan(tors_cnf_save_path)
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                tors_names = spc_dct_i['tors_names']
                ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)
                ts_bnd = frozenset(ts_bnd)
            pot = []
            if 'hind_inc' in spc_dct_i:
                scan_increment = spc_dct_i['hind_inc']
            else:
                scan_increment = 30. * phycon.DEG2RAD
            val_dct = automol.zmatrix.values(zma)
            tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                zma, tors_names, scan_increment, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)
            tors_grids = [numpy.linspace(*linspace) + val_dct[name]
                          for name, linspace in zip(tors_names, tors_linspaces)]
            tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
                zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key))
            # print('tors_names:', tors_names)
            if tors_names:
                for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                    locs_lst = []
                    enes = []
                    for grid_val in tors_grid:
                        locs_lst.append([[tors_name], [grid_val]])
                    for locs in locs_lst:
                        if scn_save_fs.leaf.exists(locs):
                            enes.append(scn_save_fs.leaf.file.energy.read(locs))
                        else:
                            enes.append(10.)
                            print('ERROR: missing grid value for torsional potential of {}'.
                                  format(spc_info[0]))
                    enes = numpy.subtract(enes, min_ene)
                    pot = list(enes*phycon.EH2KCAL)

                    # Build a potential list from only successful calculations
                    pot = hrpot_spline_fitter(pot)

                    axis = coo_dct[tors_name][1:3]
                    atm_key = axis[1]
                    if ts_bnd:
                        for atm in axis:
                            if atm in ts_bnd:
                                atm_key = atm
                                break
                    group = list(
                        automol.graph.branch_atom_keys(gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) -
                        set(axis))
                    if not group:
                        for atm in axis:
                            if atm != atm_key:
                                atm_key = atm
                        group = list(
                            automol.graph.branch_atom_keys(gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) -
                            set(axis))
                    if saddle:
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
                        # check to see if symmetry of XH3 rotor was missed
                        if sym_num == 1:
                            group2 = []
                            for idx in range(n_atm):
                                if idx not in group and idx not in axis:
                                    group2.append(idx)
                            all_H = True
                            symbols = automol.zmatrix.symbols(zma)
                            H_count = 0
                            for idx in group2:
                                if symbols[idx] != 'H' and symbols[idx] != 'X':
                                    all_H = False
                                    break
                                else:
                                    if symbols[idx] == 'H':
                                        H_count += 1
                            if all_H and H_count == 3:
                                sym_num = 3
                                lpot = int(len(pot)/3)
                                potp = []
                                potp[0:lpot] = pot[0:lpot]
                                pot = potp
                    group = list(numpy.add(group, 1))
                    axis = list(numpy.add(axis, 1))
                    if (atm_key+1) != axis[1]:
                        axis.reverse()

                    #check for dummy transformations
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

                    hind_rot_str += mess_io.writer.rotor_hindered(
                        group, axis, sym_num, pot, remdummy=remdummy)

                    proj_rotors_str += projrot_io.writer.rotors(
                        axis, group, remdummy=remdummy)
                    sym_factor /= sym_num

                # Write the string for the ProjRot input
                coord_proj = 'cartesian'
                grad = ''
                projrot_inp_str = projrot_io.writer.rpht_input(
                    tors_geo, grad, hess, rotors_str=proj_rotors_str,
                    coord_proj=coord_proj)

                bld_locs = ['PROJROT', 0]
                bld_save_fs = autofile.fs.build(tors_save_path)
                bld_save_fs.leaf.create(bld_locs)
                path = bld_save_fs.leaf.path(bld_locs)
                print('Build Path for Partition Functions')
                print(path)
                proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                with open(proj_file_path, 'w') as proj_file:
                    proj_file.write(projrot_inp_str)

                run_script(projrot_script_str, path)

                zpe_har_no_tors = har_zpe
                if pot:
                    rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                        path+'/hrproj_freq.dat')
                    freqs = rthrproj_freqs
                    zpe_har_no_tors = sum(freqs)*phycon.WAVEN2KCAL/2.

                # now try again with the other projrot parameters
                projrot_script_str2 = ("#!/usr/bin/env bash\n"
                "RPHt.exe >& /dev/null")
                run_script(projrot_script_str2, path)
                zpe_har_no_tors_2 = har_zpe
                freqs_2 = []
                if pot:
                    rthrproj_freqs_2, _ = projrot_io.reader.rpht_output(
                        path+'/hrproj_freq.dat')
                    freqs_2 = rthrproj_freqs_2
                    zpe_har_no_tors_2 = sum(freqs_2)*phycon.WAVEN2KCAL/2.

                dummy_freqs = [1000.]
                dummy_zpe = 0.0
                core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, dummy_freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )

                # create a messpf input file
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
                run_script(pf_script_str, pf_path)

                with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
                    output_string = mess_file.read()

                # Read the freqs and zpes
                tors_freqs = mess_io.reader.tors.freqs(output_string)
                tors_zpes = mess_io.reader.tors.zpves(output_string)
                tors_zpe_cor = 0.0
                tors_zpe = 0.0
                for (tors_freq, tors_1dhr_zpe) in zip(tors_freqs, tors_zpes):
                    tors_zpe_cor += tors_1dhr_zpe - tors_freq*phycon.WAVEN2KCAL/2
                    tors_zpe += tors_1dhr_zpe

                har_tors_zpe = har_zpe - zpe_har_no_tors
                har_tors_zpe_2 = har_zpe - zpe_har_no_tors_2
                del_tors_zpe = har_tors_zpe - tors_zpe
                del_tors_zpe_2 = har_tors_zpe_2 - tors_zpe
                if del_tors_zpe <= del_tors_zpe_2:
                    zpe = zpe_har_no_tors + tors_zpe
                else:
                    zpe = zpe_har_no_tors_2 + tors_zpe
                if abs(del_tors_zpe) > 0.2 and abs(del_tors_zpe_2) > 0.2:
                    print('Warning: There is a difference of {0:.2f} and {1:.2f} kcal/mol '.format(
                        del_tors_zpe, del_tors_zpe_2),
                        'between the harmonic and hindered torsional zero-point energies')
                # read torsional harmonic zpe and actual zpe
                print('zpe test in get_zpe:',zpe_har_no_tors, zpe_har_no_tors_2, 
                       tors_zpe, har_zpe, zpe)

        ret = zpe

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

    return ret, is_atom
