""" set of functions to make pf.py less terrible
"""

def species_block():
    """ Prepare the Species Block for MESS files
    """

    # Unpack the models and levels
    har_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model
    
    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_prefix)
    
    # Set the filesystem objects for various species models 
    harm = set_model_fs(thy_save_fs, spc_info, harm_level, ts=('ts_' in spc))
    harm_cnf_save_fs, harm_cnf_save_path, harm_min_cnf_locs, harm_save_path = harmfs
    if sym_level:
        sym = set_model_fs(thy_save_fs, spc_info, sym_level, ts=('ts_' in spc))
        sym_cnf_save_fs, sym_cnf_save_path, sym_min_cnf_locs, sym_save_path = symfs
    if tors_level and not rad_rad_ts:
        tors = set_model_fs(thy_save_fs, spc_info, tors_level, ts=('ts_' in spc))
        tors_cnf_save_fs, tors_cnf_save_path, tors_min_cnf_locs, tors_save_path = torsfs
    if vpt2_level:
        vpt2 = set_model_fs(thy_save_fs, spc_info, vpt2_level, ts=('ts_' in spc))
        vpt2_cnf_save_fs, vpt2_cnf_save_path, vpt2_min_cnf_locs, vpt2_save_path = vpt2fs

    # Set additional info for a saddle point
    if 'ts_' in spc:
        saddle = True
        tors_names = spc_dct_i['tors_names']
        if 'migration' in spc_dct_i['class'] or 'elimination' in spc_dct_i['class']:
            dist_names.append(spc_dct_i['dist_info'][0])
            dist_names.append(spc_dct_i['dist_info'][3])

    # Set TS information
    frm_bnd_key, brk_bnd_key = get_bnd_keys(spc_dct, saddle)

    # Set boolean to account for a radical radical reaction (not supported by vtst)
    rad_rad_ts = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True

    # Initialize electronic energy levels
    elec_levels = init_elec_levels(spc_dct, spc_info)

    # Determine the species symmetry factor using the given model
    sym_factor = get_symmetry_factor(spc_dct_i, sym_min_cnf_locs, tors_min_cnf_locs, zma, dist_names, ...):
    
    # Initialize imag frequency in case it is not set in the following if-block
    imag_freq = 0.

    # Build the MESS string for the species based on the model chosen
    spr_str = ''
    if (vib_model == 'HARM' and tors_model == 'RIGID') or rad_rad_ts:
        model_harm_vib_tors_rigid()
    elif vib_model == 'HARM' and tors_model == '1DHR':
        model_harm_vib_tors_1dhr()
    elif vib_model == 'HARM' and tors_model == 'MDHR':
        print('HARM and MDHR combination is not yet implemented')
    elif vib_model == 'HARM' and tors_model == 'TAU':
        # model_harm_vib_tors_tau()
        print('HARM and TAU combination is not yet implemented')
    elif vib_model == 'VPT2' and tors_model == 'RIGID':
        # model_vib_vpt2_tors_rigid()
        print('VPT2 and RIGID combination is not yet properly implemented')
    elif vib_model == 'VPT2' and tors_model == '1DHR':
        print('VPT2 and 1DHR combination is not yet implemented')
    elif vib_model == 'VPT2' and tors_model == 'TAU':
        print('VPT2 and TAU combination is not yet implemented')
    else:
        raise NotImplementedError

    return spc_str, imag_freq


def pst_block():
    """ prepare a Phase Space Theory species block
    """

    # Unpack the models and levels
    har_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # Set theory filesystem used throughout
    thy_save_fs_i = autofile.fs.theory(save_prefix_i)
    thy_save_fs_j = autofile.fs.theory(save_prefix_i)

    # Set the filesystem objects for the two species
    harm_i = set_model_fs(thy_save_fs_i, spc_info_i, harm_level, ts=False)
    harm_cnf_save_fs_i, harm_cnf_save_path_i, harm_min_cnf_locs_i, harm_save_path_i = harmfs_i
    harm_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, harm_level, ts=False)
    harm_cnf_save_fs_j, harm_cnf_save_path_j, harm_min_cnf_locs_j, harm_save_path_j = harmfs_j

    if sym_level:
        symfs_i = set_model_fs(thy_save_fs_i, spc_info_i, sym_level, ts=False)
        symfs_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, sym_level, ts=False)
        sym_cnf_save_fs_i, sym_cnf_save_path_i, sym_min_cnf_locs_i, sym_save_path_i = symfs_i
        sym_cnf_save_fs_j, sym_cnf_save_path_j, sym_min_cnf_locs_j, sym_save_path_j = symfs_j

    if tors_level:
        torsfs_i = set_model_fs(thy_save_fs_i, spc_info_i, tors_level, ts=False)
        torsfs_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, tors_level, ts=False)
        tors_cnf_save_fs_i, tors_cnf_save_path_i, tors_min_cnf_locs_i, tors_save_path_i = torsfs_i
        tors_cnf_save_fs_j, tors_cnf_save_path_j, tors_min_cnf_locs_j, tors_save_path_j = torsfs_j

    if vpt2_level:
        vpt2fs_i = set_model_fs(thy_save_fs_i, spc_info_i, vpt2_level, ts=False)
        vpt2fs_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, vpt2_level, ts=False)
        vpt2_cnf_save_fs_i, vpt2_cnf_save_path_i, vpt2_min_cnf_locs_i, vpt2_save_path_i = vpt2fs_i
        vpt2_cnf_save_fs_j, vpt2_cnf_save_path_j, vpt2_min_cnf_locs_j, vpt2_save_path_j = vpt2fs_j
    
    # Get the combined electronic energy levels
    elec_levels = combine_elec_levels(spc_dct_i, spc_dct_j)

    # Determine the species symmetry factor using the given model
    sym_factor_i = get_symmetry_factor(spc_dct_i, sym_min_cnf_locs, tors_min_cnf_locs, zma, dist_names, ...):
    sym_factor_j = get_symmetry_factor(spc_dct_j, sym_min_cnf_locs, tors_min_cnf_locs, zma, dist_names, ...):

    # Build the MESS string for the species based on the model chosen
    spr_str = ''
    if (vib_model == 'HARM' and tors_model == 'RIGID') or rad_rad_ts:
    elif vib_model == 'HARM' and tors_model == '1DHR':
    else:
        raise NotImplementedError

    return spc_str, imag_freq


def fake_species_block():
    """ prepare a Phase Space Theory species block
    """
    
    # Unpack the models and levels
    har_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # Set theory filesystem used throughout
    thy_save_fs_i = autofile.fs.theory(save_prefix_i)
    thy_save_fs_j = autofile.fs.theory(save_prefix_i)

    # Set the filesystem objects for the two species
    harm_i = set_model_fs(thy_save_fs_i, spc_info_i, harm_level, ts=False)
    harm_cnf_save_fs_i, harm_cnf_save_path_i, harm_min_cnf_locs_i, harm_save_path_i = harmfs_i
    harm_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, harm_level, ts=False)
    harm_cnf_save_fs_j, harm_cnf_save_path_j, harm_min_cnf_locs_j, harm_save_path_j = harmfs_j

    if sym_level:
        symfs_i = set_model_fs(thy_save_fs_i, spc_info_i, sym_level, ts=False)
        symfs_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, sym_level, ts=False)
        sym_cnf_save_fs_i, sym_cnf_save_path_i, sym_min_cnf_locs_i, sym_save_path_i = symfs_i
        sym_cnf_save_fs_j, sym_cnf_save_path_j, sym_min_cnf_locs_j, sym_save_path_j = symfs_j

    if tors_level:
        torsfs_i = set_model_fs(thy_save_fs_i, spc_info_i, tors_level, ts=False)
        torsfs_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, tors_level, ts=False)
        tors_cnf_save_fs_i, tors_cnf_save_path_i, tors_min_cnf_locs_i, tors_save_path_i = torsfs_i
        tors_cnf_save_fs_j, tors_cnf_save_path_j, tors_min_cnf_locs_j, tors_save_path_j = torsfs_j

    if vpt2_level:
        vpt2fs_i = set_model_fs(thy_save_fs_i, spc_info_i, vpt2_level, ts=False)
        vpt2fs_j = set_model_fs(thy_save_fs_j, spc_jnfo_j, vpt2_level, ts=False)
        vpt2_cnf_save_fs_i, vpt2_cnf_save_path_i, vpt2_min_cnf_locs_i, vpt2_save_path_i = vpt2fs_i
        vpt2_cnf_save_fs_j, vpt2_cnf_save_path_j, vpt2_min_cnf_locs_j, vpt2_save_path_j = vpt2fs_j
    
    # Get the combined electronic energy levels
    elec_levels = combine_elec_levels(spc_dct_i, spc_dct_j)

    # Determine the species symmetry factor using the given model
    sym_factor_i = get_symmetry_factor(spc_dct_i, sym_min_cnf_locs, tors_min_cnf_locs, zma, dist_names, ...):
    sym_factor_j = get_symmetry_factor(spc_dct_j, sym_min_cnf_locs, tors_min_cnf_locs, zma, dist_names, ...):
    sym_factor = sym_factor_i * sym_factor_j

    # Build the MESS string for the species based on the model chosen
    spr_str = ''
    if (vib_model == 'HARM' and tors_model == 'RIGID') or rad_rad_ts:
    elif vib_model == 'HARM' and tors_model == '1DHR':
    else:
        raise NotImplementedError

    return None


def get_symmetry_factor():
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
                print('ERROR: Reference geometry is missing for symmetry for species {}'.format(spc_info[0]))
                return '', 0.
            sym_geo = sym_cnf_save_fs.leaf.file.geometry.read(sym_min_cnf_locs)
            sym_ene = sym_cnf_save_fs.leaf.file.energy.read(sym_min_cnf_locs)
            if dist_names:
                zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)
                form_coords = list(automol.zmatrix.bond_idxs(zma, dist_names[0]))
                form_coords.extend(list(dist_names[1]))
            sym_factor = moldr.conformer.symmetry_factor(
                sym_geo, sym_ene, sym_cnf_save_fs, saddle, frm_bnd_key, brk_bnd_key, form_coords, tors_names)
            print('sym_factor from conformer sampling:', sym_factor)
        elif sym_model == '1DHR':
            print('Warning: the 1DHR based symmetry number has not yet been set up')
            sym_factor = 1
        else:
            print('Warning: no symmetry model requested, setting symmetry factor to 1.0')
            sym_factor = 1

    return sym_factor


def model_harm_vib_tors_rigid():
    """ Build the species string for a model of 
        harmonic vibrational frequencies and rigid torsions
    """
    # Do the freqs obtain for two species for fake and pst
    if har_min_cnf_locs is not None:
        har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
        min_ene = har_cnf_save_fs.leaf.file.energy.read(har_min_cnf_locs)

        # generate fake intermolec freqs for fake and pst
        # freqs = [30, 50, 70, 100, 200]
        # ntrans = 5
        # is_atom_i = automol.geom.is_atom(har_geo_i)
        # is_linear_i = automol.geom.is_linear(har_geo_i)
        # is_atom_j = automol.geom.is_atom(har_geo_j)
        # is_linear_j = automol.geom.is_linear(har_geo_i)
        # if is_atom_i:
        #     ntrans = ntrans - 3
        # if is_atom_j:
        #     ntrans = ntrans - 3
        # if is_linear_i:
        #     ntrans = ntrans - 2
        # if is_linear_j:
        #     ntrans = ntrans - 2
        # if is_atom_i and is_atom_j:
        #     ntrans = 0
        # freqs = freqs[0:ntrans]

        if automol.geom.is_atom(har_geo):
            # Different for pst and fake
            print('This is an atom')
            mass = ptab.to_mass(har_geo[0][0])
            spc_str = mess_io.writer.atom(
                mass, elec_levels)
        else:
            hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
            freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)
            mode_start = 6
            if 'ts_' in spc:
                mode_start = mode_start + 1
                imag_freq = freqs[0]
            if automol.geom.is_linear(har_geo):
                mode_start = mode_start - 1
            freqs = freqs[mode_start:]

            hind_rot_str = ""

            core = mess_io.writer.core_rigidrotor(har_geo, sym_factor)
            spc_str = mess_io.writer.molecule(
                core, freqs, elec_levels,
                hind_rot=hind_rot_str,
                )

            # Difference when doing the PST block
            # form_i = automol.geom.formula(har_geo_i)
            # form_j = automol.geom.formula(har_geo_j)
            # form = automol.formula.join(form_i, form_j)
            # stoich = ''
            # for key, val in form.items():
            #     stoich += key + str(val)
            # core = mess_io.writer.core_phasespace(
            #     har_geo_i, har_geo_j, sym_factor, stoich,
            #     pot_prefactor=pst_params[0], pot_power_exp=pst_params[1])
            # spc_str = mess_io.writer.molecule(
            #     core, freqs, elec_levels,
            #     hind_rot=hind_rot_str,
            #     )

    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies for species {}'.format(spc_info[0]))
        spc_str = ''
        imag_freq = 0.
    
    return spc_str, imag_freq


def model_harm_vib_tors_1dhr():
    """ Build the species string for a model of 
        harmonic vibrational frequencies and rigid torsions
    """
    if har_min_cnf_locs is not None:
        har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
        min_ene = har_cnf_save_fs.leaf.file.energy.read(har_min_cnf_locs)
        if automol.geom.is_atom(har_geo):
            mass = ptab.to_mass(har_geo[0][0])
            spc_str = mess_io.writer.atom(
                mass, elec_levels)
        else:
            hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
            freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)
            hind_rot_str = ""
            proj_rotors_str = ""

            if tors_min_cnf_locs is not None:
                zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)
                tors_names, tors_grids, tors_sym_nums = tors_params()
                idx = 0
                for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):

                    # Read the HR potential
                    enes, pot = read_hr_pot()

                    # Build a potential list from only successful calculations
                    pot = _hrpot_spline_fitter(pot)

                    # Get the HR groups using the graphs
                    group = set_groups_ini()
                    
                    # check to see if fragment group was neglected
                    if saddle:
                        pot = check_saddle_groups()
                            
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
                    idx += 1

                # Create a messpf input file and run messpf to get tors_freqs and tors_zpes
                if saddle and tors_names is not None:
                    tors_zpe = calc_tors_freqs_zpe()

                # run one version of ProjRot to get the projected frequencies for that version
                freqs = projrot_freqs_1()

                # now run the other version of ProjRot
                freqs = projrot_freqs_2()

                har_tors_zpe = har_zpe - zpe_har_no_tors
                har_tors_zpe_2 = har_zpe - zpe_har_no_tors_2
                del_tors_zpe = har_tors_zpe - tors_zpe
                del_tors_zpe_2 = har_tors_zpe_2 - tors_zpe
                if del_tors_zpe <= del_tors_zpe_2:
                    zpe = zpe_har_no_tors + tors_zpe
                else:
                    zpe = zpe_har_no_tors_2 + tors_zpe
                    freqs = freqs_2
                    imag_freq = imag_freq_2
                core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str
                    )

    return spc_str


def model_vib_harm_tors_tau():
    """ Write string for model
    """
    moldr.driver.tau_pf_write(
        name=name,
        save_prefix=thy_save_path,
        run_grad=run_grad_pf,
        run_hess=run_hess_pf,
    )

    return spc_str, imag_freq


def model_vib_vpt2_tors_rigid():
    """ Write string for model
    """
    if anh_min_cnf_locs is not None:
        anh_geo = anh_cnf_save_fs.leaf.file.geometry.read(anh_min_cnf_locs)
        min_ene = anh_cnf_save_fs.leaf.file.energy.read(anh_min_cnf_locs)
        if automol.geom.is_atom(anh_geo):
            mass = ptab.to_mass(anh_geo[0][0])
            spc_str = mess_io.writer.atom(
                mass, elec_levels)
        else:
            hess = anh_cnf_save_fs.leaf.file.hessian.read(anh_min_cnf_locs)
            freqs = elstruct.util.harmonic_frequencies(anh_geo, hess, project=True)
            mode_start = 6
            if 'ts_' in spc:
                mode_start = mode_start + 1
                imag_freq = freqs[0]
            if automol.geom.is_linear(anh_geo):
                mode_start = mode_start - 1
            freqs = freqs[mode_start:]
            zpe = sum(freqs)*phycon.WAVEN2KCAL/2.
            hind_rot_str = ""

            core = mess_io.writer.core_rigidrotor(anh_geo, sym_factor)
            spc_str = mess_io.writer.molecule(
                core, freqs, elec_levels,
                hind_rot=hind_rot_str,
                )
    else:
        print('ERROR: Reference geometry is missing for anharmonic analysis for species {}'.format(spc_info[0]))
        spc_str = ''
        imag_freq = 0.

    return spc_str, imag_freq


def projrot_freqs_1():
    """ Get frequencies from one version of ProjRot
    """
    coord_proj = 'cartesian'
    grad = ''
    projrot_inp_str = projrot_io.writer.rpht_input(
        tors_geo, grad, hess, rotors_str=proj_rotors_str,
        coord_proj=coord_proj)

    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(tors_save_path)
    bld_save_fs.leaf.create(bld_locs)
    path = bld_save_fs.leaf.path(bld_locs)
    print('Build Path for Partition Functions in species block')
    print(path)
    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_inp_str)

    moldr.util.run_script(projrot_script_str, path)

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
    if 'ts_' in spc:
        if imag_freq:
            imag_freq = imag_freq[0]
        else:
            imag_freq = freqs[-1]
            freqs = freqs[:-1]

    return freqs


def projrot_freqs_2():
    """ Get ProjRot frequencies via ProjRot 2
    """
    projrot_script_str2 = ("#!/usr/bin/env bash\n"
    "RPHt.exe >& /dev/null")
    moldr.util.run_script(projrot_script_str2, path)
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

    return freqs_2


def calc_tors_freqs_zpe():
    """ Calculate the frequencies and ZPVES of the hindered rotors
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

    moldr.util.run_script(pf_script_str, pf_path)

    with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
        output_string = mess_file.read()

    # Read the freqs and zpes
    tors_freqs = mess_io.reader.tors.freqs(output_string)
    tors_zpes = mess_io.reader.tors.zpves(output_string)

    # Calculate the torsional zpe
    tors_zpe = sum(tors_zpes) if tors_zpes else 0.0

    return tors_freqs, tors_zpes


def set_groups_ini():
    """ Set the initial set of groups
    """
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

    return group


def check_saddle_groups():
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

    return 


def fake_well_geometry(har_geo_i, har_geo_j):
    """ Put two isolated geoms in a fake well 
    """

    max_z_i = max(atom[1][2] for atom in har_geo_i)
    min_z_j = min(atom[1][2] for atom in har_geo_j)
    har_geo = har_geo_i
    har_geo_j = automol.geom.translated(har_geo_j, [0., 0., max_z_i + 6. - min_z_j])
    har_geo += har_geo_j

    return har_geo


def fake_well_frequencies(har_geo_i, har_geo_j):
    """ Frequencies for a fake well 
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
        hess_i = har_cnf_save_fs_i.leaf.file.hessian.read(har_min_cnf_locs_i)
        freqs_i = elstruct.util.harmonic_frequencies(har_geo_i, hess_i, project=False)
        mode_start = 6
        if automol.geom.is_linear(har_geo_i):
            mode_start = mode_start - 1
        freqs += freqs_i[mode_start:]
    if not is_atom_j:
        hess_j = har_cnf_save_fs_j.leaf.file.hessian.read(har_min_cnf_locs_j)
        freqs_j = elstruct.util.harmonic_frequencies(har_geo_j, hess_j, project=False)
        mode_start = 6
        if automol.geom.is_linear(har_geo_j):
            mode_start = mode_start - 1
        freqs += freqs_j[mode_start:]

    return freqs


def init_elec_levels(spc_dct, spc_info):
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


def set_model_fs(thy_save_fs, spc_info, level, ts=False):
    """ Gets filesystem objects for torsional calculations
    """
    # Set the level for the model
    levelp = level[0:3]
    levelp.append(moldr.util.orbital_restriction(spc_info, level))
    
    # Get the save fileystem path
    save_path = thy_save_fs.leaf.path(levelp[1:4])
    if ts:
        save_fs = autofile.fs.ts(save_path)
        save_fs.trunk.create()
        save_path = save_fs.trunk.path()
    
    # Get the fs object and the locs
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)

    # Get the save path for the conformers
    if min_cnf_locs:
        cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
    else:
        cnf_save_path = ''

    return conf_save_fs, conf_save_path, min_cnf_locs, save_path


def tors_params():
    """ get tors parameters
    """
    if tors_cnf_save_fs.trunk.file.info.exists():
        inf_obj_s = tors_cnf_save_fs.trunk.file.info.read()
        tors_ranges = inf_obj_s.tors_ranges
        tors_ranges = autofile.info.dict_(tors_ranges)
        tors_names = list(tors_ranges.keys())
    else:
        print('No inf obj to identify torsional angles')
        tors_names = []
    zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)

    tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(tors_min_cnf_locs)
    gra = automol.zmatrix.graph(zma, remove_stereo=True)
    coo_dct = automol.zmatrix.coordinates(zma, multi=False)

    # prepare axis, group, and projection info
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
        zma, tors_names, scan_increment,
        frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)
    tors_grids = [
        numpy.linspace(*linspace) + val_dct[name]
        for name, linspace in zip(tors_names, tors_linspaces)]
    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
        zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key))

    return tors_names, tors_grids, tors_sym_nums


def read_hr_pot():
    """ Get the potential for a hindered rotor 
    """
    locs_lst = []
    enes = []
    for grid_val in tors_grid:
        locs_lst.append([[tors_name], [grid_val]])
    for locs in locs_lst:
        if scn_save_fs.leaf.exists(locs):
            enes.append(scn_save_fs.leaf.file.energy.read(locs))
        else:
            enes.append(10.)
            print('ERROR: missing grid value for torsional potential of {}'.format(spc_info[0]))

    enes = numpy.subtract(enes, min_ene)
    pot = list(enes*phycon.EH2KCAL)

    return enes, pot
