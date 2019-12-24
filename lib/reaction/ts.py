"""
Find a TS from the grid as well as associated vdW wells
"""

import automol
import elstruct
import moldr


def ts_class(rct_zmas, prd_zmas, rad_rad, ts_mul, low_mul, high_mul, rct_cnf_save_fs_lst, prd_cnf_save_fs_lst, given_class):
    """ determine type of reaction and related ts info from the reactant and product z-matrices.
    Returns the type, the transition state z-matrix, the name of the coordinate to optimize,
    the grid of values for the initial grid search, the torsion names and symmetries, and
    whether or not to update the guess on successive steps.
    These parameters are set for both the initial and a backup evaluation for if the initial ts
    search fails.
    """

    # Convert termolecular reactions to bimolecular reactions
    rct_tors_names = []
    if len(rct_zmas) > 2 or len(prd_zmas) > 2:
        rct_zmas, prd_zmas, rct_tors_names = rxnid.conv_termol_to_bimol(
            rct_zmas, prd_zmas)

    # Determine the reaction types
    ret = rxnid.determine_reaction_type(
        rct_zmas, prd_zmas,
        ts_mul, high_mul, low_mul,
        rct_cnf_save_fs_lst, prd_cnf_save_fs_lst,
        rct_tors_names,
        given_class, rad_rad)
    [typ, bkp_typ,
     ts_zma, bkp_ts_zma,
     tors_names, bkp_tors_names,
     dist_name, bkp_dist_name, brk_name,
     frm_bnd_key, brk_bnd_key] = ret

    # Determine grid for preliminary search for all different reaction types
    dist_coo, = automol.zmatrix.coordinates(ts_zma)[dist_name]
    syms = automol.zmatrix.symbols(ts_zma)
    ts_bnd_len = tuple(sorted(map(syms.__getitem__, dist_coo)))
    grid, update_guess, bkp_grid, bkp_update_guess = rxngrid.build_grid(
        typ, bkp_typ, ts_bnd_len, ts_zma, dist_name, npoints=None)

    # Build class data lists to return from the function
    if typ:
        ts_class_data = [
            typ, ts_zma, dist_name, brk_name,
            grid, frm_bnd_key, brk_bnd_key,
            tors_names, update_guess]
    else:
        ts_class_data = []
    if bkp_typ:
        bkp_ts_class_data = [
            bkp_typ, bkp_ts_zma, bkp_dist_name,
            bkp_grid, bkp_tors_names, bkp_update_guess]
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

    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(
        *thy_info[0:2], saddle=True)

    orb_restr = moldr.util.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[0:3]
    ref_level.append(orb_restr)

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

    dist_name = dist_info[0]
    print('dist_name')
    print(dist_name)
    update_guess = dist_info[2]
    brk_name = dist_info[3]

    # Check if TS already is found, and determine if it fits original guess
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    # check to see if rxn class for already found ts is of expected class
    # do this by comparing names
    if min_cnf_locs and not overwrite:
        cnf_path = cnf_save_fs.trunk.path()
        print('Found TS at {}'.format(cnf_path))
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
        chk_bkp = False
        if automol.zmatrix.names(zma) == automol.zmatrix.names(ts_zma):
            if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1, True):
                if 'babs1' in automol.zmatrix.names(ts_zma):
                    babs1 = 170. * phycon.DEG2RAD
                    if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
                        babs1 = 85. * phycon.DEG2RAD
                    ts_zma = automol.zmatrix.set_valuess(
                        ts_zma, {'babs1': babs1})
                    ts_dct['original_zma'] = ts_zma
                    if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1):
                        chk_bkp = True
                else:
                    chk_bkp = True
        else:
            chk_bkp = True

        is_bkp = False
        if chk_bkp and bkp_ts_class_data:
            [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names,
             bkp_update_guess] = bkp_ts_class_data
            if automol.zmatrix.names(zma) == automol.zmatrix.names(bkp_ts_zma):
                if automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1, True):
                    is_bkp = True
                elif 'babs1' in automol.zmatrix.names(bkp_ts_zma):
                    babs1 = 170. * phycon.DEG2RAD
                    if automol.zmatrix.values(bkp_ts_zma)['babs1'] == babs1:
                        babs1 = 85. * phycon.DEG2RAD
                    bkp_ts_zma = automol.zmatrix.set_valuess(bkp_ts_zma, {'babs1': babs1})
                    if not automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1):
                        is_bkp = True
        if not chk_bkp:
            print("TS is type {}".format(typ))
        elif is_bkp:
            print('updating reaction class to {}'.format(bkp_typ))
            ts_dct['class'] = bkp_typ
            ts_dct['original_zma'] = bkp_ts_zma
            ts_dct['dist_info'] = bkp_dist_info
            ts_dct['tors_names'] = bkp_tors_names
            print("TS is backup type {}".format(bkp_typ))
        else:
            print("TS may not be original type or backup type")
            print("Some part of the z-matrices have changed")
        print('class test:', ts_dct['class'])
        vals = automol.zmatrix.values(zma)
        final_dist = vals[dist_name]
        dist_info[1] = final_dist
        print('dist_info is being set at end of backup checking', dist_info[1], final_dist)

    # Find TS
    else:
        fs = [None, None, ts_run_fs, ts_save_fs,
              cnf_run_fs, cnf_save_fs, None, None,
              scn_run_fs, scn_save_fs, run_fs]

        print('running ts scan:')
        rad_rad = ('radical radical' in typ)
        low_spin = ('low spin' in typ)
        if rad_rad and low_spin and 'elimination' not in ts_dct['class']:
            # run mep scan
            # multi_info = ['molpro2015', 'caspt2', 'cc-pvtz', 'RR']
            multi_info = ['molpro2015', 'caspt2', 'cc-pvdz', 'RR']

            orb_restr = moldr.util.orbital_restriction(ts_info, multi_info)
            multi_level = multi_info[0:3]
            multi_level.append(orb_restr)

            thy_run_fs = autofile.fs.theory(rxn_run_path)
            thy_run_fs.leaf.create(multi_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(multi_level[1:4])

            thy_save_fs = autofile.fs.theory(rxn_save_path)
            thy_save_fs.leaf.create(multi_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(multi_level[1:4])

            scn_run_fs = autofile.fs.scan(thy_run_path)
            scn_save_fs = autofile.fs.scan(thy_save_path)

            ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
            grid1 = grid[0]
            grid2 = grid[1]
            grid = numpy.append(grid[0], grid[1])
            high_mul = ts_dct['high_mul']
            print('starting multiref scan:', scn_run_fs.trunk.path())

            # Set a special active space for O2 otherwise, handle it below
            rcts = ts_dct['reacs']
            if 'InChI=1S/O2/c1-2' in (spc_dct[rcts[0]]['ich'], spc_dct[rcts[1]]['ich']):
                num_act_orb = 5
                num_act_elc = 7
            else:
                num_act_orb = None
                num_act_elc = None

            # run PST, VTST, VRC-TST based on RAD_RAD_TS model
            if rad_rad_ts.lower() == 'pst':
                pass
            elif rad_rad_ts.lower() == 'vtst':
                geo, zma, final_dist = vtst.run_vtst_scan(
                    ts_zma, ts_formula, ts_info, ts_dct, spc_dct,
                    high_mul, grid1, grid2, dist_name,
                    multi_level, num_act_orb, num_act_elc,
                    multi_info, ini_thy_info, thy_info,
                    run_prefix, save_prefix, scn_run_fs, scn_save_fs,
                    opt_script_str, overwrite, update_guess, **opt_kwargs)
            elif rad_rad_ts.lower() == 'vrctst':
                # vrctst.calc_vrctst_rates()
                pass

        else:
            if 'elimination' in typ:
                grid1, grid2 = grid
                grid_dct = {dist_name: grid1, brk_name: grid2}
            else:
                grid_dct = {dist_name: grid}
            print('grid_dct')
            print(grid_dct)
            moldr.scan.run_scan(
                zma=ts_zma,
                spc_info=ts_info,
                thy_level=ref_level,
                grid_dct=grid_dct,
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                script_str=opt_script_str,
                saddle=False,
                overwrite=overwrite,
                update_guess=update_guess,
                reverse_sweep=False,
                fix_failures=False,
                **opt_kwargs,
                )
            if 'elimination' in typ:
                moldr.scan.save_scan(
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    coo_names=[dist_name, brk_name],
                    )
            else:
                moldr.scan.save_scan(
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    coo_names=[dist_name],
                    )

            # Find the structure at the maximum on the grid opt scan
            if 'elimination' in typ:
                max_zma, max_ene = rxngrid.find_max_2D(
                    grid1, grid2, dist_name, brk_name, scn_save_fs)
            else:
                max_zma, max_ene = rxngrid.find_max_1D(
                    typ, grid, ts_zma, dist_name, scn_save_fs)
            print('geometry for maximum along scan:', max_zma)
            print('energy for maximum along scan:', max_ene)

            # Optimize the saddle point using the max structure
            print('optimizing ts')
            driver.run_job(
                job='optimization',
                script_str=opt_script_str,
                run_fs=run_fs,
                geom=max_zma,
                spc_info=ts_info,
                thy_level=ref_level,
                saddle=True,
                overwrite=overwrite,
                **opt_kwargs,
                )

            opt_ret = moldr.driver.read_job(
                job='optimization',
                run_fs=run_fs,
            )
            if opt_ret is not None:
                inf_obj, _, out_str = opt_ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

                print(" - Saving...")
                print(" - Save path: {}".format(ts_save_path))

                ts_save_fs.trunk.file.energy.write(ene)
                ts_save_fs.trunk.file.geometry.write(geo)
                ts_save_fs.trunk.file.zmatrix.write(zma)

                vals = automol.zmatrix.values(zma)
                final_dist = vals[dist_name]
                dist_info[1] = final_dist
                print('Test final distance for reactant coordinate', final_dist)
                run_single_conformer(ts_info, ref_level, fs, overwrite, saddle=True, dist_info=dist_info)

            elif 'addition' in typ and bkp_ts_class_data and attempt > 2:
                bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess = bkp_ts_class_data
                print('TS find failed. Attempting to find with new reaction class: {}'.format(bkp_typ))
                bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
                ts_dct['class'] = bkp_typ
                ts_dct['original_zma'] = bkp_ts_zma
                ts_dct['dist_info'] = bkp_dist_info
                ts_dct['tors_names'] = bkp_tors_names
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
                    bkp_grid, None, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            elif ('addition ' in typ or 'abstraction' in typ) and attempt < 3:
                babs1 = 170. * phycon.DEG2RAD
                if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
                    babs1 = 85. * phycon.DEG2RAD
                print('TS find failed. Attempting to find with new angle of attack: {:.1f}'.format(babs1))
                ts_zma = automol.zmatrix.set_values(ts_zma, {'babs1': babs1})
                ts_dct['original_zma'] = ts_zma
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
                    bkp_ts_class_data, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            elif 'beta scission' in typ and bkp_ts_class_data and attempt < 2:
                [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess] = bkp_ts_class_data
                print('TS find failed. Attempting to find with new reaction class: {}'.format(bkp_typ))
                bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
                ts_dct['class'] = bkp_typ
                ts_dct['original_zma'] = bkp_ts_zma
                ts_dct['dist_info'] = bkp_dist_info
                ts_dct['tors_names'] = bkp_tors_names
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
                    bkp_grid, None, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            else:
                geo = 'failed'
                zma = 'failed'
                final_dist = 0.
    return geo, zma, final_dist


def opt_ts():
    """ Optimize the transition state structure obtained from the grid search
    """

    # Run the transition state optimization
    moldr.driver.run_job(
        job='optimization',
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=max_zma,
        spc_info=ts_info,
        thy_level=ref_level,
        saddle=True,
        overwrite=overwrite,
        **opt_kwargs,
        )

    # Read the contents of the optimization
    opt_ret = moldr.driver.read_job(
        job='optimization',
        run_fs=run_fs,
    )
    if opt_ret is not None:

        # If successful, Read the geom and energy from the optimization
        inf_obj, _, out_str = opt_ret
        prog = inf_obj.prog
        method = inf_obj.method
        ene = elstruct.reader.energy(prog, method, out_str)
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = elstruct.reader.opt_zmatrix(prog, out_str)

        # Save the information into the filesystem
        print(" - Saving...")
        print(" - Save path: {}".format(ts_save_path))

        ts_save_fs.trunk.file.energy.write(ene)
        ts_save_fs.trunk.file.geometry.write(geo)
        ts_save_fs.trunk.file.zmatrix.write(zma)

        # Run single conformer to get intitial conformer in filesystem
        vals = automol.zmatrix.values(zma)
        final_dist = vals[dist_name]
        dist_info[1] = final_dist
        print('Test final distance for reactant coordinate', final_dist)
        run_single_conformer(ts_info, ref_level, fs, overwrite, saddle=True, dist_info=dist_info)

    elif 'addition' in typ and bkp_ts_class_data and attempt > 2:
        # Try to find addition rxn TS in reverse direction
        bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess = bkp_ts_class_data
        print('TS find failed. Attempting to find with new reaction class: {}'.format(bkp_typ))
        bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
        ts_dct['class'] = bkp_typ
        ts_dct['original_zma'] = bkp_ts_zma
        ts_dct['dist_info'] = bkp_dist_info
        ts_dct['tors_names'] = bkp_tors_names
        attempt += 1
        geo, zma, final_dist = find_ts(
            spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
            bkp_grid, None, ini_thy_info, thy_info, run_prefix,
            save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
            attempt=attempt)
    elif ('addition ' in typ or 'abstraction' in typ) and attempt < 3:
        # Try to find addition rxn TS in reverse direction with additional tricks
        babs1 = 170. * phycon.DEG2RAD
        if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
            babs1 = 85. * phycon.DEG2RAD
        print('TS find failed. Attempting to find with new angle of attack: {:.1f}'.format(babs1))
        ts_zma = automol.zmatrix.set_values(ts_zma, {'babs1': babs1})
        ts_dct['original_zma'] = ts_zma
        attempt += 1
        geo, zma, final_dist = find_ts(
            spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
            bkp_ts_class_data, ini_thy_info, thy_info, run_prefix,
            save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
            attempt=attempt)
    elif 'beta scission' in typ and bkp_ts_class_data and attempt < 2:
        # Run reverse for a beta scission reaction
        [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names, bkp_update_guess] = bkp_ts_class_data
        print('TS find failed. Attempting to find with new reaction class: {}'.format(bkp_typ))
        bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
        ts_dct['class'] = bkp_typ
        ts_dct['original_zma'] = bkp_ts_zma
        ts_dct['dist_info'] = bkp_dist_info
        ts_dct['tors_names'] = bkp_tors_names
        attempt += 1
        geo, zma, final_dist = find_ts(
            spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
            bkp_grid, None, ini_thy_info, thy_info, run_prefix,
            save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
            attempt=attempt)
    else:
        # Give up and return failed information
        geo = 'failed'
        zma = 'failed'
        final_dist = 0.
    return geo, zma, final_dist


def check_ts_in_fs():
    """ Check to see if the ts for the class is in the filesystem
        do this by comparing names
    """

    cnf_path = cnf_save_fs.trunk.path()
    print('Found TS at {}'.format(cnf_path))
    geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
    zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
    chk_bkp = False
    if automol.zmatrix.names(zma) == automol.zmatrix.names(ts_zma):
        if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1, True):
            if 'babs1' in automol.zmatrix.names(ts_zma):
                babs1 = 170. * phycon.DEG2RAD
                if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
                    babs1 = 85. * phycon.DEG2RAD
                ts_zma = automol.zmatrix.set_values(
                    ts_zma, {'babs1': babs1})
                ts_dct['original_zma'] = ts_zma
                if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1):
                    chk_bkp = True
            else:
                chk_bkp = True
    else:
        chk_bkp = True

    is_bkp = False
    if chk_bkp and bkp_ts_class_data:
        [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names,
         bkp_update_guess] = bkp_ts_class_data
        if automol.zmatrix.names(zma) == automol.zmatrix.names(bkp_ts_zma):
            if automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1, True):
                is_bkp = True
            elif 'babs1' in automol.zmatrix.names(bkp_ts_zma):
                babs1 = 170. * phycon.DEG2RAD
                if automol.zmatrix.values(bkp_ts_zma)['babs1'] == babs1:
                    babs1 = 85. * phycon.DEG2RAD
                bkp_ts_zma = automol.zmatrix.set_values(bkp_ts_zma, {'babs1': babs1})
                if not automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1):
                    is_bkp = True
    if not chk_bkp:
        print("TS is type {}".format(typ))
    elif is_bkp:
        print('updating reaction class to {}'.format(bkp_typ))
        ts_dct['class'] = bkp_typ
        ts_dct['original_zma'] = bkp_ts_zma
        ts_dct['dist_info'] = bkp_dist_info
        ts_dct['tors_names'] = bkp_tors_names
        print("TS is backup type {}".format(bkp_typ))
    else:
        print("TS may not be original type or backup type")
        print("Some part of the z-matrices have changed")
    print('class test:', ts_dct['class'])
    vals = automol.zmatrix.values(zma)
    final_dist = vals[dist_name]
    dist_info[1] = final_dist
    print('dist_info is being set at end of backup checking', dist_info[1], final_dist)

    return final_dst
