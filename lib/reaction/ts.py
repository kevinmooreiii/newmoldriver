"""
Find a TS from the grid as well as associated vdW wells
"""

import automol
import elstruct
import moldr


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
        ts_zma = automol.zmatrix.set_value(ts_zma, {'babs1': babs1})
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
                ts_zma = automol.zmatrix.set_value(
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
                bkp_ts_zma = automol.zmatrix.set_value(bkp_ts_zma, {'babs1': babs1})
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
