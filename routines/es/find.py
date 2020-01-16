"""
Find a TS from the grid as well as associated vdW wells
"""

import numpy
import autofile
import automol
import elstruct

# New
from lib.reaction import grid as rxngrid
from lib.phydat import phycon
from lib.runner import driver
from lib.runner import par as runpar
from lib.filesystem import orb as fsorb
from lib.filesystem import minc as fsmin
from routines.es import variational
from routines.es import conformer
from routines.es import scan


def find_ts(
        spc_dct, ts_dct, ts_zma, typ, dist_info, grid,
        bkp_ts_class_data, ini_thy_info, thy_info, run_prefix, save_prefix,
        overwrite, attempt=1,
        rad_rad_ts='vtst'):
    """ find the ts geometry
    """
    print('prepping ts scan for:', typ)

    # spc_dct[sadpt]['original_zma'],
    # spc_dct[sadpt]['class'],
    # spc_dct[sadpt]['dist_info'],
    # spc_dct[sadpt]['grid'],
    # spc_dct[sadpt]['bkp_data'],

    [_, _,
     rxn_run_path, rxn_save_path] = ts_dct[sadpt]['rxn_fs']
    ts_info = (ts_dct['ich'],
               ts_dct['chg'],
               ts_dct['mul'])

    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
        *thy_info[0:2])

    orb_restr = fsorb.orbital_restriction(ts_info, thy_info)
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
    min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)
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
                    bkp_ts_zma = automol.zmatrix.set_valuess(
                        bkp_ts_zma, {'babs1': babs1})
                    if not automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1):
                        is_bkp = True
        if not chk_bkp:
            print("TS is type {}".format(typ))
        elif is_bkp:
            print('updating reaction class to {}'.format(bkp_typ))
            ts_dct['class'] = bkp_typ
            ts_dct['original_zma'] = bkp_ts_zma
            # ts_dct['dist_info'] = bkp_dist_info
            ts_dct['tors_names'] = bkp_tors_names
            print("TS is backup type {}".format(bkp_typ))
        else:
            print("TS may not be original type or backup type")
            print("Some part of the z-matrices have changed")
        print('class test:', ts_dct['class'])
        vals = automol.zmatrix.values(zma)
        final_dist = vals[dist_name]
        dist_info[1] = final_dist
        print('dist_info is being set at end of backup checking',
              dist_info[1], final_dist)
        # Add an angle check which is added to spc dct for TS
        angle = lts.check_angle(
            ts_dct['original_zma'],
            ts_dct['dist_info'],
            ts_dct['class'])
        ts_dct['dist_info'][1] = final_dist
        ts_dct['dist_info'].append(angle)

    # Find TS
    else:
        filesys = [None, None, ts_run_fs, ts_save_fs,
                   cnf_run_fs, cnf_save_fs, None, None,
                   scn_run_fs, scn_save_fs, run_fs]

        print('running ts scan:')
        rad_rad = ('radical radical' in typ)
        low_spin = ('low spin' in typ)
        if rad_rad and low_spin and 'elimination' not in ts_dct['class']:
            # run mep scan
            # multi_info = ['molpro2015', 'caspt2', 'cc-pvtz', 'RR']
            multi_info = ['molpro2015', 'caspt2', 'cc-pvdz', 'RR']

            orb_restr = fsorb.orbital_restriction(ts_info, multi_info)
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
            rlst = (spc_dct[rcts[0]]['ich'], spc_dct[rcts[1]]['ich'])
            if 'InChI=1S/O2/c1-2' in rlst:
                num_act_orb = 5
                num_act_elc = 7
            else:
                num_act_orb = None
                num_act_elc = None

            # run PST, VTST, VRC-TST based on RAD_RAD_TS model
            if rad_rad_ts.lower() == 'pst':
                pass
            elif rad_rad_ts.lower() == 'vtst':
                geo, zma, final_dist = variational.vtst.run_vtst_scan(
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
            scan.run_scan(
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
                scan.save_scan(
                    scn_run_fs=scn_run_fs,
                    scn_save_fs=scn_save_fs,
                    coo_names=[dist_name, brk_name],
                    )
            else:
                scan.save_scan(
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

            opt_ret = driver.read_job(
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
                print('Test final distance for reactant coord', final_dist)
                conformer.single_conformer(
                    ts_info, ref_level, filesys, overwrite,
                    saddle=True, dist_info=dist_info)

            elif 'addition' in typ and bkp_ts_class_data and attempt > 2:
                [bkp_typ, bkp_ts_zma, bkp_dist_name,
                 bkp_grid, bkp_tors_names,
                 bkp_update_guess] = bkp_ts_class_data
                print('TS find failed. Attempting to find with',
                      'new reaction class: {}'.format(bkp_typ))
                bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
                ts_dct['class'] = bkp_typ
                ts_dct['original_zma'] = bkp_ts_zma
                ts_dct['dist_info'] = bkp_dist_info
                ts_dct['tors_names'] = bkp_tors_names
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, bkp_ts_zma,
                    bkp_typ, bkp_dist_info,
                    bkp_grid, None, ini_thy_info, thy_info,
                     rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            elif ('addition ' in typ or 'abstraction' in typ) and attempt < 3:
                babs1 = 170. * phycon.DEG2RAD
                if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
                    babs1 = 85. * phycon.DEG2RAD
                print('TS find failed. Attempting to find with',
                      'new angle of attack: {:.1f}'.format(babs1))
                ts_zma = automol.zmatrix.set_values(ts_zma, {'babs1': babs1})
                ts_dct['original_zma'] = ts_zma
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
                    bkp_ts_class_data, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            elif 'beta scission' in typ and bkp_ts_class_data and attempt < 2:
                [bkp_typ, bkp_ts_zma, bkp_dist_name,
                 bkp_grid, bkp_tors_names,
                 bkp_update_guess] = bkp_ts_class_data
                print('TS find failed. Attempting to find with',
                      'new reaction class: {}'.format(bkp_typ))
                bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
                ts_dct['class'] = bkp_typ
                ts_dct['original_zma'] = bkp_ts_zma
                ts_dct['dist_info'] = bkp_dist_info
                ts_dct['tors_names'] = bkp_tors_names
                attempt += 1
                geo, zma, final_dist = find_ts(
                    spc_dct, ts_dct, ts_info, bkp_ts_zma,
                    bkp_typ, bkp_dist_info,
                    bkp_grid, None, ini_thy_info, thy_info, run_prefix,
                    save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
                    attempt=attempt)
            else:
                geo = 'failed'
                zma = 'failed'
                final_dist = 0.
    return geo, zma, final_dist
# def opt_ts():
#     """ Optimize the transition state structure obtained from the grid search
#     """
#
#     # Run the transition state optimization
#     driver.run_job(
#         job='optimization',
#         script_str=opt_script_str,
#         run_fs=run_fs,
#         geom=max_zma,
#         spc_info=ts_info,
#         thy_level=ref_level,
#         saddle=True,
#         overwrite=overwrite,
#         **opt_kwargs,
#         )
#
#     # Read the contents of the optimization
#     opt_ret = driver.read_job(
#         job='optimization',
#         run_fs=run_fs,
#     )
#     if opt_ret is not None:
#
#         # If successful, Read the geom and energy from the optimization
#         inf_obj, _, out_str = opt_ret
#         prog = inf_obj.prog
#         method = inf_obj.method
#         ene = elstruct.reader.energy(prog, method, out_str)
#         geo = elstruct.reader.opt_geometry(prog, out_str)
#         zma = elstruct.reader.opt_zmatrix(prog, out_str)
#
#         # Save the information into the filesystem
#         print(" - Saving...")
#         print(" - Save path: {}".format(ts_save_path))
#
#         ts_save_fs.trunk.file.energy.write(ene)
#         ts_save_fs.trunk.file.geometry.write(geo)
#         ts_save_fs.trunk.file.zmatrix.write(zma)
#
#         # Run single conformer to get intitial conformer in filesystem
#         vals = automol.zmatrix.values(zma)
#         final_dist = vals[dist_name]
#         dist_info[1] = final_dist
#         print('Test final distance for reactant coordinate', final_dist)
#         run_single_conformer(ts_info, ref_level, fs, overwrite,
#                              saddle=True, dist_info=dist_info)
#
#     elif 'addition' in typ and bkp_ts_class_data and attempt > 2:
#         # Try to find addition rxn TS in reverse direction
#         bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names,
#         bkp_update_guess = bkp_ts_class_data
#         print('TS find failed. Attempting to find with new',
#               'reaction class: {}'.format(bkp_typ))
#         bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
#         ts_dct['class'] = bkp_typ
#         ts_dct['original_zma'] = bkp_ts_zma
#         ts_dct['dist_info'] = bkp_dist_info
#         ts_dct['tors_names'] = bkp_tors_names
#         attempt += 1
#         geo, zma, final_dist = find_ts(
#             spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
#             bkp_grid, None, ini_thy_info, thy_info, run_prefix,
#             save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
#             attempt=attempt)
#     elif ('addition ' in typ or 'abstraction' in typ) and attempt < 3:
#         # Try to find addition rxn TS in reverse direction with addn tricks
#         babs1 = 170. * phycon.DEG2RAD
#         if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
#             babs1 = 85. * phycon.DEG2RAD
#         print('TS find failed. Attempting to find with '
#               'new angle of attack: {:.1f}'.format(babs1))
#         ts_zma = automol.zmatrix.set_values(ts_zma, {'babs1': babs1})
#         ts_dct['original_zma'] = ts_zma
#         attempt += 1
#         geo, zma, final_dist = find_ts(
#             spc_dct, ts_dct, ts_info, ts_zma, typ, dist_info, grid,
#             bkp_ts_class_data, ini_thy_info, thy_info, run_prefix,
#             save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
#             attempt=attempt)
#     elif 'beta scission' in typ and bkp_ts_class_data and attempt < 2:
#         # Run reverse for a beta scission reaction
#         [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid,
#          bkp_tors_names, bkp_update_guess] = bkp_ts_class_data
#         print('TS find failed. Attempting to find with'
#               'new reaction class: {}'.format(bkp_typ))
#         bkp_dist_info = [bkp_dist_name, 0., bkp_update_guess]
#         ts_dct['class'] = bkp_typ
#         ts_dct['original_zma'] = bkp_ts_zma
#         ts_dct['dist_info'] = bkp_dist_info
#         ts_dct['tors_names'] = bkp_tors_names
#         attempt += 1
#         geo, zma, final_dist = find_ts(
#             spc_dct, ts_dct, ts_info, bkp_ts_zma, bkp_typ, bkp_dist_info,
#             bkp_grid, None, ini_thy_info, thy_info, run_prefix,
#             save_prefix, rxn_run_path, rxn_save_path, overwrite=True,
#             attempt=attempt)
#     else:
#         # Give up and return failed information
#         geo = 'failed'
#         zma = 'failed'
#         final_dist = 0.
#     return geo, zma, final_dist
#
# def check_ts_in_fs():
#     """ Check to see if the ts for the class is in the filesystem
#         do this by comparing names
#     """
#
#     cnf_path = cnf_save_fs.trunk.path()
#     print('Found TS at {}'.format(cnf_path))
#     geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
#     zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
#     chk_bkp = False
#     if automol.zmatrix.names(zma) == automol.zmatrix.names(ts_zma):
#         if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1, True):
#             if 'babs1' in automol.zmatrix.names(ts_zma):
#                 babs1 = 170. * phycon.DEG2RAD
#                 if automol.zmatrix.values(ts_zma)['babs1'] == babs1:
#                     babs1 = 85. * phycon.DEG2RAD
#                 ts_zma = automol.zmatrix.set_values(
#                     ts_zma, {'babs1': babs1})
#                 ts_dct['original_zma'] = ts_zma
#                 if not automol.zmatrix.almost_equal(zma, ts_zma, 4e-1):
#                     chk_bkp = True
#             else:
#                 chk_bkp = True
#     else:
#         chk_bkp = True
#
#     is_bkp = False
#     if chk_bkp and bkp_ts_class_data:
#         [bkp_typ, bkp_ts_zma, bkp_dist_name, bkp_grid, bkp_tors_names,
#          bkp_update_guess] = bkp_ts_class_data
#         if automol.zmatrix.names(zma) == automol.zmatrix.names(bkp_ts_zma):
#             if automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1, True):
#                 is_bkp = True
#             elif 'babs1' in automol.zmatrix.names(bkp_ts_zma):
#                 babs1 = 170. * phycon.DEG2RAD
#                 if automol.zmatrix.values(bkp_ts_zma)['babs1'] == babs1:
#                     babs1 = 85. * phycon.DEG2RAD
#                 bkp_ts_zma = automol.zmatrix.set_values(
#                  bkp_ts_zma, {'babs1': babs1})
#                 if not automol.zmatrix.almost_equal(zma, bkp_ts_zma, 4e-1):
#                     is_bkp = True
#     if not chk_bkp:
#         print("TS is type {}".format(typ))
#     elif is_bkp:
#         print('updating reaction class to {}'.format(bkp_typ))
#         ts_dct['class'] = bkp_typ
#         ts_dct['original_zma'] = bkp_ts_zma
#         ts_dct['dist_info'] = bkp_dist_info
#         ts_dct['tors_names'] = bkp_tors_names
#         print("TS is backup type {}".format(bkp_typ))
#     else:
#         print("TS may not be original type or backup type")
#         print("Some part of the z-matrices have changed")
#     print('class test:', ts_dct['class'])
#     vals = automol.zmatrix.values(zma)
#     final_dist = vals[dist_name]
#     dist_info[1] = final_dist
#     print('dist_info is being set at end of backup checking',
#           dist_info[1], final_dist)
#
#     return final_dst
