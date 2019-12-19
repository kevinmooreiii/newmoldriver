"""
Find a TS from the grid as well as associated vdW wells
"""

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


def find_vdw(ts_name, spc_dct, thy_info, ini_thy_info, vdw_params,
             nsamp_par, run_prefix, save_prefix, kickoff_size, kickoff_backward,
             projrot_script_str, overwrite):
    """ find van der Waals structures for all the pairs of species in a reaction list
    """
    new_vdws = []
    _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_info[:2])
    mul = spc_dct[ts_name]['low_mul']
    vdw_names_lst = []
    if vdw_params[0]:
        vdw_names_lst.append([sorted(spc_dct[ts_name]['reacs']), mul, 'r'])
    if vdw_params[1]:
        vdw_names_lst.append([sorted(spc_dct[ts_name]['prods']), mul, 'p'])

    for names, ts_mul, label in vdw_names_lst:
        if len(names) < 2:
            print("Cannot find van der Waals well for unimolecular reactant or product")
        ichs = list(map(lambda name: spc_dct[name]['ich'], names))
        chgs = list(map(lambda name: spc_dct[name]['chg'], names))
        muls = list(map(lambda name: spc_dct[name]['mul'], names))

        # theory
        prog = thy_info[0]
        method = thy_info[1]
        _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(prog, method)

        geos = []
        ntaudof = 0.
        for name, ich, chg, mul in zip(names, ichs, chgs, muls):
            spc_info = [ich, chg, mul]
            orb_restr = moldr.util.orbital_restriction(spc_info, ini_thy_info)
            ini_thy_level = ini_thy_info[0:3]
            ini_thy_level.append(orb_restr)
            orb_restr = moldr.util.orbital_restriction(spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info)
            spc_run_path = spc_run_fs.leaf.path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
            run_fs = autofile.fs.run(thy_run_path)

            ini_thy_save_fs = autofile.fs.theory(spc_save_path)
            ini_thy_save_fs.leaf.create(ini_thy_level[1:4])

            cnf_run_fs = autofile.fs.conformer(thy_run_path)
            cnf_save_fs = autofile.fs.conformer(thy_save_path)

            ini_fs = [None, ini_thy_save_fs]
            fs = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                  cnf_run_fs, cnf_save_fs, None, None,
                  None, None, run_fs]
    # fs = [None, None, thy_run_fs, thy_save_fs,
          # cnf_run_fs, cnf_save_fs, None, None,
            geo = moldr.geom.reference_geometry(
                spc_dct[name], thy_level, ini_thy_level, fs, ini_fs,
                kickoff_size=kickoff_size,
                kickoff_backward=kickoff_backward,
                projrot_script_str=projrot_script_str,
                overwrite=overwrite)
            geos.append(geo)
            gra = automol.geom.graph(geo)
            ntaudof += len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)
        geo1, geo2 = geos
        geo1 = automol.geom.mass_centered(geo1)
        geo2 = automol.geom.mass_centered(geo2)
        min_ene = 0.
        for idx in range(int(nsamp)):
            print('Optimizing vdw geometry {}/{}'.format(idx+1, nsamp))
            angs1 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs2 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs12 = numpy.multiply(
                numpy.random.rand(2), [1*numpy.pi, 2*numpy.pi])
            geo1 = automol.geom.euler_rotated(geo1, *angs1)
            geo2 = automol.geom.euler_rotated(geo2, *angs2)
            dist_cutoff = 3.0 * phycon.ANG2BOHR

            geo = automol.geom.join(geo1, geo2, dist_cutoff, *angs12)
            print("Species: {}".format('+'.join(names)))
            print('vdw starting geometry')
            print(automol.geom.xyz_string(geo))

   #  set up the filesystem
            ich = automol.inchi.recalculate(automol.inchi.join(ichs))
            chg = sum(chgs)
            mul = ts_mul
            spc_info = (ich, chg, mul)
            #orb_restr = moldr.util.orbital_restriction(mul, thy_info[0:3] restrict_open_shell)
            #orb_restr = restrict_open_shell
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info)
            spc_run_path = spc_run_fs.leaf.path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            orb_restr = moldr.util.orbital_restriction(spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)
            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
            run_fs = autofile.fs.run(thy_run_path)
   #  generate reference geometry
   #  generate the z-matrix and sampling ranges

            moldr.driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                geom=geo,
                spc_info=spc_info,
                thy_level=thy_level,
                run_fs=run_fs,
                script_str=opt_script_str,
                overwrite=overwrite,
                **opt_kwargs,
            )

   #  save info for the initial geometry (from inchi or from save directory)
            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                print('Saving reference geometry')
                print(" - Save path: {}".format(thy_save_path))

                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                geo = elstruct.reader.opt_geometry(prog, out_str)
                print('vdw ending geometry')
                print(automol.geom.xyz_string(geo))
                thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
                ene = elstruct.reader.energy(prog, method, out_str)
                if ene < min_ene:
                    min_ene = ene
                    print('ene test in vdw')
                    print(ene)
                    thy_save_fs.leaf.file.energy.write(ene, thy_level[1:4])
                    print('Saving reference geometry')
                    print(" - Save path: {}".format(thy_save_path))
                    vdw_name = label + ts_name.replace('ts', 'vdw')
                    spc_dct[vdw_name] = spc_dct[ts_name].copy()
                    spc_dct[vdw_name]['ich'] = ich
                    spc_dct[vdw_name]['mul'] = mul
                    spc_dct[vdw_name]['chg'] = chg
                    spc_dct[vdw_name]['dist_info'][1] = dist_cutoff
                    fs = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                          cnf_run_fs, cnf_save_fs, None, None,
                          None, None, run_fs]
                    #Make a fake conformer
                    cnf_save_fs = autofile.fs.conformer(thy_save_path)
                    cnf_run_fs = autofile.fs.conformer(thy_run_path)
                    cnf_save_fs.trunk.create()
                    cnf_run_fs.trunk.create()
                    tors_range_dct = {}
                    cinf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
                    cinf_obj.nsamp = 1
                    cnf_save_fs.trunk.file.info.write(cinf_obj)
                    locs_lst = cnf_save_fs.leaf.existing()
                    if not locs_lst:
                        cid = autofile.system.generate_new_conformer_id()
                        locs = [cid]
                    else:
                        locs = locs_lst[0]
                    cnf_save_fs.leaf.create(locs)
                    cnf_run_fs.leaf.create(locs)
                    cnf_save_fs.leaf.file.geometry_info.write(
                        inf_obj, locs)
                    cnf_save_fs.leaf.file.geometry_input.write(
                        inp_str, locs)
                    cnf_save_fs.leaf.file.energy.write(ene, locs)
                    cnf_save_fs.leaf.file.geometry.write(geo, locs)
        if min_ene:
            new_vdws.append(vdw_name)

    return new_vdws


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