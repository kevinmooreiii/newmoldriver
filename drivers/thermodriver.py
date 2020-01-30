""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs
"""

import os
import autofile.fs
import routines
from lib.load import mechanism as loadmech
from lib.filesystem import inf as finf
from lib.runner import therm as thmrunner
from lib.outpt import chemkin as cout
from lib import printmsg


def run(spc_dct,
        model_dct,
        thy_dct,
        rxn_lst,
        run_inp_dct,
        ref_scheme='basic',
        run_pf=True,
        run_thermo=True):
    """ main driver for thermo run
    """

    # Print the header message for the driver
    printmsg.program_header('thermo')

    # Pull stuff from dcts for now
    model = rxn_lst['model']
    save_prefix = run_inp_dct['save_prefix']
    ene_coeff = model_dct[model]['options']['ene_coeff']
    ene_idx = 0

    # Add reference molecules
    spc_queue = []
    for rxn in rxn_lst:
        spc_queue.extend(rxn['species'])
    refs, msg = routines.pf.therm.prepare_refs(ref_scheme, spc_dct, spc_queue)
    print('refs from thermo', refs)
    # full_queue = spc_queue + refs
    # full_queue = list(dict.fromkeys(full_queue))
    # FIX THE QUEUE BUILD
    full_queue = spc_queue + ['C1H4']
    print('full_queue', full_queue)
    print(msg)

    # Gather PF model and theory level info
    pf_levels = loadmech.set_es_model_info(model_dct['es'], thy_dct)
    pf_model = loadmech.set_pf_model_info(model_dct['pf'])

    # Write and Run MESSPF inputs to generate the partition functions
    if run_pf:

        # Get PF input header
        temp_step = 100.
        ntemps = 30
        global_pf_str = routines.pf.messf.pfblock.get_pf_header(
            temp_step, ntemps)

        # Write and Run the MESSPF file to get the partition function
        for spc in full_queue:

            # Set up the species filesystem
            spc_info = finf.get_spc_info(spc_dct[spc])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            print('TEST: spc_save_path')
            print(spc_save_path)

            # Read the ZPVE from the filesystem
            _, zpe_str = routines.pf.get_zpe(
                spc, spc_dct[spc], spc_save_path, pf_levels, pf_model)

            # Generate the partition function
            spc_block = routines.pf.messf.blocks.species_block(
                spc=spc,
                spc_dct_i=spc_dct[spc],
                spc_info=spc_info,
                spc_model=pf_model,
                pf_levels=pf_levels,
                save_prefix=spc_save_path,
                )
            spc_str = spc_block[0]

            # Write the MESSPF input file
            harm_thy_info = pf_levels[2]
            pf_input = routines.pf.messf.pfblock.get_pf_input(
                spc, spc_str, global_pf_str, zpe_str)
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            routines.pf.messf.pfblock.write_pf_input(pf_input, pf_path)

            # Run the MESSPF File
            thmrunner.run_pf(pf_path)

    # Use MESS partition functions to compute thermo quantities
    if run_thermo:

        # Setup the CHEMKIN level comment string
        # ene_str = cout.get_ckin_ene_lvl_str(pf_levels, ene_coeff)

        # Read the high-level energy
        for spc in full_queue:
            spc_info = finf.get_spc_info(spc_dct[spc])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            ene = routines.pf.get_high_level_energy(
                spc_info=spc_info,
                thy_low_level=pf_levels[0],
                thy_high_level=pf_levels[1],
                save_prefix=save_prefix,
                saddle=False)
            zpe = routines.pf.get_zero_point_energy(
                spc, spc_dct[spc], pf_levels, pf_model,
                elec_levels=((0., 1)), sym_factor=1.0,
                save_prefix=spc_save_path)
            print('therm ene test')
            print(ene)
            print(zpe)
            spc_dct[spc]['ene'] = ene*ene_coeff[ene_idx]
            spc_dct[spc]['zpe'] = zpe[0]  # returns zpe, isatom

        # Need to get zpe for reference molecules too
        calc_bas = True
        if isinstance(ref_scheme, list):
            spc_bas = ref_scheme
            calc_bas = False
        elif routines.pf.therm.is_scheme(ref_scheme) or not ref_scheme:
            reference_function = routines.pf.therm.get_function_call(
                ref_scheme)
        for spc in spc_queue:
            if calc_bas:
                spc_bas, clist = routines.pf.therm.get_ref(
                    spc, spc_dct, reference_function)
            hf0k = routines.pf.therm.get_hf0k(spc, spc_dct, spc_bas, clist)
            spc_dct[spc]['Hfs'] = [hf0k]

        # Write the NASA polynomials in CHEMKIN format
        chemkin_header_str = cout.run_ckin_header(
            pf_levels, pf_model)
        chemkin_set_str = chemkin_header_str
        for spc in spc_queue:

            # Set up the paths for running jobs
            harm_thy_info = pf_levels[2]
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            starting_path = thmrunner.go_to_path(nasa_path)

            # Write and run the thermp file to get Hf0k and ...
            thmrunner.write_thermp_inp(spc_dct[spc])
            # not getting spc str, so this isnt working, fix this
            # if spc_dct[spc]['ene'] == 0.0 or spc_dct[spc]['spc_str'] == '':
            #     print('Cannot generate thermo for species',
            #           '{} '.format(spc_dct[spc]['ich']),
            #           'because information is still missing:')
            #     continue
            hf298k = thmrunner.run_thermp(pf_path, nasa_path)
            spc_dct[spc]['Hfs'].append(hf298k)

            # Run PAC99 to get a NASA polynomial string in its format
            pac99_poly_str = thmrunner.run_pac(spc_dct[spc], nasa_path)

            # Convert the polynomial from PAC99 to CHEMKIN
            chemkin_poly_str = cout.run_ckin_poly(
                spc, spc_dct[spc], pac99_poly_str)

            # Write a string for a single spc file to a set
            chemkin_spc_str = chemkin_header_str + chemkin_poly_str
            chemkin_set_str += chemkin_poly_str

            # Write the CHEMKIN string to a file
            thmrunner.go_to_path(starting_path)
            ckin_path = thmrunner.prepare_path(starting_path, 'ckin')
            if not os.path.exists(ckin_path):
                os.makedirs(ckin_path)
            cout.write_nasa_file(
                spc_dct[spc], ckin_path, nasa_path, chemkin_spc_str)

        # Write the NASA polynomial strings to a CHEMKIN-formatted file
        with open(os.path.join(ckin_path, 'automech.ckin'), 'w') as nasa_file:
            nasa_file.write(chemkin_set_str)
