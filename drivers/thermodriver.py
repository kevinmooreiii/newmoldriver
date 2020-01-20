""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs
"""

import os
import autofile.fs
import routines
from lib.runner import therm as thmrunner
from lib.outpt import chemkin as cout
from lib import msg


def run(spc_dct, model_dct,
        run_inp_dct,
        run_jobs_dct,
        run_pf=True,
        run_thermo=True):
    """ main driver for thermo run
    """

    # Print the header message for the driver
    msg.program_header('thermo')

    # Add reference molecules
    refs, msg = routines.pf.therm.prepare_refs(ref, spc_dct, spc_queue)
    full_queue = spc_queue + refs
    full_queue = list(dict.fromkeys(full_queue))
    print(msg)

    # Write and Run MESSPF inputs to generate the partition functions
    if run_pf:

        # Get PF input header
        temp_step = 100.
        ntemps = 30
        global_pf_str = messpf.get_pf_header(temp_step, ntemps)

        # Gather PF model and theory level info
        pf_levels, ref_levels, ts_model = lmech.set_model_info(
            tsk_info_lst)

        # Write and Run the MESSPF file to get the partition function
        for spc in full_queue:

            # Set up the species filesystem
            spc_info = scripts.es.get_spc_info(spc_dct[spc])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)

            # Read the ZPVE from the filesystem
            _, zpe_str = routines.pf.get_zpe(
                spc, spc_dct[spc], spc_save_path, pf_levels, spc_model)

            # Generate the partition function
            spc_block = routines.pf.messf.blocks.species_block(
                spc=spc,
                spc_dct_i=spc_dct[spc],
                spc_info=spc_info,
                spc_model=spc_model,
                pf_levels=pf_levels,
                save_prefix=spc_save_path,
                )
            spc_str = spc_block[0]

            # Write the MESSPF input file
            pf_input = messpf.get_pf_input(
                spc, spc_str, global_pf_str, zpe_str)
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            messpf.write_pf_input(pf_input, pf_path)

            # Run the MESSPF File
            thmrunner.run_pf(pf_path)

    # Use MESS partition functions to compute thermo quantities
    if run_thermo:

        # Setup the CHEMKIN level comment string
        ene_str = cout.get_ckin_ene_lvl_str(
            ts_tsk_lst, thy_dct, ene_coeff)

        # Read the high-level energy
        for spc in full_queue:
            spc_info = spc_dct[spc]['spc_info']
            ene = moldr.pf.get_high_level_energy(
                spc_info=spc_info,
                thy_low_level=geo_thy_info,
                thy_high_level=sp_thy_info,
                save_prefix=save_prefix,
                saddle=False)
            print('ene test:', ene, ene_coeff[ene_idx], ene_idx)
            spc_dct[spc]['ene'] += ene*ene_coeff[ene_idx]

        # Need to get zpe for reference molecules too
        calc_bas = True
        if isinstance(ref, list):
            spc_bas = ref
            calc_bas = False
        elif routines.pf.therm.is_scheme(ref) or not ref:
            reference_function = routines.pf.therm.get_function_call(ref)
        for spc in spc_queue:
            if calc_bas:
                spc_bas, clist = routines.pf.therm.get_ref(
                    spc, spc_dct, reference_function)
            hf0k = routines.of.therm.get_hf0k(spc, spc_dct, spc_bas, clist)
            spc_dct[spc]['Hfs'] = [hf0k]

        # Write the NASA polynomials in CHEMKIN format
        chemkin_header_str = cout.run_ckin_header(
            pf_levels, ref_levels, spc_model)
        chemkin_set_str = chemkin_header_str
        for spc in spc_queue:

            # Set up the paths for running jobs
            pf_path = spc_dct[spc]['pf_path']
            nasa_path = spc_dct[spc]['nasa_path']
            starting_path = thmrunner.go_to_path(nasa_path)

            # Write and run the thermp file to get Hf0k and ...
            thmrunner.write_thermp_inp(spc_dct[spc])
            if spc_dct[spc]['ene'] == 0.0 or spc_dct[spc]['spc_str'] == '':
                print('Cannot generate thermo for species',
                      '{} '.format(spc_dct[spc]['ich']),
                      'because information is still missing:')
                continue
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
