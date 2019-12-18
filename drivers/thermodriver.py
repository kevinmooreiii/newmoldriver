""" drivers for thermochemistry evaluations
"""
import os
import automol.inchi
import automol.geom
import scripts.es
import thermo.heatform
import esdriver
import autofile.fs

# Calling the new libs
from lib.phydat import phycon
from lib.filesystem import build as fbuild


REF_CALLS = {"basic": "get_basic",
             "cbh0": "get_cbhzed",
             "cbh1": "get_cbhone",
             "cbh2": "get_cbhtwo"}


def run(tsk_info_lst, es_dct, spcdct, spc_queue, ref,
        run_prefix, save_prefix, ene_coeff=[1.],
        options=[True, True, True, False]):
    """ main driver for thermo run
    """

    # Determine options
    runes = options[0]      # Run electronic structure theory
    runmess = options[1]    # Run mess/just make the mess input file
    runthermo = options[2]
    if not runmess:
        runthermo = False

    # Fix any issues in tsk_list
    tsk_info_lst = fix(tsk_info_lst)
    # Add reference molecules
    for spc in spc_queue:
        if 'ich' not in spcdct[spc]:
            spcdct[spc]['ich'] = automol.smiles.inchi(spcdct[spc]['smiles'])
    refs, msg = prepare_refs(ref, spcdct, spc_queue)
    full_queue = spc_queue + refs
    full_queue = list(dict.fromkeys(full_queue))
    print(msg)

    # Prepare prefix filesystem
    fbuild.prefix_filesystem(run_prefix, save_prefix)

    # Run ESDriver
    if runes:
        runspecies = [{'species': full_queue, 'reacs': [], 'prods': []}]
        esdriver.run(
            tsk_info_lst, es_dct, runspecies, spcdct, run_prefix, save_prefix)

    if runmess:
        geo_lvl = ''
        harm_lvl = ''
        anharm_lvl = ''
        tors_lvl = ''
        sym_lvl = ''
        harm_lvl_ref = ''
        anharm_lvl_ref = ''
        tors_lvl_ref = ''
        sym_lvl_ref = ''

        # Get PF input header
        temp_step = 100.
        ntemps = 30
        global_pf_str = scripts.thermo.get_pf_header(temp_step, ntemps)

        # Gather PF model and theory level info
        spc_model = ['RIGID', 'HARM', '']
        geom = False
        hess = False
        for tsk in tsk_info_lst:
            if 'samp' in tsk[0] or 'geom' in tsk[0]:
                geo_lvl = tsk[1]
                geom = True
            if 'grad' in tsk[0] or 'hess' in tsk[0]:
                harm_lvl = tsk[1]
                harm_lvl_ref = tsk[2]
                # if 'grad' in tsk[0]:
                #     grad = True
                if 'hess' in tsk[0]:
                    hess = True
                if not geom:
                    ene_lvl = tsk[1]
                    geo_lvl = tsk[1]
            if 'hr' in tsk[0] or 'tau' in tsk[0]:
                tors_lvl = tsk[1]
                tors_lvl_ref = tsk[2]
                if 'md' in tsk[0]:
                    spc_model[0] = 'MDHR'
                if 'tau' in tsk[0]:
                    spc_model[0] = 'TAU'
                else:
                    spc_model[0] = '1DHR'
            if 'anharm' in tsk[0] or 'vpt2' in tsk[0]:
                anharm_lvl = tsk[1]
                anharm_lvl_ref = tsk[2]
                spc_model[1] = 'ANHARM'
                if not hess:
                    geo_lvl = tsk[1]
            if 'sym' in tsk[0]:
                sym_lvl = tsk[1]
                sym_lvl_ref = tsk[2]
                if 'samp' in tsk[0]:
                    spc_model[2] = 'SAMPLING'
                if '1DHR' in tsk[0]:
                    spc_model[2] = '1DHR'
        geo_thy_info = get_thy_info(es_dct, geo_lvl)
        harm_thy_info = get_thy_info(es_dct, harm_lvl)
        tors_thy_info = get_thy_info(es_dct, tors_lvl)
        anharm_thy_info = get_thy_info(es_dct, anharm_lvl)
        sym_thy_info = get_thy_info(es_dct, sym_lvl)
        pf_levels = [harm_thy_info, tors_thy_info,
                     anharm_thy_info, sym_thy_info]

        harm_ref_thy_info = get_thy_info(es_dct, harm_lvl_ref)
        tors_ref_thy_info = get_thy_info(es_dct, tors_lvl_ref)
        anharm_ref_thy_info = get_thy_info(es_dct, anharm_lvl_ref)
        sym_ref_thy_info = get_thy_info(es_dct, sym_lvl_ref)
        ref_levels = [harm_ref_thy_info, tors_ref_thy_info,
                      anharm_ref_thy_info, sym_ref_thy_info]

        # Collect the PF input for each species
        # Initialize the ene for each of the species
        for spc in full_queue:
            spc_info = scripts.es.get_spc_info(spcdct[spc])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)

            zpe, zpe_str = scripts.thermo.get_zpe(
                spc, spcdct[spc], spc_save_path, pf_levels, spc_model)
            spc_str = scripts.thermo.get_spc_input(
                spc, spcdct[spc], spc_info, spc_save_path,
                pf_levels, spc_model)

            spcdct[spc]['spc_info'] = spc_info
            spcdct[spc]['spc_save_path'] = spc_save_path
            spcdct[spc]['zpe'] = zpe
            spcdct[spc]['zpe_str'] = zpe_str
            spcdct[spc]['spc_str'] = spc_str
            spcdct[spc]['ene'] = 0

        # Make and Run the PF file
    if runthermo:
        for spc in spc_queue:
            spc_save_path = spcdct[spc]['spc_save_path']
            spc_str = spcdct[spc]['spc_str']
            zpe_str = spcdct[spc]['zpe_str']
            spc_info = spcdct[spc]['spc_info']
            pf_input = scripts.thermo.get_pf_input(
                spc, spc_str, global_pf_str, zpe_str)

            pf_path, nasa_path = scripts.thermo.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            spcdct[spc]['pf_path'] = pf_path
            spcdct[spc]['nasa_path'] = nasa_path

            scripts.thermo.write_pf_input(pf_input, pf_path)
            scripts.thermo.run_pf(pf_path)

        # Compute Hf0K
        ene_strl = []
        ene_str = '! energy level:'
        ene_lvl = ''
        ene_idx = 0
        for tsk in tsk_info_lst:
            if 'ene' in tsk[0]:
                if ene_idx > len(ene_coeff)-1:
                    print('WARNING:',
                          'insufficient energy coefficient list was provided')
                    break
                ene_lvl = tsk[1]
                geo_lvl = tsk[2]
                geo_thy_info = scripts.es.get_thy_info(es_dct[geo_lvl])
                ene_thy_info = scripts.es.get_thy_info(es_dct[ene_lvl])
                ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
                    ene_coeff[ene_idx], ene_thy_info[3],
                    ene_thy_info[1], ene_thy_info[2],
                    geo_thy_info[3], geo_thy_info[1], geo_thy_info[2]))

                for spc in full_queue:
                    spc_info = spcdct[spc]['spc_info']
                    ene = scripts.thermo.get_electronic_energy(
                        spc_info, geo_thy_info, ene_thy_info, save_prefix)
                    print('ene test:', ene, ene_coeff[ene_idx], ene_idx)
                    spcdct[spc]['ene'] += ene*ene_coeff[ene_idx]
                ene_idx += 1
        ene_str += '!               '.join(ene_strl)
        pf_levels.append(ene_str)
        # Need to get zpe for reference molecules too
        calc_bas = True
        if isinstance(ref, list):
            spc_bas = ref
            calc_bas = False
        elif is_scheme(ref) or not ref:
            reference_function = get_function_call(ref)
        for spc in spc_queue:
            if calc_bas:
                spc_bas, clist = get_ref(spc, spcdct, reference_function)
            hf0k = scripts.thermo.get_hf0k(spc, spcdct, spc_bas, clist)
            spcdct[spc]['Hfs'] = [hf0k]

        chemkin_header_str = scripts.thermo.run_ckin_header(
            pf_levels, ref_levels, spc_model)
        chemkin_set_str = chemkin_header_str
        for spc in spc_queue:
            pf_path = spcdct[spc]['pf_path']
            nasa_path = spcdct[spc]['nasa_path']
            starting_path = scripts.thermo.go_to_path(nasa_path)
            # Change back to starting dir after running thermp and pac99
            # or rest of code is confused
            scripts.thermo.write_thermp_inp(spcdct[spc])
            # run thermp creats thermo and also passed back the 298 K Hf
            if spcdct[spc]['ene'] == 0.0 or spcdct[spc]['spc_str'] == '':
                print('Cannot generate thermo for species',
                      '{} '.format(spcdct[spc]['ich']),
                      'because information is still missing:')
                continue
            hf298k = scripts.thermo.run_thermp(pf_path, nasa_path)
            spcdct[spc]['Hfs'].append(hf298k)
            pac99_poly_str = scripts.thermo.run_pac(spcdct[spc], nasa_path)
            chemkin_poly_str = scripts.thermo.run_ckin_poly(
                spc, spcdct[spc], pac99_poly_str)
            chemkin_spc_str = chemkin_header_str + chemkin_poly_str
            chemkin_set_str += chemkin_poly_str
            scripts.thermo.go_to_path(starting_path)
            ckin_path = scripts.thermo.prepare_path(starting_path, 'ckin')
            if not os.path.exists(ckin_path):
                os.makedirs(ckin_path)
            scripts.thermo.write_nasa_file(
                spcdct[spc], ckin_path, nasa_path, chemkin_spc_str)

        with open(os.path.join(ckin_path, 'automech.ckin'), 'w') as nasa_file:
            nasa_file.write(chemkin_set_str)


def is_scheme(entry):
    """ Check whether this is a basis set scheme
    """
    calls = REF_CALLS
    return entry in calls.keys()


def get_function_call(scheme):
    """ get function call
    """
    calls = REF_CALLS
    if not scheme:
        scheme = 'basic'
    call = getattr(thermo.heatform, calls[scheme])
    return call


def get_ref(spc, spcdct, call):
    """ references for a single species
    """
    ref, clist = call(spcdct[spc]['ich'])
    return ref, clist


def get_refs(species, spcs, scheme):
    """ references for a set of species
    """
    call = get_function_call(scheme)
    ref = []
    msg = ''
    for spc in species:
        spc_ref, _ = call(spcs[spc]['ich'])
        msg += 'Species {} with basis {}\n'.format(spc, ', '.join(spc_ref))
        ref.extend(spc_ref)
    return list(dict.fromkeys(ref)), msg


def prepare_refs(refscheme, spcdct, spc_queue):
    """ add refs to species list as necessary
    """
    if isinstance(refscheme, str):
        msg = 'Determining {} reference molecules for: \n'.format(refscheme)
        if is_scheme(refscheme):
            refs, newmsg = get_refs(spc_queue, spcdct, refscheme)
            msg += newmsg
    else:
        msg = 'Reference set = {}: \n'.format(', '.join(refscheme))
        refs = refscheme
    unique_refs = []
    for ref in refs:
        needtoadd = True
        for spc in spcdct:
            if spcdct[spc]['ich'] == ref:
                needtoadd = False
                unique_refs.append(spc)
                break
        if needtoadd:
            msg += 'Adding reference species ref_{}\n'.format(ref)
            spcdct['ref_' + ref] = create_spec(ref)
            unique_refs.append('ref_' + ref)
    return unique_refs, msg


def create_spec(val, charge=0,
                mc_nsamp=[True, 3, 1, 3, 100, 12], hind_inc=360.):
    """ add a species to the species dictionary
    """
    spec = {}
    if isinstance(val, str):
        ich = val
        print('ich test in create_spec:', ich)
        geo = automol.inchi.geometry(ich)
        zma = automol.geom.zmatrix(geo)
        spec['zmatrix'] = zma
    else:
        geo = val
        ich = automol.geom.inchi(geo)
    form = automol.inchi.formula_dct(ich)
    rad = automol.formula.electron_count(form) % 2
    if rad:
        mult = 2
    else:
        mult = 1
    spec['ich'] = ich
    spec['chg'] = charge
    spec['mul'] = mult
    spec['mc_nsamp'] = mc_nsamp
    spec['hind_inc'] = hind_inc * phycon.DEG2RAD
    return spec


def get_thy_info(es_dct, key):
    """ setup theory info file from es dictionary
    """
    ret = []
    if key:
        ret = scripts.es.get_thy_info(es_dct[key])
    return ret


def fix(tsk_info_lst):
    """ fix the tsk_info_list to make sure the necessary jobs precede it
    """
    has_sp = False
    last_geom = []
    for tsk in tsk_info_lst:
        if tsk[0] == 'conf_energy':
            has_sp = True
        if 'find' in tsk[0] or tsk[0] == 'opt':
            last_geom = ['conf_energy']
            last_geom.extend(tsk[1:-1])
            last_geom.append(False)
    if not has_sp:
        tsk_info_lst.append(last_geom)
    return tsk_info_lst
