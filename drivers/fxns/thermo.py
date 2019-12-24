""" drivers for thermochemistry evaluations
"""

import automol
import scripts
import autofile

# new libs
from lib.phydat import phycon, eleclvl, symm
import lib.filesystem.inf.get_thy_info


REF_CALLS = {"basic": "get_basic",
             "cbh0": "get_cbhzed",
             "cbh1": "get_cbhone",
             "cbh2": "get_cbhtwo"}


def fix_tsk_lst(tsk_info_lst):
    """ Fix any issues in tsk_list
    """
    return fix(tsk_info_lst)


def add_spc_refs_to_lst(spc_queue, spc_dct, ref):
    """ Add reference molecules to the species list
    """
    for spc in spc_queue:
        if 'ich' not in spc_dct[spc]:
            spc_dct[spc]['ich'] = automol.smiles.inchi(spc_dct[spc]['smiles'])
    refs, msg = prepare_refs(ref, spc_dct, spc_queue)
    full_queue = spc_queue + refs
    full_queue = list(dict.fromkeys(full_queue))
    print(msg)

    return full_queue


def run_es_for_thermo(full_queue, tsk_info_lst, spc_dct,
                      run_prefix, save_prefix):
    """  Run the ESDriver code to calculate thermo data
    """
    run_species = [
        {'species': full_queue,
         'reacs': [],
         'prods': []}
    ]
    esdriver.driver.run(
        tsk_info_lst, eleclvl.DCT, run_species, spc_dct,
        run_prefix, save_prefix)


def set_pf_theory_levels():
    """ something
    """
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
            if 'grad' in tsk[0]:
                grad = True
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
    geo_thy_info = get_thy_info(eleclvl.DCT, geo_lvl)
    harm_thy_info = get_thy_info(eleclvl.DCT, harm_lvl)
    tors_thy_info = get_thy_info(eleclvl.DCT, tors_lvl)
    anharm_thy_info = get_thy_info(eleclvl.DCT, anharm_lvl)
    sym_thy_info = get_thy_info(eleclvl.DCT, sym_lvl)
    pf_levels = [
        harm_thy_info,
        tors_thy_info,
        anharm_thy_info,
        sym_thy_info]

    harm_ref_thy_info = get_thy_info(eleclvl.DCT, harm_lvl_ref)
    tors_ref_thy_info = get_thy_info(eleclvl.DCT, tors_lvl_ref)
    anharm_ref_thy_info = get_thy_info(eleclvl.DCT, anharm_lvl_ref)
    sym_ref_thy_info = get_thy_info(eleclvl.DCT, sym_lvl_ref)
    ref_levels = [
        harm_ref_thy_info,
        tors_ref_thy_info,
        anharm_ref_thy_info,
        sym_ref_thy_info]

    return spc_model, pf_levels, ref_levels, global_pf_str


def get_pf_inp_and_ene(full_queue, spc_dct,
                       spc_model, pf_levels,
                       save_prefix):
    """ Collect the PF input for each species
        Initialize the ene for each of the species
    """
    for spc in full_queue:
        spc_info = scripts.es.get_spc_info(spc_dct[spc])
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

        zpe, zpe_str = scripts.thermo.get_zpe(
            spc, spc_dct[spc], spc_save_path, pf_levels, spc_model)
        spc_str = scripts.thermo.get_spc_input(
            spc, spc_dct[spc], spc_info, spc_save_path, pf_levels, spc_model)

        spc_dct[spc]['spc_info'] = spc_info
        spc_dct[spc]['spc_save_path'] = spc_save_path
        spc_dct[spc]['zpe'] = zpe
        spc_dct[spc]['zpe_str'] = zpe_str
        spc_dct[spc]['spc_str'] = spc_str
        spc_dct[spc]['ene'] = 0

    return spc_dct


def run_messpf(spc_queue, spc_dct,
               harm_thy_info, global_pf_str):
    """ Write and run the MESSPF files
    """
    for spc in spc_queue:
        spc_save_path = spc_dct[spc]['spc_save_path']
        spc_str = spc_dct[spc]['spc_str']
        zpe_str = spc_dct[spc]['zpe_str']
        spc_info = spc_dct[spc]['spc_info']
        pf_input = scripts.thermo.get_pf_input(
            spc, spc_str, global_pf_str, zpe_str)

        pf_path, nasa_path = scripts.thermo.get_thermo_paths(
            spc_save_path, spc_info, harm_thy_info)
        spc_dct[spc]['pf_path'] = pf_path
        spc_dct[spc]['nasa_path'] = nasa_path

        scripts.thermo.write_pf_input(pf_input, pf_path)
        scripts.thermo.run_pf(pf_path)


def set_pf_levels(tsk_info_lst, spc_dct, full_queue, save_prefix, ene_coeff):
    """ Build list of theo levels later used to write methods in CHEMKIN file
    """
    ene_strl = []
    ene_str = '! energy level:'
    ene_lvl = ''
    ene_idx = 0
    for tsk in tsk_info_lst:
        if 'ene' in tsk[0]:
            if ene_idx > len(ene_coeff)-1:
                print('WARNING:',
                      'an insufficient energy coefficient list was provided')
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
                spc_info = spc_dct[spc]['spc_info']
                ene = scripts.thermo.get_electronic_energy(
                    spc_info, geo_thy_info, ene_thy_info, save_prefix)
                print('ene test:', ene, ene_coeff[ene_idx], ene_idx)
                spc_dct[spc]['ene'] += ene*ene_coeff[ene_idx]
            ene_idx += 1
    ene_str += '!               '.join(ene_strl)
    pf_levels.append(ene_str)

    return pf_levels

def get_hf0k(spc_dct, spc_queue, ref):
    """ Get the Hf_0K  value
    """
    # Need to get zpe for reference molecules too
    calc_bas = True
    if isinstance(ref, list):
        spc_bas = ref
        calc_bas = False
    elif is_scheme(ref) or not ref:
        reference_function = get_function_call(ref)
    
    # Calculate Hf_0K for each species
    for spc in spc_queue:
        if calc_bas:
            spc_bas, clist = get_ref(spc, spc_dct, reference_function)
        hf0k = scripts.thermo.get_hf0k(spc, spc_dct, spc_bas, clist)
        spc_dct[spc]['Hfs'] = [hf0k]


def get_nasa_polynomial(spc_dct, spc_model, spc_queue,
                        pf_levels, ref_levels):
    """ Get the NASA polynomials
    """

    # Write the CHEMKIN-formatted NASA polynomials for each species
    chemkin_header_str = scripts.thermo.run_ckin_header(
        pf_levels, ref_levels, spc_model)
    chemkin_set_str = chemkin_header_str
    for spc in spc_queue:
        # Set paths for changing dirs, needed to run thermo codes
        pf_path = spc_dct[spc]['pf_path']
        nasa_path = spc_dct[spc]['nasa_path']
        starting_path = scripts.thermo.go_to_path(nasa_path)

        # Write the ThermP input file
        scripts.thermo.write_thermp_inp(spc_dct[spc])
        
        # Check if info is available to calculate NASA polynomial
        if spc_dct[spc]['ene'] == 0.0 or spc_dct[spc]['spc_str'] == '':
            print('Cannot generate thermo for species {} because information is still missing:'.format(spc_dct[spc]['ich']))
            continue

        # Run ThermP and PAC99 to calculate NASA polynomial
        _ = scripts.thermo.run_thermp(pf_path, nasa_path)
        pac99_poly_str = scripts.thermo.run_pac(spc_dct[spc], nasa_path)

        # Get the NASA polynomial in CHEMKIN format
        chemkin_poly_str = scripts.thermo.run_ckin_poly(
            spc, spc_dct[spc], pac99_poly_str)

        # Write the NASA polynomial to a CHEMKIN file
        chemkin_spc_str = chemkin_header_str + chemkin_poly_str
        chemkin_set_str += chemkin_poly_str
        scripts.thermo.go_to_path(starting_path)
        ckin_path = scripts.thermo.prepare_path(starting_path, 'ckin')
        if not os.path.exists(ckin_path):
            os.makedirs(ckin_path)
        scripts.thermo.write_nasa_file(
            spc_dct[spc], ckin_path, nasa_path, chemkin_spc_str)

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


def create_spec(val, charge=0, mc_nsamp=[True, 3, 1, 3, 100, 12], hind_inc=360.):
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
