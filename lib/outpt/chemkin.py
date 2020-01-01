"""
  Library for writing output in CHEMKIN formatted files
"""

import os
import automol
import thermo


def run_ckin_header(pf_info, ref_info, spc_model):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """
    tors_model, vib_model, _ = spc_model
    har_info, tors_info, vpt2_info, _, sp_str, = pf_info
    har_ref_info, tors_ref_info, vpt2_ref_info, _, _ = ref_info

    # Convert the pac99 polynomial to chemkin polynomial
    chemkin_header_str = '! tors model: {0}\n'.format(tors_model)
    chemkin_header_str += '! vib model: {0}\n'.format(vib_model)
    if har_info:
        chemkin_header_str += '! har level: {}{}/{}//{}{}/{}\n'.format(
            har_info[3], har_info[1], har_info[2],
            har_ref_info[3], har_ref_info[1], har_ref_info[2])
    if tors_info:
        chemkin_header_str += '! tors level: {}{}/{}//{}{}/{}\n'.format(
            tors_info[3], tors_info[1], tors_info[2],
            tors_ref_info[3], tors_ref_info[1], tors_ref_info[2])
    if vpt2_info:
        chemkin_header_str += '! vpt2 level: {}{}/{}//{}{}/{}\n'.format(
            vpt2_info[3], vpt2_info[1], vpt2_info[2],
            vpt2_ref_info[3], vpt2_ref_info[1], vpt2_ref_info[2])
    if sp_str:
        chemkin_header_str += sp_str

    return chemkin_header_str


def run_ckin_poly(spc, spc_dct_i, pac99_poly_str):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """

    hf_str = '! Hf(0 K) = {:.2f}, Hf(298 K) = {:.2f} kcal/mol\n'.format(
        float(spc_dct_i['Hfs'][0]), float(spc_dct_i['Hfs'][1]))
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)
    atom_dict = thermo.util.get_atom_counts_dict(formula)
    chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(
        spc, atom_dict, hf_str, pac99_poly_str)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)
    return chemkin_poly_str


def write_nasa_file(spc_dct_i, ckin_path, nasa_path, chemkin_poly_str):
    """ write out the nasa polynomials
    """
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)
    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
