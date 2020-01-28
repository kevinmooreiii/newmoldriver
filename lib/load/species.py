""" read species
"""

import os
import chemkin_io
import automol
import autoparse.find as apf
import autofile
from lib.load import ptt
from lib.phydat import symm, eleclvl
from lib.reaction import rxnid


CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'
CLA_INP = 'inp/class.csv'


def build_spc_dct(job_path, spc_type):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """
    spc_csv_str = ptt.read_inp_str(job_path, CSV_INP)
    if spc_type == 'csv':
        spc_dct = csv_dct(spc_csv_str, check_stereo=False)
    else:
        raise NotImplementedError

    # Modify spc dct with params from the AMech file
    mod_spc_dct = modify_spc_dct(job_path, spc_dct)

    return mod_spc_dct


def csv_dct(spc_str, check_stereo):
    """ read the species file in a .csv format
    """
    # Read in the initial CSV information
    smi_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'smiles')
    ich_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'inchi')
    mul_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'mult')
    chg_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'charge')
    sens_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'sens')

    # Rebuild with stereochemistry if possible
    print('check_stereo test:', check_stereo, type(check_stereo))
    if check_stereo:
        spc_str = 'name,SMILES,InChI,mult,charge,sens \n'
        for name in ich_dct:
            ich = ich_dct[name]
            smi = smi_dct[name]
            mul = mul_dct[name]
            chg = chg_dct[name]
            sens = sens_dct[name]
            if not automol.inchi.is_complete(ich):
                print('adding stereochemistry for {0}, {1}, {2}'.format(
                    name, smi, ich))
                # this returns a list of ichs w/ different possible stereo vals
                # for now just taking the first of these
                ich = automol.inchi.add_stereo(ich)
                print('new ich possibilities:', ich)
                ich = ich[-1]
                print('new ich:', ich)
                ich_dct[name] = ich
            spc_str += '{0},\'{1}\',\'{2}\',{3},{4},{5} \n'.format(
                name, smi, ich, mul, chg, sens)

        stereo_path = 'species_stereo.csv'
        with open(stereo_path, 'w') as stereo_csv_file:
            stereo_csv_file.write(spc_str)

    # Build the final dictionary
    spc_names = []
    spc_dct = {}
    for name in mul_dct:
        spc_dct[name] = {}
        spc_dct[name]['smi'] = smi_dct[name]
        spc_dct[name]['ich'] = ich_dct[name]
        spc_dct[name]['chg'] = chg_dct[name]
        spc_dct[name]['mul'] = mul_dct[name]
        spc_names.append(name)

    return spc_dct


def read_spc_amech(job_path):
    """ Read an amech style input file for the species
    """

    spc_amech_str = ptt.read_inp_str(job_path, DAT_INP)

    spc_dct = {}
    if spc_amech_str:
        spc_sections = apf.all_captures(
            ptt.end_section_wname2('spc'), spc_amech_str)
        if spc_sections:
            for section in spc_sections:
                name = section[0]
                keyword_dct = ptt.build_keyword_dct(section[1])
                spc_dct[name] = keyword_dct

    return spc_dct


def modify_spc_dct(job_path, spc_dct):
    """ Modify the species dct using input from the additional AMech file
    """

    # Read in other dcts
    amech_dct = read_spc_amech(job_path)
    geom_dct = geometry_dictionary(job_path)

    mod_spc_dct = {}
    for spc in spc_dct:
        # Set the ich and mult
        ich = spc_dct[spc]['ich']
        mul = spc_dct[spc]['mul']
        chg = spc_dct[spc]['chg']
        mod_spc_dct[spc] = {}
        mod_spc_dct[spc]['ich'] = ich
        mod_spc_dct[spc]['mul'] = mul
        mod_spc_dct[spc]['chg'] = chg

        # Add the parameters from amech file
        if spc in amech_dct:
            if 'hind_inc' in amech_dct[spc]:
                mod_spc_dct[spc]['hind_inc'] = amech_dct[spc]['hind_inc']
            if 'hind_def' in amech_dct[spc]:
                mod_spc_dct[spc]['hind_def'] = amech_dct[spc]['hind_def']
            if 'elec_levs' in amech_dct[spc]:
                mod_spc_dct[spc]['elec_levs'] = amech_dct[spc]['elec_levs']
            if 'sym' in amech_dct[spc]:
                mod_spc_dct[spc]['sym'] = amech_dct[spc]['sym']
            if 'mc_nsamp' in amech_dct[spc]:
                mod_spc_dct[spc]['mc_nsamp'] = amech_dct[spc]['mc_nsamp']

        # Add the parameters from std lib if needed
        if 'elec_levs' not in mod_spc_dct[spc]:
            if (ich, mul) in eleclvl.DCT:
                mod_spc_dct[spc]['elec_levs'] = eleclvl.DCT[(ich, mul)]
            else:
                mod_spc_dct[spc]['elec_levs'] = [[0.0, mul]]
        if 'sym' not in mod_spc_dct[spc]:
            if (ich, mul) in symm.DCT:
                mod_spc_dct[spc]['sym'] = symm.DCT[(ich, mul)]

        # Add geoms from geo dct (prob switch to amech file)
        if ich in geom_dct:
            mod_spc_dct[spc]['geo_obj'] = geom_dct[ich]

    return mod_spc_dct


def read_class_dct(job_path):
    """ Read the class dictionary
    """
    cla_str = ptt.read_inp_str(job_path, CLA_INP)
    cla_str = cla_str if cla_str else 'REACTION,RCLASS\nEND'
    cla_dct = chemkin_io.parser.mechanism.reac_class_dct(cla_str, 'class')

    return cla_dct


def build_spc_dct_for_sadpts(spc_dct, rxn_lst, rxn_name_lst,
                             rct_names_lst, prd_names_lst, cla_dct):
    """ build dct
    """
    ts_dct = {}
    ts_idx = 0
    for idx, rxn in enumerate(rxn_lst):
        reacs = rxn['reacs']
        prods = rxn['prods']
        tsname = 'ts_{:g}'.format(ts_idx)
        ts_dct[tsname] = {}
        rname = rxn_name_lst[idx]
        rname_eq = '='.join(rname.split('=')[::-1])
        if rname in cla_dct:
            ts_dct[tsname]['given_class'] = cla_dct[rname]
        elif rname_eq in cla_dct:
            ts_dct[tsname]['given_class'] = cla_dct[rname_eq]
            reacs = rxn['prods']
            prods = rxn['reacs']
        else:
            ts_dct[tsname]['given_class'] = None
        if reacs and prods:
            ts_dct[tsname]['reacs'] = reacs
            ts_dct[tsname]['prods'] = prods
        ts_dct[tsname]['ich'] = ''
        ts_chg = 0
        for rct in rct_names_lst[idx]:
            print(spc_dct[rct])
            ts_chg += spc_dct[rct]['chg']
        ts_dct[tsname]['chg'] = ts_chg
        mul_low, _, rad_rad = rxnid.ts_mul_from_reaction_muls(
            rct_names_lst[idx], prd_names_lst[idx], spc_dct)
        ts_dct[tsname]['mul'] = mul_low
        ts_dct[tsname]['rad_rad'] = rad_rad
        ts_idx += 1

    return ts_dct


def geometry_dictionary(job_path):
    """ read in dictionary of saved geometries
    """
    geom_path = os.path.join(job_path, 'data', 'geoms')
    geom_dct = {}
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_path = os.path.join(dir_path, file_name)
            if file_path.endswith('.xyz'):
                xyz_str = autofile.file.read_file(file_path)
                geo = automol.geom.from_xyz_string(xyz_str)
                ich = automol.geom.inchi(geo)
                if ich in geom_dct:
                    print('Warning: Dupilicate xyz geometry for ', ich)
                geom_dct[ich] = geo
    return geom_dct
