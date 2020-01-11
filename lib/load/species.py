""" read species
"""

import os
import chemkin_io
import automol
import autofile
from lib.load import ptt
from lib.phydat import symm, eleclvl


SPC_INP = 'inp/species.dat'


def build_spc_dct(job_path, spc_type):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """
    spc_str = ptt.read_inp_str(job_path, SPC_INP)
    if spc_type == 'csv':
        spc_dct = csv_dct(spc_str, check_stereo=False)
    else:
        raise NotImplementedError

    # Add other things to the species dct
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


def read_spc_amech():
    """ Read an amech style input file for the species
    """
    return None


def modify_spc_dct(job_path, spc_dct):
    """ Modify the species dct
    """
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

        # Add the optional things
        if (ich, mul) in eleclvl.DCT:
            mod_spc_dct[spc]['elec_levs'] = eleclvl.DCT[(ich, mul)]
        if (ich, mul) in symm.DCT:
            mod_spc_dct[spc]['sym'] = symm.DCT[(ich, mul)]
        if ich in geom_dct:
            mod_spc_dct[spc]['geo_obj'] = geom_dct[ich]
        # if 'hind_inc' not in spc_dct:
        #     mod_spc_dct[spc]['hind_inc'] = 30.0 * phycon.DEG2RAD
        # else:
        #     mod_spc_dct[spc]['hind_inc'] = spc_dct[spc]['hind_inc']

    return mod_spc_dct


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
