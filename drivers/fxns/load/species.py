""" read species
"""

import chemkin_io
import automol
from ptt import read_inp_str


SPC_INP = 'inp/species.dat'


def read_species_file():
    """ Read the species input file into a string
    """
    return read_inp_str(SPC_INP)


def build_spc_dct(spc_str, spc_type):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """
    if spc_type == 'csv':
        spc_dct = csv_dct(spc_str, check_stereo=False)
    else:
        raise NotImplementedError

    return spc_dct


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


def modify_spc_dct():
    """ Modify the species dct
    """
    for spc in SPC_DCT:
        ich = SPC_DCT[spc]['ich']
        mul = SPC_DCT[spc]['mul']
        if (ich, mul) in eleclvl.DCT:
            SPC_DCT[spc]['elec_levs'] = eleclvl.DCT[(ich, mul)]
        if (ich, mul) in symm.DCT:
            SPC_DCT[spc]['sym'] = symm.DCT[(ich, mul)]
        if ich in GEOM_DCT:
            SPC_DCT[spc]['geo_obj'] = GEOM_DCT[ich]
        SPC_DCT[spc]['hind_inc'] = PARAMS.HIND_INC * phycon.DEG2RAD

    return spc_dct