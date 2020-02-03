"""
  Therm Calculations
"""

import automol.inchi
import automol.geom
import thermo
from lib.phydat import phycon


# FUNCTIONS TO PREPARE THE LIST OF REFERENCE SPECIES NEEDED FOR THERM CALCS #
REF_CALLS = {"basic": "get_basic",
             "cbh0": "get_cbhzed",
             "cbh1": "get_cbhone",
             "cbh2": "get_cbhtwo"}


def prepare_refs(ref_scheme, spc_dct, spc_queue):
    """ add refs to species list as necessary
    """
    # Determine the function to be used to get the thermochemistry ref species
    if ref_scheme in REF_CALLS:
        get_ref_fxn = getattr(thermo.heatform, REF_CALLS[ref_scheme])
    else:
        raise NotImplementedError

    # Determine the reference species, list of inchis
    msg = 'Determining {} reference molecules for: \n'.format(ref_scheme)
    refs = []
    for spc in spc_queue:
        spc_ref, _ = get_ref_fxn(spc_dct[spc]['ich'])
        refs.extend(spc_ref)
        msg += 'Species {} with basis {}\n'.format(spc, ', '.join(spc_ref))
    refs = list(dict.fromkeys(refs))

    # Get a list of the reference species in a list
    unique_refs_dct = {}
    for ref in refs:
        if ref not in spc_dct:
            msg += 'Adding reference species ref_{}\n'.format(ref)
            unique_refs_dct[ref] = create_spec(ref)

    return unique_refs_dct, msg


def create_spec(ich, charge=0,
                mc_nsamp=(True, 3, 1, 3, 100, 12),
                hind_inc=30.):
    """ add a species to the species dictionary
    """
    spec = {}
    rad = automol.formula.electron_count(automol.inchi.formula(ich)) % 2
    mult = 1 if not rad else 2
    spec['zmatrix'] = automol.geom.zmatrix(automol.inchi.geometry(ich))
    spec['ich'] = ich
    spec['chg'] = charge
    spec['mul'] = mult
    spec['mc_nsamp'] = mc_nsamp
    spec['hind_inc'] = hind_inc * phycon.DEG2RAD
    return spec


# FUNCTIONS TO CALCULATE ENERGIES FOR THERMOCHEMICAL PARAMETERS #
def basis_energy(spc_bas, spc_dct):
    """ Return the electronic + zero point energies for a set of species.
    """
    h_basis = []
    for ich in spc_bas:
        for key in spc_dct:
            if ich == spc_dct[key]['ich']:
                ene = spc_dct[key]['ene'] + spc_dct[key]['zpe']/phycon.EH2KCAL
                h_basis.append(ene)
                break
    return h_basis


def get_hf0k(spc, spc_dct, spc_bas, coeff, ref_set='ANL0'):
    """ Determine the 0 K heat of formation from the
        species dictionary and a set of references species.
    """
    spc_ene = spc_dct[spc]['ene'] + spc_dct[spc]['zpe']/phycon.EH2KCAL
    h_basis = basis_energy(spc_bas, spc_dct)
    print('hf0k test:', spc, spc_dct[spc]['ene'], spc_dct[spc]['zpe'], spc_ene,
          spc_bas, h_basis, coeff)

    # Get the 0 K heat of formation
    # ref_set should be a parameter for this routine
    h0form = thermo.heatform.calc_hform_0k(
        spc_ene, h_basis, spc_bas, coeff, ref_set=ref_set)
    return h0form
