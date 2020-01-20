"""
  Therm Calculations
"""

import automol.inchi
import automol.geom
import thermo
from lib.phydat import phycon


REF_CALLS = {"basic": "get_basic",
             "cbh0": "get_cbhzed",
             "cbh1": "get_cbhone",
             "cbh2": "get_cbhtwo"}


# Put in load in the keyword parsing
def is_scheme(entry):
    """ Check whether this is a basis set scheme
    """
    calls = REF_CALLS
    return entry in calls.keys()


def get_function_call(scheme):
    """ get function call
    """
    return getattr(thermo.heatform, REF_CALLS[scheme])


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


def prepare_refs(ref_scheme, spc_dct, spc_queue):
    """ add refs to species list as necessary
    """
    # Determine the reference species, lsit of inchis
    if ref_scheme in REF_CALLS:
        msg = 'Determining {} reference molecules for: \n'.format(ref_scheme)
        refs, newmsg = get_refs(spc_queue, spc_dct, ref_scheme)
        msg += newmsg
    else:
        msg = 'Reference set = {}: \n'.format(', '.join(ref_scheme))
        refs = ref_scheme

    # Get a list of the reference species in a list
    unique_refs = []
    spc_ichs = [spc_dct[spc]['ich'] for spc in spc_dct]
    for ref in refs:
        if ref in spc_ichs:
            msg += 'Adding reference species ref_{}\n'.format(ref)
            spc_dct[ref] = create_spec(ref)
            unique_refs.append(ref)
            # spc_dct['ref_' + ref] = create_spec(ref)
            # unique_refs.append('ref_' + ref)
    return unique_refs, msg


def create_spec(val, charge=0,
                mc_nsamp=(True, 3, 1, 3, 100, 12),
                hind_inc=30.):
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


def basis_energy(spc_bas, spc_dct):
    """ return the electronic + zero point energies for a set of species
    """
    h_basis = []
    for ich in spc_bas:
        for key in spc_dct:
            if ich == spc_dct[key]['ich']:
                ene = spc_dct[key]['ene'] + spc_dct[key]['zpe']/phycon.EH2KCAL
                h_basis.append(ene)
                break
    return h_basis


def get_hf0k(spc, spc_dct, spc_bas, coeff):
    """ determine the 0 K heat of formation from the
    species dictionary and a set of references species
    """
    spc_ene = spc_dct[spc]['ene'] + spc_dct[spc]['zpe']/phycon.EH2KCAL
    h_basis = basis_energy(spc_bas, spc_dct)
    print('hf0k test:', spc, spc_dct[spc]['ene'], spc_dct[spc]['zpe'], spc_ene,
          spc_bas, h_basis, coeff)

    # Get the 0 K heat of formation
    # ref_set should be a parameter for this routine
    h0form = thermo.heatform.calc_hform_0k(
        spc_ene, h_basis, spc_bas, coeff, ref_set='ANL0')
    return h0form
