"""
  Therm Calculations
"""

import thermo

# New
from lib.phydat import phycon


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
