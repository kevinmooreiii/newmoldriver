""" test
"""

def calc_shift_ene(rxn_spc, rxn, pes_ref_spc, spc_dct):
    """ """
    # Calculate
    ref_ene1, ref_ene2 = 0.0, 0.0
    for spc in pre_ref_spc:
        ref_zero_ene1 += get_high_energy(thy1, mod1) + get_zpe(thy1, mod1)
        ref_zero_ene2 += get_high_energy(thy2, mod2) + get_zpe(thy2, mod2)

    # Determine energies at reference lvl
    rxn_enes_1, rxn_enes_2 = {}, {}
    fix_lst = []
    for spc in rxn:
        # Calculate for level 1
        ene1 = get_high_energy(thy1, mod1)
        zpe1 = get_zpe(thy1, mod1)
        if ene1 and zpe1:
            rxn_enes_1[spc] = ene1 + zpe1
        else:
            fix_lst.append(spc)
        # Calculate for level 2
        ene2 = get_high_energy(thy2, mod2)
        zpe2 = get_zpe(thy2, mod2)
        if ene2 and zpe2:
            rxn_enes_2[spc] = ene2 + zpe2
        else:
            fix_lst.append(spc)

    # Calculate relative ene for rxn_spc at level
    for spc, ene in rxn_enes_2.items():
        if spc != rxn_spc:
            dene_2 = rxn_enes_2[rxn_spc] - rxn_enes_2[spc]
            comp_spc = spc
            break

    # Calculate relative ene at ref level to ref species
    rel_ene_1 = rxn_enes_1[comp_spc] - ref_zero_ene1
    
    # Calculate shifted energy
    ene = dene_2 - rel_ene_1

    return ene


def shifted_sadpt_ene():
    """ calculate the energy of a TS using two methods
    """
    dene_mol1_0 = ene_mol1_0 - ene_ref_0
    dene_mol2_0 = ene_mol2_0 - ene_ref_0
    barrier_mol1_1 = ene_sadpt_1 - e_mol1_1 
    barrier_mol1_2 = ene_sadpt_1 - e_mol1_2
    
    dene_sadpt = (ene_sadpt_1 - e_mol1) + dene_mol1_0

    return dene_sadpt
