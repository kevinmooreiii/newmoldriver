"""
Read the mechanism file
"""

import chemkin_io
import automol
from ptt import read_inp_str


MECH_INP = 'inp/mechanism.dat'


def read_theory_file():
    """ Read the mechanism input file into a string
    """
    return read_inp_str(MECH_INP)


def parse_mechanism_file(mech_str, mech_type,
                         check_stereo=False, sort_rxns=False,
                         rad_rad_sort=False):
    """ Get the reactions and species from the mechanism input
    """
    # parsing moved to the input parsing module I am writing
    if mech_type.lower() == 'chemkin':
        spc_dct, rct_names, prd_names, rxn_name, form_strs = _parse_chemkin(
            mech_str, sort_rxns)
    # elif mech_type.lower() == 'json':
    #     spc_dct, rct_names, prd_names, rxn_name, form_strs = _parse_json(
    #         mech_path, mech_file, check_stereo, rad_rad_sort)
    else:
        raise NotImplementedError

    return spc_dct, rct_names, prd_names, rxn_name, form_strs


def _parse_chemkin(mech_str, sort_rxns):
    """ parse a chemkin formatted mechanism file
    """

    # Read the reactions and participaring species from the mech file
    rxn_block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    rxn_strs = chemkin_io.parser.reaction.data_strings(rxn_block_str)
    rct_names_lst = list(
        map(chemkin_io.parser.reaction.reactant_names, rxn_strs))
    prd_names_lst = list(
        map(chemkin_io.parser.reaction.product_names, rxn_strs))

    # Sort reactant and product name lists by formula to facilitate
    # multichannel, multiwell rate evaluations
    formula_str = ''
    rxn_name_lst = []
    formula_str_lst = []
    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        rxn_name_lst.append(rxn_name)
        rct_ichs = list(map(ich_dct.__getitem__, rct_names))
        formula_dct = ''
        for rct_ich in rct_ichs:
            formula_i_dct = automol.inchi.formula_dct(rct_ich)
            formula_dct = automol.formula.join(formula_dct, formula_i_dct)
        formula_str = automol.formula.string(formula_dct)
        formula_str_lst.append(formula_str)

    rxn_info_lst = list(
        zip(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst))

    # Sort the reactions if desired
    print('sort_rxns test:', sort_rxns)
    if sort_rxns:
        rxn_info_lst.sort()
        formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst = zip(
            *rxn_info_lst)

    # Build the total PES dct
    pes_dct = build_pes_dct(
        formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst)

    return pes_dct


def build_pes_dct(formula_str_lst, rct_names_lst,
                  prd_names_lst, rxn_name_lst):
    """ Build a dictionary of the PESs
    """
    pes_dct = {}
    current_formula = ''
    for fidx, formula in enumerate(formula_str_lst):
        if current_formula == formula:
            pes_dct[formula]['rct_names_lst'].append(rct_names_lst[fidx])
            pes_dct[formula]['prd_names_lst'].append(prd_names_lst[fidx])
            pes_dct[formula]['rxn_name_lst'].append(rxn_name_lst[fidx])
        else:
            current_formula = formula
            pes_dct[formula] = {}
            pes_dct[formula]['rct_names_lst'] = [rct_names_lst[fidx]]
            pes_dct[formula]['prd_names_lst'] = [prd_names_lst[fidx]]
            pes_dct[formula]['rxn_name_lst'] = [rxn_name_lst[fidx]]

    return pes_dct


def print_pes(pes_dct):
    """ Print the PES
    """
    for pes_idx, pes in enumerate(pes_dct, start=1):
        print('pes test:', pes_idx, pes)
        pes_rxn_name_lst = pes_dct[pes]['rxn_name_lst']
        pes_rct_names_lst = pes_dct[pes]['rct_names_lst']
        pes_prd_names_lst = pes_dct[pes]['prd_names_lst']
        for chn_idx, _ in enumerate(pes_rxn_name_lst):
            print('channel {}: {} = {}'.format(
                chn_idx+1, ' + '.join(pes_rct_names_lst[chn_idx]),
                ' + '.join(pes_prd_names_lst[chn_idx])))


def determine_connected_pes_channels(pes_dct):
    """ Sort the PES lst and set the connected channels
    """
    for pes_idx, pes in enumerate(pes_dct, start=1):
        if pes_idx in pesnums:
            pes_rct_names_lst = pes_dct[pes]['rct_names_lst']
            pes_prd_names_lst = pes_dct[pes]['prd_names_lst']
            pes_rxn_name_lst = pes_dct[pes]['rxn_name_lst']
            if isinstance(channels, str):
                if channels == 'all':
                    print(len(pes_rxn_name_lst))
                    pes_chns = numpy.arange(len(pes_rxn_name_lst)+1)
                elif '-' in channels:
                    start, end = channels.split('-')
                    pes_chns = numpy.arange(int(start), int(end)+1)
                elif '[' in channels:
                    nums = channels.replace('[', '').replace(']', '').split(',')
                    pes_chns = [int(num) for num in nums]
            print('for pes:', pes_idx)

            # Split up sub-pes within a formula
            subpes_idx = 0
            conndct = {}
            connchnls = {}
            for chnl_idx, _ in enumerate(pes_rxn_name_lst):
                connected_to = []
                chnl_species = [list(pes_rct_names_lst[chnl_idx]),
                                list(pes_prd_names_lst[chnl_idx])]
                for conn_chnls_idx in conndct:
                    for spc_pair in chnl_species:
                        if len(spc_pair) == 1:
                            if spc_pair in conndct[conn_chnls_idx]:
                                if conn_chnls_idx not in connected_to:
                                    connected_to.append(conn_chnls_idx)
                            elif spc_pair[::-1] in conndct[conn_chnls_idx]:
                                if conn_chnls_idx not in connected_to:
                                    connected_to.append(conn_chnls_idx)
                if not connected_to:
                    conndct[subpes_idx] = chnl_species
                    connchnls[subpes_idx] = [chnl_idx]
                    subpes_idx += 1
                else:
                    conndct[connected_to[0]].extend(chnl_species)
                    connchnls[connected_to[0]].append(chnl_idx)
                    if len(connected_to) > 1:
                        for cidx, cval in enumerate(connected_to):
                            if cidx > 0:
                                conn_specs = conndct.pop(cval, None)
                                conn_chnls = connchnls.pop(cval, None)
                                conndct[connected_to[0]].extend(conn_specs)
                                connchnls[connected_to[0]].extend(conn_chnls)
                    for cidx in conndct:
                        conndct[cidx].sort()
                        conndct[cidx] = [
                            conndct[cidx][i] for i in
                            range(len(conndct[cidx])) if i == 0 or
                            conndct[cidx][i] != conndct[cidx][i-1]]

            for cidx, cvals in enumerate(connchnls.values()):
                print('PES{}_{}: {} surface'.format(
                    str(pes_idx), str(cidx+1), pes))
                for cval in cvals:
                    print(pes_rxn_name_lst[cval])

            return connchnls


# def parse_json():
#     """ parse a json file mechanism file
#     """
#
#     with open(os.path.join(mech_path, mech_file)) as mechfile:
#         inp_data = json.load(
#             mechfile, object_pairs_hook=collections.OrderedDict)
#         mech_data = []
#     for reaction in inp_data:
#         if isinstance(reaction, dict):
#             mech_data = inp_data
#             break
#         else:
#             for entry in inp_data[reaction]:
#                 mech_data.append(entry)
#
#     # Convert essential pieces of json file to chemkin formatted data so
#     # (i) can easily remove species that don't really exist
#     # (ii) revise products of reactions for species that don't exist
#     # (iii) do the stereochemistry generation only one
#
#     formula_str = ''
#     formula_str_lst = []
#     rxn_name_lst = []
#     rct_names_lst = []
#     prd_names_lst = []
#     rct_smis_lst = []
#     rct_ichs_lst = []
#     rct_muls_lst = []
#     prd_smis_lst = []
#     prd_ichs_lst = []
#     prd_muls_lst = []
#     prd_names_lst = []
#     rxn_sens = []
#     rxn_unc = []
#     rxn_val = []
#     rxn_fam = []
#     unq_rxn_lst = []
#     fll_rxn_lst = []
#     # idxp = 0
#     for idx, reaction in enumerate(mech_data):
#         if 'Reactants' in reaction and 'Products' in reaction:
#             print(idx, reaction['name'])
#             if reaction['name'] in fll_rxn_lst:
#                 print('duplicate reaction found:', reaction['name'], idx)
#             else:
#                 unq_rxn_lst.append(reaction['name'])
#             fll_rxn_lst.append(reaction['name'])
#     print('reaction duplicate test:', len(unq_rxn_lst), len(fll_rxn_lst))
#
#     for _, reaction in enumerate(mech_data):
#         # set up reaction info
#         rct_smis = []
#         rct_ichs = []
#         rct_muls = []
#         rct_names = []
#         prd_smis = []
#         prd_ichs = []
#         prd_muls = []
#         prd_names = []
#         if 'Reactants' in reaction and 'Products' in reaction:
#             for rct in reaction['Reactants']:
#                 rct_names.append(rct['name'])
#                 rct_smis.append(rct['SMILES'][0])
#                 ich = rct['InChi']
#                 if check_stereo:
#                     if not automol.inchi.is_complete(ich):
#                         print('adding stereochemsiry for {}'.format(ich))
#                         ich = automol.inchi.add_stereo(rct['InChi'])[0]
#                 rct_ichs.append(ich)
#                 rct_muls.append(rct['multiplicity'])
#             rad_rad_reac = True
#             if len(rct_ichs) == 1:
#                 rad_rad_reac = False
#             else:
#                 if min(rct_muls) == 1:
#                     rad_rad_reac = False
#             for prd in reaction['Products']:
#                 prd_names.append(prd['name'])
#                 prd_smis.append(prd['SMILES'][0])
#                 ich = prd['InChi']
#                 if check_stereo:
#                     if not automol.inchi.is_complete(ich):
#                         print('adding stereochemsiry for {}'.format(ich))
#                         ich = automol.inchi.add_stereo(prd['InChi'])[0]
#                 prd_ichs.append(ich)
#                 prd_muls.append(prd['multiplicity'])
#             rad_rad_prod = True
#             if len(prd_ichs) == 1:
#                 rad_rad_prod = False
#             else:
#                 if min(prd_muls) == 1:
#                     rad_rad_prod = False
#
#             print('rad_rad_sort test:', rad_rad_sort)
#             if rad_rad_sort and not rad_rad_reac and not rad_rad_prod:
#                 continue
#             rct_smis_lst.append(rct_smis)
#             rct_ichs_lst.append(rct_ichs)
#             rct_muls_lst.append(rct_muls)
#             rct_names_lst.append(rct_names)
#             prd_smis_lst.append(prd_smis)
#             prd_ichs_lst.append(prd_ichs)
#             prd_muls_lst.append(prd_muls)
#             prd_names_lst.append(prd_names)
#         rxn_name_lst.append(reaction['name'])
#         if 'Sensitivity' in reaction:
#             rxn_sens.append(reaction['sensitivity'])
#         else:
#             rxn_sens.append('')
#         if 'Uncertainty' in reaction:
#             rxn_unc.append(reaction['uncertainty'])
#         else:
#             rxn_unc.append('')
#         if 'value' in reaction:
#             rxn_val.append(reaction['value'])
#         else:
#             rxn_val.append('')
#         if 'family' in reaction:
#             rxn_fam.append(reaction['family'])
#         else:
#             rxn_fam.append('')
#
#         formula_dct = ''
#         for rct_ich in rct_ichs:
#             formula_i_dct = automol.inchi.formula_dct(rct_ich)
#             formula_dct = automol.formula.join(
#                 formula_dct, formula_i_dct)
#         formula_str = automol.formula.string(formula_dct)
#         formula_str_lst.append(formula_str)
#
#     print('rct_muls_lst:', rct_muls_lst)
#     unq_ich_lst = []
#     unq_mul_lst = []
#     unq_smi_lst = []
#     unq_lab_lst = []
#     unq_lab_idx_lst = []
#     csv_str = 'name,smiles,mult'
#     csv_str += '\n'
#     spc_str = 'species'
#     spc_str += '\n'
#     for ichs, muls, smis in zip(rct_ichs_lst, rct_muls_lst, rct_smis_lst):
#         for ich, mul, smi in zip(ichs, muls, smis):
#             unique = True
#             for unq_ich, unq_mul in zip(unq_ich_lst, unq_mul_lst):
#                 if ich == unq_ich and mul == unq_mul:
#                     unique = False
#             if unique:
#                 unq_ich_lst.append(ich)
#                 unq_mul_lst.append(mul)
#                 unq_smi_lst.append(smi)
#
#                 formula_dct = automol.inchi.formula_dct(ich)
#                 lab = automol.formula.string(formula_dct)
#
#                 unq_lab_lst.append(lab)
#                 lab_idx = -1
#                 for lab_i in unq_lab_lst:
#                     if lab == lab_i:
#                         lab_idx += 1
#                 unq_lab_idx_lst.append(lab_idx)
#                 if lab_idx == 0:
#                     label = lab
#                 else:
#                     label = lab + '(' + str(lab_idx) + ')'
#                 csv_str += ','.join([label, smi, str(mul)])
#                 csv_str += '\n'
#                 spc_str += label
#                 spc_str += '\n'
#     for ichs, muls, smis in zip(prd_ichs_lst, prd_muls_lst, prd_smis_lst):
#         for ich, mul, smi in zip(ichs, muls, smis):
#             unique = True
#             for unq_ich, unq_mul in zip(unq_ich_lst, unq_mul_lst):
#                 if ich == unq_ich and mul == unq_mul:
#                     unique = False
#             if unique:
#                 unq_ich_lst.append(ich)
#                 unq_mul_lst.append(mul)
#                 unq_smi_lst.append(smi)
#
#                 formula_dct = automol.inchi.formula_dct(ich)
#                 lab = automol.formula.string(formula_dct)
#
#                 unq_lab_lst.append(lab)
#                 lab_idx = -1
#                 for lab_i in unq_lab_lst:
#                     if lab == lab_i:
#                         lab_idx += 1
#                 unq_lab_idx_lst.append(lab_idx)
#                 if lab_idx == 0:
#                     label = lab
#                 else:
#                     label = lab + '(' + str(lab_idx) + ')'
#                 csv_str += ','.join([label, smi, str(mul)])
#                 csv_str += '\n'
#                 spc_str += label
#                 spc_str += '\n'
#
#     spc_str += 'END'
#     spc_str += '\n'
#     spc_str += '\n'
#
#     sort_smiles_path = os.path.join(mech_path, 'smiles_sort.csv')
#     with open(sort_smiles_path, 'w') as sorted_csv_file:
#         sorted_csv_file.write(csv_str)
#
#     rxn_info_lst = list(zip(
#        formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst, rxn_sens,
#         rxn_unc, rxn_val, rxn_fam, rct_smis_lst, rct_ichs_lst, rct_muls_lst,
#         prd_smis_lst, prd_ichs_lst, prd_muls_lst))
#     rxn_info_lst = sorted(rxn_info_lst, key=lambda x: (x[0]))
#     old_formula = rxn_info_lst[0][0]
#     sens = rxn_info_lst[0][4]
#     ordered_formula = []
#     ordered_sens = []
#     for entry in rxn_info_lst:
#         formula = entry[0]
#         if formula == old_formula:
#             sens = max(sens, entry[4])
#         else:
#             ordered_sens.append(sens)
#             ordered_formula.append(old_formula)
#             sens = entry[4]
#             old_formula = formula
#     ordered_sens.append(sens)
#     ordered_formula.append(old_formula)
#     sens_dct = {}
#     for i, sens in enumerate(ordered_sens):
#         sens_dct[ordered_formula[i]] = sens
#     rxn_info_lst = sorted(
#         rxn_info_lst, key=lambda x: (sens_dct[x[0]], x[4]), reverse=True)
#
#     formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst, rxn_sens,
#     rxn_unc, rxn_val, rxn_fam, rct_smis_lst, rct_ichs_lst, rct_muls_lst,
#     prd_smis_lst, prd_ichs_lst, prd_muls_lst = zip(*rxn_info_lst)
#
#     rxn_namep_lst = []
#     rxn_namep_str = 'REACTIONS   KCAL/MOLE   MOLES'
#     rxn_namep_str += '\n'
#     for i, _ in enumerate(rxn_name_lst):
#         rxn_namep = []
#         rct_labs = []
#         rct_dat = zip(rct_smis_lst[i], rct_ichs_lst[i], rct_muls_lst[i])
#         spc_dat = zip(unq_ich_lst, unq_mul_lst, unq_lab_lst, unq_lab_idx_lst)
#         for _, rct_ich, rct_mul in rct_dat:
#             for ich, mul, lab, lab_idx in spc_dat:
#                 if rct_ich == ich and rct_mul == mul:
#                     if lab_idx == 0:
#                         rct_lab = lab
#                     else:
#                         rct_lab = lab + '(' + str(lab_idx) + ')'
#                     break
#             rct_labs.append(rct_lab)
#         rct_label = '+'.join(rct_labs)
#         prd_labs = []
#         prd_dat = zip(prd_smis_lst[i], prd_ichs_lst[i], prd_muls_lst[i])
#         spc_dat = zip(unq_ich_lst, unq_mul_lst, unq_lab_lst, unq_lab_idx_lst)
#         for _, prd_ich, prd_mul in prd_dat:
#             for ich, mul, lab, lab_idx in spc_dat:
#                 if prd_ich == ich and prd_mul == mul:
#                     if lab_idx == 0:
#                         prd_lab = lab
#                     else:
#                         prd_lab = lab + '(' + str(lab_idx) + ')'
#                     break
#             prd_labs.append(prd_lab)
#         prd_label = '+'.join(prd_labs)
#         rate_str = str('  1.e10   1.0   10000.  ! Sens = ')
#         rxn_namep = (
#             rct_label + ' <=> ' + prd_label + rate_str + str(rxn_sens[i]))
#         rxn_namep_str += rxn_namep
#         rxn_namep_str += '\n'
#         rxn_namep_lst.append(rxn_namep)
#
#     mech_str = spc_str + rxn_namep_str
#     mech_str += 'END'
#     mech_str += '\n'
#
#   with open(os.path.join(mech_path, 'mech_sort.txt'), 'w') as sort_mech_file:
#         sort_mech_file.write(mech_str)
#
#     # set up species info
#     spc_names = []
#     chg_dct = {}
#     mul_dct = {}
#     spc_dct = {}
#     tot_lst = rct_names_lst + prd_names_lst
#     for i, spc_names_lst in enumerate(tot_lst):
#         for j, spc_name in enumerate(spc_names_lst):
#             chg = 0
#             if spc_name not in spc_names:
#                 spc_names.append(spc_name)
#                 chg_dct[spc_name] = chg
#                 print('rct_muls test1:', spc_name, i, j)
#                 print('rct_muls test2:', rct_muls_lst[i])
#                 mul_dct[spc_name] = rct_muls_lst[i][j]
#                 spc_dct[spc_name] = {}
#                 spc_dct[spc_name]['chg'] = chg
#                 spc_dct[spc_name]['ich'] = rct_ichs_lst[i][j]
#                 spc_dct[spc_name]['mul'] = rct_muls_lst[i][j]
#     for i, spc_names_lst in enumerate(prd_names_lst):
#         for j, spc_name in enumerate(spc_names_lst):
#             chg = 0
#             if spc_name not in spc_names:
#                 spc_names.append(spc_name)
#                 chg_dct[spc_name] = chg
#                 mul_dct[spc_name] = prd_muls_lst[i][j]
#                 spc_dct[spc_name] = {}
#                 spc_dct[spc_name]['chg'] = chg
#                 spc_dct[spc_name]['ich'] = prd_ichs_lst[i][j]
#                 spc_dct[spc_name]['mul'] = prd_muls_lst[i][j]
#     rxn_info_lst = list(
#         zip(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst))
#
#   return spc_dct, rct_names_lst, prd_names_lst, rxn_name_lst, formula_str_lst