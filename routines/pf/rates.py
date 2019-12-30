"""
Write and Read MESS files for Rates
"""

import copy
import numpy
import automol
import mess_io
import ratefit
import autofile

# New libs
from lib.phydat import phycon
from routines.pf.messf import blocks


# Writer
def rate_headers(
        rct_ichs, temps, press, exp_factor, exp_power, exp_cutoff, eps1, eps2,
        sig1, sig2, mass1):
    """ makes the standard header and energy transfer sections for MESS input file
    """
    # header section
    header_str = mess_io.writer.global_reaction(temps, press)
    print(header_str)
    tot_mass = 0.
    for rct_ich in rct_ichs:
        geo = automol.convert.inchi.geometry(rct_ich)
        masses = automol.geom.masses(geo)
        for mass in masses:
            tot_mass += mass

    # energy transfer section
    energy_trans_str = mess_io.writer.energy_transfer(
        exp_factor, exp_power, exp_cutoff,
        eps1, eps2, sig1, sig2, mass1, tot_mass)

    return header_str, energy_trans_str


def make_all_species_data(rxn_lst, spc_dct, save_prefix, model_info,
                          pf_info, ts_found, projrot_script_str):
    """ generate the MESS species blocks for all the species
    """
    species = {}
    spc_save_fs = autofile.fs.species(save_prefix)
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        if tsname in ts_found:
            specieslist = rxn['reacs'] + rxn['prods']
            for name in specieslist:
                if name not in species:
                    species[name], _ = make_species_data(
                        name, spc_dct[name], spc_save_fs,
                        model_info, pf_info, projrot_script_str)
            if 'radical radical addition' not in spc_dct[tsname]['class']:
                ret1, ret2 = make_species_data(
                    tsname, spc_dct[tsname], save_prefix,
                    model_info, pf_info, projrot_script_str)
                species[tsname] = ret1
                spc_dct[tsname]['imag_freq'] = ret2
    return species


def make_species_data(spc, spc_dct_i, spc_save_fs, spc_model,
                      pf_levels, projrot_script_str):
    """ makes the main part of the MESS species block for a given species
    """
    spc_info = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    if 'ts_' in spc:
        save_path = spc_dct_i['rxn_fs'][3]
    else:
        spc_save_fs.leaf.create(spc_info)
        save_path = spc_save_fs.leaf.path(spc_info)
    species_data = blocks.species_block(
        spc=spc,
        spc_dct_i=spc_dct_i,
        spc_info=spc_info,
        spc_model=spc_model,
        pf_levels=pf_levels,
        projrot_script_str=projrot_script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix=save_path,
        )
    return species_data


def make_fake_species_data(spc_dct_i, spc_dct_j, spc_save_fs,
                           spc_model, pf_levels, projrot_script_str):
    """ make a fake MESS species block to represent the van der Waals well
    arising from the combination of two fragment species
    """
    spc_info_i = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    spc_info_j = (spc_dct_j['ich'], spc_dct_j['chg'], spc_dct_j['mul'])
    spc_save_fs.leaf.create(spc_info_i)
    spc_save_fs.leaf.create(spc_info_j)
    save_path_i = spc_save_fs.leaf.path(spc_info_i)
    save_path_j = spc_save_fs.leaf.path(spc_info_j)
    print('spc_model test:', spc_model)
    species_data = blocks.fake_species_block(
        spc_dct_i=spc_dct_i,
        spc_dct_j=spc_dct_j,
        spc_info_i=spc_info_i,
        spc_info_j=spc_info_j,
        spc_model=spc_model,
        pf_levels=pf_levels,
        projrot_script_str=projrot_script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix_i=save_path_i,
        save_prefix_j=save_path_j
        )
    return species_data


def make_channel_pfs(
        tsname, rxn, species_data, spc_dct, idx_dct, strs, first_ground_ene,
        spc_save_fs, spc_model, pf_levels, multi_info, projrot_script_str,
        pst_params=(1.0, 6),
        rad_rad_ts='pst'):
    """ make the partition function strings for each of the channels
    includes strings for each of the unimolecular wells, bimolecular fragments,
    and transition states connecting them.
    It also includes a special treatment for abstraction to include phase space
    blocks and coupling bimolecular fragments to fake van der Waals wells
    """
    bim_str, well_str, ts_str = strs

    # Find the number of uni and bimolecular wells already in the dictionary
    pidx = 1
    widx = 1
    fidx = 1
    for val in idx_dct.values():
        if 'P' in val:
            pidx += 1
        elif 'W' in val:
            widx += 1
        elif 'F' in val:
            fidx += 1

    # Set up new well for the reactants if that combo isn't already in the dct
    reac_label = ''
    reac_ene = 0.
    bimol = False
    if len(rxn['reacs']) > 1:
        bimol = True
    well_data = []
    spc_label = []
    for reac in rxn['reacs']:
        spc_label.append(automol.inchi.smiles(spc_dct[reac]['ich']))
        well_data.append(species_data[reac])
        reac_ene += (
            (spc_dct[reac]['ene'] + spc_dct[reac]['zpe']/phycon.EH2KCAL) *
            phycon.EH2KCAL)
    well_dct_key1 = '+'.join(rxn['reacs'])
    well_dct_key2 = '+'.join(rxn['reacs'][::-1])
    if well_dct_key1 not in idx_dct:
        if well_dct_key2 in idx_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                reac_label = 'P' + str(pidx)
                pidx += 1
                if not first_ground_ene:
                    first_ground_ene = reac_ene
                ground_energy = reac_ene - first_ground_ene
                bim_str += ' \t ! {} + {} \n'.format(
                    rxn['reacs'][0], rxn['reacs'][1])
                bim_str += mess_io.writer.bimolecular(
                    reac_label, spc_label[0], well_data[0], spc_label[1],
                    well_data[1], ground_energy)
                idx_dct[well_dct_key1] = reac_label
            else:
                if not first_ground_ene:
                    first_ground_ene = reac_ene
                reac_label = 'W' + str(widx)
                widx += 1
                zero_energy = reac_ene - first_ground_ene
                well_str += ' \t ! {} \n'.format(rxn['reacs'][0])
                well_str += mess_io.writer.well(
                    reac_label, well_data[0], zero_energy)
                idx_dct[well_dct_key1] = reac_label
    if not reac_label:
        reac_label = idx_dct[well_dct_key1]

    # Set up a new well for the products if that combo isn't already in the dct
    prod_label = ''
    prod_ene = 0.
    bimol = False
    if len(rxn['prods']) > 1:
        bimol = True
    well_data = []
    spc_label = []
    for prod in rxn['prods']:
        spc_label.append(automol.inchi.smiles(spc_dct[prod]['ich']))
        well_data.append(species_data[prod])
        prod_ene += (
            (spc_dct[prod]['ene'] + spc_dct[prod]['zpe']/phycon.EH2KCAL) *
            phycon.EH2KCAL)
    zero_energy = prod_ene - reac_ene
    well_dct_key1 = '+'.join(rxn['prods'])
    well_dct_key2 = '+'.join(rxn['prods'][::-1])
    if well_dct_key1 not in idx_dct:
        if well_dct_key2 in idx_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                prod_label = 'P' + str(pidx)
                ground_energy = prod_ene - first_ground_ene
                bim_str += ' \t ! {} + {} \n'.format(
                    rxn['prods'][0], rxn['prods'][1])
                bim_str += mess_io.writer.bimolecular(
                    prod_label, spc_label[0], well_data[0],
                    spc_label[1], well_data[1],
                    ground_energy)
                idx_dct[well_dct_key1] = prod_label
            else:
                prod_label = 'W' + str(widx)
                zero_energy = prod_ene - first_ground_ene
                well_str += ' \t ! {} \n'.format(rxn['prods'][0])
                well_str += mess_io.writer.well(
                    prod_label, well_data[0], zero_energy)
                idx_dct[well_dct_key1] = prod_label
    if not prod_label:
        prod_label = idx_dct[well_dct_key1]

    # For abstraction first make fake wells and PST TSs
    # Then for tight TS's make ts_str
    ts_ene = (
        (spc_dct[tsname]['ene'] + spc_dct[tsname]['zpe']/phycon.EH2KCAL) *
        phycon.EH2KCAL)
    zero_energy = ts_ene - first_ground_ene
    ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
    imag_freq = 0
    if 'imag_freq' in spc_dct[tsname]:
        imag_freq = abs(spc_dct[tsname]['imag_freq'])
    if not imag_freq:
        print('No imaginary freq for ts: {}'.format(tsname))

    fake_wellr_label = ''
    fake_wellp_label = ''
    print('class test:', spc_dct[tsname]['class'])
    rad_rad = 'radical radical' in spc_dct[tsname]['class']
    low_spin = 'high spin' not in spc_dct[tsname]['class']
    abst_rxn = 'abstraction' in spc_dct[tsname]['class']
    addn_rxn = 'addition' in spc_dct[tsname]['class']
    subs_rxn = 'substitution' in spc_dct[tsname]['class']
    if abst_rxn or subs_rxn:
        # Make fake wells and PST TSs as needed
        well_dct_key1 = 'F' + '+'.join(rxn['reacs'])
        well_dct_key2 = 'F' + '+'.join(rxn['reacs'][::-1])
        pst_r_ts_str = ''
        print('well_dct_key1 test:', well_dct_key1)
        print('well_dct_key2 test:', well_dct_key2)
        print('idx_dct:', idx_dct)
        if well_dct_key1 not in idx_dct:
            if well_dct_key2 in idx_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellr_label = 'F' + str(fidx)
                fidx += 1
                vdwr_ene = reac_ene - 1.0
                zero_energy = vdwr_ene - first_ground_ene
                well_str += ' \t ! Fake Well for {}\n'.format(
                    '+'.join(rxn['reacs']))
                fake_wellr = make_fake_species_data(
                    spc_dct[rxn['reacs'][0]], spc_dct[rxn['reacs'][1]],
                    spc_save_fs, spc_model, pf_levels, projrot_script_str)
                well_str += mess_io.writer.well(
                    fake_wellr_label, fake_wellr, zero_energy)
                idx_dct[well_dct_key1] = fake_wellr_label

                pst_r_label = 'FRB' + str(int(tsname.replace('ts_', ''))+1)
                idx_dct[well_dct_key1.replace('F', 'FRB')] = pst_r_label
                spc_dct_i = spc_dct[rxn['reacs'][0]]
                spc_dct_j = spc_dct[rxn['reacs'][1]]
                pst_r_ts_str = blocks.pst_block(
                    spc_dct_i, spc_dct_j, spc_model=spc_model,
                    pf_levels=pf_levels, projrot_script_str=projrot_script_str,
                    spc_save_fs=spc_save_fs,
                    pst_params=pst_params)
            print('fake_wellr_label test:', fake_wellr_label)
            if not fake_wellr_label:
                print('well_dct_key1 test:', well_dct_key1)
                fake_wellr_label = idx_dct[well_dct_key1]
                pst_r_label = idx_dct[well_dct_key1.replace('F', 'FRB')]
            zero_energy = reac_ene - first_ground_ene
            tunnel_str = ''
            print('ts_str input test:', pst_r_label, reac_label,
                  fake_wellr_label, pst_r_ts_str, zero_energy, tunnel_str)
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                pst_r_label, reac_label, fake_wellr_label, pst_r_ts_str,
                zero_energy, tunnel_str)
            print('ts_str test:', ts_str)
        else:
            fake_wellr_label = idx_dct[well_dct_key1]
        well_dct_key1 = 'F' + '+'.join(rxn['prods'])
        well_dct_key2 = 'F' + '+'.join(rxn['prods'][::-1])
        if well_dct_key1 not in idx_dct:
            if well_dct_key2 in idx_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellp_label = 'F' + str(fidx)
                fidx += 1
                vdwp_ene = prod_ene - 1.0
                zero_energy = vdwp_ene - first_ground_ene
                well_str += ' \t ! Fake Well for {}\n'.format(
                    '+'.join(rxn['prods']))
                fake_wellp = make_fake_species_data(
                    spc_dct[rxn['prods'][0]], spc_dct[rxn['prods'][1]],
                    spc_save_fs, spc_model, pf_levels, projrot_script_str)
                well_str += mess_io.writer.well(
                    fake_wellp_label, fake_wellp, zero_energy)
                idx_dct[well_dct_key1] = fake_wellp_label

                pst_p_label = 'FPB' + str(int(tsname.replace('ts_', ''))+1)
                idx_dct[well_dct_key1.replace('F', 'FPB')] = pst_p_label
                spc_dct_i = spc_dct[rxn['prods'][0]]
                spc_dct_j = spc_dct[rxn['prods'][1]]
                pst_p_ts_str = blocks.pst_block(
                    spc_dct_i, spc_dct_j, spc_model=spc_model,
                    pf_levels=pf_levels, projrot_script_str=projrot_script_str,
                    spc_save_fs=spc_save_fs,
                    pst_params=pst_params)
            if not fake_wellp_label:
                fake_wellp_label = idx_dct[well_dct_key1]
                pst_p_label = idx_dct[well_dct_key1.replace('F', 'FPB')]
            zero_energy = prod_ene - first_ground_ene
            tunnel_str = ''
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                pst_p_label, prod_label, fake_wellp_label, pst_p_ts_str,
                zero_energy, tunnel_str)
        else:
            fake_wellp_label = idx_dct[well_dct_key1]

        # Print inner TS data for radical radical call vtst or vrctst
        if rad_rad and low_spin and rad_rad_ts == 'pst':
            # zero_energy = SOMETHING
            pst_str = blocks.pst_block(
                spc_dct_i, spc_dct_j, spc_model=spc_model,
                pf_levels=pf_levels, projrot_script_str=projrot_script_str,
                spc_save_fs=spc_save_fs,
                pst_params=pst_params)
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                ts_label, reac_label, prod_label, pst_str, zero_energy)
        elif rad_rad and low_spin:
            ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
            if 'P' in reac_label:
                spc_ene = reac_ene - first_ground_ene
                spc_zpe = (
                    spc_dct[rxn['reacs'][0]]['zpe'] +
                    spc_dct[rxn['reacs'][1]]['zpe'])
            else:
                spc_ene = prod_ene - first_ground_ene
                spc_zpe = (
                    spc_dct[rxn['prods'][0]]['zpe'] +
                    spc_dct[rxn['prods'][1]]['zpe'])
            ts_str += '\n' + blocks.vtst_with_no_saddle_block(
                spc_dct[tsname], ts_label, fake_wellr_label, fake_wellp_label,
                spc_ene, spc_zpe, projrot_script_str,
                multi_info, elec_levels=[[0., 1]], sym_factor=1.
                )
        else:
            vdwr_ene = reac_ene - 1.0
            vdwp_ene = prod_ene - 1.0
            zero_energy = ts_ene - first_ground_ene
            ts_reac_barr = ts_ene - vdwr_ene
            ts_prod_barr = ts_ene - vdwp_ene
            if ts_reac_barr < 0.:
                ts_reac_barr = 0.1
            if ts_prod_barr < 0.:
                ts_prod_barr = 0.1
            tunnel_str = mess_io.writer.tunnel_eckart(
                imag_freq, ts_reac_barr, ts_prod_barr)
            ts_str += '\n' + mess_io.writer.ts_sadpt(
                ts_label, fake_wellr_label, fake_wellp_label,
                species_data[tsname], zero_energy, tunnel_str)
    elif rad_rad and addn_rxn and low_spin:
        # For radical radical addition call vtst or vrctst
        ts_label = 'B' + str(int(tsname.replace('ts_', ''))+1)
        if 'P' in reac_label:
            spc_ene = reac_ene - first_ground_ene
            spc_zpe = (
                spc_dct[rxn['reacs'][0]]['zpe'] +
                spc_dct[rxn['reacs'][1]]['zpe'])
        else:
            spc_ene = prod_ene - first_ground_ene
            spc_zpe = (
                spc_dct[rxn['prods'][0]]['zpe'] +
                spc_dct[rxn['prods'][1]]['zpe'])
        ts_str += '\n' + blocks.vtst_with_no_saddle_block(
            spc_dct[tsname], ts_label, reac_label, prod_label,
            spc_ene, spc_zpe, projrot_script_str,
            multi_info, elec_levels=[[0., 1]], sym_factor=1.)
    elif rad_rad and addn_rxn and low_spin and rad_rad_ts == 'pst':
        # zero_energy = SOMETHING
        pst_str = blocks.pst_block(
            spc_dct_i, spc_dct_j, spc_model=spc_model,
            pf_levels=pf_levels, projrot_script_str=projrot_script_str,
            spc_save_fs=spc_save_fs,
            pst_params=pst_params)
        ts_str += '\n' + mess_io.writer.ts_sadpt(
            ts_label, reac_label, prod_label, pst_str, zero_energy)
    else:
        ts_reac_barr = ts_ene - reac_ene
        ts_prod_barr = ts_ene - prod_ene
        if ts_reac_barr < 0.:
            ts_reac_barr = 0.1
        if ts_prod_barr < 0.:
            ts_prod_barr = 0.1
        tunnel_str = mess_io.writer.tunnel_eckart(
            imag_freq, ts_reac_barr, ts_prod_barr)
        ts_str += '\n' + mess_io.writer.ts_sadpt(
            ts_label, reac_label, prod_label,
            species_data[tsname], zero_energy, tunnel_str)

    return [well_str, bim_str, ts_str], first_ground_ene


# Readers
def read_rates(rct_lab, prd_lab, mess_path, assess_pdep_temps,
               pdep_tolerance=20.0, no_pdep_pval=1.0,
               pdep_low=None, pdep_high=None, bimol=False):
    """ Read the rate constants from the MESS output and
        (1) filter out the invalid rates that are negative or undefined
        and obtain the pressure dependent values
    """

    # Dictionaries to store info; indexed by pressure (given in fit_ps)
    calc_k_dct = {}
    valid_calc_tk_dct = {}
    ktp_dct = {}

    # Read the MESS output file into a string
    with open(mess_path+'/rate.out', 'r') as mess_file:
        output_string = mess_file.read()
    # with open(mess_path+'/mess.inp', 'r') as mess_file:
    #     input_string = mess_file.read()

    # Read the temperatures and pressures out of the MESS output
    mess_temps, _ = mess_io.reader.rates.get_temperatures(
        output_string)
    mess_pressures, punit = mess_io.reader.rates.get_pressures(
        output_string)
    # mess_pressures, punit = mess_io.reader.rates.get_pressures_input(
    #     input_string)

    # Loop over the pressures obtained from the MESS output
    for pressure in mess_pressures:

        # Read the rate constants
        if pressure == 'high':
            rate_ks = mess_io.reader.highp_ks(
                output_string, rct_lab, prd_lab)
        else:
            rate_ks = mess_io.reader.pdep_ks(
                output_string, rct_lab, prd_lab, pressure, punit)

        # Store in a dictionary
        calc_k_dct[pressure] = rate_ks

    # Remove k(T) vals at each P where where k is negative or undefined
    # If ANY valid k(T,P) vals at given pressure, store in dct
    for pressure, calc_ks in calc_k_dct.items():
        filtered_temps, filtered_ks = ratefit.fit.util.get_valid_tk(
            mess_temps, calc_ks, bimol)
        if filtered_ks.size > 0:
            valid_calc_tk_dct[pressure] = numpy.concatenate(
                (filtered_temps, filtered_ks))

    # Filter the ktp dictionary by assessing the presure dependence
    if list(valid_calc_tk_dct.keys()) == ['high']:
        ktp_dct['high'] = valid_calc_tk_dct['high']
    else:
        rxn_is_pdependent = ratefit.err.assess_pressure_dependence(
            valid_calc_tk_dct, assess_pdep_temps,
            tolerance=pdep_tolerance, plow=pdep_low, phigh=pdep_high)
        if rxn_is_pdependent:
            # Set dct to fit as copy of dct to do PLOG fits at all pressures
            ktp_dct = copy.deepcopy(valid_calc_tk_dct)
        else:
            # Set dct to have single set of k(T, P) vals: P is desired pressure
            if no_pdep_pval in valid_calc_tk_dct:
                ktp_dct['high'] = valid_calc_tk_dct[no_pdep_pval]
    if 'high' not in ktp_dct and 'high' in valid_calc_tk_dct.keys():
        ktp_dct['high'] = valid_calc_tk_dct['high']

    return ktp_dct
