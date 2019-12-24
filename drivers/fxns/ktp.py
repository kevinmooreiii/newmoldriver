""" driver for running the ktp
"""

# new libs
import lib.filesystem.inf.get_thy_info


def build_spc_list():
    """ set species list
    """

    spc_queue = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[rxn])
        rxn_spc.extend(list(prd_names_lst[rxn]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)

    spc_tsk_lst = []
    ts_tsk_lst = []
    ts_tsk = False
    for tsk in tsk_info_lst:
        if 'find_ts' in tsk[0]:
            ts_tsk = True
        if ts_tsk:
            ts_tsk_lst.append(tsk)
        else:
            spc_tsk_lst.append(tsk)

    return spc_tsk_lst, ts_tsk_lst


def init_es_run():
    """ First run ESDriver for the species on the PES so that
        exothermicity info is available to sort reactions by exothermicity
    """

    runspecies = [{'species': spc_queue, 'reacs': [], 'prods': []}]
    esdriver.driver.run(
        spc_tsk_lst, es_dct, runspecies, spc_dct,
        run_prefix, save_prefix, vdw_params,
        pst_params=pst_params,
        rad_rad_ts=rad_rad_ts)


def ts_tsks():
    """ set ts task lists
    """

    # Form the reaction list
    rxn_lst = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_lst.append(
            {'species': [],
             'reacs': list(rct_names_lst[rxn]),
             'prods': list(prd_names_lst[rxn])})

    # Add addtional dictionary items for all the TSs
    # This presumes that es has been run previously for species list
    # to produce energies in save file system
    if ts_tsk_lst:
        print('\nBegin transition state prep')
        for ts in spc_dct:
            if 'ts_' in ts:
                # Figure out how to pass ts_info or have it recalculated in esdriver
                # Exothermicity reordering requires elec energy 
                # which requires theory level and tsk_info
                es_ini_key = ts_tsk_lst[0][2]
                es_run_key = ts_tsk_lst[0][1]
                ini_thy_info = esdriver.driver.get_es_info(es_dct, es_ini_key)
                thy_info = esdriver.driver.get_es_info(es_dct, es_run_key)

                # Generate rxn data, reorder if necessary, put in spc_dct for ts
                rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul = scripts.es.rxn_info(
                    run_prefix, save_prefix, ts, spc_dct, thy_info, ini_thy_info)
                spc_dct[ts]['rxn_ichs'] = rxn_ichs
                spc_dct[ts]['rxn_chgs'] = rxn_chgs
                spc_dct[ts]['rxn_muls'] = rxn_muls
                spc_dct[ts]['low_mul'] = low_mul
                spc_dct[ts]['high_mul'] = high_mul

                # Generate rxn_fs from rxn_info stored in spc_dct
                rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path = scripts.es.get_rxn_fs(
                    run_prefix, save_prefix, spc_dct[ts])
                spc_dct[ts]['rxn_fs'] = [rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path]

                rct_zmas, prd_zmas, rct_cnf_save_fs, prd_cnf_save_fs = scripts.es.get_zmas(
                    spc_dct[ts]['reacs'], spc_dct[ts]['prods'], spc_dct,
                    ini_thy_info, save_prefix, run_prefix, KICKOFF_SIZE,
                    KICKOFF_BACKWARD, substr.PROJROT)
                ret = scripts.es.ts_class(
                    rct_zmas, prd_zmas, spc_dct[ts]['rad_rad'],
                    spc_dct[ts]['mul'], low_mul, high_mul,
                    rct_cnf_save_fs, prd_cnf_save_fs, spc_dct[ts]['given_class'])
                ret1, ret2 = ret

                # Do something with spc_dct
                if ret1:
                    rxn_class, ts_zma, dist_name, brk_name, grid, frm_bnd_key, brk_bnd_key, tors_names, update_guess = ret1
                    spc_dct[ts]['class'] = rxn_class
                    spc_dct[ts]['grid'] = grid
                    spc_dct[ts]['tors_names'] = tors_names
                    spc_dct[ts]['original_zma'] = ts_zma
                    dist_info = [dist_name, 0., update_guess, brk_name]
                    spc_dct[ts]['dist_info'] = dist_info
                    spc_dct[ts]['frm_bnd_key'] = frm_bnd_key
                    spc_dct[ts]['brk_bnd_key'] = brk_bnd_key

                    # Adding in the rct and prd zmas for vrctst
                    spc_dct[ts]['rct_zmas'] = rct_zmas
                    spc_dct[ts]['prd_zmas'] = prd_zmas

                    # Add the backup data to the spc_dct
                    if ret2:
                        spc_dct[ts]['bkp_data'] = ret2
                    else:
                        spc_dct[ts]['bkp_data'] = None
                else:
                    spc_dct[ts]['class'] = None
                    spc_dct[ts]['bkp_data'] = None

        print('End transition state prep\n')

        #Run ESDriver to get the transition state?
        if runes:
            ts_found = esdriver.driver.run(
                ts_tsk_lst, es_dct, rxn_lst, spc_dct,
                run_prefix, save_prefix, vdw_params)


def determine_mess_models():
    """ Figure out the model and theory levels for the MESS files
    """
    for ts in spc_dct:
        if 'original_zma' in spc_dct[ts]:
            pes_formula = automol.geom.formula(
                automol.zmatrix.geometry(spc_dct[ts]['original_zma']))
            print('Starting mess file preparation for {}:'.format(pes_formula))
            break

    # Figure out the model and theory levels for the MESS files
    geo_lvl = ''
    harm_lvl = ''
    anharm_lvl = ''
    tors_lvl = ''
    sym_lvl = ''
    geo_lvl_ref = ''
    harm_lvl_ref = ''
    anharm_lvl_ref = ''
    tors_lvl_ref = ''
    sym_lvl_ref = ''

    ts_model = ['RIGID', 'HARM', '']
    for tsk in ts_tsk_lst:
        if 'samp' in tsk[0] or 'find' in tsk[0]:
            geo_lvl = tsk[1]
            geom = True
            if 'find' in tsk[0]:
                geo_lvl_ref = geo_lvl
        if 'grad' in tsk[0] or 'hess' in tsk[0]:
            harm_lvl = tsk[1]
            harm_lvl_ref = tsk[2]
            if 'grad' in tsk[0]:
                grad = True
            if 'hess' in tsk[0]:
                hess = True
            if not geom:
                ene_lvl = tsk[1]
                geo_lvl = tsk[1]
        if 'hr' in tsk[0] or 'tau' in tsk[0]:
            # print('found')
            tors_lvl = tsk[1]
            tors_lvl_ref = tsk[2]
            if 'md' in tsk[0]:
                ts_model[0] = 'MDHR'
            if 'tau' in tsk[0]:
                ts_model[0] = 'TAU'
            else:
                ts_model[0] = '1DHR'
        if 'anharm' in tsk[0] or 'vpt2' in tsk[0]:
            anharm_lvl = tsk[1]
            anharm_lvl_ref = tsk[2]
            ts_model[1] = 'ANHARM'
            if not hess:
                geo_lvl = tsk[1]
        if 'sym' in tsk[0]:
            sym_lvl = tsk[1]
            sym_lvl_ref = tsk[2]
            if 'samp' in tsk[0]:
                ts_model[2] = 'SAMPLING'
            if '1DHR' in tsk[0]:
                ts_model[2] = '1DHR'


def determine_model_thy_levels():
    """ Figure out the model and theory levels for the MESS files
    """
    def get_thy_info(es_dct, key):
        """ setup theory info file from es dictionary
        """
        ret = []
        if key:
            ret = scripts.es.get_thy_info(es_dct[key])
        return ret

    geo_thy_info = get_thy_info(es_dct, geo_lvl)
    geo_thy_info_ref = get_thy_info(es_dct, geo_lvl_ref)
    harm_thy_info = get_thy_info(es_dct, harm_lvl)
    tors_thy_info = None
    anharm_thy_info = None
    sym_thy_info = None
    harm_ref_thy_info = None
    tors_ref_thy_info = None
    anharm_ref_thy_info = None
    sym_ref_thy_info = None
    if tors_lvl:
        tors_thy_info = get_thy_info(es_dct, tors_lvl)
    if anharm_lvl:
        anharm_thy_info = get_thy_info(es_dct, anharm_lvl)
    if sym_lvl:
        sym_thy_info = get_thy_info(es_dct, sym_lvl)
    if harm_lvl_ref:
        harm_ref_thy_info = get_thy_info(es_dct, harm_lvl_ref)
    if tors_lvl_ref:
        tors_ref_thy_info = get_thy_info(es_dct, tors_lvl_ref)
    if anharm_lvl_ref:
        anharm_ref_thy_info = get_thy_info(es_dct, anharm_lvl_ref)
    if sym_lvl_ref:
        sym_ref_thy_info = get_thy_info(es_dct, sym_lvl_ref)
    pf_levels = [harm_thy_info, tors_thy_info, anharm_thy_info, sym_thy_info]
    multi_levels = []
    ref_levels = [
        harm_ref_thy_info, tors_ref_thy_info, anharm_ref_thy_info, sym_ref_thy_info]

    # Collect ground energies and zero-point energies
    spc_save_fs = autofile.fs.species(save_prefix)
    ts_queue = []
    for spc in spc_dct:
        if spc in ts_found:
            ts_queue.append(spc)
    print('getting ready for zpe:')
    for spc in spc_queue + ts_queue:
        spc_info = (spc_dct[spc]['ich'], spc_dct[spc]['chg'], spc_dct[spc]['mul'])
        if 'ts_' in spc:
            spc_save_path = spc_dct[spc]['rxn_fs'][3]
            saddle = True
            save_path = spc_save_path
        else:
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            saddle = False
            save_path = save_prefix
        zpe, _ = scripts.thermo.get_zpe(
            spc, spc_dct[spc], spc_save_path, pf_levels, ts_model)
        spc_dct[spc]['zpe'] = zpe
        ene_strl = []
        ene_lvl = ''
        ene_lvl_ref = ''
        ene_idx = 0
        spc_dct[spc]['ene'] = 0.
        ene_str = '! energy level:'
        # print('looking at ts tasks')
        for tsk in ts_tsk_lst:
            if 'ene' in tsk[0]:
                if ene_idx > len(ene_coeff)-1:
                    print('Warning - an insufficient energy coefficient list was provided')
                    break
                ene_lvl = tsk[1]
                ene_lvl_ref = tsk[2]
                ene_ref_thy_info = scripts.es.get_thy_info(es_dct[ene_lvl_ref])
                ene_thy_info = scripts.es.get_thy_info(es_dct[ene_lvl])
                ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
                    ene_coeff[ene_idx], ene_thy_info[3], ene_thy_info[1], ene_thy_info[2],
                    ene_ref_thy_info[3], ene_ref_thy_info[1], ene_ref_thy_info[2]))
                ene = scripts.thermo.get_electronic_energy(
                    spc_info, ene_ref_thy_info, ene_thy_info, save_path, saddle)
                # print('ene test:', ene_idx, ene_coeff[ene_idx], ene)
                spc_dct[spc]['ene'] += ene*ene_coeff[ene_idx]
                ene_idx += 1
    ene_str += '!               '.join(ene_strl)


def write_mess_reaction_files()
    """ Make the channel strings for MESS
    """

    # Set multireference info lists for temporary writing
    multi_info = ['molpro2015', 'caspt2', 'cc-pVDZ', 'RR']
    # multi_info = ['molpro2015', 'caspt2', 'cc-pVTZ', 'RR']

    # Collect formula and header string for the PES
    tsname_0 = 'ts_0'
    rct_ichs = spc_dct[tsname_0]['rxn_ichs'][0]

    # Write the MESS global keyword section and energy trans section
    header_str, energy_trans_str = scripts.ktp.pf_headers(
        rct_ichs, TEMPS, PRESS, *etrans)

    # Get all of the data for MESS species blocks
    species = scripts.ktp.make_all_species_data(
        rxn_lst, spc_dct, save_prefix, ts_model, pf_levels, ts_found, substr.PROJROT)

    # Write the MESS species blocks
    mess_strs = ['', '', '']
    idx_dct = {}
    first_ground_ene = 0.
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        if tsname in ts_found:
            tsform = automol.geom.formula(
                automol.zmatrix.geometry(spc_dct[tsname]['original_zma']))
            if tsform != pes_formula:
                print('Reaction list contains reactions on different potential energy',
                      'surfaces: {} and {}'.format(tsform, pes_formula))
                print('Will proceed to construct only {}'.format(pes_formula))
                continue
            mess_strs, first_ground_ene = scripts.ktp.make_channel_pfs(
                tsname, rxn, species, spc_dct, idx_dct, mess_strs,
                first_ground_ene, spc_save_fs, ts_model, pf_levels,
                multi_info, substr.PROJROT,
                pst_params=pst_params)

    # Collate the MESS strings
    well_str, bim_str, ts_str = mess_strs
    ts_str += '\nEnd\n'
    print(well_str)
    print(bim_str)
    print(ts_str)


def run_mess_rates():
    """ call lib to run the rates
    """
    mess_path = scripts.ktp.run_rates(
        header_str, energy_trans_str, well_str, bim_str, ts_str,
        spc_dct[tsname_0], geo_thy_info_ref, spc_dct[tsname_0]['rxn_fs'][3])

# fit rate output to modified Arrhenius forms and print in ChemKin format
# lib.fit.fit_rates()
