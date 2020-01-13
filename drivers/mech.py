""" estoktp driver functions with no
better home
"""

import automol.geom
from drivers import esdriver
from drivers import ktpdriver
from routines.pf import get_high_level_energy
from routines.pf import rates as messrates
from lib.filesystem import inf as finf
from lib.filesystem import path as fpath
from lib.filesystem import read as fread
from lib.reaction import rxnid


def run_driver(pes_dct, conn_chnls_dct,
               spc_dct, cla_dct,
               thy_dct,
               model_dct,
               es_tsk_lst,
               run_jobs_lst,
               run_inp_dct,
               run_options_dct,
               driver='es_spc'):
    """ Run the ktp driver for the PESs
    """
    for formula in pes_dct:

        # Build the names list
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        # Select names from the names list corresponding to chnls to run
        # conn_chnls_dct[formula] = {sub_pes_idx: [channel_idxs]}
        for cvals in conn_chnls_dct[formula].values():
            run_pes = False
            rct_names_lst = []
            prd_names_lst = []
            rxn_name_lst = []
            for chn_idx, _ in enumerate(pes_rxn_name_lst):
                if chn_idx in cvals:
                    run_pes = True
                    rct_names_lst.append(pes_rct_names_lst[chn_idx])
                    prd_names_lst.append(pes_prd_names_lst[chn_idx])
                    rxn_name_lst.append(pes_rxn_name_lst[chn_idx])
                    print('running channel {}: {} = {}'.format(
                        str(chn_idx+1),
                        ' + '.join(pes_rct_names_lst[chn_idx]),
                        ' + '.join(pes_prd_names_lst[chn_idx])))
            if run_pes:
                rxn_lst = []
                for rxn, _ in enumerate(rct_names_lst):
                    rxn_lst.append(
                        {'species': [],
                         'reacs': list(rct_names_lst[rxn]),
                         'prods': list(prd_names_lst[rxn])})

                # Add the stationary points to the spc dcts
                ts_dct = build_spc_dct_for_sadpts(
                    spc_dct, rxn_lst, rxn_name_lst,
                    rct_names_lst, prd_names_lst, cla_dct)
                spc_dct.update(ts_dct)

                # Run the appropriate driver
                if driver == 'es_spc':
                    spc_esdriver(
                        rct_names_lst, prd_names_lst, es_tsk_lst,
                        spc_dct, thy_dct,
                        model_dct,
                        run_options_dct,
                        run_inp_dct)
                if driver == 'es_rxn':
                    rxn_esdriver(
                        rct_names_lst, prd_names_lst, es_tsk_lst,
                        spc_dct, thy_dct,
                        model_dct,
                        run_options_dct,
                        run_inp_dct)
                elif driver == 'ktp':
                    rate_jobs = [bool('rates' in run_jobs_lst),
                                 bool('params' in run_jobs_lst)]
                    ktpdriver.run(
                        spc_dct,
                        thy_dct,
                        es_tsk_lst,
                        pes_rct_names_lst,
                        pes_prd_names_lst,
                        model_dct,
                        run_inp_dct,
                        rate_jobs=rate_jobs)


def spc_esdriver(rct_names_lst, prd_names_lst, es_tsk_lst,
                 spc_dct, thy_dct,
                 model_dct,
                 run_options_dct,
                 run_inp_dct):
    """ Call the ESDriver Routines for species in the spc dct
    """

    # Get the reactant and product speceies in a list to be run
    spc_queue = form_spc_queue(rct_names_lst, prd_names_lst)
    spc_run_lst = format_run_spc_lst(spc_queue)

    # Format the task info list
    spc_tsk_lst, _ = format_tsk_lst(es_tsk_lst)

    print('in spc_es_driver')
    print(spc_tsk_lst)
    # Execute ESDriver
    esdriver.run(spc_run_lst, spc_dct,
                 spc_tsk_lst,
                 thy_dct, model_dct,
                 run_options_dct, run_inp_dct)


def rxn_esdriver(rct_names_lst, prd_names_lst, es_tsk_lst,
                 spc_dct, thy_dct,
                 model_dct,
                 run_options_dct,
                 run_inp_dct):
    """ Call the ESDriver Routines for reactions in the spc dct
    """
    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']
    kickoff = run_options_dct['kickoff']

    # Format the task info list
    _, ts_tsk_lst = format_tsk_lst(es_tsk_lst)

    # Form the reaction list
    rxn_lst = format_run_rxn_lst(rct_names_lst, prd_names_lst)

    # Add addtional dictionary items for all the TSs
    # This presumes that es has been run previously for species list
    # to produce energies in save file system
    if ts_tsk_lst:
        print('\nBegin transition state prep')
        for spc in spc_dct:
            if 'ts_' in spc:
                spc_dct[spc] = set_sadpt_info(
                    ts_tsk_lst, spc_dct, spc, thy_dct,
                    run_prefix, save_prefix,
                    kickoff)
        print('End transition state prep\n')

        # Execute ESDriver
        esdriver.run(rxn_lst, spc_dct,
                     ts_tsk_lst,
                     thy_dct, model_dct,
                     run_options_dct, run_inp_dct)


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


def form_spc_queue(rct_names_lst=(), prd_names_lst=()):
    """ form the species queue from tht elist
    """
    spc_queue = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[rxn])
        rxn_spc.extend(list(prd_names_lst[rxn]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)

    return spc_queue


def form_spc_queue_2(rxn_lst):
    """ second spc queue; redundant to one above
    """
    spc_queue = []
    for _, rxn in enumerate(rxn_lst):
        reacs = rxn['reacs']
        prods = rxn['prods']
        spc_queue.extend(rxn['species'])
        spc_queue.extend(reacs)
        spc_queue.extend(prods)
    spc_queue = list(dict.fromkeys(spc_queue))
    return spc_queue


def format_tsk_lst(tsk_info_lst):
    """ Format the input task list appropriate for spc or rxns
    """
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


def format_run_spc_lst(spc_queue):
    """ get the list of species to run
    """
    return [{'species': spc_queue, 'reacs': [], 'prods': []}]


def format_run_rxn_lst(rct_names_lst, prd_names_lst):
    """ Get the lst of reactions to be run
    """
    run_lst = []
    for rxn, _ in enumerate(rct_names_lst):
        run_lst.append(
            {'species': [], 'reacs': list(rct_names_lst[rxn]), 'prods':
             list(prd_names_lst[rxn])})

    return run_lst


def set_sadpt_info(ts_tsk_lst, spc_dct, spc, thy_dct,
                   run_prefix, save_prefix, kickoff):
    """ set the saddle point dct with info
    """
    es_ini_key = ts_tsk_lst[0][2]
    es_run_key = ts_tsk_lst[0][1]
    ini_thy_info = finf.get_es_info(es_ini_key, thy_dct)
    thy_info = finf.get_es_info(es_run_key, thy_dct)

    # Generate rxn data, reorder if necessary, and put in spc_dct for given ts
    rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul = finf.rxn_info(
        save_prefix, spc, spc_dct, thy_info, ini_thy_info)
    spc_dct[spc]['rxn_ichs'] = rxn_ichs
    spc_dct[spc]['rxn_chgs'] = rxn_chgs
    spc_dct[spc]['rxn_muls'] = rxn_muls
    spc_dct[spc]['low_mul'] = low_mul
    spc_dct[spc]['high_mul'] = high_mul

    # Generate rxn_fs from rxn_info stored in spc_dct
    [kickoff_size, kickoff_backward] = kickoff
    rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path = fpath.get_rxn_fs(
        run_prefix, save_prefix, spc_dct[spc])
    spc_dct[spc]['rxn_fs'] = [
        rxn_run_fs,
        rxn_save_fs,
        rxn_run_path,
        rxn_save_path]
    rct_zmas, prd_zmas, rct_cnf_save_fs, prd_cnf_save_fs = fread.get_zmas(
        spc_dct[spc]['reacs'], spc_dct[spc]['prods'], spc_dct,
        ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward)
    ret = rxnid.ts_class(
        rct_zmas, prd_zmas, spc_dct[spc]['rad_rad'],
        spc_dct[spc]['mul'], low_mul, high_mul,
        rct_cnf_save_fs, prd_cnf_save_fs, spc_dct[spc]['given_class'])
    ret1, ret2 = ret
    if ret1:
        [rxn_class, spc_zma,
         dist_name, brk_name, grid,
         frm_bnd_key, brk_bnd_key,
         tors_names, update_guess] = ret1
        spc_dct[spc]['class'] = rxn_class
        spc_dct[spc]['grid'] = grid
        spc_dct[spc]['tors_names'] = tors_names
        spc_dct[spc]['original_zma'] = spc_zma
        dist_info = [dist_name, 0., update_guess, brk_name]
        spc_dct[spc]['dist_info'] = dist_info
        spc_dct[spc]['frm_bnd_key'] = frm_bnd_key
        spc_dct[spc]['brk_bnd_key'] = brk_bnd_key
        # Adding in the rct and prd zmas for vrctst
        spc_dct[spc]['rct_zmas'] = rct_zmas
        spc_dct[spc]['prd_zmas'] = prd_zmas
        if ret2:
            spc_dct[spc]['bkp_data'] = ret2
        else:
            spc_dct[spc]['bkp_data'] = None
    else:
        spc_dct[spc]['class'] = None
        spc_dct[spc]['bkp_data'] = None

    return spc_dct[spc]


def set_pf_model_info(pf_model):
    """ Set the PF model list based on the input
    """
    tors_model = pf_model['tors'] if pf_model['tors'] else 'RIGID'
    vib_model = pf_model['vib'] if pf_model['vib'] else 'HARM'
    sym_model = pf_model['sym'] if pf_model['sym'] else ''

    pf_models = [tors_model, vib_model, sym_model]
    return pf_models


def set_es_model_info(es_model, thy_dct):
    """ Set the model info
    """
    # Read the ES models from the model dictionary
    geo_lvl = es_model['geo'] if es_model['geo'] else None
    ene_lvl = es_model['ene'] if es_model['ene'] else None
    harm_lvl = es_model['harm'] if es_model['harm'] else None
    anharm_lvl = es_model['anharm'] if es_model['anharm'] else None
    sym_lvl = es_model['sym'] if es_model['sym'] else None

    # Torsional Scan which needs a reference for itself
    tors_lvl_sp = es_model['tors'][0] if es_model['tors'] else None
    tors_lvl_scn = es_model['tors'][1] if es_model['tors'] else None

    # Set the theory info objects
    geo_thy_info = finf.get_thy_info(geo_lvl, thy_dct)
    ene_thy_info = finf.get_thy_info(ene_lvl, thy_dct)
    harm_thy_info = finf.get_thy_info(harm_lvl, thy_dct)
    anharm_thy_info = (finf.get_thy_info(anharm_lvl, thy_dct)
                       if anharm_lvl else None)
    sym_thy_info = (finf.get_thy_info(sym_lvl, thy_dct)
                    if sym_lvl else None)
    tors_sp_thy_info = (finf.get_thy_info(tors_lvl_sp, thy_dct)
                        if tors_lvl_sp else None)
    tors_scn_thy_info = (finf.get_thy_info(tors_lvl_scn, thy_dct)
                         if tors_lvl_scn else None)

    # Combine levels into a list
    es_levels = [
        geo_thy_info,
        ene_thy_info,
        harm_thy_info,
        anharm_thy_info,
        sym_thy_info,
        [tors_sp_thy_info, tors_scn_thy_info]
    ]

    return es_levels

    # # Read the ES models from the dictionary
    # # geo_lvl = es_model['geo'][0] if es_model['geo'] else None
    # geo_lvl_ref = es_model['geo'][1] if es_model['geo'] else None
    # harm_lvl = es_model['harm'][0] if es_model['harm'] else None
    # harm_lvl_ref = es_model['harm'][1] if es_model['harm'] else None
    # anharm_lvl = es_model['anharm'][0] if es_model['anharm'] else None
    # anharm_lvl_ref = es_model['anharm'][1] if es_model['anharm'] else None
    # tors_lvl = es_model['tors'][0] if es_model['tors'] else None
    # tors_lvl_ref = es_model['tors'][1] if es_model['tors'] else None
    # sym_lvl = es_model['sym'][0] if es_model['sym'] else None
    # sym_lvl_ref = es_model['sym'][1] if es_model['sym'] else None

    # # Set the theory info objects
    # harm_thy_info = finf.get_thy_info(harm_lvl, thy_dct)
    # geo_ref_thy_info = finf.get_thy_info(geo_lvl_ref, thy_dct)
    # tors_thy_info = (finf.get_thy_info(tors_lvl, thy_dct)
    #                  if tors_lvl else None)
    # anharm_thy_info = (finf.get_thy_info(anharm_lvl, thy_dct)
    #                    if anharm_lvl else None)
    # sym_thy_info = (finf.get_thy_info(sym_lvl, thy_dct)
    #                 if sym_lvl else None)
    # harm_ref_thy_info = (finf.get_thy_info(harm_lvl_ref, thy_dct)
    #                      if harm_lvl_ref else None)
    # tors_ref_thy_info = (finf.get_thy_info(tors_lvl_ref, thy_dct)
    #                      if tors_lvl_ref else None)
    # anharm_ref_thy_info = (finf.get_thy_info(anharm_lvl_ref, thy_dct)
    #                        if anharm_lvl_ref else None)
    # sym_ref_thy_info = (finf.get_thy_info(sym_lvl_ref, thy_dct)
    #                     if sym_lvl_ref else None)

    # # Combine levels into a list
    # pf_levels = [
    #     harm_thy_info, tors_thy_info,
    #     anharm_thy_info, sym_thy_info]
    # ref_levels = [
    #     harm_ref_thy_info, tors_ref_thy_info,
    #     anharm_ref_thy_info, sym_ref_thy_info,
    #     geo_ref_thy_info]

    # return pf_levels, ref_levels


def set_pes_formula(spc_dct):
    """ Set pes formula using zma
    """
    for spc_2 in spc_dct:
        if 'original_zma' in spc_dct[spc_2]:
            pes_formula = automol.geom.formula(
                automol.zmatrix.geometry(spc_dct[spc_2]['original_zma']))
            print('Starting mess file preparation for {}:'.format(
                pes_formula))
            break
    return pes_formula


def get_high_energy(ts_tsk_lst, thy_dct, spc_info,
                    save_path, saddle, ene_coeff):
    """ calc high energy
    """
    spc_ene = 0.0
    ene_idx = 0
    for tsk in ts_tsk_lst:
        if 'ene' in tsk[0]:
            if ene_idx > len(ene_coeff)-1:
                print('Warning - an insufficient ',
                      'energy coefficient list was provided')
                break
            ene_lvl = tsk[1]
            ene_lvl_ref = tsk[2]
            ene_ref_thy_info = finf.get_thy_info(ene_lvl_ref, thy_dct)
            ene_thy_info = finf.get_thy_info(ene_lvl, thy_dct)
            ene = get_high_level_energy(
                spc_info=spc_info,
                thy_low_level=ene_ref_thy_info,
                thy_high_level=ene_thy_info,
                save_prefix=save_path,
                saddle=saddle)
            spc_ene += ene*ene_coeff[ene_idx]
            ene_idx += 1

    return spc_ene


def get_ckin_ene_lvl_str(ts_tsk_lst, thy_dct, ene_coeff):
    """ Write the comment lines for the enrgy lvls for ckin
    """
    ene_strl = []
    ene_idx = 0
    ene_str = '! energy level:'
    for tsk in ts_tsk_lst:
        if 'ene' in tsk[0]:
            if ene_idx > len(ene_coeff)-1:
                print('Warning - an insufficient ',
                      'energy coefficient list was provided')
                break
            ene_lvl = tsk[1]
            ene_lvl_ref = tsk[2]
            ene_ref_thy_info = finf.get_thy_info(ene_lvl_ref, thy_dct)
            ene_thy_info = finf.get_thy_info(ene_lvl, thy_dct)
            ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
                ene_coeff[ene_idx],
                ene_thy_info[3],
                ene_thy_info[1],
                ene_thy_info[2],
                ene_ref_thy_info[3],
                ene_ref_thy_info[1],
                ene_ref_thy_info[2]))
            ene_idx += 1
    ene_str += '!               '.join(ene_strl)

    return ene_str


def write_channel_mess_strs(spc_dct, rxn_lst, pes_formula,
                            ts_model, pf_levels, multi_info, pst_params,
                            spc_save_fs, save_prefix,
                            idx_dct, mess_strs):
    """ Write all the MESS input file strings for the reaction channels
    """
    first_ground_ene = 0.
    species = messrates.make_all_species_data(
        rxn_lst, spc_dct, save_prefix, ts_model, pf_levels)
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        tsform = automol.geom.formula(
            automol.zmatrix.geometry(spc_dct[tsname]['original_zma']))
        if tsform != pes_formula:
            print('Reaction ist contains reactions on different potential',
                  'energy surfaces: {} and {}'.format(tsform, pes_formula))
            print('Will proceed to construct only {}'.format(pes_formula))
            continue
        mess_strs, first_ground_ene = messrates.make_channel_pfs(
            tsname, rxn, species, spc_dct, idx_dct, mess_strs,
            first_ground_ene, spc_save_fs, ts_model, pf_levels,
            multi_info,
            pst_params=pst_params)
    well_str, bim_str, ts_str = mess_strs
    ts_str += '\nEnd\n'
    print(well_str)
    print(bim_str)
    print(ts_str)

    return well_str, bim_str, ts_str
