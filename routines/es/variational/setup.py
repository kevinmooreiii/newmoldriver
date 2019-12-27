""" 
 Setup the multireference stuff for the scan
"""

def setup_multiref_mep():
    """ Set-up information for the multireference MEP scan
    """
    # multi_info = ['molpro2015', 'caspt2', 'cc-pvtz', 'RR']
    multi_info = ['molpro2015', 'caspt2', 'cc-pvdz', 'RR']

    orb_restr = moldr.util.orbital_restriction(ts_info, multi_info)
    multi_level = multi_info[0:3]
    multi_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(rxn_run_path)
    thy_run_fs.leaf.create(multi_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(multi_level[1:4])

    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(multi_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(multi_level[1:4])

    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)

    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
    grid1 = grid[0]
    grid2 = grid[1]
    grid = numpy.append(grid[0], grid[1])
    high_mul = ts_dct['high_mul']
    print('starting multiref scan:', scn_run_fs.trunk.path())
            
    return None

def active_space(rcts)
    """ Determine the active space for the multireference MEP scan
    """
    rcts = ts_dct['reacs']
    if 'InChI=1S/O2/c1-2' in (spc_dct[rcts[0]]['ich'], spc_dct[rcts[1]]['ich']):
        num_act_orb = 5
        num_act_elc = 7
    else:
        num_act_orb = None
        num_act_elc = None

    return num_act_orb, num_act_elc
