""" pf stuff
"""

import os
import mess_io


def get_pf_header(temp_step, ntemps):
    """ prepare partition function header string
    """
    global_pf_str = mess_io.writer.global_pf(
        [], temp_step, ntemps, rel_temp_inc=0.001,
        atom_dist_min=0.6)
    return global_pf_str


def get_pf_input(spc, spc_str, global_pf_str, zpe_str):
    """ prepare the full pf input string for running messpf
    """

    # create a messpf input file
    spc_head_str = 'Species ' + spc
    print('pf string test:', global_pf_str, spc_head_str, spc_str, zpe_str)
    pf_inp_str = '\n'.join(
        [global_pf_str, spc_head_str,
         spc_str, zpe_str, '\n'])
    return pf_inp_str


def write_pf_input(pf_inp_str, pf_path):
    """ write the pf.inp file
    """
    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)
