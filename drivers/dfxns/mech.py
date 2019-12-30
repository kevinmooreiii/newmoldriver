""" estoktp driver functions with no
better home
"""

import os
import sys

from routines import util


def get_user_input():
    """ Read the options from the command line
    """
    # Set mechanism and type to be read based on user input
    data_path = sys.argv[1]
    mechanism_name = sys.argv[2]
    mech_type = sys.argv[3]
    mech_path = os.path.join(data_path, 'data', mechanism_name)
    mech_file = 'mech.json'

    # Set further parameters for what reactions and PESs to be run
    if len(sys.argv) > 4:
        pesnums = sys.argv[4]
    if len(sys.argv) > 5:
        channels = sys.argv[5]
    print('pesnums and params.channels:', pesnums, channels)

    return data_path, mech_path, mech_type, mech_file, pesnums, channels


def build_geom_dct(data_path):
    """ Obtain dct containing geometries to use as input
    """
    geom_path = os.path.join(data_path, 'data', 'geoms')
    geom_dct = util.geometry_dictionary(geom_path)
    return geom_dct


def etrans_lst(params):
    """ set the etrans list
    """
    return [params.EXP_FACTOR,
            params.EXP_POWER,
            params.EXP_CUTOFF,
            params.EPS1,
            params.EPS2,
            params.SIG1,
            params.SIG2,
            params.MASS1]
