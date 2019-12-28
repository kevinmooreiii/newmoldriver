"""
Library of functions to parse the mechanism files
"""

import os

# new libs
from lib import moldr


def build_geom_dct(data_path):
    """ Obtain dct containing geometries to use as input
    """
    geom_path = os.path.join(data_path, 'data', 'geoms')
    geom_dct = moldr.util.geometry_dictionary(geom_path)
    return geom_dct
