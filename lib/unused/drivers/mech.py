"""
Library of functions to parse the mechanism files
"""

import sys
import os
import collections
import json
import numpy
import chemkin_io
import automol
import moldr
import ktpdriver
from submission import read_dat

# new libs
from lib.phydat import phycon, eleclvl, symm


def build_geom_dct(data_path):
    """ Obtain dct containing geometries to use as input
    """
    geom_path = os.path.join(data_path, 'data', 'geoms')
    geom_dct = moldr.util.geometry_dictionary(geom_path)
    return geom_dct
