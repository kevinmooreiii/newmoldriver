"""
utility functions
"""

import automol


def get_stoich(harm_min_cnf_locs_i, harm_min_cnf_locs_j,
               harm_cnf_save_fs_i, harm_cnf_save_fs_j):
    """ get the overall combined stoichiometry
    """
    if harm_min_cnf_locs_i is not None:
        harm_geo_i = harm_cnf_save_fs_i.leaf.file.geometry.read(
            harm_min_cnf_locs_i)
        if harm_min_cnf_locs_j is not None:
            harm_geo_j = harm_cnf_save_fs_j.leaf.file.geometry.read(
                harm_min_cnf_locs_j)

    form_i = automol.geom.formula(harm_geo_i)
    form_j = automol.geom.formula(harm_geo_j)
    form = automol.formula.join(form_i, form_j)
    stoich = ''
    for key, val in form.items():
        stoich += key + str(val)

    return stoich
