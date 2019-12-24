""" libraries to build wavefunctions for Molpro multireference calculations
"""

import automol
import elstruct
import autofile
import moldr
import scripts.es


def cas_options_1(spc_info, formula, num_act_elc, num_act_orb, high_mul):
    """ prepare cas options for multireference wavefunctions
    """

    elec_count = automol.formula.electron_count(formula)
    closed_orb = (elec_count-num_act_elc)//2
    occ_orb = closed_orb + num_act_orb
    closed_orb = (elec_count-num_act_elc)//2 - 2
    two_spin = spc_info[2]-1
    high_two_spin = high_mul - 1
    chg = spc_info[1]
    cas_opt = [
        elstruct.option.specify(
            elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(
            elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.WFN_, elec_count, 1, two_spin, chg)
        ]

    wfn_str = (
        "{{uhf,maxit=300;wf,{0},1,{1}}}\n"
        "{{multi,maxit=40;closed,{2};occ,{3};wf,{0},1,{4};orbprint,3}}"
    ).format(elec_count, high_two_spin, closed_orb, occ_orb, two_spin)

    return cas_opt, wfn_str


def cas_options_2(spc_info, formula, num_act_elc, num_act_orb, high_mul):
    """ prepare cas options for multireference wavefunctions
    """

    elec_count = automol.formula.electron_count(formula)
    closed_orb = (elec_count-num_act_elc)//2
    occ_orb = closed_orb + num_act_orb
    two_spin = spc_info[2]-1
    high_two_spin = high_mul - 1
    chg = spc_info[1]
    cas_opt = [
        elstruct.option.specify(
            elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(
            elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.WFN_, elec_count, 1, two_spin, chg)
        ]

    wfn_str = (
        "{{uhf,maxit=300;wf,{0},1,{1}}}\n"
        "{{multi,maxit=40;closed,{2};occ,{3};wf,{0},1,{4};orbprint,3}}"
    ).format(elec_count, high_two_spin, closed_orb, occ_orb, two_spin)

    return cas_opt, wfn_str


def multiref_wavefunction_guess(high_mul, zma,
                                spc_info, thy_level,
                                casscf_options):
    """ prepare wavefunction template for multireference electronic structure calcs
    """

    charge = spc_info[1]
    mul = spc_info[2]
    basis = thy_level[2]
    prog = thy_level[0]
    prog = 'molpro2015'

    guess_str1 = elstruct.writer.energy(
        geom=zma,
        charge=charge,
        mult=high_mul,
        method='hf',
        basis=basis,
        prog=prog,
        orb_restricted=False,
        mol_options=['nosym'],
        )
    guess_str1 += '\n\n'
    guess_str1 = '\n'.join(guess_str1.splitlines()[2:-6])

    guess_str2 = elstruct.writer.energy(
        geom=zma,
        charge=charge,
        mult=mul,
        method='casscf',
        basis=basis,
        prog=prog,
        orb_restricted=True,
        casscf_options=casscf_options[0],
        mol_options=['nosym'],
        )
    guess_str2 += '\n\n'
    guess_str2 = '\n'.join(guess_str2.splitlines()[2:-6])

    guess_str = guess_str1 + guess_str2 

    if len(casscf_options) > 1:

        guess_str3 = elstruct.writer.energy(
            geom=zma,
            charge=charge,
            mult=mul,
            method='casscf',
            basis=basis,
            prog=prog,
            orb_restricted=True,
            casscf_options=casscf_options[1],
            mol_options=['nosym'],
            )
        guess_str3 += '\n\n'
        guess_str3 = '\n'.join(guess_str3.splitlines()[2:])

        guess_str += guess_str3

    return guess_str
