"""
  Library of Runtime Messages
"""


def tsk_msg(tsk, spc, thy_info, ini_thy_info):
    """ print a task message
    """
    print('\nTask {} \t {}//{} \t Species {}'.format(
        tsk, '/'.join(thy_info), '/'.join(ini_thy_info), spc))


def sadpt_tsk_msg(tsk, ts, spc_dct, thy_info, ini_thy_info):
    """ print a task message for a TS
    """
    print('Task {} \t for {} \t {}//{} \t {} = {}'.format(
        tsk, ts, '/'.join(thy_info), '/'.join(ini_thy_info),
        '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))


def ini_info_noavail_msg(tsk):
    """ print a message saying initial input info not available for task
    """
    print(
        'Initial level of theory for conformers must be ',
        'run before {} '.format(tsk))


def run_tsk_msg(tsk):
    """ print a message saying a task is running
    """
    print('running task {}'.format(tsk))


KTPMSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
          ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
          ================================================================
          ==                         KTPDRIVER                          ==
          ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
          ===      Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,  ===
          ==            Ahren Jasper, Stephen Klippenstein              ==
          ================================================================\n"""

THMMSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
          ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
          ================================================================
          ==                        THERMODRIVER                        ==
          ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
          ===    Murat Keceli, Yuri Georgievski, Stephen Klippenstein   ==
          ================================================================\n"""

ESMSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
          ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
          ================================================================
          ==                          ESDRIVER                          ==
          ====        Sarah Elliott, Andreas Copan, Kevin Moore,      ====
          ==            Carlo Cavolotti, Stephen Klippenstein           ==
          ================================================================\n"""
