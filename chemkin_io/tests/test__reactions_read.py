""" Tests the parsing of reaction strings
"""

import os
import numpy as np
import ioformat
from chemkin_io.parser.reaction import get_rxn_param_dct as parser


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


def test_arr():
    """ Tests the Arrhenius reader
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2+M=OH+O+M     1.000E+15     0.000    25000\n'
        '     N2/1.400/   AR/1.000/   \n\n\n'  # note three spaces before \n
        'END\n\n')
    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2=OH+O     1.000E+15     0.000    25000\n'
        'DUP\n'
        'H+O2=OH+O     1.000E+15     0.000    25000\n'
        'DUP\n\n\n'
        'END\n\n')

    rxn_param_dct1 = parser(ckin_str1, 'cal/mole', 'moles')
    rxn_param_dct2 = parser(ckin_str2, 'cal/mole', 'moles')

    for params in rxn_param_dct1.values():  # should only be one rxn
        arr_tuples = params.arr
        arr_collid = params.arr_collid
        assert len(arr_tuples) == 1
        for arr_tuple in arr_tuples:
            assert np.allclose(arr_tuple, [1e15, 0, 25000])
        for spc, eff in arr_collid.items():
            assert spc in ('N2', 'AR')
            assert np.isclose(1.4, eff) or np.isclose(1.0, eff)

    for params in rxn_param_dct2.values():
        arr_tuples = params.arr
        assert len(arr_tuples) == 2  # should be a duplicate
        for arr_tuple in arr_tuples:
            assert np.allclose(arr_tuple, [1e15, 0, 25000])


def test_plog():
    """ Tests the PLOG reader and writer
    """

    ckin_str1 = 'REACTIONS     CAL/MOLE     MOLES\n\n' \
        'H+O2=OH+O     1.000E+15     0.000    25000   ! Duplicates exist at ' \
        '1 atm (see below); only single 1-atm fit is written\n' \
        '    PLOG /1.000E-01   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+02   1.000E+15     0.000    25000 /\n\n\nEND\n\n'

    # This one is for checking the incorrect way of writing duplicate PLOGS
    ckin_str2 = 'REACTIONS     CAL/MOLE     MOLES\n\n' \
        'H+O2=OH+O     1.000E+15     0.000    25000   ! Duplicates exist at ' \
        '1 atm (see below); only single 1-atm fit is written\n' \
        '    PLOG /1.000E-01   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+02   1.000E+15     0.000    25000 /\nDUP\n' \
        'H+O2=OH+O     1.000E+15     0.000    25000   ! Duplicates exist at ' \
        '1 atm (see below); only single 1-atm fit is written\n' \
        '    PLOG /1.000E-01   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+02   1.000E+15     0.000    25000 /\nDUP\n\n\nEND' \
        '\n\n'

    rxn_param_dct1 = parser(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        plog_dct = params.plog
        for pressure, arr_tuples in plog_dct.items():
            for arr_tuple in arr_tuples:
                assert np.allclose(arr_tuple, [1e15, 0, 25000])
            if pressure == 1:  # this pressure should have duplicate Arrhenius
                assert len(arr_tuples) == 2

    # Check the duplicate case
    rxn_param_dct2 = parser(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        plog_dct = params.plog
        plog_dups = params.plog_dups
        # Check the plog_dct
        for pressure, arr_tuples in plog_dct.items():
            for arr_tuple in arr_tuples:
                assert np.allclose(arr_tuple, [1e15, 0, 25000])
            if pressure == 1:  # this pressure should have duplicate Arrhenius
                assert len(arr_tuples) == 2
        # Check the plog_dups
        assert len(plog_dups) == 1
        dup_dct = plog_dups[0]
        for pressure, arr_tuples in dup_dct.items():
            for arr_tuple in arr_tuples:
                assert np.allclose(arr_tuple, [1e15, 0, 25000])
            if pressure == 1:  # this pressure should have duplicate Arrhenius
                assert len(arr_tuples) == 2


def test_cheb():
    """ Tests the Chebyshev reader and writer
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+00     0.000        0\n'
        '    TCHEB/      500.00     2000.00 /\n'
        '    PCHEB/        0.03      100.00 /\n'
        '    CHEB /           6           4 /\n'
        '    CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '    CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '    CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '    CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '    CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '    CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n\n\n'
        'END\n\n')

    # For duplicates
    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+00     0.000        0\n'
        '    TCHEB/      500.00     2000.00 /\n'
        '    PCHEB/        0.03      100.00 /\n'
        '    CHEB /           6           4 /\n'
        '    CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '    CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '    CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '    CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '    CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '    CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n'
        'DUP\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+00     0.000        0\n'
        '    TCHEB/      500.00     2000.00 /\n'
        '    PCHEB/        0.03      100.00 /\n'
        '    CHEB /           6           4 /\n'
        '    CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '    CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '    CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '    CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '    CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '    CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n'
        'DUP\n\n\nEND\n\n')

    rxn_param_dct1 = parser(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        cheb_dct = params.cheb
        tlim = cheb_dct['tlim']
        plim = cheb_dct['plim']
        alpha = cheb_dct['alpha']
        assert tlim == (500, 2000)
        assert plim == (0.03, 100)
        assert np.shape(alpha) == (6, 4)

    # Do the duplicates check
    rxn_param_dct2 = parser(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        cheb_dct = params.cheb
        tlim = cheb_dct['tlim']
        plim = cheb_dct['plim']
        alpha = cheb_dct['alpha']
        cheb_dups = params.cheb_dups
        assert tlim == (500, 2000)
        assert plim == (0.03, 100)
        assert np.shape(alpha) == (6, 4)
        assert len(cheb_dups) == 1


def test_troe():
    """ Tests the Troe reader and writer
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '    TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '     AR/1.400/   N2/1.700/   \n\n\nEND\n\n')

    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '    TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '    TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n\n\nEND\n\n')

    rxn_param_dct1 = parser(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        troe_dct = params.troe
        highp_arr = troe_dct['highp_arr']
        lowp_arr = troe_dct['lowp_arr']
        troe_params = troe_dct['troe_params']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        assert np.allclose(troe_params, [1.5, 8e3, 1e2, 1e3])

    # Do the duplicates check
    rxn_param_dct2 = parser(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        troe_dct = params.troe
        highp_arr = troe_dct['highp_arr']
        lowp_arr = troe_dct['lowp_arr']
        troe_params = troe_dct['troe_params']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        assert np.allclose(troe_params, [1.5, 8e3, 1e2, 1e3])
        assert len(params.troe_dups) == 1


def test_lind():
    """ Tests the Lindemann reader and writer
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '     AR/1.400/   N2/1.700/   \n\n\nEND\n\n')

    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n\n\nEND\n\n')

    rxn_param_dct1 = parser(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        lind_dct = params.lind
        highp_arr = lind_dct['highp_arr']
        lowp_arr = lind_dct['lowp_arr']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])

    # Do the duplicates check
    rxn_param_dct2 = parser(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        lind_dct = params.lind
        highp_arr = lind_dct['highp_arr']
        lowp_arr = lind_dct['lowp_arr']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        assert len(params.lind_dups) == 1


def test_rxn_names():
    """ test mechanalyzer.parser.reaction
    """
    ckin_str = ioformat.pathtools.read_file(DAT_PATH, 'rxn_block.dat')
    rxn_param_dct = parser(ckin_str, 'cal/mole', 'moles')
    print(rxn_param_dct)


if __name__ == '__main__':
    test_arr()
    test_plog()
    test_cheb()
    test_troe()
    test_lind()
    test_rxn_names()
