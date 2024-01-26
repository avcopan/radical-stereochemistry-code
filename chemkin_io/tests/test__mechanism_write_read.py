""" Tests parsing of mechanism blocks
"""

from chemkin_io.parser import mechanism

SPECIES_BLOCK_GOOD = 'SPECIES\nCH4 H2 O2 N2\nEND'
SPECIES_BLOCK_BAD = 'SPECIES\nCH4 H2 O2 N2\n'
REACTION_BLOCK_GOOD = 'REACTIONS\nCH4+H=CH3+H2\nEND'
REACTION_BLOCK_BAD = 'REACTIONS\nCH4+H=CH3+H2\n'
THERMO_BLOCK_GOOD = 'THERMO\nC7H15OOH-1\nEND'
THERMO_BLOCK_BAD = 'THERMO\nC7H15OOH-1\n'
ELEMENT_BLOCK_GOOD = 'ELEMENTS\nC H O N\nEND'
ELEMENT_BLOCK_BAD = 'ELEMENTS\nC H O N\n'


def test_species_block():
    """ Tests the species_block function
    """
    species_block_good = mechanism.species_block(SPECIES_BLOCK_GOOD)
    species_block_bad = mechanism.species_block(SPECIES_BLOCK_BAD)
    assert species_block_good == '\nCH4 H2 O2 N2\n'
    assert species_block_bad is None


def test_reaction_block():
    """ Tests the reaction_block function
    """
    reaction_block_good = mechanism.reaction_block(REACTION_BLOCK_GOOD)
    reaction_block_bad = mechanism.reaction_block(REACTION_BLOCK_BAD)
    assert reaction_block_good == '\nCH4+H=CH3+H2\n'
    assert reaction_block_bad is None


def test_thermo_block():
    """ Tests the thermo_block function
    """
    thermo_block_good = mechanism.thermo_block(THERMO_BLOCK_GOOD)
    thermo_block_bad = mechanism.thermo_block(THERMO_BLOCK_BAD)
    assert thermo_block_good == '\nC7H15OOH-1\n'
    assert thermo_block_bad is None


def test_element_block():
    """ Tests the element_block function
    """
    element_block_good = mechanism.element_block(ELEMENT_BLOCK_GOOD)
    element_block_bad = mechanism.element_block(ELEMENT_BLOCK_BAD)
    assert element_block_good == '\nC H O N\n'
    assert element_block_bad is None


if __name__ == '__main__':
    test_species_block()
    test_reaction_block()
    test_thermo_block()
    test_element_block()
