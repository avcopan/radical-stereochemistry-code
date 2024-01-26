""" Writes the parts of a Chemkin file pertaining to species, i.e., the
    species block and the elements block
"""


def write_species(mech_spc_dct):
    """ Writes the contents of a mech_spc_dct in Chemkin format, with notes

        :param mech_spc_dct: species data for a mechanism
        :type mech_spc_dct: dct {spc_name:data}
        :return spc_str: str containing the species information
        :rtype: str 
    """ 

    # Get the max species name length
    max_spc_len = 0 
    for spc in mech_spc_dct.keys():
        if len(spc) > max_spc_len:
            max_spc_len = len(spc)

    # Get the max SMILES name length
    max_smiles_len = 0 
    for spc, ident_dct in mech_spc_dct.items():
        if len(ident_dct['smiles']) > max_smiles_len:
            max_smiles_len = len(ident_dct['smiles'])

    buffer = 5 

    # Write the spc_str
    spc_str = ''
    for spc, ident_dct in mech_spc_dct.items():
        spc_str += (
            '{0:<'+str(max_spc_len+buffer)+'s}{1:>9s}{2:<' +
            str(max_smiles_len+buffer)+'s}{3:>9s}{4:<9s}\n').format(
                spc, '! SMILES: ',
                ident_dct['smiles'], 'ChI: ', ident_dct['inchi'])

    return spc_str


def write_elements(mech_spc_dct):
    """ Writes the unique elements in a mech_spc_dct to a string

        :param mech_spc_dct: species data for a mechanism
        :type mech_spc_dct: dct {spc_name:data}
        :return elem_str: a str containing the list of unique elements
        :rtype: str
    """
    def _unique_elems(mech_spc_dct):
        """ Gets the unique elements in a mech_spc_dct 
        """
    
        # Build a non-unique list of all elements in the set
        all_elems = []
        for spc, spc_dct in mech_spc_dct.items():
            fml = spc_dct['fml']  # assumes that the fml exists!
            all_elems.extend(list(fml.keys()))
    
        # Get the unique elements via a set
        unique_elems = list(set(all_elems))
    
        return unique_elems

    unique_elems = _unique_elems(mech_spc_dct)
    elem_str = ''
    for unique_elem in unique_elems:
        elem_str += unique_elem + '\n'

    return elem_str

