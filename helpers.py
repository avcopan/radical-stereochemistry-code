import functools
import re
from typing import Optional

import automol
import autoreact
import chemkin_io
import pandas


# Generic
def lookup_rows(
    df: pandas.DataFrame, col_val_dct: dict[str, object]
) -> pandas.DataFrame:
    """Look up rows in a pandas DataFrame

    :param df: A DataFrame
    :param col_val_dct: Key-value pairs for selecting rows by column value
    """
    sels = [df[k] == v for k, v in col_val_dct.items()]
    return df[functools.reduce(lambda x, y: x & y, sels)]


def lookup_row(df: pandas.DataFrame, col_val_dct: dict[str, object]) -> pandas.Series:
    """Look up one row in a pandas DataFrame

    Finds the first matching row

    :param df: A DataFrame
    :param col_val_dct: Key-value pairs for selecting rows by column value
    """
    sel_df = lookup_rows(df, col_val_dct)
    return sel_df.iloc[0]


def lookup_value(
    df: pandas.DataFrame, key: str, col_val_dct: dict[str, object]
) -> object:
    """Look up one row value in a pandas DataFrame

    Finds the first matching row

    :param df: A DataFrame
    :param key: The key value to look up
    :param col_val_dct: Key-value pairs for selecting rows by column value
    """
    row = lookup_row(df, col_val_dct)
    return row[key]


def lookup_value_series(
    df: pandas.DataFrame, key: str, col_val_dcts: list[dict[str, object]]
) -> object:
    """Look up one row values for a series in a pandas DataFrame

    Finds the first matching row for each member of the series

    :param df: A DataFrame
    :param key: The key value to look up
    :param col_val_dcts: Key-value pairs for selecting rows by column value
    """
    rows = [lookup_row(df, c) for c in col_val_dcts]
    return tuple(r[key] for r in rows)


def with_columns_from_other(
    df: pandas.DataFrame, other_df: pandas.DataFrame, merge_on: list[str]
) -> pandas.DataFrame:
    """Add missing columns from another DataFrame

    No new rows are added -- performs a non-redundant left merge

    :param df: A DataFrame
    :param other_df: A DataFrame whose columns are to be added
    :param merge_on: The columns to merge on
    :return: The DataFrame with added columns from `other_df`
    """
    other_df = other_df.drop_duplicates(subset=merge_on)
    df = pandas.merge(df, other_df, how="left", on=merge_on, suffixes=(None, "_y"))
    df.drop(df.filter(regex="_y$").columns, axis=1, inplace=True)
    return df


# Species
def racemic_species_name(name: str, chi: str) -> str:
    """Get the racemic form of a species name (has no effect on achiral species names)

    :param name: A CHEMKIN species name
    :param chi: An InChI or AMChI string
    :return: The CHEMKIN name of the racemized species
    """
    if not automol.amchi.is_enantiomer(chi):
        return name
    assert name.endswith("0") or name.endswith("1")
    return re.sub("[01]$", "r", name)


def reflect_species_name(name: str, chi: str) -> str:
    """Get the name of the mirror image of a species, which is either the species itself
    or its enantiomer

    :param name: A CHEMKIN species name
    :param chi: An InChI or AMChI string
    :return: The CHEMKIN name of the mirror image
    """
    if not automol.amchi.is_enantiomer(chi):
        return name
    assert name.endswith("0") or name.endswith("1")
    return re.sub("1$", "0", name) if name.endswith("1") else re.sub("0$", "1", name)


# Reaction
def chiral_reactant_count(name: str, chi_dct: dict[str, str]) -> int:
    """Get the number of chiral reactants from a CHEMKIN name

    :param name: A CHEMKIN reaction name
    :param chi_dct: A dictionary mapping species names onto InChI or AMChI strings
    :return: The number of chiral reactants
    """
    rct_names, _, _ = chemkin_io.parser.reaction.get_rxn_name(name)
    return sum(map(automol.amchi.is_enantiomer, map(chi_dct.get, rct_names)))


def racemic_reaction_name(name: str, chi_dct: dict[str, str]) -> str:
    """Get the racemic form of a reaction name (has no effect on achiral species names)

    :param name: A CHEMKIN reaction name
    :param chi_dct: A dictionary mapping species names onto InChI or AMChI strings
    :return: The CHEMKIN name of the racemized reaction
    """
    rct_name0s, prd_name0s, third_body = chemkin_io.parser.reaction.get_rxn_name(name)
    rct_names = [racemic_species_name(n, chi_dct[n]) for n in rct_name0s]
    prd_names = [racemic_species_name(n, chi_dct[n]) for n in prd_name0s]
    return chemkin_io.writer.format_rxn_name((rct_names, prd_names, third_body))


def reflect_reaction_name(name: str, chi_dct: dict[str, str]) -> Optional[str]:
    """Get the name of the mirror image of a reaction (returns None if the reaction is achiral)

    :param name: A CHEMKIN reaction name
    :param chi_dct: A dictionary mapping species names onto InChI or AMChI strings
    :return: The CHEMKIN name of the mirror image, or None if the reaction is achiral
    """
    rct_name0s, prd_name0s, third_body = chemkin_io.parser.reaction.get_rxn_name(name)
    rct_chis = list(map(chi_dct.get, rct_name0s))
    prd_chis = list(map(chi_dct.get, prd_name0s))
    if not automol.amchi.is_enantiomer_reaction(rct_chis, prd_chis):
        return None

    rct_names = [reflect_species_name(n, c) for n, c in zip(rct_name0s, rct_chis)]
    prd_names = [reflect_species_name(n, c) for n, c in zip(prd_name0s, prd_chis)]
    return chemkin_io.writer.format_rxn_name((rct_names, prd_names, third_body))


# Mechanism
def reactions_with_params_objects(rxn_df: pandas.DataFrame) -> pandas.DataFrame:
    """Ensure that the reaction params are stored as objects, not strings

    Acts on the 'params' column

    :param rxn_df: A reaction DataFrame
    :return: A reaction DataFrame with reaction objects
    """
    if pandas.api.types.infer_dtype(rxn_df["params"]) == "string":
        rxn_df["params"] = rxn_df["params"].apply(
            autoreact.params.RxnParams.from_string
        )
    return rxn_df


def mechanism_dict(rxn_df: pandas.DataFrame) -> dict[tuple, object]:
    """Get the mechanism dictionary (reaction parameters by reagent names)

    :param rxn_df: A reaction DataFrame
    :return: The mechanism dictionary
    """
    rxn_df = reactions_with_params_objects(rxn_df)
    rxns = rxn_df["name"].map(chemkin_io.parser.reaction.get_rxn_name)
    rxn_dct = dict(zip(rxns, rxn_df["params"]))
    return rxn_dct


def species_for_mechanism(
    spc_df: pandas.DataFrame, rxn_dct: dict[tuple, object]
) -> pandas.DataFrame:
    """Get the subset of species for a mechanism dictionary

    :param spc_df: A species DataFrame
    :param rxn_dct: A mechanism dictionary
    :return: The species DataFrame containing only the species for this mechanism
    """
    spc_set = {s for r, p, _ in rxn_dct for s in r + p}
    spc_df = spc_df[spc_df["name"].isin(spc_set)]
    assert (
        set(spc_df["name"]) == spc_set
    ), f'Missing species: {set(spc_df["name"]) - spc_set}'
    return spc_df


def mechanism_string(rxn_dct: dict[tuple, object], spc_df: pandas.DataFrame) -> str:
    """Get a CHEMKIN-formatted mechanism string

    :param rxn_dct: A mechanism dictionary
    :param spc_df: A species DataFrame
    :return: The CHEMKIN-formatted mechanism string
    """
    spc_df["fml"] = spc_df["inchi"].map(automol.chi.formula)
    spc_df = spc_df[~spc_df["name"].duplicated(keep="first")]
    spc_df = spc_df.set_index("name")
    spc_dct = spc_df.to_dict("index")
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        rxn_param_dct=rxn_dct, mech_spc_dct=spc_dct
    )
    return mech_str
