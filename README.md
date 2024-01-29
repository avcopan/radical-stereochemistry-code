# Code for "Radical Stereochemistry: Accounting for Diastereomers in Kinetic Mechanism Development"

## Description

This repository contains the code for the paper "Radical Stereochemistry: Accounting
for Diastereomers in Kinetic Mechanism Development" by Andreas V. Copan, Kevin B. Moore
III, Sarah N. Elliott, Clayton R.  Mulvihill, Luna Pratalli Maffei, and Stephen J.
Klippenstein.
Both the source code and the Jupyter notebooks used to generate the results are included.

Main algorithms referenced in the paper:

- Calculate stereo parities (Scheme 1): `automol/graph/base/_08canon.py::_calculate_stereo_core()`
- Enumerate stereoisomers (Scheme 2): `automol/graph/base/_11canon.py::_expand_stereo_core()`
- Calculate TS stereo parities (Scheme 3): `automol/graph/base/_08canon.py::_calculate_ts_stereo()`
- Generate AMChI string (Figure 11): `automol/graph/base/_09amchi.py::_amchi_with_numbers()`

The procedure for generating a non-redundantly stereoexpanded mechanism is implemented
in the notebook `06_generate-non-redundantly-expanded-submechanism.ipynb`.

The procedure for generating a SMILES string with resonance bond stereochemistry is
shown in `automol/graph/base/_10smiles.py::_connected_smiles()`.

## Run the code

1. [Install pixi](https://pixi.sh/).
2. `pixi run jupyter notebook`
3. Run the Jupyter notebooks in numerical order.

## License

[MIT](https://choosealicense.com/licenses/mit/)