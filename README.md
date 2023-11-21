# conformer
Conformation generator with RDKit

This is a modified code of  https://github.com/iwatobipen/rdk_confgen.git (https://iwatobipen.wordpress.com/)

## Description
 - This is simple script for generating conformers and tautomer from single molfiles
 - This script generates conformers with MMFF94s/UFF FF which is implemented in rdkit

## Requirements
 - RDKit
 - os
 - datetime
 - argparse

## Install
 
 ```
 $ git clone https://github.com/spectroscop/conformer.git
 $ cd conformer
 $ pip install -e .
 ```
## Basic usage
 - User need to prepare input molfile as a template.
 - Then run conformer command (with deafault settings) 

 ```
 $ python conformer.py -i ./sample/ligand.mol -o ./sample
## UFF
 $ python conformer.py -i ./sample/ligand.mol -o ./sample -f UFF
## tautomer 
 $ python conformer.py -i ./sample/ligand.mol -o ./sample -t true
## stereomer 
 $ python conformer.py -i ./sample/ligand.mol -o ./sample/ -s true 
## RMSD
 $ python conformer.py -i ./sample/ligand.mol -o ./sample -r 0.9
## macrocycle
 $ python conformer.py -i ./sample/ligand.mol -o ./sample -c 30
## number of conformers per structure
 $ python conformer.py -i ./sample/ligand.mol -o ./sample -m true
## not write FF energy
 $ python conformer.py -i ./sample/ligand.mol -o ./sample -a false
 ``` 

## Useful informations
 - [Coformer Generation using RDKit](https://www.rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf)
 - [Freely Available Conformer Generation Methods: How Good Are They?](https://pubs.acs.org/doi/abs/10.1021/ci2004658)