#This file was created on 08/12/2018

The purpose of Zybel (ZStruct + Pybel) is to allow easy prototyping of reaction
coordinates for the Growing String Method.


There are two examples in this folder: 
1) Bimolecular zstruct

To run bimolecular zstruct for butadiene dimer
 ./bimolecular_example.py --substrate1="C=CC=C" --substrate2="C=CC=C" 
or 
 ./bimolecular_example.py --substrate1="C=CC=C" --substrate2="C=CC=C"  --read_xyz
which requires react1.xyz and react2.xyz

The reactive atoms are determined by smarts patterns.


2) Tethered zstruct
To run tethered zstruct run the script.sh. This creates driving coordinates and conformers 
for a tethered penta-enyl-benzene. The conformers are created in-situ using the smiles 
strings. 




