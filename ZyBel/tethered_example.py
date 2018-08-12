#!/usr/bin/env python
import sys
import os
import pybel 
import numpy as np
import openbabel as ob
import itertools
import argparse
import iupac as ip
import confab 
import zstruct2


def main(substrate,group,sub_id,group_id,position,ref):

    # ========================> Generate Smiles < ========================= #
    if not group:
        smiles=substrate
    else:
        if position==1:
            # ============= Put group at the ortho positions ============== #
            smiles=substrate[0:2] +"(" + group +")" + substrate[2:]
        elif position==2:
            # ============= Put group at the meta positions =============== #
            smiles=substrate[0:3] +"(" + group +")" + substrate[3:]
        elif position==3:
            # ============= Put group at the para positions =============== #
            smiles=substrate[0:4] +"(" + group +")" + substrate[4:]
    print(smiles)

    # ============== Generate Folder Name ========== #
    folder=ip.iupac_name(smiles,substrate,group,sub_id,group_id,position)
    print folder
    os.system('mkdir %s' % folder)
    os.chdir(folder)
    os.system('mkdir scratch')

    # ============== Copy data to folder =========== #
    #TODO

    # ============== Generate Molecule ========== #
    mol=confab.gen3d(smiles)

    # ============ Driving Coordinate Idx Generation ==============#
    smarts1 = pybel.Smarts(substrate)
    smarts2 = pybel.Smarts("cccccc")
    r1=smarts1.findall(mol)
    r2=smarts2.findall(mol)

    #============ confab and align to reference  ==================#
    r1=[r1[0][-2], r1[0][-1]]
    r2=list(r2[0])
    mol=confab.confab(mol,r1,r2)
    # => align to substructure <=#
    path=os.getcwd()+"/tmp.xyz"
    if ref is not None:
        cmd_str="obabel %s %s -O %s -s %s --align" % (ref,path,path,substrate)
        print cmd_str
        os.system(cmd_str)
        #read aligned geom
        mol = confab.read_molecules( path, single = False )
        mol=mol[1]
    
    #============ do zstruct ===================================== #
    zstruct2.zstruct1(mol,r1,r2,None,doOne=False,doTwo=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse Bool")   
    parser.add_argument('--substrate', help='substrate smile', type=str, required=True)
    parser.add_argument('--group', help='ewg or edg group smile', type=str, required=True)
    parser.add_argument('--group_id', help='backup id for group in case iupac fail', type=int, required=False)
    parser.add_argument('--sub_id', help='backup id for sub in case iupac fail', type=int, required=True)
    parser.add_argument('--position', help='ortho', type=int, required=True)
    parser.add_argument('--ref', help='reference geometry', type=str, required=False)
    args = parser.parse_args()
    if args.ref==None:
        main(args.substrate,args.group,args.sub_id,args.group_id,args.position,None) 
    else:
        main(args.substrate,args.group,args.sub_id,args.group_id,args.position,args.ref) 

