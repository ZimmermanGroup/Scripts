#!/usr/bin/env python
import sys
import os
import pybel 
import numpy as np
import openbabel as ob
import itertools
import iupac
import align

# => Please see license.txt for licensing and copyright information <= 
#       =>    Zimmerman Group, University of Michigan <= #

def zstruct1(omol,r1,r2,b1,doOne=True,doTwo=False):
    i=0
    output_format = 'xyz'
    substrate_path=os.getcwd()
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)


    # generate all one-permutations of r1 and r2 
    if (doOne==True):
        z=[]
        for atom in r2:
           for atom2 in r1:
               z.append([(atom,atom2)])
        #print(z)
        for isomer in z:
            output_path=substrate_path+"/scratch/initial%04d.xyz" % i
            fname=str("ISOMERS%04d" % i)
            s1=" ".join([str(j) for j in isomer[0]])
            i+=1
            with open(fname, 'w') as f:
                f.write("ADD " + s1 + "\n")
            with open(output_path, 'w') as f:
                f.write(obconversion.WriteString(omol))

    # ====== generate all one-permutations of add and break ===== #
    if not (b1 is None):
        zb=[[x[0],b1] for x in z]
        #print(zb)
        for isomer in zb:
            output_path=substrate_path+"/scratch/initial%04d.xyz" % i
            fname=str("ISOMERS%04d" % i)
            s1=" ".join([str(j) for j in isomer[0]])
            s2=" ".join([str(j) for j in isomer[1]])
            i+=1
            with open(fname, 'w') as f:
                f.write("ADD " + s1 + "\n")
                f.write("BREAK " + s2 + "\n")
            with open(output_path, 'w') as f:
                f.write(obconversion.WriteString(omol))
    # ========generate all two-permutations of r1 and r2 ====== #
    if (doTwo==True):
        z2=[zip(x,r1) for x in itertools.permutations(r2,len(r1))]

        for isomer in z2:
            output_path=substrate_path+"/scratch/initial%04d.xyz" % i
            fname=str("ISOMERS%04d" % i)
            s1=" ".join([str(j) for j in isomer[0]])
            s2=" ".join([str(j) for j in isomer[1]])
            i+=1
            r1c=[isomer[0][0],isomer[1][0]]
            r2c=[isomer[0][1],isomer[1][1]]
            with open(fname, 'w') as f:
                f.write("ADD " + s1 + "\n")
                f.write("ADD " + s2 + "\n")
            with open(output_path, 'w') as f:
                f.write(obconversion.WriteString(omol))


def zstruct2(mol1,mol2,r1,r2,b1,doOne=True,doTwo=False):
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat("xyz")
    i=0
    # generate all one-permutations of r1 and r2 
    if (doOne==True):
        z=[]
        for atom in r2:
           for atom2 in r1:
               z.append([(atom,atom2)])
        #print(z)
        for isomer in z:
            output_path=os.getcwd()+"/initial%04d.xyz" % i
            omol=align.align(mol1,mol2,1,list(isomer[0]))
            fname=str("ISOMERS%04d" % i)
            s1=" ".join([str(j) for j in isomer[0]])
            i+=1
            with open(fname, 'w') as f:
                f.write("ADD " + s1 + "\n")
            with open(output_path, 'w') as f:
                f.write(obconversion.WriteString(omol))
    # ====== generate all one-permutations of add and break ===== #
    if not (b1 is None):
        zb=[[x[0],b1] for x in z]
        #print(zb)
        for isomer in zb:
            output_path=os.getcwd()+"/initial%04d.xyz" % i
            fname=str("ISOMERS%04d" % i)
            s1=" ".join([str(j) for j in isomer[0]])
            s2=" ".join([str(j) for j in isomer[1]])
            i+=1
            omol=align.align(mol1,mol2,1,list(isomer[0]))
            with open(fname, 'w') as f:
                f.write("ADD " + s1 + "\n")
                f.write("BREAK " + s2 + "\n")
            with open(output_path, 'w') as f:
                f.write(obconversion.WriteString(omol))
    # ========generate all two-permutations of r1 and r2 ====== #
    if (doTwo==True):
        z2=[zip(x,r1) for x in itertools.permutations(r2,len(r1))]
        for isomer in z2:
            output_path=os.getcwd()+"/initial%04d.xyz" % i
            fname=str("ISOMERS%04d" % i)
            s1=" ".join([str(j) for j in isomer[0]])
            s2=" ".join([str(j) for j in isomer[1]])
            isom_list = list(isomer[0]) + list(isomer[1])
            i+=1
            omol=align.align(mol1,mol2,2,isom_list)
            with open(fname, 'w') as f:
                f.write("ADD " + s1 + "\n")
                f.write("ADD " + s2 + "\n")
            with open(output_path, 'w') as f:
                f.write(obconversion.WriteString(omol))

