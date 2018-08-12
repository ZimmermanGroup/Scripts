#!/usr/bin/env python
import sys
import os
import pybel 
import numpy as np
import openbabel as ob
import itertools
import requests
from urllib2 import urlopen
import re
import argparse
import mm_params


def vectorAngle(v1,v2):
    dp = np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    if dp< -0.99999:
        dp=-0.99999
    if dp> 0.99999:
        dp=0.99999
    return np.arccos(dp)*180/np.pi 

def vdw_energy_1(mol,idx1,idx2):
    a1=mol.OBMol.GetAtom(idx1+1)
    a2=mol.OBMol.GetAtom(idx2+1)
    r=a2.GetDistance(a1)
    R = mm_params.ffR_table[a1.GetAtomicNum()]
    R += mm_params.ffR_table[a2.GetAtomicNum()]
    eps = mm_params.eps_table[a1.GetAtomicNum()]
    eps *= mm_params.eps_table[a2.GetAtomicNum()]
    eps = np.sqrt(eps)
    Rr = R/r
    Rr6 = Rr*Rr*Rr
    E = eps * (Rr6*Rr6 - 2 * Rr6 ) 
    return E

def vdw_energy(mol,natoms1,natoms2):
    r1=np.array(range(natoms1))
    r2=np.array(range(natoms1,natoms2+natoms1))
    E = 0.0
    for i in r1: 
        for j in r2:
            E += vdw_energy_1(mol,i,j)
     
    return E

def vdw_vector_opt(mol,natoms1,natoms2,add_idx,v):
    E=0.
    pE=0.
    Em=0.
    THRESH=0.001
    for i in range(200):
        pE=E
        E=vdw_energy(mol,natoms1,natoms2)
        if (E>pE+THRESH):
            break
        for j in range(natoms2):
            tmpatom=mol.OBMol.GetAtom(j+natoms1+1)
            tmpvec=np.array([tmpatom.GetX(),tmpatom.GetY(),tmpatom.GetZ()])
            tmpvec= tmpvec-0.1*v
            mol.OBMol.GetAtom(j+natoms1+1).SetVector(tmpvec[0],tmpvec[1],tmpvec[2])
    return mol


def cross_with_axis(atom):
    #x1=v1
    for nbr in ob.OBAtomAtomIter(atom):
       x1=np.array([atom.GetX()-nbr.GetX(),atom.GetY()-nbr.GetY(),atom.GetZ()-nbr.GetZ()])
    x2=np.array([0,0,1])
    
    v1=np.cross(x1,x2)
    norm=np.linalg.norm(v1)
    if norm<0.001:
        x2=np.array([0,1,0])
        v1=np.cross(x1,x2)
    return v1

def planar_cross(atom,nbrs):
    xs=[]
    for nbr in ob.OBAtomAtomIter(atom):
       xs.append([-atom.GetX()+nbr.GetX(),-atom.GetY()+nbr.GetY(),-atom.GetZ()+nbr.GetZ()])
    
    xs[0]=xs[0]/np.linalg.norm(xs[0])
    xs[1]=xs[1]/np.linalg.norm(xs[1])

    v1=np.cross(xs[0],xs[1])
    return v1

def GetNewBondVector(mol,idx,nvpf):
    atom=mol.OBMol.GetAtom(idx)
    nbrvecs=[]
    newbond=np.zeros((3))
    nbrs=[]
    for nbr in ob.OBAtomAtomIter(atom):
       # shouldn't really be += ...
       newbond+=np.array([atom.GetX()-nbr.GetX(),atom.GetY()-nbr.GetY(),atom.GetZ()-nbr.GetZ()])
       nbrvecs.append([atom.GetX()-nbr.GetX(),atom.GetY()-nbr.GetY(),atom.GetZ()-nbr.GetZ()])
       nbrs.append(nbr)
    
    if atom.GetValence()==2:
        ang=mol.OBMol.GetAngle(nbrs[0],atom,nbrs[1])
        #print ang
        if ang> 175.:
            newbond=cross_with_axis(atom)

    if atom.GetValence()==3:
        imptorv=mol.OBMol.GetTorsion(nbrs[0],atom,nbrs[1],nbrs[2])
        if abs(imptorv):
            newbond=planar_cross(atom,nbrs)
            nvpf=nvpf+1
    
    newbond=newbond/np.linalg.norm(newbond)
    return newbond,nvpf

def RotAboutAxisByAngle(v,angle):
    angle=angle*np.pi/180.
    s=np.sin(angle)
    c=np.cos(angle)
    t=1.-c
    vtmp=v
    vtmp=vtmp/np.linalg.norm(vtmp)
    x=vtmp[0]
    y=vtmp[1]
    z=vtmp[2]
    rotm =np.empty((3,3),dtype=float)
    rotm[0][0] = t*x*x + c
    rotm[0][1] = t*x*y - s*z
    rotm[0][2] = t*x*z + s*y
                             
    rotm[1][0] = t*x*y + s*z
    rotm[1][1] = t*y*y + c
    rotm[1][2] = t*y*z - s*x
                             
    rotm[2][0] = t*x*z - s*y
    rotm[2][1] = t*y*z + s*x
    rotm[2][2] = t*z*z + c
    return rotm

def align_to_x(mol,v2):
    axis=np.zeros((3))
    axis[0]=1.; 
    taxis=np.zeros((3))
    taxis[1]=1.; 
    taxis[2]=1.; 
    taxis=taxis/np.linalg.norm(taxis)
    v2=np.asarray(v2)
    v2=v2/np.linalg.norm(v2)

    # => Calculate the angle of the vector from the axis <= #
    angle=vectorAngle(v2,axis)

    if angle>0:
        # => Get axis for rotation <= #
        axisofrot=np.cross(axis,v2)
        axisofrot=axisofrot/np.linalg.norm(axisofrot)

        # => Now rotate v2 by angle in the axis of rotation <= #
        rotm=RotAboutAxisByAngle(axisofrot,-angle)
        for i in range(len(mol.atoms)):
            tmpatom=mol.OBMol.GetAtom(i+1)
            tmpvec=np.array([tmpatom.GetX(),tmpatom.GetY(),tmpatom.GetZ()])
            tmpvec=rotm.dot(tmpvec)
            mol.OBMol.GetAtom(i+1).SetVector(tmpvec[0],tmpvec[1],tmpvec[2])
            
    return mol


# => align the reactants by add indices <= #
def align(mol1,mol2,nadds,add_idx):
    nvpf1 = 0
    nvpf2 = 0
    v1a=[]
    v1b=[]
    print add_idx
    # =>  adjust for the 
    for i in range(nadds):
        add_idx[(i*2)+1] -= len(mol1.atoms)
    # => Get the New bond Vectors <= #
    for i in range(nadds):   
        tmpvec1,nvpf1=GetNewBondVector(mol1,add_idx[(i*2)+0],nvpf1)
        tmpvec2,nvpf2=GetNewBondVector(mol2,add_idx[(i*2)+1],nvpf2)
        v1a.append(tmpvec1)
        v1b.append(tmpvec2)
        if nvpf1>1:
            print "TODO"
        #TODO
        # if nvpf1 > 1:
        # v1a[i*3,(i+1)*3= align_v()
        # if nvpf2 > 1:
        # v1b[i*3,(i+1)*3= align_v()

    # => Averaging nadds v1 to get v2 <= #
    v2=	 A=np.empty((2,3),dtype=float)
    v2[0]=v1a[0]
    v2[1]=v1b[0]
    if nadds==2:
        v2[0]+=v1a[1]
        v2[1]+=v1b[1]
    n1=np.linalg.norm(v2[0])
    n2=np.linalg.norm(v2[1])
    if (n1<0.00001):
       #TODO
       atom=mol1.OBMol.GetAtom(add_idx[(nadds-1)+0])
       v1a=cross_with_axis(atom)
       n1=np.linalg.norm(v1a)
       v2[0]=v1a[0]
    if n2<0.000001:
       #TODO
       atom=mol2.OBMol.GetAtom(add_idx[(nadds-1)+1])
       v1b=cross_with_axis(atom)
       n2=np.linalg.norm(v1b)
       v2[1]=v1b[0]
    v2[0]=v2[0]/n1
    v2[1] =v2[1]/n2

    # => Get the central location for each index <= #
    c=	 A=np.empty((2,3),dtype=float)
    cfound=0
    for i in range(nadds):
        id1=add_idx[(i*2)+0]
        id2=add_idx[(i*2)+1]
        atom1=mol1.OBMol.GetAtom(id1)
        atom2=mol2.OBMol.GetAtom(id2)
        c[0] += np.array([atom1.GetX(),atom1.GetY(),atom1.GetZ()])
        c[1] += np.array([atom2.GetX(),atom2.GetY(),atom2.GetZ()])
        cfound+=1
    if (cfound>1):
        c/=cfound

    # => Move geom to center <= #
    for i in range(len(mol1.atoms)):
        tmpatom=mol1.OBMol.GetAtom(i+1)
        tmpvec=np.array([tmpatom.GetX(),tmpatom.GetY(),tmpatom.GetZ()])
        tmpvec = tmpvec-c[0]
        mol1.OBMol.GetAtom(i+1).SetVector(tmpvec[0],tmpvec[1],tmpvec[2])
    # => Move geom to center <= #
    for i in range(len(mol2.atoms)):
        tmpatom=mol2.OBMol.GetAtom(i+1)
        tmpvec=np.array([tmpatom.GetX(),tmpatom.GetY(),tmpatom.GetZ()])
        tmpvec = tmpvec-c[1]
        mol2.OBMol.GetAtom(i+1).SetVector(tmpvec[0],tmpvec[1],tmpvec[2])
    #mol1.write("xyz","cmol1.xyz",overwrite=True)
    #mol2.write("xyz","cmol2.xyz",overwrite=True)
   
    # => algin v2 along x/-x direction
    mol1=align_to_x(mol1,v2[0])
    mol2=align_to_x(mol2,v2[1])

    # => Face eachother at X Angstrom <= #
    for i in range(len(mol2.atoms)):
        tmpatom=mol2.OBMol.GetAtom(i+1)
        tmpvec=np.array([tmpatom.GetX(),tmpatom.GetY(),tmpatom.GetZ()])+ np.array([10.0,0.0,0.0])
        mol2.OBMol.GetAtom(i+1).SetVector(tmpvec[0],tmpvec[1],tmpvec[2])
    #mol1.write("xyz","amol1.xyz",overwrite=True)
    #mol2.write("xyz","amol2.xyz",overwrite=True)
     
    # => Combine mols <= #
    mol = ob.OBMol()
    for i in range(len(mol1.atoms)):
        tmpatom=mol1.OBMol.GetAtom(i+1)
        a=mol.NewAtom()
        a.SetAtomicNum(tmpatom.GetAtomicNum())
        tmpvec=np.array([tmpatom.GetX(),tmpatom.GetY(),tmpatom.GetZ()])
        mol.GetAtom(i+1).SetVector(tmpvec[0],tmpvec[1],tmpvec[2])
    j=len(mol1.atoms)
    for i in range(len(mol2.atoms)):
        tmpatom=mol2.OBMol.GetAtom(i+1)
        a=mol.NewAtom()
        a.SetAtomicNum(tmpatom.GetAtomicNum())
        tmpvec=np.array([tmpatom.GetX(),tmpatom.GetY(),tmpatom.GetZ()])
        mol.GetAtom(j+i+1).SetVector(tmpvec[0],tmpvec[1],tmpvec[2])


    mol = pybel.Molecule(mol)
    #mol.write("xyz", "beforeopt.xyz",overwrite=True)
    for i in range(nadds):
        add_idx[(i*2)+1] = add_idx[(i*2)+1] + len(mol1.atoms)

    # => Vdw opt <= #
    vx=np.zeros((3))
    vx[0]=1.0
    mol=vdw_vector_opt(mol,len(mol1.atoms),len(mol2.atoms),add_idx,vx)
    #mol = pybel.Molecule(mol)
    #mol.write("xyz", "afteropt.xyz",overwrite=True)
    return mol.OBMol


