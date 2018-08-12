import pybel
import openbabel as ob
import numpy as np
import os 

def read_molecules(filepath, single = False):
    in_format = filepath.strip().split( '.' )[-1]
    obconversion = ob.OBConversion()
    obconversion.SetInFormat( in_format )
    obmol = ob.OBMol()

    molecules = []
    notatend = obconversion.ReadFile( obmol, filepath )
    while notatend:
        molecules.append( obmol )
        obmol = ob.OBMol()
        notatend = obconversion.Read( obmol )

    if single:
        assert( len(molecules) == 1 )
        return molecules[0]
    else:
        return molecules

def gen3d(smi):
    #gen3d
    mol = pybel.readstring("smi", smi)
    mol.addh()
    mol.make3D()
    mol.localopt()
    return mol

def RMSD(conf_strings,r1,r2):
    #calculate RMSD between ethene and bezene for all conformations
    dist=[]
    pmols=[pybel.readstring("xyz",x) for x in conf_strings ] 
    for pmol in pmols:
        molcoords = [atom.coords for atom in pmol]
        r1_coords= [molcoords[x-1] for x in r1]
        #print(r1_coords)
        r2_coords= [molcoords[x-1] for x in r2]
        #print(r2_coords)
        center_r2=np.zeros(3)
        for atom in r2_coords:
            center_r2= center_r2+np.asarray(atom)
        center_r2=center_r2/float(len(r2))
        center_r1=np.zeros(3)
        for atom in r1_coords:
            center_r1= center_r1+np.asarray(atom)
        center_r1=center_r1/float(len(r1))
        dist.append(np.linalg.norm(np.subtract(center_r1,center_r2)))
        #dist.append(np.linalg.norm(np.subtract(r1_coords,r2_coords)))
    return dist

def confab(mol,r1,r2):
    rmsd_cutoff = 0.5
    conf_cutoff = 400000
    energy_cutoff = 100.0
    confab_verbose = False
    output_format = 'xyz'
    # Run Confab conformer generation
    conf_strings=[]
    pff = ob.OBForceField_FindType( "mmff94" )
    omol=mol.OBMol
    assert( pff.Setup(omol) ) # Make sure setup works OK
    pff.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, confab_verbose)
    pff.GetConformers(omol);
    confs_to_write = omol.NumConformers()
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)
    for conf_num in xrange(confs_to_write):
    	omol.SetConformer(conf_num);
    	conf_strings.append( obconversion.WriteString(omol))
    dist=RMSD(conf_strings,r1,r2)
    #get the minimum
    minidx=min(xrange(len(dist)), key=dist.__getitem__)
    omol.SetConformer(minidx)
    #write file 
    output_path=os.getcwd()+"/tmp.xyz"
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)
    with open(output_path, 'w') as f:
        f.write(obconversion.WriteString(omol))
    return omol

#def frozen_confab(input_path,atoms_to_unfreeze):
def frozen_confab(omol,atoms_to_unfreeze):
    rmsd_cutoff = 0.5
    conf_cutoff = 400000
    energy_cutoff = 100.0
    confab_verbose = False
    output_format = 'xyz'
    # Run Confab conformer generation
    #nmol = read_molecules( input_path, single = True )
    nmol=omol
    pff2 = ob.OBForceField_FindType( "mmff94" )
    assert( pff2.Setup(nmol) ) # Make sure setup works OK
    conf_strings=[]
    if len(atoms_to_unfreeze) > 0:
        print '%d atoms will be allowed to move; freezing others' % len( atoms_to_unfreeze )
        constraints = ob.OBFFConstraints()
        for atom in ob.OBMolAtomIter(nmol):
            atom_id = atom.GetIndex()+1
            if atom_id not in atoms_to_unfreeze:
                constraints.AddAtomConstraint(atom_id)
            else:
                print(atom_id)
        pff2.SetConstraints( constraints )
    pff2.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, confab_verbose)
    pff2.GetConformers(nmol);
    confs_to_write = nmol.NumConformers()-1
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)
    #output_path=os.getcwd()+"/test.xyz"
    #for conf_num in xrange(confs_to_write):
    #    nmol.SetConformer(conf_num);
    #    conf_strings.append( obconversion.WriteString(nmol))
    #with open(output_path, 'w') as f:
    #    for conf_num in xrange(confs_to_write):
    #        nmol.SetConformer(conf_num);
    #        f.write(obconversion.WriteString(nmol))
    nmol.SetConformer(confs_to_write-1)
    return nmol
