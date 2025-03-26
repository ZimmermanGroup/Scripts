#!/usr/bin/env python

import re
import sys

typesFile = sys.argv[1] ##File atomTypes.txt with atom types to match in format "charmm charge newtype"
paramFile = sys.argv[2] ##Files with bonds, angles and torsion from cgenff and formated with "bond" "angle" "torsion" on each line


inFile = open(typesFile, 'r')

atomsType = (inFile.read()).splitlines()

inFile.close()

atomType = []
charmType = []

for line in atomsType:
		singleLine = line.split()
		atomType.append(singleLine[2])
		charmType.append(singleLine[0])

matches = zip(atomType, charmType)
matchDict = dict(matches)

#reading params
inFile = open(paramFile, 'r')

content = (inFile.read()).splitlines()

inFile.close()
#making list of charmtypes and bonds values
bondsPar = []

for line in content:
		if re.search('bond', line):
				x = line.split()
				params = x[1:]
				bondsPar.append(params)
#creating dict for angles charmtypes with values
anglesPar = {}
ureybradPar = {}

for line in content:
		if re.search('angle', line):
				x = line.split()
				params = tuple(x[1:4])
				values = x[4] + " " + x[5]
				anglesPar[params] = values
				if (len(x) == 8):
						values = x[6] + " " + x[7]
						ureybradPar[params] = values

torsPar = {}

for line in content:
		if re.search('torsion', line):
				x = line.split()
				params = tuple(x[1:5])
				values = x[5] + " " + x[6] + " " + x[7]
				if params not in torsPar:
						torsPar[params] = values
				else:
						valAtKey = torsPar[params]
						extVal = valAtKey + " " + values
						torsPar[params] = extVal

#reading required values and atom types
inFile = open('out.log', 'r')
content = (inFile.read()).splitlines()
inFile.close()

bonds = []
for line in content:
    if re.search('Bond', line):
        x = line.split()
        params = x[3:]
        bonds.append(params)

angles = []
for line in content:
		if re.search('Angle', line):
				x = line.split()
				params = x[4:]
				angles.append(params)

torsions = []
for line in content:
    if re.search('Torsion', line):
        x = line.split()
        params = x[5:]
        torsions.append(params)

outFile = open('bonds.prm', 'w')

for line in bonds:
		atom1 = matchDict[line[0]]
		atom2 = matchDict[line[1]]
		#print(atom1 + " " + atom2 + " " + line[0] + " " + line[1])
		for param in bondsPar:
				param1 = param[0]
				param2 = param[1]
				if (atom1 == param1 and atom2 == param2) or (atom1 == param2 and atom2 == param1):
						paramLine = "bond " + line[0] + " " + line[1] + " " + param[2] + " " + param[3] 
						outFile.write(paramLine + '\n')

outFile.close()

outFile = open('angles.prm', 'w')

ureybradList = []
for line in angles:
		atoms = []
		for i in range(3):
				atom = matchDict[line[i]]
				atoms.append(atom)
		atoms = tuple(atoms)
		if atoms not in anglesPar:
				atoms = tuple(reversed(atoms))	
		params = anglesPar[atoms]
		paramLine = "angle " + " " + line[0] + " " + line[1] + " " + line[2] + " " + params
		if atoms in ureybradPar:
				parUB = ureybradPar[atoms]
				paramLineUB = "ureybrad" + " " + line[0] + " " + line[1] + " " + line[2] + " " + parUB
				ureybradList.append(paramLineUB)
		outFile.write(paramLine + '\n')

outFile.close()

outFile = open('ureybrad.prm', 'w')
for item in ureybradList:
    outFile.write(item + '\n')
outFile.close()

outFile = open('torsions.prm', 'w')

for line in torsions:
		atoms = []
		for i in range(4):
				atom = matchDict[line[i]]
				atoms.append(atom)
		atoms = tuple(atoms)
		if atoms not in torsPar:
				atoms = tuple(reversed(atoms))
		params = torsPar[atoms]
		atomTypes = line[0] + " " + line[1] + " " + line[2] + " " + line[3] + " "
		paramLine = "torsion " + atomTypes + params
		outFile.write(paramLine + '\n')

outFile.close()

