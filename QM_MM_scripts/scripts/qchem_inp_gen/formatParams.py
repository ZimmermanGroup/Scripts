#!/usr/bin/env python

import re
from sys import argv
paramFile = argv[1] #select file with params from tinker without forcefield section

inFile = open(paramFile, 'r')

content = inFile.readlines()

inFile.close()

multiplicity = 2 ## multiplicty to add to Improper for tinker format (without multiplicity)

atomTypes = []
charges = []
vdw = []
bonds = []
angles = []
torsions = []
ureybrad = []
improper = []

def formatType(string, num):
		for i in range(1,num+1):
				string[i] = "{}".format("-"+x[i])
				
for line in content:
		if re.search('atom', line):
				x = line.split()
				formatType(x, 1)
				atomTypes.append(x[1:2])
		if re.search('vdw', line):
				x = line.split()
				formatType(x, 1)
				vdw.append(x[2:4])
		if re.search('charge', line):
				x = line.split()
				formatType(x, 1)
				charges.append(x[2:3])
		if re.search('bond', line):
				x = line.split()
				formatType(x, 2)
				bonds.append(x[1:])
		if re.search('angle', line):
				x = line.split()
				formatType(x, 3)
				angles.append(x[1:])
		if re.search('torsion', line):
				x = line.split()
				formatType(x, 4)
				torsions.append(x[1:])
		if re.search('ureybrad', line):
				x = line.split()
				formatType(x, 3)
				ureybrad.append(x[1:])
		if re.search('improper', line):
				x = line.split()
				formatType(x, 4)
				improper.append(x[1:])

newTypes = []
for i in range(len(atomTypes)):
		new = atomTypes[i] + charges[i] + vdw[i]
		newTypes.append(new)

outFile = open('qchem_formated.prm', 'w')

for line in newTypes:
    outFile.write("AtomType  ")
    for item in line:
        outFile.write(item + "  ")
    outFile.write('\n')
outFile.write('\n')
for line in bonds:
		outFile.write("Bond  ")
		for item in line:
				outFile.write(item + "  ")
		outFile.write('\n')
outFile.write('\n')
for line in angles:
    outFile.write("Angle  ")
    for item in line:
        outFile.write(item + "  ")
    outFile.write('\n')
outFile.write('\n')
for line in torsions:
    outFile.write("Torsion  ")
    for item in line:
        outFile.write(item + "  ")
    outFile.write('\n')
outFile.write('\n')
for line in ureybrad:
		outFile.write("UreyBrad ")
		for item in line:
			outFile.write(item + "  ")
		outFile.write('\n')
for line in improper:
		outFile.write("Improper ")
		newOrder = [1, 2, 0, 3, 4, 5]
		line = [line[i] for i in newOrder]
		line.append(multiplicity)
		print(line)
		for item in line:
			outFile.write(item + "  ")
		outFile.write('\n')

outFile.write("$end")
outFile.write('\n')

outFile.close()
