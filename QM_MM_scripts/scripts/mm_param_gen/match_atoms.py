#!/usr/bin/env python

import sys

matchesFile = sys.argv[1] # atomsType.txt

inFile = open(matchesFile, 'r')
atomTypes = (inFile.read()).splitlines()
inFile.close()

atomsDict = {}

for line in atomTypes:
		x = line.split()
		atomType = tuple(x[0:2])
		#print(atomType)
		newType = x[2]
		atomsDict[atomType] = newType

inFile = open('atoms.prm', 'r')
atoms = (inFile.read()).splitlines()
inFile.close()

outFile = open('matched_atoms.prm', 'w')
counter = 1
for line in atoms:
		x = line.split()
		index = x[0]
		charmmType = x[1]
		atom = tuple(x[1:3])
		newType = atomsDict[atom]
		matchLine = str(counter) + "   " + index + "   " + newType + "   " + charmmType
		outFile.write(matchLine + '\n')
		counter += 1

outFile.close()
