#!/bin/python
from sys import argv


sampleFile=argv[1] #sample.qmmm.inp with MM params
infile=argv[2]	#xyz with joined connectivity
outfile=argv[3]	#name of outFile

i = open(infile, 'r')
o = open(outfile, 'w+')

#Insert the "head" portion of the input file
c=open(sampleFile,'r')

clines = c.readlines()
for k in range(0,len(clines)):
	tmp = clines[k]
	o.write(tmp)


#this sets up the line for each atom
lines = i.readlines()
del lines[0:2]

def formatLine(inString):
  prt = "{0: <5}".format(inString[1])
  prt = prt +"{0: <15}".format(inString[2])
  prt = prt +"{0: <15}".format(inString[3])
  prt = prt +"{0: <15}".format(inString[4])
  prt = prt +"{0: <9}".format("-" + inString[5])
  connectNum = len(inString) - 6
  zeroNum = 10 - len(inString)
  for i in range(connectNum):
    prt = prt +"{0: <6}".format(inString[6+i])
  for i in range(zeroNum):
    prt = prt +"{0: <6}".format("0")
  return prt


for line in lines:
  string = line.split()
  prt = formatLine(string)
  prt = prt+"\n"
  o.write(prt)	
o.write("$end" + '\n')
