#usr/bin/env python


### Written by: ####################################################
#  _____              _____               _ _      _               # 
# |_   _|            |  __ \             | | |    | |              # 
#   | |  __ _ _ __   | |__) |__ _ __   __| | | ___| |_ ___  _ __   # 
#   | | / _` | '_ \  |  ___/ _ \ '_ \ / _` | |/ _ \ __/ _ \| '_ \  #
#  _| || (_| | | | | | |  |  __/ | | | (_| | |  __/ || (_) | | | | #
# |_____\__,_|_| |_| |_|   \___|_| |_|\__,_|_|\___|\__\___/|_| |_| #
#                                                                  #
### contact info at www.pendletonian.com ###########################

### NOTES ##############################################################################
########################################################################################
## ONLY WORKS FOR QCHEM ################################################################
## Want to learn more about why we do this? see http://gaussian.com/thermo/           ##
## Other references: http://pubs.acs.org/doi/suppl/10.1021/acs.macromol.6b01648       ## 
########################################################################################

########################################################################################
# USER INPUT VARIABLES # 
Freq = 50.0   #### Set the minimum wavenumber cutoff *frequency cutoff*
Temp = 273.15 #### Temperature in K for the entropy calculation - should match isotopes temp in qchem
imagfreq = 20  #### For freq calculations on TSs this will set the number of imaginary frequencies to leave in file

########################################################################################
##module list###
import sys
import os
import re
import math
import numpy as np
from scipy import constants
########################################################################################
##### Hardcode - set values - do not adjust####
directory=os.path.dirname(os.path.abspath(__file__)) 
MF = np.array([Freq])
C_h=constants.value('Planck constant')
C_c=constants.value('speed of light in vacuum') 
C_kb=constants.value('Boltzmann constant')
C_ug=1.9872036 #gas constant

########################################################################################
##################                    Code Start                 #######################
########################################################################################

########################################################################################
# Reads in isotope presence from frequency calculation more information here: 
# https://www.q-chem.com/qchem-website/manual/qchem43_manual/sect0066.html (Isotopes 
# section of qchem) 
########################################################################################
#### Section for default qchem frequency calculations (No ISOTOPES section) ############
#Get frequencies from .out file into a numpy array
def FreqArray(file1):
 with open(file1, 'r') as searchfile:
  x=[]
  for linenum, line in enumerate(searchfile):
   if "Frequency:" in line:
    freq=re.findall("[-+]?\d+[\.]?\d*", line)
    for item in freq:
     item2 = float(item)
     x.append(item2)
  array=np.asarray(x)
  return array

def OutputConfigure1(file1, fname, out, c):
 with open(file1, 'r') as searchfile:
  x = []
  for linenum, line in enumerate(searchfile, 1):
# if "Vibman isotope loop:     2" in line: 
#  for line in searchfile:
   if "Translational Entropy:" in line:
    TE_s=re.findall("[-+]?\d+[\.]?\d*", line)
    for item in TE_s:
     item2 = str(item)
     x.append(item2)
    out.write('%s -->' %fname)
#    out.write(" Vib. Entropy: %s" % c)  ## prints calculated vibration entropy
    out.write(' ')
#    out.write('%s ' % item2) ## Enable to print translational entropy to output file
   if "Rotational Entropy:" in line:
    RE_s=re.findall("[-+]?\d+[\.]?\d*", line)
    for item3 in RE_s:
     item4 = str(item3)
     x.append(item4)
#    out.write(' %s' % item4) ## Enable to print Rotational Entropy to output file
   if "Total Entropy:" in line:
    RE_s=re.findall("[-+]?\d+[\.]?\d*", line)
    for item5 in RE_s:
     item6 = str(item5)
  x.append(c) 
  array=np.asarray(x)
  array2=array.astype(np.float)
  Sum=np.sum(array2)
  out.write(' ')
  out.write('Total Entropy:')
  out.write(' ')
  out.write('%s cal/mol.K' %Sum)
  out.write(' ')
  out.write('Original Entropy:')
  out.write(' ')
  out.write('%s cal/mol.K' %item6)
  out.write('\n')

#### Section for qchem frequency calculations with ISOTOPES TRUE  #####################
### Getting frequency table read into NP for ISOTOPES LOOP 2 ##
#Get frequencies from .out file into a numpy array
def FreqArray2(file1):
 with open(file1, 'r') as searchfile:
  for linenum, line in enumerate(searchfile, 1):
   if "Vibman isotope loop:     2" in line: 
    freqarray=np.array([])
    count=1
    x = []
    for line in searchfile:
     if "Frequency:" in line:
      freq=re.findall("[-+]?\d+[\.]?\d*", line)
      for item in freq:
       item2 = float(item)
       x.append(item2)
    array=np.asarray(x)
    return array

def OutputConfigure2(file1, fname, out, c):
 with open(file1, 'r') as searchfile:
  for linenum, line in enumerate(searchfile, 1):
   if "Vibman isotope loop:     2" in line: 
    x = []
    for line in searchfile:
     if "Translational Entropy:" in line:
      TE_s=re.findall("[-+]?\d+[\.]?\d*", line)
      for item in TE_s:
       item2 = str(item)
       x.append(item2)
      out.write('%s -->' %fname)
#      out.write(" Vib. Entropy: %s" % c)  ## Prints calculated vibration entropy
      out.write(' ')
#      out.write(item2) ## Enable to print translational entropy to output file
     if "Rotational Entropy:" in line:
      RE_s=re.findall("[-+]?\d+[\.]?\d*", line)
      for item3 in RE_s:
       item4 = str(item3)
       x.append(item4)
#      out.write(' ')
#      out.write(item4) ## Enable to print Rotational Entropy to output file
     if "Total Entropy:" in line:
      RE_s=re.findall("[-+]?\d+[\.]?\d*", line)
      for item5 in RE_s:
       item6 = str(item5)
    x.append(c) 
    array=np.asarray(x)
    array2=array.astype(np.float)
    Sum=np.sum(array2)
    out.write(' ')
    out.write('Total Entropy:')
    out.write(' ')
    out.write('%s cal/mol.K' %Sum)
    out.write(' ')
    out.write('Original Entropy:')
    out.write(' ')
    out.write('%s cal/mol.K' %item6)
    out.write('\n')

########################################################################################

#########Array assembly from frequency.out file #####
def FixArray(a):
 count=1
 w=[]
 y=[]
 z=[]
 for x in np.nditer(a.T.copy(order='C')):
  if x <=0:
   if count <=imagfreq:
    count+=1
    x2=float(x)
    y.append(x2)  
   else:    
    y.append(MF) 
    w.append(MF) 
    z.append(x)
  elif x > 0 and x <=MF:
   MF2=float(MF)
   y.append(MF2)
   x2=float(x)
   w.append(MF)
   z.append(x2)
  else:
   x2=float(x)
   y.append(x2)
## Full array with all frequencyies and replacements organizing ##
 b=np.asarray(y)
 c=np.vstack([a,b])
## Final replaced values array organizing ##
 i=np.asarray(z)
 n=np.asarray(w)
 n=np.swapaxes(n,0,1)
 j=np.vstack([i,n])
## Final Arrays ## 
 l=np.swapaxes(j,0,1) ### Replaced Values output in this organized array (fortran format, no reason could make contiguous)
 k=np.swapaxes(c,0,1)  #### all value output in the organized array (fortran format) fortran looks pretty!
# a, original frequency array
# b, corrected frequency array
# l, ONLY replaced values output in two column format, ########## need to have this output ########
# k, all values output for comparison in two column format #debugging#
 return (a,b,l,k)

### Math operations on lines in corrected frequency file ###
def WkArray(a):
 y=[]
 for x in np.nditer(a.T.copy(order='C')):
  Theta=(x*C_c*100*(C_h/C_kb)) 
  n=abs(Theta/Temp) #using the absolute value of the frequency (THIS MEANS IMAGINARY FREQUENCY IS TREATED POSITIVE)
  d=math.exp(n)-1
  n2=math.exp(-1*n)
  LN=math.log((1-n2))
  END=(n/d-LN)
  y.append(END)
 ea=np.asarray(y) #Array of individual contributions from specific vibrations to the total vibrational entropy
 Sum=np.sum(ea) 
 VibEntrpy=Sum*C_ug
 return VibEntrpy


########################################################################################
######File Handling - Ensure program operates on all available XYZ files in directory###
########################################################################################
 
##checks number of isotope loops in the calculation  (ONLY WORKS ON MAX 2 Loops) ### 
def iso(file1):
 with open(file1, 'r') as readfile:
  lc=readfile.read()
  if "ISOTOPES TRUE" in lc:
   return "True"
  else: return "False"

#Reads in .out files and relays to other scripts, formats output
def config(freqfile, file):
   outfile=open('frequency_corrected.txt', "a")
   outfile2=open('frequency_SI_values.txt', "a")
   CheckIso=iso(freqfile)
   if CheckIso=="True":
    FreqObjArray=FreqArray2(freqfile) 
    (a,b,l,k)=FixArray(FreqObjArray) 
# a, original frequency array
# b, corrected frequency array
# l, ONLY replaced values output in two column format, ########## need to have this output ########
# k, all values output for comparison in two column format #debugging#
    c = str(WkArray(b))
    OutputConfigure2(freqfile, file, outfile, c)
    outfile2.write(file)
    outfile2.write('\n')
    np.savetxt(outfile2, l, fmt='%.3f')
    outfile2.write('\n')
   if CheckIso=="False":
    frequencyarray=FreqArray(freqfile) 
    (a,b,l,k)=FixArray(frequencyarray) 
    c = str(WkArray(b))
    OutputConfigure1(freqfile, file, outfile, c)
    outfile2.write(file)
    outfile2.write('\n')
    np.savetxt(outfile2, l, fmt='%.3f')
    outfile2.write('\n')
            
def failureCheck(dir):
 lst=os.listdir(dir)
 lst.sort()
 for file in lst: 
  if file.endswith(".out"):
   freqfile=(os.path.join(dir, file))
   outfile=open('frequency_corrected.txt', "a")
   with open(freqfile, 'r') as a:
    lc=a.read()
    if  "Thank you very much for using Q-Chem.  Have a nice day." in lc:
     config(freqfile, file)
    else: outfile.write('%s : Job Failed to Complete - No Data \n' %file)
     




failureCheck(directory)
