#!/bin/bash

paramFile=$1

sed -i '/^bond/d' $paramFile
sed -i '/^angle/d' $paramFile
sed -i '/^torsion/d' $paramFile
sed -i '/^ureybrad/d' $paramFile


sort -u bonds.prm >> $paramFile
#echo " " >> $paramFile
sort -u angles.prm >> $paramFile
#echo " " >> $paramFile
sort -u torsions.prm >> $paramFile
#echo " " >> $paramFile
sort -u ureybrad.prm >> $paramFile
