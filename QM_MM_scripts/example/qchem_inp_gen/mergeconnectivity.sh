#!/bin/bash

txyzFile=$1
xyzFile=$2

awk '{ print $6"\t"$7"\t"$8"\t"$9"\t"$10 }' $txyzFile >> tmp.txt
nl tmp.txt > connectivity.txt
rm tmp.txt
echo "Exctracting connectivity from: " "$txyzFile"

newName=${xyzFile/.xyz/_merged}.xyz
nl $xyzFile > coords.txt
echo "Exctracting coords from: " "$xyzFile"

join coords.txt connectivity.txt > $newName
echo "Merging xyz and txyz to: " "$newName"
