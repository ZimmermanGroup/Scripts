#!/bin/bash

rm lgatoms.rtf

filename=$1
lgFile=$2
atomsList=`cat $filename`

for atom in $atomsList

do
	echo $atom
	line=`awk '/^'$atom'/' $lgFile | awk '{print $1,$2,$4,$3}'`
	echo $line >> lgatoms.rtf
done

#rm temp.prm
