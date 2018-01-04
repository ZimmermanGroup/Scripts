#/bin/bash

# splits an stringfile.xyz into its xyz parts
# still need to make sure there are no realignments in molden



# first arg
file=$1
# second arg
#nnodes=$2
 
numlines=`cat $file | wc -l`
echo " lines: $numlines"


#echo "working on $file"

# csplit operates on file, splits into 2 files (xx00 and xx01) at line 2
csplit $file 2

# cats the first line, and stores the number of atoms
natoms=`cat xx00`

#echo " natoms: $natoms"

# natomsp2 is the length of geometry
let natomsp2=natoms+2
echo " natomsp2: $natomsp2"

let nnodes=numlines/natomsp2-2
echo " nnodes: $nnodes"

#echo " natomsp2: $natomsp2"

i=0

let a=natomsp2+1
while [ $i -le $nnodes ]
do
  echo "list:$list"
  list="$list $a"
  let a=a+natomsp2  

  let i=i+1
done

#exit

echo "list: $list"

csplit $file $list 


i=0

let nnodes_p1=nnodes+1

while [ $i -le $nnodes_p1 ]
do
  ID=`printf "%0*d\n" 2 ${i}`

  echo " working on xx$ID"
  python ~/bin/removeH.py xx$ID
  python ~/bin/distance.py xx$ID 

  let i=i+1
done

cat xx*.flat > $file.a
#rm xx*
rm noH.*
