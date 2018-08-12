#!/bin/bash
#This file was created on 08/05/2018
export OMP_NUM_THREADS=1

bases="c1ccccc1CCCC=C" # c1ccccc1COCCC=C c1ccccc1CCOC=C c1ccccc1CCCOC=C c1ccccc1OCCOC=C"
groups="C(Cl)=O" # C=O C(C)=O N#C C(Cl)(Cl)C N(=O)=O OC O OC(=O)C C N NC(OC)=O NC N(C)C C=C"

count=1
for base in $bases;do
    count2=1
     # => Do substrates with no groups <= #
    ./tethered_example.py --substrate=$base --group="" --sub_id=$count --group_id=$count2 --position=99 

    #calculate all the ISOMERS
    dir=`cat dirname`
    #TODO

    # => Get path to reference <= #
    # => Reference is unsubstituted mol <=#
    dir=`cat dirname`
    path_to_ref=`echo "$(pwd -P)/$(basename "$dir")/tmp.xyz"`
    echo $path_to_ref

     # => Do substrates with ortho groups <= #
    count2=1
    for group in $groups;do 
        ./tethered_example.py --substrate=$base --group=$group --sub_id=$count --group_id=$count2 --position=1 --ref=$path_to_ref
        #calculate all the ISOMERS
        dir=`cat dirname`
        #TODO
        
        let count2++
    done

     # => Do substrates with meta groups <= #
    count2=1
    for group in $groups;do 
        ./tethered_example.py --substrate=$base --group=$group --sub_id=$count --group_id=$count2 --position=2 --ref=$path_to_ref

        #calculate all the ISOMERS
        let count2++
    done

     # => Do substrates with para groups <= #
    count2=1
    for group in $groups;do 
        ./tethered_example.py --substrate=$base --group=$group --sub_id=$count --group_id=$count2 --position=3 --ref=$path_to_ref

        #calculate all the ISOMERS
        dir=`cat dirname`
        let count2++
    done
    cd ..
    let count++
done

