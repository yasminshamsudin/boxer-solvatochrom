#!/bin/bash

# By Yasmin 2021-07-02 
# Tested locally 2021-07-02
# Tested on Sherlock 2021-07-05

# This script prepares job files for GROMACS simulations. 
# This script submits files on the Stanford Sherlock cluster.
# The submit file also runs a script for calculating electric fields. 

# Requires access on Stanford's Sherlock cluster. Can be modified for other clusters.

############################# EDIT THESE PARAMETERS #########################################

ligandDirectory=LIGANDS                     # Path to the ligand directory
solventDirectory=SOLVENTS                   # Path to directory of solvents
FF=GAFF-Test                                # default force field: GAFF-Yasmin

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

######################## NO MORE EDITING PAST THIS LINE! ###################################

cd $ligandDirectory

# Repeat for all ligands
for ligname in *
 do
    cd $ligname
    echo -e "Solvent Average Field Standard deviation" > $ligname.header.tmp
    cp $ligname.header.tmp combinedfields.tmp

# Repeat for all solvents

        cd $workingDirectory/$solventDirectory/$FF/
        
        # Get the solvent name
        for solvent in *
	 do
           cd $workingDirectory/$ligandDirectory/$ligname/
           echo -e "$solvent" > $ligname.$solvent.name.tmp

        # Calculate the average fields
           cd ${ligname}_${solvent}
           awk '{ total += $4 } END { print total/NR }' BLCO_md.txt > ../$ligname.$solvent.aveField.tmp

        # Calculate the standard deviation
           awk '{x+=$4;y+=$4^2}END{print sqrt(y/NR-(x/NR)^2)}' BLCO_md.txt > ../$ligname.$solvent.stdDev.tmp

        # Put it all together
        cd ..
        paste $ligname.$solvent.name.tmp $ligname.$solvent.aveField.tmp $ligname.$solvent.stdDev.tmp > $ligname.$solvent.tmp
        cat combinedfields.tmp $ligname.$solvent.tmp > combined.tmp
        cp combined.tmp combinedfields.tmp
        done

        cp combinedfields.tmp $ligname.fields.log
        rm *.tmp

    cd ..
done
