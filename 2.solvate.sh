#!/bin/bash

# By Yasmin 2021-06-24 
# Tested locally 2021-06-24????

# This script solvates ligands in a chosen solvent for GROMACS simulations.

# Requires GROMACS installation.

############################# EDIT THESE PARAMETERS #########################################

workingDirectory=/mnt/c/Users/Yasmin/Desktop/Test # Path to where this script is located
ligandDirectory=LIGANDS                     # Path to the ligand directory
solventDirectory=SOLVENTS                   # Path to directory of solvents
scriptDirectory=SCRIPTS                     # Path to python scripts
pythonVersion=python3                       # Installed python version

############ Advanced modelling. Only edit if you know what you are doing! #################

FF=GAFF-Test         # default force field: GAFF-Yasmin
boxtype=cubic        # default solvent box: cubic
boxsize=2.0          # default radius (in nm): 2.0
nosolutes=1          # default number of ligands: 1

######################## NO MORE EDITING PAST THIS LINE! ###################################

cd $ligandDirectory

# Repeat for all ligands
for ligname in *
 do

# Center the solute in the box

            gmx editconf -f $ligname/$ligname"_GMX.gro" -o $ligname/$ligname.gro -bt $boxtype -d $boxsize -c
            wait

# Repeat for all solvents in the SOLVENTS folder

    cd $workingDirectory/$solventDirectory/$FF/

    for solvents in *   # Get the name of each solvent
     do
      cd $solvents
        for i in *.pdb
            do
            solvname=$( echo "$i" | sed -e 's/\.pdb//g')
            cd $workingDirectory/$ligandDirectory/$ligname/

# Solvate the box with solvent

            gmx solvate -cp $ligname.gro -cs $workingDirectory/$solventDirectory/$FF/$solvents/$solvname.pdb -o ${ligname}_${solvents}.gro
            wait

# Fetch the number of solvent molecules and store as a variable

            sed 'x;$!d' ${ligname}_${solvents}.gro > tmp.txt        # Fetch the second last row
            awk '{print $1}' tmp.txt > tmp2.txt                     # Get the first column
            totalresidues=$( grep -o -E '[0-9]+' tmp2.txt)          # Separate the residue number
            nosolvent=$(($totalresidues - $nosolutes))              # Get number of solvents
            rm tmp*.txt                                             # Clean up

# Create master topology (.top) file and parameter file for solvent (.itp) 

            cp $ligname"_GMX.top" ${ligname}_${solvents}.top        # Create a master topology file
            cp $workingDirectory/$solventDirectory/$FF/$solvents/$solvname.itp .    # Copy solvent file to ligand directory
            sed -n -r '/atomtypes/,/^\s*$/p' $solvname.itp > solvatomtypes.tmp    # Get atomtypes of solvent
            sed -n -r '/atomtypes/,/^\s*$/p' ${ligname}_${solvents}.top > ligatomtypes.tmp    # Get atomtypes of solvent
            sed -i '1,2d;$d' *.tmp                                     # Remove first and last (empty) lines
            cat solvatomtypes.tmp ligatomtypes.tmp > complex.tmp    # Put it together
            awk '!_[$0]++' complex.tmp                              # Remove duplicate lines
            wait
        done
     wait
    done
done
