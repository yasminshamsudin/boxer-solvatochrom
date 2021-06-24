#!/bin/bash

# By Yasmin 2021-06-24 
# Tested locally 2021-06-24????

# This script solvates the ligand in a chosen solvent for GROMACS simulations.

# Requires GROMACS installation.

####################### EDIT ONLY THESE PARAMETERS #########################################

workingDirectory=/mnt/c/Users/Yasmin/Desktop/Test # Path to where this script is located
ligandDirectory=LIGANDS                     # Path to the ligand directory
solventDirectory=SOLVENTS                   # Path to directory of solvents
scriptDirectory=SCRIPTS                     # Path to python scripts
pythonVersion=python3                       # Installed python version

############ Advanced modelling. Only edit if you know what you are doing! #################

FF=GAFF-Test         # default: GAFF-Yasmin
boxtype=cubic        # default: cubic
boxsize=2.0          # default (in nm): 2.0

######################## NO MORE EDITING PAST THIS LINE! ###################################

cd $ligandDirectory

# Repeat for all ligands
for ligname in *
 do
    cd $workingDirectory/$solventDirectory/$FF/

# Repeat for all solvents in the SOLVENTS folder

    for solvents in *   # Get the name of each solvent
     do
      cd $solvents
        for i in *.pdb
            do
            solvname=$( echo "$i" | sed -e 's/\.pdb//g')
            cd $workingDirectory/$ligandDirectory/$ligname/

# Center the solute in the box

            gmx editconf -f $ligname"_GMX.gro" -o $ligname.gro -bt $boxtype -d $boxsize -c

# Solvate the box with solvent

            gmx solvate -cp $ligname.gro -cs $workingDirectory/$solventDirectory/$FF/$solvents/$solvname.pdb -o ${ligname}_${solvents}.gro
            wait
        done
     wait
    done
done