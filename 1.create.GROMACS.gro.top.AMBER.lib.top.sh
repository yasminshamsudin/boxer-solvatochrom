#!/bin/bash

# By Yasmin 2021-06-24 
# Tested locally 2021-06-24

# This script creates ligand AMBER format toppology and lib files. 
# AMBER top and lib files are then converted to GROMACS gro and top files.

# Requires an AMBER antechamber installation to get AMBER files.
# Requires acpype.py for conversion to GROMACS formatted files.

####################### EDIT ONLY THESE PARAMETERS #########################################

AMBERPATH=/home/yasmin/amber20/amber.sh     # Path for the Amber installation
ligandDirectory=LIGANDS                     # Path to the ligand directory
scriptDirectory=SCRIPTS                     # Path to python scripts
pythonVersion=python3                       # Installed python version
workingDirectory=/mnt/c/Users/Yasmin/Desktop/Test

############ Advanced modelling. Only edit if you know what you are doing! #################

chargeModel=bcc         # default: bcc = AM1-BCC
verbosity=2             # default: 2
FFPATH=/home/yasmin/amber20_src/dat/leap/cmd/oldff/leaprc.ff99SBildn   # default: ~/amber20_src/dat/leap/cmd/oldff/leaprc.ff99SBildn
ligFF=leaprc.gaff       # default: leaprc.gaff

######################## NO MORE EDITING PAST THIS LINE! ###################################

# Set the AMBER path to find all modules

source $AMBERPATH

# Convert ligand pdb to mol2 using antechamber

cd $ligandDirectory

for i in *.pdb
 do
    ligname=$( echo "$i" | sed -e 's/\.pdb//g')
     antechamber -i $ligname.pdb -fi pdb -o $ligname.mol2 -fo mol2 -c $chargeModel -s $verbosity
 wait
 done

 # Check the parameters using parmchk2

for i in *.mol2
 do
    ligname=$( echo "$i" | sed -e 's/\.mol2//g')
     parmchk2 -i $ligname.mol2 -f mol2 -o $ligname.frcmod
 wait
 done

 # Create the .lib and .prmtop files using tleap

 for i in *.frcmod
 do
    ligname=$( echo "$i" | sed -e 's/\.frcmod//g')
     echo -e "source $FFPATH"  > tleap.$ligname.in
     echo -e "source $ligFF" >> tleap.$ligname.in
     echo -e "$ligname=loadmol2 $ligname.mol2" >> tleap.$ligname.in
     echo -e "loadamberparams $ligname.frcmod" >> tleap.$ligname.in
     echo -e "check $ligname" >> tleap.$ligname.in
     echo -e "saveoff $ligname $ligname.lib" >> tleap.$ligname.in
     echo -e "saveamberparm $ligname $ligname.prmtop $ligname.inpcrd" >> tleap.$ligname.in
     echo -e "Quit" >> tleap.$ligname.in

    tleap -f tleap.$ligname.in
 wait
 done

 # Convert to GROMACS files using acpype.py

 for i in *.prmtop
 do
   ligname=$( echo "$i" | sed -e 's/\.prmtop//g')
   $pythonVersion $workingDirectory/$scriptDirectory/acpype.py -p $workingDirectory/$ligandDirectory/$ligname.prmtop -x $workingDirectory/$ligandDirectory/$ligname.inpcrd
 wait
 done

 ###### DONE! ##########

