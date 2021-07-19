#!/bin/bash

# By Yasmin 2021-06-24 
# 2021-07-07 Tested locally 
# 2021-07-18 Updated to remove SCRIPTS directory

# This script creates ligand AMBER format toppology and lib files. 
# AMBER top and lib files are then converted to GROMACS gro and top files.

# Requires an AMBER antechamber installation to get AMBER files.
# Requires acpype.py for conversion to GROMACS formatted files.

####################### EDIT ONLY THESE PARAMETERS #########################################

AMBERPATH=/home/groups/sboxer/amberV100/amber.sh     # Path for the Amber installation (locally on Yasmin: /home/yasmin/amber20/amber.sh)
ligandDirectory=LIGANDS                     # Path to the ligand directory
pythonVersion=python                       # Installed python version

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

chargeModel=bcc                             # default: bcc = AM1-BCC
verbosity=2                                 # default: 2
FFPATH=/home/groups/sboxer/amber.sota/dat/leap/cmd/oldff/leaprc.ff99SBildn   # default: ~/amber20_src/dat/leap/cmd/oldff/leaprc.ff99SBildn
ligFF=leaprc.gaff                           # default: leaprc.gaff

######################## NO MORE EDITING PAST THIS LINE! ###################################

# Load modules on Sherlock
ml reset
ml gcc/8.1.0 cuda/10.2.89 openmpi/4.0.3

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

 # Convert to GROMACS files using acpype.py (uses GAFF parameters and ff99SB forcefield)

 for i in *.prmtop
 do
   ligname=$( echo "$i" | sed -e 's/\.prmtop//g')
   $pythonVersion $workingDirectory/acpype.py -p $workingDirectory/$ligandDirectory/$ligname.prmtop -x $workingDirectory/$ligandDirectory/$ligname.inpcrd
 wait

 # Clean up the filenames and ligand names in files

   mv MOL_GMX.top $ligname"_GMX.top"
   mv MOL_GMX.gro $ligname"_GMX.gro"

 # Clean up the folder

   mkdir $ligname
   mv $ligname"_GMX.gro" $ligname/
   mv $ligname"_GMX.top" $ligname/

   mkdir $ligname/AMBER
   mv $ligname.* $ligname/AMBER
   mv *.$ligname.in $ligname/AMBER
 done

   rm *.*
 
 ###### DONE! ##########

