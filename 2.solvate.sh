#!/bin/bash

# By Yasmin 2021-06-24 
# Tested locally 2021-07-07

# This script solvates ligands in a chosen solvent for GROMACS simulations.

# Requires GROMACS installation.

############################# EDIT THESE PARAMETERS #########################################

FF=GAFF-Tested                                # default force field: GAFF-Yasmin

ligandDirectory=LIGANDS                     # Path to the ligand directory
solventDirectory=SOLVENTS                   # Path to directory of solvents

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

boxtype=cubic                               # default solvent box: cubic
boxsize=2.0                                 # default radius (in nm): 2.0
noligands=1                                 # default number of ligands: 1

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
            nosolvent=$(($totalresidues - $noligands))              # Get number of solvents
            rm tmp*.txt                                             # Clean up

# Create master topology (.top) files 

            # Get the header until atom types
            sed -n -r '/MOL_GMX/,/Amb/p' $ligname"_GMX.top" > header.tmp    

            # Get the atomtypes of ligand and solvent and remove duplicates
            cp $workingDirectory/$solventDirectory/$FF/$solvents/$solvname.itp .    # Copy solvent file to ligand directory
            sed -n -r '/atomtypes/,/^\s*$/p' $solvname.itp > solvatomtypes.tmp    # Get atomtypes of solvent
            sed -n -r '/atomtypes/,/^\s*$/p' $ligname"_GMX.top" > ligatomtypes.tmp    # Get atomtypes of ligand
            sed -i '1,2d;$d' *atomtypes.tmp                                     # Remove first and last (empty) lines
            cat solvatomtypes.tmp ligatomtypes.tmp > complexatomtypes.tmp    # Put it together
            sort complexatomtypes.tmp | uniq > complex_sorted.tmp                   # Remove duplicate lines
            echo -e "" >> complex_sorted.tmp                                   # Adds blank line at the end
            
            # Get ligand parameters and include solvent parameters
            sed -n -r '/moleculetype/,/system/p' $ligname"_GMX.top" > ligparams.tmp # Get the ligand parameters
            sed -i '$ d' ligparams.tmp                                     # Remove first and last (empty) lines
            echo -e "#include \"../$solvname.itp\"\n" >> ligparams.tmp # Add solvent parameters

            sed -n -r '/system/,//p' $ligname"_GMX.top" > system.tmp # Get the system name parameters
            echo -e " SOL              $nosolvent" >> system.tmp # Add the number of solvent molecules

            # Put everything together into top files (charged and no-charge)
            cat header.tmp complex_sorted.tmp ligparams.tmp system.tmp > ${ligname}_${solvents}.top # Create charged top-file
            sed -i 's/MOL/LIG/g' ${ligname}_${solvents}.top     # Change all instances of MOL to LIG

            cp ${ligname}_${solvents}.top ${ligname}_${solvents}"_0q.top" # Make no-charge version
            sed -i "s/$solvname/$solvname"_0q"/g" ${ligname}_${solvents}"_0q.top"   # Change path to no-charge solvent itp  
            
# Create parameter files for charged and uncharged solvents (.itp)
            
            # Create a solvent itp file from file copied from source for importing
            sed -i '/moleculetype/,$!d' $solvname.itp                   # Make the charged solvent itp file

            # Create a no-charge version for the atom types
            sed -n -r '/moleculetype/,/mass/p' $solvname.itp > header_0q.tmp
            sed -n -r '/mass/,/^\s*$/p' $solvname.itp > solvatoms.tmp    # Get atoms section of solvent
            sed -i '1d;$d' solvatoms.tmp                                     # Remove first line and blank last line
            awk '{print "     "$1"         "$2"      "$3"    "$4"     "$5"      "$6"    0.00000  "$8}' solvatoms.tmp  > solvatoms_0q.tmp
            echo -e "" >> solvatoms_0q.tmp                                   # Adds blank line at the end

            sed -n -r '/bonds/,//p' $solvname.itp > end.tmp         # Get the system name parameters                      
            cat header_0q.tmp solvatoms_0q.tmp end.tmp > $solvname"_0q.itp"    # Create no-charge solvent itp file

# Clean up the folder

            rm *.tmp                                            
            wait
            cd $workingDirectory/$solventDirectory/$FF/
        done
    done
    cd $workingDirectory/$ligandDirectory/$ligname/
     mkdir GMXPREP
     mv $ligname"_GMX".* GMXPREP/
     mv $ligname.gro GMXPREP/
     cd $workingDirectory/$ligandDirectory
done
