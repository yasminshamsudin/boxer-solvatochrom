#!/bin/bash

# By Yasmin 2021-07-02 
# 2021-07-02 Tested locally 
# 2021-07-05 Tested on Sherlock 
# 2021-07-14 Added capability to process several replicates
# 2021-07-14 Added capability to merge several ligand runs into one file
# 2021-07-14 Tested on Sherlock 

# This script analyzes electric field calculations from GROMACS simulations. 

############################# EDIT THESE PARAMETERS #########################################

date=$( date +%y%m%d)                       # Today's date
identifier=$date-Triplicate                 # Unique identifiers for the group of simulations
ligandDirectory=LIGANDS                     # Path to the ligand directory

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

######################## NO MORE EDITING PAST THIS LINE! ###################################

echo -e "Results of $identifier" > $identifier.header.tmp
cp $identifier.header.tmp summary.tmp

cd $ligandDirectory

# Repeat for all ligands
for ligname in *
 do
    cd $ligname
    echo -e "Solvent <Field> Standard_deviation Removed_lines" > $ligname.header.tmp
    cp $ligname.header.tmp combinedfields.tmp
    echo -e "${ligname}_Solvent <Field> Standard_error <StdDev> Data_points" > $ligname.header.tmp
    cp $ligname.header.tmp final.tmp



# Repeat for all solvents
        
        # Get the solvent name
        for solvents in *_0q.itp
	      do
           solvent=$( echo "$solvents" | sed -e 's/_0q\.itp//g')
           
           # Initiate for calculating averages  and SEM over replicates
            echo -e "" > $solvent.replicate.tmp
            cp $solvent.replicate.tmp $solvent.stdDev.replicate.tmp
            cp $solvent.replicate.tmp $solvent.counted.tmp
            
           # Initiate by making the first column with the solvent name
           echo -e "$solvent" > $ligname.$solvent.name.tmp
           cp $ligname.$solvent.name.tmp combinedsolvents.tmp

            # Repeat for all replicates
            for solvreplica in ${ligname}_${solvent}*
            do

            cd $solvreplica

            # Remove lines where probe is at the edge of the box and count the number of removed lines
            awk '$6 <= 1' BLCO_md.txt > counted_lines.txt
            awk '$6 > 1' BLCO_md.txt > removed_lines.txt
            awk 'END{print NR}' removed_lines.txt > ../$ligname.$solvent.noremovedlines.tmp
            awk 'END{print NR}' counted_lines.txt > ../$ligname.$solvent.nocountedlines.tmp
            head -1 ../$ligname.$solvent.nocountedlines.tmp >> ../$solvent.counted.tmp # Output the value to calculate average of averages

            # Calculate the average fields
            awk '{ total += $4 } END { print total/NR }' counted_lines.txt > ../$ligname.$solvent.aveField.tmp
            head -1 ../$ligname.$solvent.aveField.tmp >> ../$solvent.replicate.tmp # Output the value to calculate average of averages

            # Calculate the standard deviation
            awk '{x+=$4;y+=$4^2}END{print sqrt(y/NR-(x/NR)^2)}' counted_lines.txt > ../$ligname.$solvent.stdDev.tmp
            head -1 ../$ligname.$solvent.stdDev.tmp >> ../$solvent.stdDev.replicate.tmp # Output the value to calculate average of stDev

            # Put all columns together
            cd ..
            paste $ligname.$solvent.aveField.tmp $ligname.$solvent.stdDev.tmp $ligname.$solvent.noremovedlines.tmp > $solvreplica.tmp

            # Put all replicates of the same solvent together
            paste combinedsolvents.tmp $solvreplica.tmp > combination.tmp
            cp combination.tmp combinedsolvents.tmp
            echo -e "$replicate" > noreplicate.tmp

            done

            cp combinedsolvents.tmp $ligname.$solvent.tmp

            # Calculate average of replicates
            sed -i '1d' $solvent.replicate.tmp          # Remove the top blank line
            awk '{ total += $1 } END { print total/NR }' $solvent.replicate.tmp > $solvent.replicate.aveField.tmp

            sed -i '1d' $solvent.stdDev.replicate.tmp          # Remove the top blank line
            awk '{ total += $1 } END { print total/NR }' $solvent.stdDev.replicate.tmp > $solvent.replicate.aveStdDev.tmp

            # Calculate the standard error of the replicates
            awk '{x+=$1;y+=$1^2}END{print (sqrt(y/NR-(x/NR)^2))/sqrt(NR)}' $solvent.replicate.tmp > $solvent.replicate.stdErr.tmp

            # Count the number of data points
            awk '{SUM+=$1}END{print SUM}' $solvent.counted.tmp > $solvent.replicate.counted.tmp

            # Paste the average and standard error into separate file
            paste $ligname.$solvent.name.tmp $solvent.replicate.aveField.tmp $solvent.replicate.stdErr.tmp $solvent.replicate.aveStdDev.tmp $solvent.replicate.counted.tmp > $ligname.stdErr.tmp
            cat final.tmp $ligname.stdErr.tmp > combined.final.tmp
            cp combined.final.tmp final.tmp
            
            # Put all solvents together
            cat combinedfields.tmp $ligname.$solvent.tmp > combined.tmp
            cp combined.tmp combinedfields.tmp

        done

            cp combinedfields.tmp $ligname.fields.log
            cp final.tmp $ligname.aveFields.log
            cat $workingDirectory/summary.tmp $ligname.aveFields.log > $workingDirectory/combinedave.tmp
            cp $workingDirectory/combinedave.tmp $workingDirectory/summary.tmp
# Tidy up 
        rm *.tmp

    cd $workingDirectory/$ligandDirectory/


done
cd $workingDirectory
mv summary.tmp $identifier.log
rm *.tmp
