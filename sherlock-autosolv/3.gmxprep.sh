#!/bin/bash

# By Yasmin 2021-07-01 
# Tested locally 2021-07-02
# Tested on Sherlock 2021-07-07

# This script prepares input files for GROMACS simulations and analysis scripts.

# Requires probe atoms to be in consecutive lines in the original pdb file.
# If ligand has more than two possible instances of probe atoms (e.g two carbonyls)
# only the first instance will be assigned. Make sure the probe atoms of interest
# appear before the second instance in the initial ligand structure file.

############################# EDIT THESE PARAMETERS #########################################

probeAtom1=C                                 # Default: C (For other constructs, manually set probeatomno1)
probeAtom2=O                                 # Default: O (Could be set to N, D, or other)
probeAtomno1=0                               # Default for automated C atomnumber allocation: 0
probeAtomno2=0                               # Default for automated atomnumber allocation: 0

ligandDirectory=LIGANDS                     # Path to the ligand directory
solventDirectory=SOLVENTS                   # Path to directory of solvents
scriptDirectory=SCRIPTS                     # Path to python scripts
pythonVersion=python3                       # Installed python version

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

nptsteps=200000                             # Default: 100000 (stepsize (1fs) * runtime (100K) = 100ps
mdsteps=500000                              # Default: 500000 (stepsize (2fs) * runtime (500K) = 1ns

######################## NO MORE EDITING PAST THIS LINE! ###################################

cd $ligandDirectory

# Repeat for all ligands
for ligname in *
 do
    cd $ligname

# Create the probe atoms list file (.ndx) and the python script for calculating electric fields
        
         cd GMXPREP/
         for topfile in *.top
            do
             # Get atom list of ligand
             sed -n -r '/atoms/,/^\s*$/p' $topfile > ligatomtypes.tmp    
                   
             # For C=O probes 
             if [[ $probeAtom1=C && $probeAtom2=O ]]
                then
                   # For automated assignment of C=O probe atom numbers (only finds first occurrence)
                   if [[ $probeAtomno1 = 0 && $probeAtomno2 = 0 ]]
                        then
                        # Extract the lines of the probe atoms
                        sed '/ c /,/ o /!d;/ o /q' ligatomtypes.tmp > probes.tmp
                        
                        # Get the probe atom numbers
                        awk '{print $1}' probes.tmp > probeno.tmp
                        carbonno=$( head -1 probeno.tmp )
                        oxygenno=$( tail -1 probeno.tmp )

                        # Get the probe atom charges
                        awk '{print $7}' probes.tmp > charges.tmp
                        carboncharge=$( head -1 charges.tmp )
                        oxygencharge=$( tail -1 charges.tmp )
                        
                        else # If user has specified the probe atoms
                        carbonno=$(($probeAtomno1))
                        oxygenno=$(($probeAtomno2))

                        # Get the probe atom charges
                        grep "\s$carbono\s" ligatomtypes.tmp > cprobe.tmp
                        grep "\s$oxygenno\s" ligatomtypes.tmp > oprobe.tmp
                        cat cprobe.tmp oprobe.tmp > probes.tmp
                        awk '{print $7}' probes.tmp > charges.tmp
                        carboncharge=$( head -1 charges.tmp )
                        oxygencharge=$( tail -1 charges.tmp )
                        fi

                else # If not C=O probe
                carbonno=$(($probeAtomno1))
                oxygenno=$(($probeAtomno2))

                # Get the probe atom charges
                grep "\s$carbono\s" ligatomtypes.tmp > cprobe.tmp
                grep "\s$oxygenno\s" ligatomtypes.tmp > oprobe.tmp
                cat cprobe.tmp oprobe.tmp > probes.tmp
                awk '{print $7}' probes.tmp > charges.tmp
                carboncharge=$( head -1 charges.tmp )
                oxygencharge=$( tail -1 charges.tmp )

                fi
                rm *.tmp
        done
        cd ..
        
        # Make the probe index file (.ndx) (only one per ligand needed)
        echo -e "[$probeAtom1"-"$probeAtom2]"'\n' > probe.ndx
        echo -e "$carbonno $oxygenno" >> probe.ndx

        # Make the Electric field calculation script
        echo -e "# The purpose of this script is to calculate electric fields for a solute in solvent with carbonyl vibrational probe." > BLCO.py
        echo -e "# This file has been automatically generated."'\n' >> BLCO.py
        echo -e "import sys" >> BLCO.py
        echo -e "import re" >> BLCO.py
        echo -e "import itertools" >> BLCO.py
        echo -e "import os" >> BLCO.py
        echo -e "import numpy as np"'\n' >> BLCO.py
        echo -e "def isnumber(string):" >> BLCO.py
        echo -e '\t'"try:" >> BLCO.py
        echo -e '\t\t'"float(string)" >> BLCO.py
        echo -e '\t\t'"return True" >> BLCO.py
        echo -e '\t'"except ValueError:" >> BLCO.py
        echo -e '\t\t'"return False"'\n' >> BLCO.py
        echo -e "# Get the relevant forces and coords from trajectories"'\n' >> BLCO.py
        echo -e "f_x = open('co_coords.xvg' , 'r')" >> BLCO.py
        echo -e "f_f = open('co_forces.xvg' , 'r')" >> BLCO.py
        echo -e "f_f0q = open ('co_forces_0q.xvg'  , 'r')"'\n' >> BLCO.py
        echo -e "# Calculate the fields" >> BLCO.py
        echo -e "fields = open('BLCO_md.txt'  , 'w')" >> BLCO.py
        echo -e "for line_x, line_f, line_f0q in itertools.izip(f_x, f_f, f_f0q):" >> BLCO.py
        echo -e '\t'"if isnumber( line_x.split()[0] ):" >> BLCO.py
        echo -e '\t\t'"#Get the time" >> BLCO.py
        echo -e '\t\t'"time_x = float( line_x.split()[0] )" >> BLCO.py
        echo -e '\t\t'"time_f0q = float( line_f0q.split()[0] )"'\n' >> BLCO.py
        echo -e '\t\t'"#Check the info" >> BLCO.py
        echo -e '\t\t'"if time_x != time_f0q:" >> BLCO.py
        echo -e '\t\t\t'"print \"time indeces do not match. Error!\"" >> BLCO.py
        echo -e '\t\t\t'"sys.exit()"'\n' >> BLCO.py
        echo -e '\t\t'"#Get the C and O coordinates" >> BLCO.py
        echo -e '\t\t'"#Calculate the CO bond length and CO unit vector" >> BLCO.py
        echo -e '\t\t'"x_C = np.array([float(x) for x in line_x.split()[1:4]])" >> BLCO.py
        echo -e '\t\t'"x_O = np.array([float(x) for x in line_x.split()[4:7]])" >> BLCO.py
        echo -e '\t\t'"COvec = x_O - x_C" >> BLCO.py
        echo -e '\t\t'"COlen = np.sqrt((COvec**2).sum())" >> BLCO.py
        echo -e '\t\t'"COunitvec = COvec/ COlen"'\n' >> BLCO.py
        echo -e '\t\t'"#Get the C and O forces" >> BLCO.py
        echo -e '\t\t'"f_C = np.array([float(x) for x in line_f.split()[1:4]])" >> BLCO.py
        echo -e '\t\t'"f_O = np.array([float(x) for x in line_f.split()[4:7]])"'\n' >> BLCO.py
        echo -e '\t\t'"#Get the C and O discharged forces" >> BLCO.py
        echo -e '\t\t'"f0q_C = np.array([float(x) for x in line_f0q.split()[1:4]])" >> BLCO.py
        echo -e '\t\t'"f0q_O = np.array([float(x) for x in line_f0q.split()[4:7]])"'\n' >> BLCO.py
        echo -e '\t\t'"#Calculate the electrostatic forces" >> BLCO.py
        echo -e '\t\t'"fe_C = f_C - f0q_C" >> BLCO.py
        echo -e '\t\t'"fe_O = f_O - f0q_O"'\n' >> BLCO.py
        echo -e '\t\t'"#Calculate the electrostatic force projection onto CO" >> BLCO.py
        echo -e '\t\t'"feproj_C = np.dot( fe_C, COunitvec )" >> BLCO.py
        echo -e '\t\t'"feproj_O = np.dot( fe_O, COunitvec )"'\n' >> BLCO.py
        echo -e '\t\t'"#Calculate the electric field and convert" >> BLCO.py
        echo -e '\t\t'"Eproj_C = (feproj_C /$carboncharge) * 0.1036427" >> BLCO.py
        echo -e '\t\t'"Eproj_O = (feproj_O /$oxygencharge) * 0.1036427" >> BLCO.py
        echo -e '\t\t'"Eproj = (Eproj_C + Eproj_O)/2" >> BLCO.py
        echo -e '\t\t'"EprojDrop = (Eproj_O-Eproj_C)"'\n' >> BLCO.py
        echo -e '\t\t'"#Print to file" >> BLCO.py
        echo -e '\t\t'"fieldInfo = str( time_x ) + '\\\t' +str( Eproj_C ) + '\\\t' + str( Eproj_O ) + '\\\t' + str( Eproj )+ '\\\t' + str( EprojDrop ) + '\\\t' + str( COlen ) + '\\\n'" >> BLCO.py
        echo -e '\t\t'"fields.write ( fieldInfo )"'\n' >> BLCO.py
        echo -e "f_x.close()" >> BLCO.py
        echo -e "f_f.close()" >> BLCO.py
        echo -e "f_f0q.close()" >> BLCO.py
        echo -e "fields.close()" >> BLCO.py

# Create minimization, heating, equilibration, and md run input files (.mdp)

        # Create minimization input file
        echo -e "; min.mdp - used as input into grompp to generate em.tpr" > min.mdp
        echo -e "integrator	= steep		; Algorithm (steep = steepest descent minimization)" >> min.mdp
        echo -e "emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm" >> min.mdp
        echo -e "emstep      = 0.01      ; Energy step size" >> min.mdp
        echo -e "nsteps		= 50000	  	; Maximum number of (minimization) steps to perform"'\n' >> min.mdp
        echo -e "; Parameters describing how to find the neighbors of each atom and how to calculate the interactions" >> min.mdp
        echo -e "nstlist		    = 10	    ; Frequency to update the neighbor list and long range forces" >> min.mdp
        echo -e "cutoff-scheme   = Verlet" >> min.mdp
        echo -e "ns_type		    = grid		; Method to determine neighbor list (simple, grid)" >> min.mdp
        echo -e "coulombtype	    = PME		; Treatment of long range electrostatic interactions" >> min.mdp
        echo -e "rcoulomb	    = 1.0		; Short-range electrostatic cut-off" >> min.mdp
        echo -e "rvdw		    = 1.0		; Short-range Van der Waals cut-off" >> min.mdp
        echo -e "pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)" >> min.mdp

        # Create nvt (heating) input file
        echo -e "title		= Electric Fields solvatochromism NVT equilibration" > nvt.mdp
        echo -e "define		= -DPOSRES	; position restrain the protein"'\n' >> nvt.mdp
        echo -e "; Run parameters" >> nvt.mdp
        echo -e "integrator	= sd		; leap-frog integrator" >> nvt.mdp
        echo -e "nsteps		= 100000		; 1 * 100000 = 100 ps" >> nvt.mdp
        echo -e "dt		    = 0.001		; Step-size 1 fs"'\n' >> nvt.mdp
        echo -e "; Output control" >> nvt.mdp
        echo -e "nstxout		= 1000		; save coordinates every 1.0 ps" >> nvt.mdp
        echo -e "nstvout		= 1000		; save velocities every 1.0 ps" >> nvt.mdp
        echo -e "nstenergy	= 1000		; save energies every 1.0 ps" >> nvt.mdp
        echo -e "nstlog		= 1000		; update log file every 1.0 ps"'\n' >> nvt.mdp
        echo -e "; Bond parameters" >> nvt.mdp
        echo -e "continuation	        = no		; first dynamics run" >> nvt.mdp
        echo -e "constraint_algorithm    = lincs	    ; holonomic constraints" >> nvt.mdp
        echo -e "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained" >> nvt.mdp
        echo -e "lincs_iter	            = 1		    ; accuracy of LINCS" >> nvt.mdp
        echo -e "lincs_order	            = 4		    ; also related to accuracy"'\n' >> nvt.mdp
        echo -e "; Neighborsearching" >> nvt.mdp
        echo -e "cutoff-scheme   = Verlet" >> nvt.mdp
        echo -e "ns_type		    = grid		; search neighboring grid cells" >> nvt.mdp
        echo -e "nstlist		    = 10		; 20 fs, largely irrelevant with Verlet" >> nvt.mdp
        echo -e "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)" >> nvt.mdp
        echo -e "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)" >> nvt.mdp
        echo -e "; Electrostatics" >> nvt.mdp
        echo -e "coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics" >> nvt.mdp
        echo -e "pme_order	    = 4		; cubic interpolation" >> nvt.mdp
        echo -e "fourierspacing	= 0.16	; grid spacing for FFT"'\n'  >> nvt.mdp
        echo -e "; Temperature coupling is on" >> nvt.mdp
        echo -e "tcoupl		= V-rescale	            ; modified Berendsen thermostat" >> nvt.mdp
        echo -e "tc-grps		= LIG    SOL	; two coupling groups - more accurate" >> nvt.mdp
        echo -e "tau_t		= 2	  2           ; time constant, in ps" >> nvt.mdp
        echo -e "ref_t		= 300 	  300           ; reference temperature, one for each group, in K"'\n' >> nvt.mdp
        echo -e "; Pressure coupling is off" >> nvt.mdp
        echo -e "pcoupl		= no 		; no pressure coupling in NVT"'\n' >> nvt.mdp
        echo -e "; Periodic boundary conditions" >> nvt.mdp
        echo -e "pbc		        = xyz 		; 3-D Periodic Boundary Conditions (yes/no)"'\n' >> nvt.mdp
        echo -e "; Dispersion correction" >> nvt.mdp
        echo -e "DispCorr	= EnerPres	; account for cut-off vdW scheme"'\n' >> nvt.mdp
        echo -e "; Velocity generation" >> nvt.mdp
        echo -e "gen_vel		= yes		; assign velocities from Maxwell distribution" >> nvt.mdp
        echo -e "gen_temp	= 300		; temperature for Maxwell distribution" >> nvt.mdp
        echo -e "gen_seed	= -1		; generate a random seed" >> nvt.mdp

        # Create npt (equilibration) input file
        echo -e "title		= Electric Fields solvatochromism NPT equilibration" > npt.mdp
        echo -e "define		= -DPOSRES	; position restrain the protein"'\n' >> npt.mdp
        echo -e "; Run parameters" >> npt.mdp
        echo -e "integrator	= sd		; leap-frog integrator" >> npt.mdp
        echo -e "nsteps		= $nptsteps		; 2 * 500000 = 1000 ps = 1 ns" >> npt.mdp
        echo -e "dt		    = 0.001		; Step-size 1 fs"'\n' >> npt.mdp
        echo -e "; Output control" >> npt.mdp
        echo -e "nstxout		= 1000		; save coordinates every 1.0 ps" >> npt.mdp
        echo -e "nstvout		= 1000		; save velocities every 1.0 ps" >> npt.mdp
        echo -e "nstenergy	= 1000		; save energies every 1.0 ps" >> npt.mdp
        echo -e "nstlog		= 1000		; update log file every 1.0 ps"'\n' >> npt.mdp
        echo -e "; Bond parameters" >> npt.mdp
        echo -e "continuation	        = yes		; Restarting after NVT run" >> npt.mdp
        echo -e "constraint_algorithm    = lincs	    ; holonomic constraints" >> npt.mdp
        echo -e "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained" >> npt.mdp
        echo -e "lincs_iter	            = 1		    ; accuracy of LINCS" >> npt.mdp
        echo -e "lincs_order	            = 4		    ; also related to accuracy"'\n' >> npt.mdp
        echo -e "; Neighborsearching" >> npt.mdp
        echo -e "cutoff-scheme   = Verlet" >> npt.mdp
        echo -e "ns_type		    = grid		; search neighboring grid cells" >> npt.mdp
        echo -e "nstlist		    = 10	; 20 fs, largely irrelevant with Verlet" >> npt.mdp
        echo -e "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)" >> npt.mdp
        echo -e "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)" >> npt.mdp
        echo -e "; Electrostatics" >> npt.mdp
        echo -e "coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics" >> npt.mdp
        echo -e "pme_order	    = 4		; cubic interpolation" >> npt.mdp
        echo -e "fourierspacing	= 0.16	; grid spacing for FFT"'\n'  >> npt.mdp
        echo -e "; Temperature coupling is on" >> npt.mdp
        echo -e "tcoupl		= V-rescale	            ; modified Berendsen thermostat" >> npt.mdp
        echo -e "tc-grps		= LIG    SOL	; two coupling groups - more accurate" >> npt.mdp
        echo -e "tau_t		= 2	  2           ; time constant, in ps" >> npt.mdp
        echo -e "ref_t		= 300 	  300           ; reference temperature, one for each group, in K"'\n' >> npt.mdp
        echo -e "; Pressure coupling is on" >> npt.mdp
        echo -e "pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT" >> npt.mdp
        echo -e "pcoupltype	        = isotropic	            ; uniform scaling of box vectors" >> npt.mdp
        echo -e "tau_p		        = 2.0		            ; time constant, in ps" >> npt.mdp
        echo -e "ref_p		        = 1.0		            ; reference pressure, in bar" >> npt.mdp
        echo -e "compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1" >> npt.mdp
        echo -e "refcoord_scaling    = com"'\n' >> npt.mdp
        echo -e "; Periodic boundary conditions" >> npt.mdp
        echo -e "pbc		        = xyz 		; 3-D Periodic Boundary Conditions (yes/no)"'\n' >> npt.mdp
        echo -e "; Dispersion correction" >> npt.mdp
        echo -e "DispCorr	= EnerPres	; account for cut-off vdW scheme"'\n' >> npt.mdp
        echo -e "; Velocity generation" >> npt.mdp
        echo -e "gen_vel		= no		; Velocity generation is off " >> npt.mdp

        # Create md (data collection run) input file
        echo -e "title		= Electric Fields solvatochromism MD simulations"'\n' > md.mdp
        echo -e "; Run parameters" >> md.mdp
        echo -e "integrator	= sd		; leap-frog integrator" >> md.mdp
        echo -e "nsteps		= $mdsteps		; 2 * 100000 = 100 ps" >> md.mdp
        echo -e "dt		    = 0.002		; Step-size 2 fs"'\n' >> md.mdp
        echo -e "; Output control" >> md.mdp
        echo -e "nstxout		= 100		; save coordinates every 0.2 ps" >> md.mdp
        echo -e "nstfout		= 100		; save forces every 0.2 ps" >> md.mdp
        echo -e "nstvout		= 0		; do not save velocities" >> md.mdp
        echo -e "nstenergy	= 5000		; save energies every 10.0 ps" >> md.mdp
        echo -e "nstlog		= 5000 		; update log file every 10.0 ps" >> md.mdp
        echo -e "nstxout-compressed  = 5000      ; save compressed coordinates every 10.0 ps (replaces nstxtcout)" >> md.mdp
        echo -e "compressed-x-grps   = System    ; replaces xtc-grps"'\n' >> md.mdp
        echo -e "; Bond parameters" >> md.mdp
        echo -e "continuation	        = yes		; Restarting after md run" >> md.mdp
        echo -e "constraint_algorithm    = lincs	    ; holonomic constraints" >> md.mdp
        echo -e "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained" >> md.mdp
        echo -e "lincs_iter	            = 1		    ; accuracy of LINCS" >> md.mdp
        echo -e "lincs_order	            = 4		    ; also related to accuracy"'\n' >> md.mdp
        echo -e "; Neighborsearching" >> md.mdp
        echo -e "cutoff-scheme   = Verlet" >> md.mdp
        echo -e "ns_type		    = grid		; search neighboring grid cells" >> md.mdp
        echo -e "nstlist		    = 10		; 20 fs, largely irrelevant with Verlet" >> md.mdp
        echo -e "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)" >> md.mdp
        echo -e "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)" >> md.mdp
        echo -e "; Electrostatics" >> md.mdp
        echo -e "coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics" >> md.mdp
        echo -e "pme_order	    = 4		; cubic interpolation" >> md.mdp
        echo -e "fourierspacing	= 0.16	; grid spacing for FFT"'\n'  >> md.mdp
        echo -e "; Temperature coupling is on" >> md.mdp
        echo -e "tcoupl		= V-rescale	            ; modified Berendsen thermostat" >> md.mdp
        echo -e "tc-grps		= LIG    SOL	; two coupling groups - more accurate" >> md.mdp
        echo -e "tau_t		= 2	  2           ; time constant, in ps" >> md.mdp
        echo -e "ref_t		= 300 	  300           ; reference temperature, one for each group, in K"'\n' >> md.mdp
        echo -e "; Pressure coupling is on" >> md.mdp
        echo -e "pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in md" >> md.mdp
        echo -e "pcoupltype	        = isotropic	            ; uniform scaling of box vectors" >> md.mdp
        echo -e "tau_p		        = 2.0		            ; time constant, in ps" >> md.mdp
        echo -e "ref_p		        = 1.0		            ; reference pressure, in bar" >> md.mdp
        echo -e "compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1"'\n' >> md.mdp
        echo -e "; Periodic boundary conditions" >> md.mdp
        echo -e "pbc		        = xyz 		; 3-D Periodic Boundary Conditions (yes/no)"'\n' >> md.mdp
        echo -e "; Dispersion correction" >> md.mdp
        echo -e "DispCorr	= EnerPres	; account for cut-off vdW scheme"'\n' >> md.mdp
        echo -e "; Velocity generation" >> md.mdp
        echo -e "gen_vel		= no		; Velocity generation is off " >> md.mdp

    cd ..
done