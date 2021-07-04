#!/bin/bash

# By Yasmin 2021-07-02 
# Tested locally 2021-07-02
# Tested on Sherlock 2021-07-02

# This script prepares job files for GROMACS simulations. 
# This script submits files on the Stanford Sherlock cluster.
# The submit file also runs a script for calculating electric fields. 

# Requires access on Stanford's Sherlock cluster. Can be modified for other clusters.

############################# EDIT THESE PARAMETERS #########################################

walltime=3:00:00                            # Default for solvatochromism: 3:00:00

ligandDirectory=LIGANDS                     # Path to the ligand directory
scriptDirectory=SCRIPTS                     # Path to python scripts
pythonVersion=python3                       # Installed python version

partitions=hns,iric,owners                  # Default for Boxers: hns

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

nodes=1                                     # Default: 1
cores=16                                    # Default: 20

######################## NO MORE EDITING PAST THIS LINE! ###################################

cd $ligandDirectory

# Repeat for all ligands
for ligname in *
 do
    cd $ligname

# Create generic job-files

        # Assign Sherlock flags for running
        echo -e "#!/bin/bash" > COMPLEX.job
        echo -e "#SBATCH --mail-type=ALL" >> COMPLEX.job
        echo -e "#SBATCH --job-name=COMPLEX" >> COMPLEX.job
        echo -e "#SBATCH -p $partitions" >> COMPLEX.job
        echo -e "#SBATCH -N $nodes" >> COMPLEX.job
        echo -e "#SBATCH -n $cores" >> COMPLEX.job
        echo -e "#SBATCH --time=$walltime"'\n' >> COMPLEX.job

        # Load gromacs and python modules
        echo -e "module load chemistry gromacs" >> COMPLEX.job
        echo -e "module load py-numpy/1.14.3_py27"'\n' >> COMPLEX.job

        echo -e "mkdir GMXRUN"'\n' >> COMPLEX.job

        echo -e "gmx grompp -f min.mdp -c COMPLEX.gro -p COMPLEX.top -o COMPLEX_min.tpr" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_min"'\n' >> COMPLEX.job
        echo -e "gmx grompp -f nvt.mdp -c COMPLEX_min.gro -p COMPLEX.top -o COMPLEX_nvt.tpr" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_nvt"'\n' >> COMPLEX.job
        echo -e "gmx grompp -f npt.mdp -c COMPLEX_nvt.gro -p COMPLEX.top -o COMPLEX_npt.tpr" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_npt"'\n' >> COMPLEX.job
        echo -e "gmx grompp -f md.mdp -c COMPLEX_npt.gro -p COMPLEX.top -o COMPLEX_md.tpr" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_md"'\n' >> COMPLEX.job

        echo -e "gmx traj -s COMPLEX_md.tpr -f COMPLEX_md.trr -n probe.ndx -ox co_coords.xvg" >> COMPLEX.job
        echo -e "gmx traj -s COMPLEX_md.tpr -f COMPLEX_md.trr -n probe.ndx -of co_forces.xvg" >> COMPLEX.job

        echo -e "gmx grompp -f md.mdp -c COMPLEX_npt.gro -p COMPLEX_0q.top -o COMPLEX_md_0q.tpr" >> COMPLEX.job
        echo -e "gmx mdrun -rerun COMPLEX_md.trr -v -deffnm COMPLEX_md_0q"'\n' >> COMPLEX.job

        echo -e "gmx traj -s COMPLEX_md_0q.tpr -f COMPLEX_md_0q.trr -n probe.ndx -ox co_coords_0q.xvg" >> COMPLEX.job
        echo -e "gmx traj -s COMPLEX_md_0q.tpr -f COMPLEX_md_0q.trr -n probe.ndx -of co_forces_0q.xvg"'\n' >> COMPLEX.job

        echo -e "python BLCO.py" >> COMPLEX.job
        echo -e "wait"'\n' >> COMPLEX.job
        echo -e "mv *COMPLEX* GMXRUN/" >> COMPLEX.job

# Customize job-files for all solvents

    for i in *0q.top   # Get the name of each solvent
     do
        complex=$( echo "$i" | sed -e 's/\_0q.top//g')
        sed "s/COMPLEX/$complex/g" COMPLEX.job > $complex.job
        sbatch $complex.job 
     wait
    done
    cd ..
done