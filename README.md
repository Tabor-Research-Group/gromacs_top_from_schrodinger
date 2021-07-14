### Code Overview

This code is used to convert Schrodinger force field output to GROMACS format and also generate other input files used to initiate GROMACS simulation.

### Required Resources

1. Schrodinger 

2. You will need to have the Atomic Simulation Environment (ASE) (https://wiki.fysik.dtu.dk/ase/install.html) installed in your python environment. You can install it with:

$ pip install --update --user ase

### Execution of Code

1. Place the genTOP.py and run.sh in your directory that contains the geometry of your monomer that you are interested. The format of the geometry can be either .xyz/.gro/.pdb that supported by the ASE.

2. Change the filename in the run.sh file to match with your monomer, do not include the extension (.xyz/.pdb/.gro etc.)

MOL_GEOM=YOUR_FILE_NAME

3. The default run.sh file is built for a SLURM cluster at Texas A&M. Modify and submit the job for your cluster. 

$ sbatch run.sh

4. The code will generate several output files:

  *.top # Main topology file in gromacs format with FF included.

  ffnonbonded.itp # The GROMACS force field file contains the non-bonding parameters

  topol.top # Default GROMACS topology

  mol.gro # GROMACS input geometry file with updated atomic type match with the topology. If your input geometry is in .gro format, the output will keep the same box, otherwise, the output gro has no box, be sure to use the 'gmx editconf' command to add the box before you use it for simulation. 

