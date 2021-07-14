# genTOP

# Currently works on Terra

This code is used to convert Schrodinger force field output to GROMACS format and also generate other input files used to initiate GROMACS simulation.

1. You will need to have the Atomic Simulation Environment (ASE) (https://wiki.fysik.dtu.dk/ase/install.html) installed in your python environment. You can install it with:

$ pip install --update --user ase

2. Then, put the genTOP.py and run.sh in your directory along with the geometry of your monomer that you are interested. The format of the geometry can be either .xyz/.gro/.pdb that supported by the ASE.

3. Change the filename in the run.sh file to match with your monomer, do not include the extension (.xyz/.pdb/.gro etc.)

MOL_GEOM=YOUR_FILE_NAME

4. Submit the job script to Terra:

$ sbatch run.sh

5. The code will generate several output files:

1) *.top # this is the main topology file in gromacs format with FF included.
2) ffnonbonded.itp # This is the GROMACS force field file contains the non-bonding parameters
3) topol.top # Default GROMACS topology
4) mol.gro # GROMACS input geometry file with updated atomic type match with the topology. The file has no initial box, be sure to use the 'gmx editconf' command to add the box before you use it for simulation. 

