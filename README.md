### genTOP

### Currently works on Terra

### Keep in mind that the FF included here are extracted from OPLS_2005. To obtain accurate simulation results, it is alway recommanded to parameterized the distinct dihedrals.  

This code is used to convert Schrodinger force field output to GROMACS format and also generate other input files used to initiate GROMACS simulation.

0. You will need to have the Atomic Simulation Environment (ASE) (https://wiki.fysik.dtu.dk/ase/install.html) installed in your python environment. You can install it with:

$ pip install --update --user ase

1. Then, put the genTOP.py and run.sh in your directory along with the geometry of your monomer that you are interested. The format of the geometry can be either .xyz/.gro/.pdb that supported by the ASE.

2. Change the filename in the run.sh file to match with your monomer, do not include the extension (.xyz/.pdb/.gro etc.)

MOL_GEOM=YOUR_FILE_NAME

3. Submit the job script to Terra:

$ sbatch run.sh

4. The code will generate several output files:

  *.top # Main topology file in gromacs format with FF included.

  ffnonbonded.itp # The GROMACS force field file contains the non-bonding parameters

  topol.top # Default GROMACS topology

  mol.gro # GROMACS input geometry file with updated atomic type match with the topology. The file has no initial box, be sure to use the 'gmx editconf' command to add the box before you use it for simulation. 

