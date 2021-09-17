#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=makeTOP       #Set the job name 
#SBATCH --time=10:00 #Set the wall clock limit 
#SBATCH --ntasks=1                   #Request tasks
#SBATCH --mem=10GB                  #Request Memory in MB per node
#SBATCH --output=job.%j      #Send stdout/err 

export PATH=/sw/group/lms/sw/schrodinger/2020-1/utilities/:$PATH

MOL_GEOM=YOUR_FILE_NAME

#Convert input geometry to sdf file
obabel $MOL_GEOM.xyz -O $MOL_GEOM.sdf

#Generate schrodinger force field file
ffld_server -isdf $MOL_GEOM.sdf -print_parameters -out_file $MOL_GEOM.log -version 14

#Convert schrodinger FF file to GROMACS format
python genTOP.py -f $MOL_GEOM.xyz -l $MOL_GEOM.log -o mol.gro -t $MOL_GEOM.top

#Clean up
rm $MOL_GEOM.log $MOL_GEOM.sdf

