#!/bin/bash
#SBATCH --job-name=profumo_HC+OCD_job              # Name of the job
#SBATCH --output=/gpfs/gibbs/pi/n3/Studies/Pushkarskaya.PVDM/analysis/full_analysis_unharm/profumo_HC+OCD_output.log # File to write standard output
#SBATCH --error=/gpfs/gibbs/pi/n3/Studies/Pushkarskaya.PVDM/analysis/full_analysis_unharm/profumo_HC+OCD_error.log   # File to write standard error
#SBATCH --partition=week                        # Partition/queue name
#SBATCH --time=3-00:00:00                       # Time limit (6 days)
#SBATCH --mem=819200                          # Memory request (in MB)
#SBATCH --cpus-per-task=24                      # Number of CPU cores
#SBATCH --ntasks=1                              # Number of tasks

cd /gpfs/gibbs/pi/n3/software/profumo
 
singularity exec --cleanenv /gpfs/gibbs/pi/n3/software/profumo/profumo.sif \
  /bin/bash -c "source /etc/profile; PROFUMO /gpfs/gibbs/pi/n3/Studies/Pushkarskaya.PVDM/analysis/full_analysis_unharm/HC+OCD.json 50 /gpfs/gibbs/pi/n3/Studies/Pushkarskaya.PVDM/analysis/full_analysis_unharm_HC+OCD.pfm"

