#!/bin/bash
#SBATCH --job-name=run_slidetag_processing   # Name of your job
#SBATCH --account=def-shreejoy                # Compute Canada account
#SBATCH --time=16:00:00                        # Job time (HH:MM:SS)
#SBATCH --nodes=1                              # Number of nodes
#SBATCH --ntasks=1                             # Number of tasks
#SBATCH --cpus-per-task=4                      # Request 4 CPU cores per task
#SBATCH --mem=32G                              # Request 64 GB of memory
#SBATCH --output=/home/hudsonhu/scratch/NF_MouseDataFirstRun.log
#SBATCH --error=/home/hudsonhu/scratch/NF_MouseDataFirstRun.err   # Error file name
#SBATCH --mail-user=hugo.hudson@mail.utoronto.ca                    
#SBATCH --mail-type=BEGIN,END,FAIL             

# Load the Apptainer module
module load StdEnv/2023                     # Load default environment
module load apptainer/1.3.4                    # Adjust version if needed
module load nextflow

# Run Cell Bender
nextflow run /home/hudsonhu/scratch/SlideTagNextflow/main.nf --samplesheet=/home/hudsonhu/scratch/SlideTagNextflow/samplesheetFirst.csv -profile standard