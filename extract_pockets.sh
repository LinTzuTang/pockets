#!/bin/bash
#SBATCH --job-name=extract_pockets    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --partition=bigmem
#SBATCH --qos=yanjun.li-b
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=extract_pockets.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate hariboss

python extract_pockets.py

date