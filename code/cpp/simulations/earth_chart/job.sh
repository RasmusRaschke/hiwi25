#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --partition=obelix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --mem-per-cpu=2gb
#SBATCH --output=%j_%x.log
#SBATCH --mail-user=rasmus.raschke@uni-hamburg.de
#SBATCH --mail-type=END,FAIL

############################
### USER COMMANDS TO EXECUTE
############################

### Purge all modules first, just to be safe
# module purge
# module load gcc/13.1
# module load python/3.12.9
module purge
module load anaconda3
module load texlive
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/.conda/envs/visual
python3 earth_grid.py
