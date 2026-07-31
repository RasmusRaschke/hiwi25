#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --partition=obelix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=%j_%x.log

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
manim -pqh main.py Anim
