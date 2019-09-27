#!/bin/bash
#SBATCH --job-name=Minus_recombinant_regions_constraints_3
#SBATCH --time=740:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=10240

module load beagle-lib/3.0.2-fosscuda-2018b

java -jar beast_1104_prerelease_051118.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite Minus_recombinant_regions_constraints_3.xml
