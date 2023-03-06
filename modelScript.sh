#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=PEIXOTT@ccf.org
#SBATCH -n %%%CORES%%%
#SBATCH --mem=%%%MEM%%%
#SBATCH --job-name=%%%NAME%%%
#SBATCH -o %%%NAME%%%.out        # Standard output
#SBATCH -e %%%NAME%%%.err        # Standard error

module load bcftools/1.9
module load tabix
module load plink/1.90
module load Eagle/2.4.1
module loaf python/3.8.6
module load RFMix/1.5.4
module load python/2.7.12

