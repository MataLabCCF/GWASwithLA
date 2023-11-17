# GWASwithLA

In construction


## Parameters 

![Parameters](./Figures/Parameters.png?style=centerme)

## Preparing the data

![Common](./Figures/mergePhaseSplit.png?style=centerme)

## Local Ancestry with RFMix

### Sequential

![Sequential](./Figures/Sequential.png?style=centerme)

### Bot (HPC)

![Sequential](./Figures/Bot.png?style=centerme)

## Command Line example

```
python3.8 main3.py \
-i /home/peixott/beegfs/Analysis/DataClean/CleanData/LARGE_hg38/LARGE_chr\* -I AAAAAA \
-o LARGE_PD -R /cm/shared/apps/RFMix/1.5.4/RunRFMix.py -Y python3.8 -y python2.7 \
-O /home/peixott/beegfs/Analysis/LocalAncestry/TesteParallel_12_22/ -t 48 \
-V /home/peixott/beegfs/Analysis/LocalAncestry/VCF2RFMix1.py -r RFMix1ToRFMix2.py \
-l /home/peixott/beegfs/References/1KGP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr\*.filtered.shapeit2-duohmm-phased.vcf.gz \
-p /home/peixott/beegfs/References/1KGP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr\*.filtered.shapeit2-duohmm-phased.vcf.gz \
-b 12 -e 22 -s 100 -g /home/peixott/beegfs/References/geneticMap/genetic_map_hg38_withX.txt.gz -c correspondenceT.txt \
-f /home/peixott/beegfs/References/Fasta/hg38/Homo_sapiens_assembly38.fasta -n 10 \
-m /home/peixott/beegfs/Analysis/LocalAncestry/modelScript.sh \
-j 20 -q squeue -Q sbatch -M 90G -C 24
```
## Contact

If you chat about this on going project, please sent an email to me

Thiago Peixoto Leal, PhD

Email: PEIXOTT@ccf.org

