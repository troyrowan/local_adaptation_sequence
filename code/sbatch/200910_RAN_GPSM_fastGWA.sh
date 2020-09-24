#!/bin/sh
#SBATCH --partition=hpc5,BioCompute
#SBATCH --account=animalsci
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH --qos=biolong
#SBATCH --time=7-00:00:00
#SBATCH --job-name=fastGWA
#SBATCH --output=./output/200910_RAN/log/slurm_out/200910_RAN.fastGWA.log

psrecord "code/gcta_1.93.2beta/gcta64 --bfile output/200910_RAN/plink_convert/200910_RAN.850K --fastGWA-mlm --model-only --pheno output/200910_RAN/phenotypes/200910_RAN.age.txt --grm-sparse output/200910_RAN/grm/200910_RAN.850K.sparse --autosome-num 29 --threads 10 --out output/200910_RAN/gpsm/gwas/200910_RAN.gpsm.fastGWA.850K" --log output/200910_RAN/log/psrecord/gwas_gpsm/200910_RAN.fastGWA.log --interval 30
