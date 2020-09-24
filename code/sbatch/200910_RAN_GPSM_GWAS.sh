#!/bin/sh
#SBATCH --partition=hpc5,BioCompute
#SBATCH --account=animalsci
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --qos=biolong
#SBATCH --time=7-00:00:00
#SBATCH --job-name=RedAngus_GPSM_GWAS_realpriors
#SBATCH --output=./output/200910_RAN/log/slurm_out/sexgwas.log

psrecord "code/gcta_1.93.2beta/gcta64 --bfile output/200910_RAN/plink_convert/200910_RAN.850K --mlma --reml-priors 0.5401733 0.5401733 --pheno output/200910_RAN/phenotypes/200910_RAN.age.txt --grm output/200910_RAN/grm/200910_RAN.850K --autosome-num 29 --thread-num 10 --out output/200910_RAN/gpsm/gwas/200910_RAN.gpsm.850K" --log output/200907_SIM/log/psrecord/gwas_gpsm/200910_RAN.gpsm_gwas.log --interval 30
	
