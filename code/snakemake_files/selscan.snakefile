#Pulls out subsetted individuals and their respectve haplotypes

rule selscan:
	input:
		expand("output/{run_name}/downsample/{run_name}.{subset}.chr{chr}.vcf.gz",
		run_name = ["200910_RAN", "200907_SIM"],
		subset = ["oldest_100_perbreeder", "oldest_100", "oldest_500_perbreeder", "youngest_500_perbreeder"],
		chr = list(range(1,30)))
rule tabix:
	input:
		vcf = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz"
	output:
		tabix = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		tabix {input.vcf}
		"""

		
rule extract_vcf:
	input:
		vcf = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz", #Only Autosomes
		tabix = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz.tbi",
		subsets = "output/{run_name}/subsets/{run_name}.{subset}.txt"
	output:
		"output/{run_name}/downsample/{run_name}.{subset}.chr{chr}.vcf.gz"
	shell:
		"""
		module load bcftools
		bcftools view -S {input.subsets} -O z -o {output}
		"""
