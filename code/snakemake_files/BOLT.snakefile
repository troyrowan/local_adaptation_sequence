import os
#Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("output/{run_name}/log/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

#include: "GRM.snakefile"

rule envgwas:
	input:
		expand("output/{run_name}/greml/{run_name}.{phenotype}.850K.hsq",
		run_name = config["run_name"],
		phenotype = config["phenotypes"]),
		expand("output/{run_name}/gwas/single_chrom/{run_name}.{phenotype}.chr{chr}.850K.mlma",
		run_name = config["run_name"],
		phenotype = config["phenotypes"],
		chr = config["chroms"])

phenotype_dict = config["phenotype_dict"]

def columnchooser(WC):
	column = phenotype_dict[WC.phenotype]
	return column

prior_dict = config["reml_priors"]

def priorchooser(WC):
	column = prior_dict[WC.phenotype]
	return column


rule BOLT_gwas:
	input:
	params:
	threads:
	output:
	shell:
		"""
		code/BOLT-LMM_v2.3.4/bolt --Nautosomes 29 --bfile {params.inprefix} --phenoFile {params.phenotypes} --phenoCol {params.phenotypes} --numLeaveOutChunks 534 --lmmInfOnly --numThreads {threads} --modelSnps {params.coreSNPs} --statsFile {params.oprefix}"""
