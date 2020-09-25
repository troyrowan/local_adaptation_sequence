import os
#Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("output/{run_name}/log/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

#include: "GRM.snakefile"

rule envgwas:
	input:
		# expand("output/{run_name}/greml/{run_name}.{phenotype}.850K.hsq",
		# run_name = config["run_name"],
		# phenotype = config["phenotypes"]),
		expand("output/{run_name}/gwas/{run_name}.{phenotype}.850K.mlma",
		run_name = config["run_name"],
		phenotype = config["phenotypes"])

phenotype_dict = config["phenotype_dict"]

def columnchooser(WC):
	column = phenotype_dict[WC.phenotype]
	return column
#This GRM creation will just be for autosomes as there are some issues with how you'd want to go about handing the X Chromosome https://cnsgenomics.com/software/gcta/#MakingaGRM

rule greml:
	input:
		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
		phenotypes = config["phenotype_file"]
	params:
		psrecord = "output/{run_name}/log/psrecord/greml/greml.{phenotype}.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/greml/{run_name}.{phenotype}.850K",
		phenotype = columnchooser,
	threads: config["greml_threads"]
	output:
		out = "output/{run_name}/greml/{run_name}.{phenotype}.850K.hsq",
		blups = "output/{run_name}/greml/{run_name}.{phenotype}.850K.indi.blp",
		log = "output/{run_name}/greml/{run_name}.{phenotype}.850K.log",
	shell:
		# "code/gcta_1.93.2beta/gcta64 --reml --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --autosome-num 29 --reml-pred-rand --thread-num {params.threads} --out {params.oprefix}"
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --reml --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --autosome-num 29 --reml-pred-rand --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""

rule snp_effects:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"]),
		blups = "output/{run_name}/greml/{run_name}.{phenotype}.850K.indi.blp",
		phenotypes = config["phenotype_file"]
	params:
		psrecord = "output/{run_name}/log/psrecord/snp_effects/snp_effects.{phenotype}.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/greml/{run_name}.{phenotype}.850K",
	threads: config["gwas_threads"]
	output:
		blups = "output/{run_name}/greml/{run_name}.{phenotype}.850K.snp.blp"
	shell:
		# "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --blup-snp {input.blups} --thread-num {params.threads} --out {params.oprefix}"
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --blup-snp {input.blups} --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""
rule gwas:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"]),
		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
		phenotypes = config["phenotype_file"]
	params:
		psrecord = "output/{run_name}/log/psrecord/gwas/gwas.{phenotype}.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/gwas/{run_name}.{phenotype}.850K",
		phenotype = columnchooser,
		maf = config["maf"]
	threads: config["gwas_threads"]
	output:
		out = "output/{run_name}/gwas/{run_name}.{phenotype}.850K.mlma",
		log = "output/{run_name}/gwas/{run_name}.{phenotype}.850K.log"
	shell:
		# "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {params.threads} --out {params.oprefix}"
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""
