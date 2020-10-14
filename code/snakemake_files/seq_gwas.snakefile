import os
#Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = config['run_name'], rule = ["convert_seq_plink", "seq_gwas", "seq_cojo"]):
	os.makedirs(x, exist_ok = True)
for x in expand("output/{run_name}/log/psrecord/{rule}", run_name = config['run_name'], rule = ["convert_seq_plink", "seq_gwas", "seq_cojo"]):
	os.makedirs(x, exist_ok = True)

#include: "GRM.snakefile"

rule envgwas:
	input:
		expand("output/{run_name}/seq_gwas/{run_name}.{phenotype}.chr{chr}.mlma.gz",
		phenotype = config["phenotypes"],
		run_name = config["run_name"],
		chr = list(range(1,30)))

phenotype_dict = config["phenotype_dict"]

def columnchooser(WC):
	column = phenotype_dict[WC.phenotype]
	return column

#prior_dict = config["reml_priors"]

# def priorchooser(WC):
# 	column = prior_dict[WC.phenotype]
# 	return column
#This GRM creation will just be for autosomes as there are some issues with how you'd want to go about handing the X Chromosome https://cnsgenomics.com/software/gcta/#MakingaGRM

rule convert_seq_plink:
	input:
		vcf = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.filtered.vcf.gz",
		#tbi = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.filtered.vcf.gz.tbi"
	params:
		psrecord = "output/{run_name}/log/psrecord/convert_seq_plink/convert_seq_plink.{run_name}.chr{chr}.log",
		oprefix = "output/{run_name}/imputed_genotypes/seq/{run_name}.chr{chr}.filtered",
		iprefix = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.filtered"
	threads: config["plink_threads"]
	output:
		out = expand("output/{{run_name}}/imputed_genotypes/seq/{{run_name}}.chr{{chr}}.filtered.{suffix}", suffix = ["bed", "bim", "fam"])
	shell:
		# "code/gcta_1.93.2beta/gcta64 --reml --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --autosome-num 29 --reml-pred-rand --thread-num {params.threads} --out {params.oprefix}"
		"""
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		module load plink
		psrecord "plink --vcf {input.vcf} --cow --real-ref-alleles --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
		"""

# rule snp_effects:
# 	input:
# 		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
# 		suffix = ["bed", "bim", "fam"]),
# 		blups = "output/{run_name}/greml/{run_name}.{phenotype}.850K.indi.blp",
# 		phenotypes = config["phenotype_file"]
# 	params:
# 		psrecord = "output/{run_name}/log/psrecord/snp_effects/snp_effects.{phenotype}.log",
# 		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
# 		grmprefix = "output/{run_name}/grm/{run_name}.850K",
# 		oprefix = "output/{run_name}/greml/{run_name}.{phenotype}.850K",
# 	threads: config["gwas_threads"]
# 	output:
# 		blups = "output/{run_name}/greml/{run_name}.{phenotype}.850K.snp.blp"
# 	shell:
# 		# "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --blup-snp {input.blups} --thread-num {params.threads} --out {params.oprefix}"
# 		"""
# 		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
# 		module load openmpi/openmpi-3.1.3-intel-16.0.2
# 		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --blup-snp {input.blups} --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
# 		"""
rule seq_gwas:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/seq/{{run_name}}.chr{{chr}}.filtered.{suffix}", suffix = ["bed", "bim", "fam"]),
		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
		phenotypes = config["phenotype_file"]
	params:
		psrecord = "output/{run_name}/log/psrecord/seq_gwas/seq_gwas.{phenotype}.chr{chr}.log",
		iprefix = "output/{run_name}/imputed_genotypes/seq/{run_name}.chr{chr}.filtered",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/seq_gwas/{run_name}.{phenotype}.chr{chr}",
		phenotype = columnchooser,
		out = "output/{run_name}/seq_gwas/{run_name}.{phenotype}.chr{chr}.mlma",
		chrom = "{chr}",
		maf = config["maf"]
	threads: config["gwas_threads"]
	output:
		log = "output/{run_name}/seq_gwas/{run_name}.{phenotype}.chr{chr}.log",
		out = "output/{run_name}/seq_gwas/{run_name}.{phenotype}.chr{chr}.mlma.gz",
	shell:
		# "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {params.threads} --out {params.oprefix}"
		"""
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
		pigz {params.out}
		"""

rule seq_cojo:
	input:
		log = "output/{run_name}/seq_gwas/{run_name}.{phenotype}.chr{chr}.log",
		out = "output/{run_name}/seq_gwas/{run_name}.{phenotype}.chr{chr}.mlma.gz",
		plink = expand("output/{{run_name}}/imputed_genotypes/seq/{{run_name}}.chr{{chr}}.filtered.{suffix}",
		suffix = ["bed", "bim", "fam"])
	params:
		n = config["n"],
		cojofile = "output/{run_name}/seq_cojo/{run_name}.{phenotype}.chr{chr}.ma",
		bfile = "output/{run_name}/imputed_genotypes/seq/{run_name}.chr{chr}.filtered",
		cojop = config["cojop"],
		oprefix = "output/{run_name}/seq_cojo/{run_name}.{phenotype}.chr{chr}",
		psrecord = "output/{run_name}/log/psrecord/seq_cojo/seq_cojo.{phenotype}.chr{chr}.log",
	threads: config["grm_threads"]
	output:
		cojo = "output/{run_name}/seq_cojo/{run_name}.{phenotype}.chr{chr}.jma"
	shell:
		"""
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		zcat {input.assoc} | awk '{{print $2, $4, $5, $6, $7, $8, $9, {params.n}}}' > {params.cojofile}
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.bfile} --autosome-num 29 --maf 0.01 --cojo-file {params.cojofile} --cojo-slct --cojo-p {params.cojop} --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
		"""
