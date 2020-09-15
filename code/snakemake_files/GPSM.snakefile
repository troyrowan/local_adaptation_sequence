import os
# Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("output/{run_name}/log/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

#include: "GRM.snakefile"
rule gpsm:
	input:
		# expand("output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K.hsq",
		# run_name = config["run_name"],
		# phenotype = config["phenotypes"]),
		expand("output/{run_name}/gpsm/gwas/{run_name}.{phenotype}.850K.mlma",
		run_name = config["run_name"],
		phenotype = config["phenotypes"],
		chr = list(range(1,30))),
		# expand("output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K.snp.blp",
		# run_name = config["run_name"],
		# phenotype = config["phenotypes"],
		# chr = list(range(1,30))),
		# expand("output/{run_name}/gpsm/gwas/{run_name}.{phenotype}.chr{chr}.850K.loco.mlma",
		# run_name = config["run_name"],
		# phenotype = config["phenotypes"],
		# chr = list(range(1,30))),

rule convert_plink:
	input:
		vcf = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz", #Only Autosomes
		tabix = "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz.tbi"
	params:
		psrecord = "output/{run_name}/log/psrecord/convert_plink/convert_plink.chr{chr}.log",
		oprefix = "output/{run_name}/plink_convert/{run_name}.chr{chr}.850K",
		threads=config["plink_threads"],
		mem=config["plink_mem"]
	output:
		plink = expand("output/{{run_name}}/plink_convert/{{run_name}}.chr{{chr}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"]),
	shell: #
		"""
		module load plink
		psrecord "plink --vcf {input.vcf} --threads {params.threads} --memory {params.mem} --cow --real-ref-alleles --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""

rule greml_gpsm:
	input:
		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
		phenotypes = "output/{run_name}/phenotypes/{run_name}.{phenotype}.txt"
	params:
		psrecord = "output/{run_name}/log/psrecord/greml_gpsm/greml_gpsm.{phenotype}.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K",
		threads=config["greml_threads"]
	output:
		out = "output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K.hsq",
		blups = "output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K.indi.blp"
	shell:
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --reml --pheno {input.phenotypes} --grm {params.grmprefix} --autosome-num 29 --reml-pred-rand --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""

rule snp_effects_gpsm:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"]),
		blups = "output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K.indi.blp",
		phenotypes = "output/{run_name}/phenotypes/{run_name}.{phenotype}.txt"
	params:
		psrecord = "output/{run_name}/log/psrecord/snp_effects_gpsm/snp_effects_gpsm.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K",
		threads=config["gwas_threads"]
	output:
		blups = "output/{run_name}/gpsm/greml/{run_name}.{phenotype}.850K.snp.blp"
	shell:
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --blup-snp {input.blups} --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""
rule gwas_gpsm:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"]),
		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
		phenotypes = "output/{run_name}/phenotypes/{run_name}.{phenotype}.txt"
	params:
		psrecord = "output/{run_name}/log/psrecord/gwas_gpsm/gwas_gpsm.log",
		iprefix = "output/{run_name}/plink_convert/{run_name}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/gpsm/gwas/{run_name}.{phenotype}.850K",
		threads=config["gwas_threads"]
	output:
		out = "output/{run_name}/gpsm/gwas/{run_name}.{phenotype}.850K.mlma"
	shell:
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --grm {params.grmprefix} --autosome-num 29 --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""

rule loco_gwas_gpsm:
	input:
		plink = expand("output/{{run_name}}/plink_convert/{{run_name}}.chr{chr}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"],
		chr = list(range(1,30))),
		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
		phenotypes = "output/{run_name}/phenotypes/{run_name}.{phenotype}.txt",
		chrgrm = expand("output/{{run_name}}/grm/{{run_name}}.chr{{chr}}.850K.grm.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"])
	params:
		psrecord = "output/{run_name}/log/psrecord/gwas_gpsm/gwas_gpsm.log",
		iprefix = "output/{run_name}/plink_convert/{run_name}.chr{chr}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		locoprefix = "output/{run_name}/grm/{run_name}.chr{chr}.850K",
		oprefix = "output/{run_name}/gpsm/gwas/{run_name}.{phenotype}.chr{chr}.850K",
		threads=config["gwas_threads"]
	output:
		out = "output/{run_name}/gpsm/gwas/{run_name}.{phenotype}.chr{chr}.850K.loco.mlma"
	shell:
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --grm {params.grmprefix} --mlma-subtract-grm {params.locoprefix} --autosome-num 29 --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""
