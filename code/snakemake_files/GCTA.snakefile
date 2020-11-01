import os
#Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = config['run_name'], rule = ["convert_plink", "concat_plink", "build_grm_chunks", "concat_grm_chunks", "greml", "gwas", "snp_effects", "gwas_single_chrom", "cojo"]):
	os.makedirs(x, exist_ok = True)
for x in expand("output/{run_name}/log/psrecord/{rule}", run_name = config['run_name'], rule = ["convert_plink", "concat_plink", "build_grm_chunks", "concat_grm_chunks", "greml", "gwas", "snp_effects", "gwas_single_chrom", "cojo"]):
	os.makedirs(x, exist_ok = True)

include: "GRM.snakefile"

rule envgwas:
	input:
		expand("output/{run_name}/greml/{run_name}.{phenotype}.850K.hsq",
		run_name = config["run_name"],
		phenotype = config["phenotypes"]),
		# expand("output/{run_name}/cojo/{run_name}.{phenotype}.850K.jma",
		# run_name = config["run_name"],
		# phenotype = config["phenotypes"])
		expand("output/{run_name}/gwas/{run_name}.{phenotype}.850K.mlma.gz",
		run_name = config["run_name"],
		phenotype = config["phenotypes"])

phenotype_dict = config["phenotype_dict"]
n_dict = config["n"]
def columnchooser(WC):
	column = phenotype_dict[WC.phenotype]
	return column

def nchooser(WC):
	n = n_dict[WC.phenotype]
	return n

#prior_dict = config["reml_priors"]

# def priorchooser(WC):
# 	column = prior_dict[WC.phenotype]
# 	return column
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
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		psrecord "code/gcta_1.93.2beta/gcta64 --reml --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --autosome-num 29 --reml-pred-rand --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
		"""

# rule greml_prior:
# 	input:
# 		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
# 		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
# 		phenotypes = config["phenotype_file"]
# 	params:
# 		psrecord = "output/{run_name}/log/psrecord/greml/greml.{phenotype}.log",
# 		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
# 		grmprefix = "output/{run_name}/grm/{run_name}.850K",
# 		oprefix = "output/{run_name}/greml/{run_name}.{phenotype}.850K",
# 		phenotype = columnchooser,
# 		prior = priorchooser,
# 		threads=config["greml_threads"]
# 	output:
# 		out = "output/{run_name}/greml/{run_name}.{phenotype}.850K.hsq",
# 		blups = "output/{run_name}/greml/{run_name}.{phenotype}.850K.indi.blp",
# 		log = "output/{run_name}/greml/{run_name}.{phenotype}.850K.log",
# 	shell:
# 		# "code/gcta_1.93.2beta/gcta64 --reml --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --autosome-num 29 --reml-pred-rand --thread-num {params.threads} --out {params.oprefix}"
# 		"""
# 		psrecord "code/gcta_1.93.2beta/gcta64 --reml --reml-priors {params.prior} --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --autosome-num 29 --reml-pred-rand --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
# 		"""

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
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --blup-snp {input.blups} --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
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
		out = "output/{run_name}/gwas/{run_name}.{phenotype}.850K.mlma",
		#prior = priorchooser, #This uses priors to try and prevent as many iterations of REML. Doesn't appear to make much of a difference in speed
		maf = config["maf"]
	threads: config["gwas_threads"]
	output:
		log = "output/{run_name}/gwas/{run_name}.{phenotype}.850K.log",
		out = "output/{run_name}/gwas/{run_name}.{phenotype}.850K.mlma.gz",
	shell:
		# "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {params.threads} --out {params.oprefix}"
		"""
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
		pigz {params.out}
		"""

rule gwas_single_chrom:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"]),
		grm = expand("output/{{run_name}}/grm/{{run_name}}.850K.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"]),
		phenotypes = config["phenotype_file"]
	params:
		psrecord = "output/{run_name}/log/psrecord/gwas_single_chrom/gwas_single_chrom.{phenotype}.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		grmprefix = "output/{run_name}/grm/{run_name}.850K",
		oprefix = "output/{run_name}/gwas/single_chrom/{run_name}.{phenotype}.chr{chr}.850K",
		chrom = "{chr}",
		phenotype = columnchooser,
		#prior = priorchooser,
		maf = config["maf"]
	threads: config["gwas_threads"],
	output:
		out = "output/{run_name}/gwas/single_chrom/{run_name}.{phenotype}.chr{chr}.850K.mlma",
		log = "output/{run_name}/gwas/single_chrom/{run_name}.{phenotype}.chr{chr}.850K.log"
	shell:
		# "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {params.threads} --out {params.oprefix}"
		"""
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --chr {params.chrom} --mlma --pheno {input.phenotypes} --mpheno {params.phenotype} --grm {params.grmprefix} --maf {params.maf} --autosome-num 29 --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
		"""

rule cojo:
	input:
		log = "output/{run_name}/gwas/{run_name}.{phenotype}.850K.log",
		assoc = "output/{run_name}/gwas/{run_name}.{phenotype}.850K.mlma.gz",
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"])
	params:
		n = nchooser,
		cojofile = "output/{run_name}/cojo/{run_name}.{phenotype}.850K.ma",
		bfile = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		cojop = config["cojop"],
		oprefix = "output/{run_name}/cojo/{run_name}.{phenotype}.850K",
		psrecord = "output/{run_name}/log/psrecord/cojo/cojo.{phenotype}.log",
	threads: config["grm_threads"]
	output:
		cojo = "output/{run_name}/cojo/{run_name}.{phenotype}.850K.jma"
	shell:
		"""
		export OMPI_MCA_btl_openib_if_include='mlx5_3:1'
		module load openmpi/openmpi-3.1.3-intel-16.0.2
		zcat {input.assoc} | awk '{{print $2, $4, $5, $6, $7, $8, $9, {params.n}}}' > {params.cojofile}
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.bfile} --autosome-num 29 --maf 0.01 --cojo-file {params.cojofile} --cojo-slct --cojo-p {params.cojop} --thread-num {threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 60
		"""
