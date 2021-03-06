import os
# Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("output/{run_name}/log/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

rule grm_target:
	input:
		expand("output/{run_name}/grm/{run_name}.850K.grm.N.bin",
		run_name = config["run_name"])
#This GRM creation will just be for autosomes as there are some issues with how you'd want to go about handing the X Chromosome https://cnsgenomics.com/software/gcta/#MakingaGRM
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
#issues with this rule making the mergelist. Did it manually to get things to run...

rule concat_plink:
	input:
		plink = expand("output/{{run_name}}/plink_convert/{{run_name}}.chr{chr}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"],
		chr = list(range(1,30))) #Only autosomes here for GRM creation
	params:
		psrecord = "output/{run_name}/log/psrecord/concat_plink/concat_plink.log",
		oprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		chromdir = "output/{run_name}/plink_convert/*bed",
		dir = "output/{run_name}/imputed_genotypes",
		threads=config["plink_threads"],
		mem=config["plink_mem"],
		#list = "output/{run_name}/plink_convert/{run_name}.850K.mergelist.txt"
	output:
		list = "output/{run_name}/plink_convert/{run_name}.850K.mergelist.txt",
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"])
	shell:
		"""
		module load plink
		ls {params.chromdir} | tr "\\t" "\\n" | sed 's/.bed//g' > {output.list}
		psrecord "plink --merge-list {output.list} --threads {params.threads} --memory {params.mem} --cow --real-ref-alleles --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""

#This rule allows us to massively speed up GRM creation
rule build_grm_chunks:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"])
	params:
		psrecord = "output/{run_name}/log/psrecord/build_grm_chunks/build_grm_chunks.part{part}.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		oprefix = "output/{run_name}/grm/{run_name}.850K",
		part = "{part}",
		maf = config["maf"],
		threads=config["grm_threads"]
	output:
			chunks = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{{part}}.{suffix}",
			suffix = ["grm.id", "grm.bin", "grm.N.bin"])
	shell:
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --maf {params.maf} --make-grm-part 50 {params.part} --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""

rule concat_grm_chunks:
	input:
		id = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{part}.grm.id",
		part = config["parts"]),
		bin = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{part}.grm.bin",
		part = config["parts"]),
		nbin = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{part}.grm.N.bin",
		part = config["parts"])
	output:
		id = "output/{run_name}/grm/{run_name}.850K.grm.id",
		bin = "output/{run_name}/grm/{run_name}.850K.grm.bin",
		nbin = "output/{run_name}/grm/{run_name}.850K.grm.N.bin"
	shell:
		"""
		cat {input.id} > {output.id}; cat {input.bin} > {output.bin}; cat {input.nbin} > {output.nbin}
		"""

rule single_chrom_grm:
	input:
		plink = expand("output/{{run_name}}/plink_convert/{{run_name}}.chr{{chr}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"])
	params:
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.chr{chr}.850K",
		oprefix = "output/{run_name}/grm/{run_name}.chr{chr}.850K",
		threads = config["grm_threads"],
		maf = config["maf"],
		psrecord = "output/{run_name}/log/psrecord/single_chrom_grm/single_chrom_grm.chr{chr}.log"
	output:
		id = "output/{run_name}/grm/{run_name}.chr{chr}.850K.grm.id",
		bin = "output/{run_name}/grm/{run_name}.chr{chr}.850K.grm.bin",
		nbin = "output/{run_name}/grm/{run_name}.chr{chr}.850K.grm.N.bin"
	shell:
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {params.iprefix} --autosome-num 29 --maf {params.maf} --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 30
		"""
rule gemma_grm:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = ["bed", "bim", "fam"])
	params:
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		odir = "output/{run_name}/grm/",
		oprefix = "{run_name}.gemma.850K",
		psrecord = "output/{run_name}/log/psrecord/gemma_grm/gemma_grm.log"
	output:
		gemma = "output/{run_name}/grm/{run_name}.gemma.850K.sXX.txt"
	shell:
		"""
		psrecord "gemma -bfile {params.iprefix} -gk 2 -outdir {params.oprefix} -o {params.oprefix}" --log {params.psrecord} --interval 30
 		"""
