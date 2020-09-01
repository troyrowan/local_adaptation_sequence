rule grm_target:
	input:
		expand{"output/{run_name}/grm/{run_name}.850K.grm.bin",
		run_name = config["run_name"]}
#This GRM creation will just be for autosomes as there are some issues with how you'd want to go about handing the X Chromosome https://cnsgenomics.com/software/gcta/#MakingaGRM
rule concat_vcf:
	input:
		vcf = expand("/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{{run_name}}/imputed_genotypes/single_chrom/{{run_name}}.chr{chr}.reordered.vcf.gz", chr = list(range(1,30))), #Only Autosomes
		tabix = expand("/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/{{run_name}}/imputed_genotypes/single_chrom/{{run_name}}.chr{chr}.reordered.vcf.gz.tbi", chr = list(range(1,30)))
	params:
		psrecord = "output/snakemake/log/{run_name}/psrecord/concat_vcf/concat_vcf.log"
	output:
		vcf = "output/{run_name}/imputed_genotypes/{run_name}.850K.vcf.gz",
		tabix = "output/{run_name}/imputed_genotypes/{run_name}.850K.vcf.gz.tbi",
	shell: #List of all single-chromosome vcf files will be put onto the command line here for the concat command
		"""
		module load bcftools
		psrecord "bcftools concat {input.vcf} -O z -o {output.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5"""

rule convert_plink:
	input:
		vcf = "output/{run_name}/imputed_genotypes/{run_name}.850K.vcf.gz",
		tabix = "output/{run_name}/imputed_genotypes/{run_name}.850K.vcf.gz.tbi"
	params:
		psrecord = "output/snakemake/log/{run_name}/psrecord/convert_plink/convert_plink.log",
		oprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		threads=config["plink_threads"],
		mem=config["plink_mem"]
	output:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = [".bed", ".bim", ".fam"])
	shell:
		"""
		module load plink
		psrecord "plink --vcf {input.vcf} --threads {params.threads} --memory {params.mem} --cow --real-ref-alleles --make-bed --out {params.oprefix} " --log {params.psrecord} --include-children --interval 5"""
#This rule allows us to massively speed up GRM creation
rule build_grm_chunks:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = [".bed", ".bim", ".fam"])
	params:
		psrecord = "output/snakemake/log/{run_name}/psrecord/build_grm_chunks/build_grm_chunks.part{part}.log",
		iprefix = "output/{run_name}/imputed_genotypes/{run_name}.850K",
		oprefix = "output/{run_name}/grm/{run_name}.850K",
		part = "{part}",
		threads=config["grm_threads"]
	output:
		chunks = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{{part}}.{suffix}",
		suffix = ["grm.id", "grm.bin", "grm.N.bin"])
	shell:
		"""
		psrecord "code/gcta_1.93.2beta/gcta64 --bfile {input.iprefix} --autosome-num 29 --make-grm-part 50 {params.part} --thread-num {params.threads} --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5"""

rule concat_grm_chunks:
	input:
		id = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{part}.grm.id",
		part = list(range(1,51))),
		bin = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{part}.grm.bin",
		part = list(range(1,51))),
		nbin = = expand("output/{{run_name}}/grm/{{run_name}}.850K.part_50_{part}.grm.N.bin",
		part = list(range(1,51)))
	params:
		psrecord = "output/snakemake/log/{run_name}/psrecord/concat_grm_chunks/concat_grm_chunks.log"
	output:
		id = "output/{run_name}/grm/{run_name}.850K.grm.id",
		bin = "output/{run_name}/grm/{run_name}.850K.grm.bin",
		nbin = "output/{run_name}/grm/{run_name}.850K.grm.N.bin"
	shell:
	"""
	psrecord "cat {input.id} > {output.id}; cat {input.bin} > {output.bin}; cat {input.nbin} > {output.nbin}""""
