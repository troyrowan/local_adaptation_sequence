
#This GRM creation will just be for autosomes as there are some issues with how you'd want to go about handing the X Chromosome https://cnsgenomics.com/software/gcta/#MakingaGRM
rule run_envGWAS:
	input:
		plink = expand("output/{{run_name}}/imputed_genotypes/{{run_name}}.850K.{suffix}",
		suffix = [".bed", ".bim", ".fam"]),
		id = "output/{run_name}/grm/{run_name}.850K.grm.id",
		bin = "output/{run_name}/grm/{run_name}.850K.grm.bin",
		nbin = "output/{run_name}/grm/{run_name}.850K.grm.N.bin"
	params:
		psrecord = "output/snakemake/log/{run_name}/psrecord/concat_vcf/concat_vcf.log",
		grm = "output/{run_name}/grm/{run_name}.850K.grm"
