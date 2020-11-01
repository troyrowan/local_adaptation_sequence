#Pulls out subsetted individuals and their respectve haplotypes
#Comment
import os
#Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = ["200910_RAN", "200907_SIM"], rule = ["tabix", "extract_vcf", "iHS", "nSL", "xpehh", "norm", "xp_norm", "xpnsl"]):
	os.makedirs(x, exist_ok = True)
#snakemake -s code/snakemake_files/selscan.snakefile --cluster-config code/cluster/selscan.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem}" --rerun-incomplete -np

rule selscan:
	input:
		expand("output/{run_name}/selscan/{test}/{run_name}.{subset}.{test}.out.100bins.norm.gz",
		test = ["nsl"],
		run_name = ["200910_RAN"],
		subset = ["all"]),
		# expand("output/{run_name}/selscan/{test}/{run_name}.{subset}.{test}.out.100bins.norm.gz",
		# test = ["nsl"],
		# run_name = ["200910_RAN", "200907_SIM"],
		# subset = ["oldest_1000", "oldest_5000", "young_1000", "young_5000","purebred"]),
		# expand("output/{run_name}/selscan/{test}/{run_name}.{test}.out.norm.gz",
		# test = ["xpnsl"],
		# run_name = ["200910_RAN", "200907_SIM"])


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
		temp("output/{run_name}/subsets/{run_name}.{subset}.chr{chr}.vcf.gz")
	shell:
		"""
		module load bcftools
		bcftools view {input.vcf} -S {input.subsets} --force-samples -O z -o {output}
		"""

rule iHS:
	input:
		"output/{run_name}/subsets/{run_name}.{subset}.chr{chr}.vcf.gz"
	params:
		oprefix = "output/{run_name}/selscan/ihs/{run_name}.{subset}.chr{chr}"
	threads: 16
	output:
		"output/{run_name}/selscan/ihs/{run_name}.{subset}.chr{chr}.ihs.out"
	shell:
		"./code/selscan/releases/selscan-linux-1.3.0/selscan --vcf {input} --pmap --maf 0.01 --threads {threads} --ihs --out {params.oprefix}"

rule nSL:
	input:
		"output/{run_name}/subsets/{run_name}.{subset}.chr{chr}.vcf.gz"
	params:
		oprefix = "output/{run_name}/selscan/nsl/{run_name}.{subset}.chr{chr}"
	threads: 16
	output:
		"output/{run_name}/selscan/nsl/{run_name}.{subset}.chr{chr}.nsl.out"
	shell:
		"./code/selscan/releases/selscan-linux-1.3.0/selscan --vcf {input} --pmap --maf 0.01 --threads {threads} --nsl --out {params.oprefix}"

rule xpehh:
	input:
		old = "output/{run_name}/subsets/{run_name}.oldest_500_perbreeder.chr{chr}.vcf.gz",
		young = "output/{run_name}/subsets/{run_name}.youngest_500_perbreeder.chr{chr}.vcf.gz"
	params:
		oprefix = "output/{run_name}/selscan/xpehh/{run_name}.chr{chr}"
	threads: 16
	output:
		"output/{run_name}/selscan/xpehh/{run_name}.chr{chr}.xpehh.out"
	shell:
		"./code/selscan/releases/selscan-linux-1.3.0/selscan --vcf {input.young} --vcf-ref {input.old} --pmap --maf 0.01 --threads {threads} --xpehh --out {params.oprefix}"

rule xpnsl:
	input:
		old = "output/{run_name}/subsets/{run_name}.oldest_5000.chr{chr}.vcf.gz",
		young = "output/{run_name}/subsets/{run_name}.young_5000.chr{chr}.vcf.gz"
	params:
		oprefix = "output/{run_name}/selscan/xpnsl/{run_name}.chr{chr}"
	threads: 16
	output:
		"output/{run_name}/selscan/xpnsl/{run_name}.chr{chr}.xpnsl.out"
	shell:
		"./code/selscan/releases/selscan-linux-1.3.0/selscan --vcf {input.young} --vcf-ref {input.old} --pmap --maf 0.01 --threads {threads} --xpnsl --out {params.oprefix}"

rule norm:
	input:
		out = expand("output/{{run_name}}/selscan/{{test}}/{{run_name}}.{{subset}}.chr{chr}.{{test}}.out",
		chr = list(range(1,30)))
	params:
		files = "output/{run_name}/selscan/{test}/{run_name}.{subset}.chr{1..29}.{test}.out",
		test = "--{test}"
	output:
		"output/{run_name}/selscan/{test}/{run_name}.{subset}.chr{chr}.{test}.out.100bins.norm"
	shell:
		"./code/selscan/releases/selscan-linux-1.3.0/norm --files {params.files} {params.test}"

rule xp_norm:
	input:
		out = expand("output/{{run_name}}/selscan/{{test}}/{{run_name}}.chr{chr}.{{test}}.out",
		chr = list(range(1,30)))
	params:
		files = "output/{run_name}/selscan/{test}/{run_name}.chr{1..29}.{test}.out",
		test = "--{test}"
	output:
		"output/{run_name}/selscan/{test}/{run_name}.chr{chr}.{test}.out.norm"
	shell:
		"./code/selscan/releases/selscan-linux-1.3.0/norm --files {params.files} {params.test}"

rule cat_norm:
	input:
		norm = expand("output/{{run_name}}/selscan/{{test}}/{{run_name}}.{{subset}}.chr{chr}.{{test}}.out.100bins.norm",
		chr = list(range(1,30)))
	params:
		files = "output/{run_name}/selscan/{test}/{run_name}.{subset}.chr{1..29}.{test}.out.100bins.norm",
		tmp = "output/{run_name}/selscan/{test}/{run_name}.{subset}.{test}.out.100bins.norm"
	output:
		"output/{run_name}/selscan/{test}/{run_name}.{subset}.{test}.out.100bins.norm.gz"
	shell:
		"cat {params.files} > {params.tmp}; pigz {params.tmp}"

rule cat_xp_norm:
	input:
		norm = expand("output/{{run_name}}/selscan/{{test}}/{{run_name}}.chr{chr}.{{test}}.out.norm",
		chr = list(range(1,30)))
	params:
		files = "output/{run_name}/selscan/{test}/{run_name}.chr{1..29}.{test}.out.norm",
		tmp = "output/{run_name}/selscan/{test}/{run_name}.{test}.out.norm"
	output:
		"output/{run_name}/selscan/{test}/{run_name}.{test}.out.norm.gz"
	shell:
		"cat {params.files} > {params.tmp}; pigz {params.tmp}"
