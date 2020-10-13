#Pulls out subsetted individuals and their respectve haplotypes
#Comment
import os
#Make log directories if they don't exist
for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = ["200910_RAN", "200907_SIM"], rule = ["filtering", "tabix", "SweeD"]):
	os.makedirs(x, exist_ok = True)
#snakemake -s code/snakemake_files/selscan.snakefile --cluster-config code/cluster/selscan.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem}" --rerun-incomplete -np

rule selscan:
	input:
		expand("output/{run_name}/sfs_selection/SweeD_Report.{run_name}.chr{chr}",
		run_name = ["200910_RAN", "200907_SIM"],
		chr = list(range(1,30)))

rule filtering:
	input:
		vcf = "/storage/htc/schnabellab/results/9913/wgs/1kbulls_ars1.2/vcf/200627_Run8_TauInd/Chr{chr}-Run8-TAUIND-raw-toDistribute.vcf.gz",
		tbi = "/storage/htc/schnabellab/results/9913/wgs/1kbulls_ars1.2/vcf/200627_Run8_TauInd/Chr{chr}-Run8-TAUIND-raw-toDistribute.vcf.gz.tbi",
		samples = "data/{run_name}/{run_name}.sequenced_ids.txt"
	output:
		vcf = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf.gz"
	threads: 8
	shell:
		"""
		module load bcftools
		bcftools view -i 'AC>=20' -f PASS,VQSRTrancheINDEL90.00to99.00,VQSRTrancheSNP90.00to99.00 --max-alleles 2 -v snps {input.vcf} -S {input.samples} --threads {threads} -O z -o {output.vcf}
		"""
rule tabix:
	input:
		vcf = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf.gz"
	output:
		tbi = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		tabix {input.vcf}
		"""

griddict = {"1":"31707", "2":"27246", "3":"24201", "4":"23983", "5":"24018", "6":"23559", "7":"22136", "8":"22656", "9":"20929", "10":"20661", "11":"21396", "12":"17442", "13":"16693", "14":" 16478", "15":"16999", "16":"16201", "17":"14633", "18":"13164", "19":"12690", "20":"14394", "21":"13971", "22":"12154", "23":"10500", "24":"12461", "25":"8470", "26":"10398", "27":"9122", "28":"9187", "29":"10220"}

def grid_chooser(WC):
	chunks = griddict[WC.chr]
	return chunks

rule SweeD:
	input:
		vcf = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf"
	params:
		name = "{run_name}.chr{chr}",
		report = "SweeD_Report.{run_name}.chr{chr}",
		info = "SweeD_Info.{run_name}.chr{chr}",
		odir = "output/{run_name}/sfs_selection/",
		chunks = grid_chooser
	threads: 24
	output:
		"output/{run_name}/sfs_selection/SweeD_Report.{run_name}.chr{chr}"
	shell:
		"""
		code/sweed/SweeD-P -name {params.name} -input {input.vcf} -threads {threads} -grid {params.chunks} -folded
		mv {params.report} {params.odir}
		mv {params.info} {params.odir}
		"""
