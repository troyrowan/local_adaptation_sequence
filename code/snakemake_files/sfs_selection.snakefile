#Pulls out subsetted individuals and their respectve haplotypes
#Comment
# import os
# #Make log directories if they don't exist
# for x in expand("output/{run_name}/log/slurm_out/{rule}", run_name = ["200910_RAN", "200907_SIM"], rule = ["filtering", "tabix", "SweeD", "RAiSD"]):
# 	os.makedirs(x, exist_ok = True)
#snakemake -s code/snakemake_files/selscan.snakefile --cluster-config code/cluster/selscan.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem}" --rerun-incomplete -np


rule selscan:
	input:
		expand("RAiSD_Report.200907_SIM.{dataset}.chip.chr{chr}",
		chr = list(range(1,30)),
		dataset = ["somesim", "purebred"])
		# expand("output/{run_name}/sfs_selection/SweeD_Report.{run_name}.chr{chr}",
		# run_name = ["200910_RAN", "200907_SIM"],
		# chr = list(range(1,30))),
# rule filtering:
# 	input:
# 		vcf = "/storage/htc/schnabellab/results/9913/wgs/1kbulls_ars1.2/vcf/200627_Run8_TauInd/Chr{chr}-Run8-TAUIND-raw-toDistribute.vcf.gz",
# 		tbi = "/storage/htc/schnabellab/results/9913/wgs/1kbulls_ars1.2/vcf/200627_Run8_TauInd/Chr{chr}-Run8-TAUIND-raw-toDistribute.vcf.gz.tbi",
# 		samples = "data/{run_name}/{run_name}.sequenced_ids.txt"
# 	output:
# 		vcf = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf"
# 	threads: 8
# 	shell:
# 		"""
# 		module load bcftools
# 		bcftools view  --max-alleles 2 -v snps {input.vcf} -S {input.samples} --threads {threads} -o {output.vcf}
# 		"""
# # rule tabix:
# # 	input:
# # 		vcf = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf.gz"
# # 	output:
# # 		tbi = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf.gz.tbi"
# # 	shell:
# # 		"""
# # 		module load bcftools
# # 		tabix {input.vcf}
# # 		"""
#
# griddict = {"1":"31707", "2":"27246", "3":"24201", "4":"23983", "5":"24018", "6":"23559", "7":"22136", "8":"22656", "9":"20929", "10":"20661", "11":"21396", "12":"17442", "13":"16693", "14":" 16478", "15":"16999", "16":"16201", "17":"14633", "18":"13164", "19":"12690", "20":"14394", "21":"13971", "22":"12154", "23":"10500", "24":"12461", "25":"8470", "26":"10398", "27":"9122", "28":"9187", "29":"10220"}
#
# def grid_chooser(WC):
# 	chunks = griddict[WC.chr]
# 	return chunks
#
# rule SweeD:
# 	input:
# 		vcf = "output/{run_name}/sfs_selection/genotypes/{run_name}.chr{chr}.1kbulls.vcf"
# 	params:
# 		name = "{run_name}.chr{chr}",
# 		report = "SweeD_Report.{run_name}.chr{chr}",
# 		info = "SweeD_Info.{run_name}.chr{chr}",
# 		odir = "output/{run_name}/sfs_selection/",
# 		chunks = grid_chooser
# 	threads: 24
# 	output:
# 		"output/{run_name}/sfs_selection/SweeD_Report.{run_name}.chr{chr}"
# 	shell:
# 		"""
# 		code/sweed/SweeD-P -name {params.name} -input {input.vcf} -threads {threads} -grid {params.chunks} -folded
# 		mv {params.report} {params.odir}
# 		mv {params.info} {params.odir}
# 		"""
# chrlen_dict = {"1":"158533849", "2":"136228794", "3":"121005064", "4":"119916527", "5":"120089179", "6":"117795627", "7":"110681041", "8":"113279970", "9":"104647085", "10":"103307330", "11":"106979715", "12":"87209924", "13":"83467450", "14":"82397735", "15":"84995892", "16":"81005832", "17":"73165306", "18":"65818438", "19":"63448886", "20":"71972210", "21":"69856109", "22":"60772208", "23":"52498378", "24":"62305915", "25":"42349150", "26":"51989682", "27":"45611753", "28":"45937039", "29":"51098434"}
#
# nsnp_dict = {"1":"8027715", "2":"6693967", "3":"5887012", "4":"6178704", "5":"5979847", "6":"6169700", "7":"5315613", "8":"5578741", "9":"5383288", "10":"5055416", "11":"5182259", "12":"4775141", "13":"3968719", "14":"4017854", "15":"4514490", "16":"3974042", "17":"3703813", "18":"3156765", "19":"2999984", "20":"3663523", "21":"3539306", "22":"2933562", "23":"2881291", "24":"3215028", "25":"2112169", "26":"2662838", "27":"2500469", "28":"2431687", "29":"2727425"}
#
# def chrlen_chooser(WC):
# 	length = chrlen_dict[WC.chr]
# 	return length
#
# def nsnp_chooser(WC):
# 	nsnps = nsnp_dict[WC.chr]
# 	return nsnps

# rule RAiSD:
# 	input:
# 		vcf = "genotypes/{run_name}.chr{chr}.1kbulls.vcf"
# 	params:
# 		name = "{run_name}",
# 		#mover = "RAiSD_Report.{run_name}.{chr}",
# 		#outdir = "output/{run_name}/sfs_selection/",
# 		chrlen = chrlen_chooser,
# 		snps = nsnp_chooser
# 	output:
# 		"RAiSD_Report.{run_name}.{chr}"
# 	shell:
# 		"""
# 		code/RAiSD/raisd-master/bin/release/RAiSD -n {params.name} -I {input.vcf} -B {params.chrlen} {params.snps} -f -R
# 		"""
chrlen_dict = {"1":"158533849", "2":"136228794", "3":"121005064", "4":"119916527", "5":"120089179", "6":"117795627", "7":"110681041", "8":"113279970", "9":"104647085", "10":"103307330", "11":"106979715", "12":"87209924", "13":"83467450", "14":"82397735", "15":"84995892", "16":"81005832", "17":"73165306", "18":"65818438", "19":"63448886", "20":"71972210", "21":"69856109", "22":"60772208", "23":"52498378", "24":"62305915", "25":"42349150", "26":"51989682", "27":"45611753", "28":"45937039", "29":"51098434"}

def chrlen_chooser(WC):
	length = chrlen_dict[WC.chr]
	return length

chip_nsnp_dict = {"1":"325791", "2":"225115", "3":"74580", "4":"37583", "5":"38692", "6":"37584", "7":"37416", "8":"36315", "9":"32873", "10":"33535", "11":"35052", "12":"27720", "13":"26168", "14":"27166", "15":"28804", "16":"26521", "17":"24317", "18":"24190", "19":"23149", "20":"22737", "21":"22853", "22":"20174", "23":"18492", "24":"19667", "25":"15804", "26":"16579", "27":"14712", "28":"14206", "29":"34270"}

def chip_nsnp_chooser(WC):
	nsnps = chip_nsnp_dict[WC.chr]
	return nsnps


rule chip_RAiSD_RAN:
	input:
		vcf = "output/200910_RAN/subsets/200910_RAN.all.chr{chr}.vcf"
	params:
		name = "200910_RAN.chip.chr{chr}",
		#mover = "RAiSD_Report.{run_name}.{chr}",
		#outdir = "output/{run_name}/sfs_selection/",
		chrlen = chrlen_chooser,
		snps = chip_nsnp_chooser
	output:
		"RAiSD_Report.200910_RAN.chip.chr{chr}"
	shell:
		"""
		code/RAiSD/raisd-master/bin/release/RAiSD -n {params.name} -I {input.vcf} -B {params.chrlen} {params.snps} -f -R
		"""

rule chip_RAiSD_SIM:
	input:
		vcf = "output/200907_SIM/subsets/200907_SIM.{dataset}.chr{chr}.vcf"
	params:
		name = "200907_SIM.{dataset}.chip.chr{chr}",
		#mover = "RAiSD_Report.{run_name}.{chr}",
		#outdir = "output/{run_name}/sfs_selection/",
		chrlen = chrlen_chooser,
		snps = chip_nsnp_chooser
	output:
		"RAiSD_Report.200907_SIM.{dataset}.chip.chr{chr}"
	shell:
		"""
		code/RAiSD/raisd-master/bin/release/RAiSD -n {params.name} -I {input.vcf} -B {params.chrlen} {params.snps} -f -R
		"""
