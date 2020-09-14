hisat_index = '~/Lab/Genomes/pombe/hisat/sp_hisatIndex_v2.30'
gatk_jar = "/opt/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar"
gatk_ref = "/home/zdwyer/Lab/Genomes/pombe/sp_genome_v2.30.fa"

rule all:
	input:
		"snps/{sample}.vcf.idx"

rule trim:
	input:
		R1="fastq/{sample}_R1.fastq.gz",
		R2="fastq/{sample}_R2.fastq.gz"
	output:
		R1="trim/{sample}_R1.fastq.gz",
		R2="trim/{sample}_R2.fastq.gz",
		HTML="trim_report/{sample}.html",
		JSON="trim_report/{sample}.json"
	log:
		"logs/fastp/{sample}.log"
	threads: 2
	benchmark:
		"benchmarks/fastp/{sample}.txt"
	shell:
		'fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --adapter_sequence CCGAGCCCACGAGAC --adapter_sequence_r2 GACGCTGCCGACGA --html {output.HTML} --json {output.JSON} --thread {threads} 2> {log}'

rule align:
	input:
		R1="trim/{sample}_R1.fastq.gz",
		R2="trim/{sample}_R2.fastq.gz",
		index=expand("/home/zdwyer/Lab/Genomes/pombe/hisat/sp_hisatIndex_v2.30.{number}.ht2", number=range(1,9))
	output:
		bam="align/{sample}.bam"
	log:
		"logs/hisat/{sample}.log"
	threads: 2
	benchmark:
		"benchmarks/hisat/{sample}.txt"
	shell:
		'hisat2 -p {threads} --maxins 3000 --phred33 --new-summary --no-spliced-alignment --rg-id {wildcards.sample} --rg ID:{wildcards.sample} --rg LB:20181202 --rg PL:illumina --rg SM:{wildcards.sample} --rg PU:{wildcards.sample} -x %s -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh - -q 5 -F 260 | samtools sort - -o {output.bam}' % (hisat_index)

rule index_bam:
	input:
		'align/{sample}.bam'
	output:
		'align/{sample}.bam.bai'
	shell:
		'samtools index {input}'

rule call_snps:
	input:
		bam='align/{sample}.bam',
		bai='align/{sample}.bam.bai',
		ref='/home/zdwyer/Lab/Genomes/pombe/sp_genome_v2.30.fa'
	output:
		'snps/{sample}.vcf.idx'
	log:
		"logs/gatk/{sample}.log"
	threads: 2
	benchmark:
		"benchmarks/gatk/{sample}.txt"
	shell:
		'java -XX:ParallelGCThreads=1 -jar %s HaplotypeCaller --native-pair-hmm-threads {threads} -R {input.ref} -ploidy 1 -I {input.bam} -O {output} 2> {log}' % (gatk_jar)