# This smk aligns reads to the reference and post-processes the output to be used as in input to varlociraptor.
# The alignment strategy is as followed:
# first, align against with bwa against the primary assembly
# second, extract all reads that belong to HLA loci, use samtools view and then fastq
# last, map only to those extracted reads

# Step 1: align reads by bwa
rule bwa_index:
    input:
        "results/refs/hs_genome.fasta"
    output:
        idx=multiext("results/bwa-index/hs_genome", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index/hs_genome.log"
    params:
        prefix="results/bwa-index/hs_genome",
        algorithm="bwtsw",
    wrapper:
        "v2.0.0/bio/bwa/index" 

rule bwa_mem:
    input:
        reads = get_fastq_input,
        #reads = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"],
        idx = multiext("results/bwa-index/hs_genome", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "results/bwa_alignment/{sample}_mapped.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:    
        "benchmarks/bwa_mem/{sample}.tsv"
    params:
        index="results/bwa-index/hs_genome",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",             
        sort_order="coordinate", 
    threads: 40
    wrapper:
        "v2.0.0/bio/bwa/mem"

rule samtools_index_bwa:
    input:
        "results/bwa_alignment/{sample}_mapped.bam"
    output:
        "results/bwa_alignment/{sample}_mapped.bai"
    log:
        "logs/samtools_index_bwa/{sample}.log"
    threads: 10
    wrapper:
        "v2.0.0/bio/samtools/index"

# Step 2: extract reads that map to HLA loci
rule samtools_view:
    input:
        "results/bwa_alignment/{sample}_mapped.bam",
        "results/bwa_alignment/{sample}_mapped.bai",
        regions="resources/HLA_regions/regions.bed" #only contains the regions for HLA-A,HLA-B, HLA-C and HLA-DQB1
    output:
        bam="results/bwa_alignment/{sample}_mapped_extracted.bam",
        idx="results/bwa_alignment/{sample}_mapped_extracted.bai"
    log:
        "logs/samtools-view/{sample}.log",
    benchmark:    
        "benchmarks/samtools_view_extract_HLA/{sample}.tsv"
    params:
        extra=lambda wc, input: "-L {}".format(input.regions)
    resources:
        mem_mb=60000
    threads: 40
    wrapper:
        "v2.0.0/bio/samtools/view"

rule samtools_fastq_separate:
    input:
        "results/bwa_alignment/{sample}_mapped_extracted.bam"
    output:
        "results/extracted_reads/{sample}.1.fq",
        "results/extracted_reads/{sample}.2.fq",
    log:
        "logs/extracted_reads/{sample}.separate.log",
    benchmark:    
        "benchmarks/bam_to_fastq/{sample}.tsv"
    params:
        sort="-m 4G",
        fastq="-n",
    threads: 40
    wrapper:
        "v2.0.0/bio/samtools/fastq/separate"

# Step 3: map extracted reads to the pangenome with vg giraffe
rule vg_giraffe:
    input:
        reads_1 = "results/extracted_reads/{sample}.1.fq",
        reads_2 = "results/extracted_reads/{sample}.2.fq",
        idx = "resources/vg-pangenome/hprc-v1.0-mc-grch38.xg"
    output:
        "results/vg/alignment/{sample}_vg.bam"
    log:
        "logs/vg/alignment/{sample}.log"
    benchmark:    
        "benchmarks/vg_giraffe/{sample}.tsv"
    conda:
        "../envs/vg.yaml"
    threads: 40
    shell:
        "vg giraffe -x {input.idx} -f {input.reads_1} -f {input.reads_2} --max-multimaps 3 --output-format BAM -t {threads}  > {output} 2> {log}"

rule samtools_sort:
    input:
        "results/vg/alignment/{sample}_vg.bam"
    output:
        bam="results/vg/alignment/{sample}_vg.sorted.bam",
        idx="results/vg/alignment/{sample}_vg.sorted.bai"
    log:
        "logs/samtools_sort/{sample}.log",
    benchmark:    
        "benchmarks/samtools_sort/{sample}.tsv" 
    # params:
    #     extra="-m 4G",
    threads: 40
    wrapper:
        "v1.22.0/bio/samtools/sort"

# modify the header for chromosome names to be compatible with the reference genome that we acquire from ensembl
rule reheader:
    input:
        "results/vg/alignment/{sample}_vg.sorted.bam"
    output:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam"
    log:
        "logs/samtools_reheader/{sample}.log"
    benchmark:    
        "benchmarks/samtools_reheader/{sample}.tsv" 
    conda:
        "../envs/samtools.yaml"
    threads: 40
    shell:
        "samtools view -H {input} | sed 's/GRCh38.chr//g' | samtools reheader - {input} > {output} 2> {log}"
    
rule samtools_index_after_reheader:
    input:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam"
    output:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam.bai"
    log:
        "logs/samtools_index_after_reheader/{sample}.log"
    benchmark:    
        "benchmarks/samtools_index_after_reheader/{sample}.tsv"
    threads: 40
    wrapper:
        "v1.22.0/bio/samtools/index"

#keep only primary chromosomes to prevent varlociraptor throwing errors
rule keep_only_primary_chr:
    input:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam",
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam.bai",
    output:
        bam="results/vg/alignment/{sample}_vg.sorted.reheadered.extracted.bam",
        idx="results/vg/alignment/{sample}_vg.sorted.reheadered.extracted.bam.bai"
    log:
        "logs/samtools-view-primary-chr/{sample}.log",
    benchmark:    
        "benchmarks/samtools_view_primary_chr/{sample}.tsv"
    params:
        region="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
    threads: 40
    wrapper:
        "v2.0.0/bio/samtools/view"
