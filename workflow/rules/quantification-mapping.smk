rule get_hs_genome:
    output:
        "results/refs/hs_genome.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="106",
    log:
        "logs/ensembl/get_genome.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.72.0/bio/reference/ensembl-sequence"


rule genome_index:
    input:
        "results/refs/hs_genome.fasta"
    output:
        "results/refs/hs_genome.fasta.fai"
    log:
        "logs/genome_index/hs_genome.log"
    threads: 4
    wrapper:
        "v1.3.2/bio/samtools/faidx"

rule kallisto_index:
    input:
        fasta = "resources/HLA-alleles/{hla}_gen.fasta"
    output:
        index = "results/kallisto-index/{hla}_gen.idx"
    params:
        extra = "--kmer-size=11"
    log:
        "logs/kallisto/index/{hla}_gen.log"
    threads: 2
    wrapper:
        "v0.86.0/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq = get_fastq_input,
        #fastq = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"],
        index = "results/kallisto-index/{hla}_gen.idx"
    output:
        directory('results/kallisto/quant_results_{sample}_{hla}')
    params:
        extra = "-b 10 --seed=42"
    log:
        "logs/kallisto/kallisto_quant_{sample}_{hla}.log"
    threads: 10
    wrapper:
        "v0.86.0/bio/kallisto/quant"

rule bwa_index:
    input:
        "results/refs/hs_genome.fasta"
    output:
        multiext("results/bwa-index/hs_genome", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index/hs_genome.log"
    params:
        prefix="results/bwa-index/hs_genome",
        algorithm="bwtsw",
    wrapper:
        "v0.86.0/bio/bwa/index"

rule bwa_mem:
    input:
        reads = get_fastq_input,
        #reads = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"],
        index = multiext("results/bwa-index/hs_genome", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "results/mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="results/bwa-index/hs_genome",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",             
        sort_order="coordinate", 
    threads: 16
    wrapper:
        "v0.86.0/bio/bwa/mem"

rule samtools_index:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    threads: 4
    wrapper:
        "v0.86.0/bio/samtools/index"
