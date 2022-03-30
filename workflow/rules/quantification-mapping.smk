rule kallisto_index:
    input:
        fasta = "resources/HLA-alleles/hla_gen.fasta"
    output:
        index = "results/kallisto-index/hla_gen.idx"
    params:
        extra = "--kmer-size=31"
    log:
        "logs/kallisto/index/hla_gen.log"
    threads: 8
    wrapper:
        "v0.86.0/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"],
        index = "results/kallisto-index/hla_gen.idx"
    output:
        directory('results/kallisto/quant_results_{sample}')
    params:
        extra = "-b 100 --seed=42 --pseudobam"
    log:
        "logs/kallisto/kallisto_quant_{sample}.log"
    threads: 10
    wrapper:
        "v0.86.0/bio/kallisto/quant"

rule bwa_index:
    input:
        "resources/hs_genome/hs_genome.fasta"
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
        reads = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"],
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
    threads: 10
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