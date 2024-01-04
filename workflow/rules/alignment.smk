rule genome_index:
    input:
        "resources/reference/NC_045512.2.fasta"
    output:
        "resources/reference/NC_045512.2.fasta.fai"
    log:
        "logs/genome_index/genome_index.log"
    threads: 4
    wrapper:
        "v1.25.0/bio/samtools/faidx" #update wrapper version

rule orthanq_candidates:
    input:
        lineage_collection = config["lineage_collection"],
        genome = "resources/reference/NC_045512.2.fasta",
    output:
        candidate_variants = "results/orthanq-candidates/candidates.vcf",
        #fasta = expand("results/orthanq-candidates/{hla}.fasta",hla=loci)
    conda:
        "../envs/orthanq.yaml"
    log:
        "logs/orthanq-candidates/candidates.log"
    shell:
        "../../orthanq/target/debug/orthanq candidates virus --alleles {input.lineage_collection} --genome {input.genome} "
        "--output results/orthanq-candidates 2> {log}"

rule bgzip:
    input:
        "results/orthanq-candidates/candidates.vcf"
    output:
        "results/orthanq-candidates/candidates.vcf.gz"
    params:
        extra="", # optional
    log:
        "logs/bgzip/candidates.log",
    wrapper:
        "v1.23.5-27-g1638662a/bio/bgzip"

rule tabix:
    input:
        "results/orthanq-candidates/candidates.vcf.gz",
    output:
        "results/orthanq-candidates/candidates.vcf.gz.tbi"
    log:
        "logs/tabix/candidates.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        "v1.23.5-27-g1638662a/bio/tabix/index"

rule vg_autoindex:
    input:
        genome="resources/reference/NC_045512.2.fasta",
        variants="results/orthanq-candidates/candidates.vcf.gz",
        variants_tbi="results/orthanq-candidates/candidates.vcf.gz.tbi",
    output:
        "results/vg/autoindex/idx.giraffe.gbz"
    log:
        "logs/vg/autoindex.log"
    conda:
        "../envs/vg.yaml"
    threads: 40
    shell:
        "vg autoindex --workflow giraffe -r {input.genome} -v {input.variants} -p results/vg/autoindex/idx -t {threads}"

rule vg_giraffe:
    input:
        reads = get_fastq_input,
        idx = "results/vg/autoindex/idx.giraffe.gbz"
    output:
        "results/vg/alignment/{sample}_vg.bam"
    log:
        "logs/vg/alignment/{sample}.log"
    conda:
        "../envs/vg.yaml"
    threads: 40
    shell:
        "vg giraffe -Z results/vg/autoindex/idx.giraffe.gbz -f {input.reads[0]} -f {input.reads[1]} --output-format BAM -t {threads}  > {output} 2> {log}"

rule samtools_sort:
    input:
        "results/vg/alignment/{sample}_vg.bam"
    output:
        "results/vg/alignment/{sample}_vg.sorted.bam"
    log:
        "logs/samtools_sort/{sample}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v1.22.0/bio/samtools/sort"

rule samtools_index:
    input:
        "results/vg/alignment/{sample}_vg.sorted.bam"
    output:
        "results/vg/alignment/{sample}_vg.sorted.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    threads: 4
    wrapper:
        "v1.22.0/bio/samtools/index"
