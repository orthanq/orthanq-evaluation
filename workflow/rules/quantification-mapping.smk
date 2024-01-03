# ruleorder: tabix_exons > bcftools_merge

# rule bwakit_index:
#     output:
#         "results/bwakit-genome/hs38DH.fa",
#         "results/bwakit-genome/hs38DH.fa.alt"
#     log:
#         "logs/bwakit/bwakit-index.log"
#     conda:
#         "../envs/bwakit.yaml"
#     threads: 20
#     shell:
#         "mkdir -p results/bwakit-genome && cd results/bwakit-genome && run-gen-ref hs38DH && bwa index hs38DH.fa"

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

##orthanq_candidates generates both candidate variants and individual genome sequences for each locus(to be used for quantification indices in the downstream processing)
#the alleles are filtered out according to the following criteria:
#if the allele has a full name equivalent (e.g. A*24:02:02:02) in the database and it's higher than 0.05, it's IN.
#if the allele doesn't have a full name equivalent (e.g. A*24:02:02), then look for the first three fields, if still not then look for first two fields (A*24:02)
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
        "/projects/koesterlab/orthanq/orthanq/target/release/orthanq candidates --alleles {input.lineage_collection} --genome {input.genome} "
        "--output results/orthanq-candidates 2> {log}"

rule bgzip:
    input:
        "results/orthanq-candidates/candidates.vcf"
    output:
        "results/orthanq-candidates/candidates.vcf.gz"
    params:
        extra="", # optional
    threads: 1
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
        variants="results/orthanq-candidates/candidates.vcf",
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
    threads: 30
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

# rule kallisto_index:
#     input:
#         fasta = "resources/IMGTHLA-3.32.0-alpha/fasta/{hla}_gen.fasta"
#     output:
#         index = "results/kallisto-index/{hla}.idx"
#     params:
#         extra = "--kmer-size=31"
#     log:
#         "logs/kallisto/index/{hla}.log"
#     threads: 20
#     wrapper:
#         "v1.25.0/bio/kallisto/index"

# rule kallisto_quant:
#     input:
#         fastq = get_fastq_input,
#         index = "results/kallisto-index/{hla}.idx"
#     output:
#         directory('results/kallisto/quant_results_{sample}_{hla}')
#     params:
#         extra = "-b 5 --seed=42 --pseudobam"
#     log:
#         "logs/kallisto/kallisto_quant_{sample}_{hla}.log"
#     threads: 20
#     wrapper:
#         "v1.25.0/bio/kallisto/quant"