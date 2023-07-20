# This smk generates candidate variants belonging to HLA alleles using orthanq candidates.

# Genome acquisition and indexing of it
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
        "v1.25.0/bio/reference/ensembl-sequence"

rule genome_index:
    input:
        "results/refs/hs_genome.fasta"
    output:
        "results/refs/hs_genome.fasta.fai"
    log:
        "logs/genome_index/bwa_hs_genome.log"
    threads: 4
    wrapper:
        "v1.25.0/bio/samtools/faidx" #update wrapper version

# Generating candidate variants with Orthanq
# Alleles are filtered out according to the following criteria:
# if the allele has a full name equivalent (e.g. A*24:02:02:02) in the database and it's higher than 0.05, it's IN.
# if the allele doesn't have a full name equivalent (e.g. A*24:02:02), then look for the first three fields, if still not then look for first two fields (A*24:02)
rule orthanq_candidates:
    input:
        alleles = config["allele_db"],
        genome = "results/refs/hs_genome.fasta",
        xml = config["allele_db_xml"],
        allele_freqs = "resources/allele_freqs/allele_frequencies.csv"
    output:
        candidate_variants = expand("results/orthanq-candidates/{hla}.vcf",hla=loci),
        #fasta = expand("results/orthanq-candidates/{hla}.fasta",hla=loci)
    conda:
        "../envs/orthanq.yaml"
    log:
        "logs/orthanq-candidates/candidates.log"
    shell:
        "../orthanq/target/release/orthanq candidates --alleles {input.alleles} --genome {input.genome} "
        "--xml {input.xml} --allele-freq {input.allele_freqs} --wes --output results/orthanq-candidates 2> {log}" # --wes option for protein level hla type variant generation, --wgs for individual types 

rule bgzip:
    input:
        "results/orthanq-candidates/{hla}.vcf",
    output:
        "results/orthanq-candidates/{hla}.vcf.gz",
    params:
        extra="", # optional
    threads: 1
    log:
        "logs/bgzip/{hla}.log",
    wrapper:
        "v1.23.5-27-g1638662a/bio/bgzip"

rule tabix:
    input:
        "results/orthanq-candidates/{hla}.vcf.gz",
    output:
        "results/orthanq-candidates/{hla}.vcf.gz.tbi",
    log:
        "logs/tabix/{hla}.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        "v1.23.5-27-g1638662a/bio/tabix/index"

# intersect the candidate variants to contain only variants that are in coding exons.
rule intersect_exons:
    input:
        left="results/orthanq-candidates/{hla}.vcf",
        right="resources/HLA_exons/HLA_{hla}_exons.bed"
    output:
        "results/orthanq-candidates-intersected/{hla}.vcf"
    params:
        extra="-header"
    log:
        "logs/bedtools/{hla}.log"
    wrapper:
        "v1.31.1/bio/bedtools/intersect"

rule bgzip_exons:
    input:
        "results/orthanq-candidates-intersected/{hla}.vcf"
    output:
        "results/orthanq-candidates-intersected/{hla}.vcf.gz"
    params:
        extra="", # optional
    threads: 1
    log:
        "logs/bgzip-merged-exons/{hla}.log",
    wrapper:
        "v1.23.5-27-g1638662a/bio/bgzip"

rule tabix_exons:
    input:
        "results/orthanq-candidates-intersected/{hla}.vcf.gz"
    output:
        "results/orthanq-candidates-intersected/{hla}.vcf.gz.tbi"
    log:
        "logs/tabix-merged-exons/{hla}.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        "v1.23.5-27-g1638662a/bio/tabix/index"
