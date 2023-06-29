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

# rule prepare_orthanq_genome:
#     input:
#         genome="results/bwakit-genome/hs38DH.fa",
#         chromosomes="resources/genome-prep/chromosomes_to_extract.txt"
#     output:
#         "results/orthanq-genome/hg38_genome_extracted.fa"
#     log:
#         "logs/seqtk/extract-chromosomes.log"
#     threads: 2
#     shell:
#         "seqtk subseq {input.genome} {input.chromosomes} > {output} 2> {log}"

##orthanq_candidates generates both candidate variants and individual genome sequences for each locus(to be used for quantification indices in the downstream processing)
#the alleles are filtered out according to the following criteria:
#if the allele has a full name equivalent (e.g. A*24:02:02:02) in the database and it's higher than 0.05, it's IN.
#if the allele doesn't have a full name equivalent (e.g. A*24:02:02), then look for the first three fields, if still not then look for first two fields (A*24:02)
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
#
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
#merge all candidate variants in the orthanq-candidates (not orthanq-candidates-intersected)
rule bcftools_merge:
    input:
        calls=expand("results/orthanq-candidates/{hla}.vcf.gz", hla=loci)
    output:
       "results/orthanq-candidates/merged.vcf",
    log:
        "logs/bcftools-merge/merge.log"
    params:
        uncompressed_bcf=False,
        extra="",  # optional parameters for bcftools concat (except -o)
    wrapper:
        "v1.23.5/bio/bcftools/merge"

#then, do the intersection for all exons this time
rule intersect_exons_all:
    input:
        left="results/orthanq-candidates/merged.vcf",
        right="resources/HLA_exons/all_exons.bed"
    output:
        "results/orthanq-candidates-intersected-all/merged_intersected_variants.vcf"
    params:
        extra="-header"
    log:
        "logs/bedtools/intersect-all-exons.log"
    wrapper:
        "v1.31.1/bio/bedtools/intersect"

rule bgzip_exons_all:
    input:
        "results/orthanq-candidates-intersected-all/merged_intersected_variants.vcf"
    output:
        "results/orthanq-candidates-intersected-all/merged_intersected_variants.vcf.gz"
    params:
        extra="", # optional
    threads: 1
    log:
        "logs/bgzip-merged-exons/all.log",
    wrapper:
        "v1.31.1/bio/bgzip"


rule tabix_exons_all:
    input:
        "results/orthanq-candidates-intersected-all/merged_intersected_variants.vcf.gz"
    output:
        "results/orthanq-candidates-intersected-all/merged_intersected_variants.vcf.gz.tbi"
    log:
        "logs/tabix-merged-exons/all.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        "v1.31.1/bio/tabix/index"

rule vg_autoindex:
    input:
        genome="results/refs/hs_genome.fasta",
        variants="results/orthanq-candidates-intersected-all/merged_intersected_variants.vcf.gz"
    output:
        "results/vg/autoindex/idx.giraffe.gbz"
    log:
        "logs/vg/autoindex.log"
    conda:
        "../envs/vg.yaml"
    threads: 40
    shell:
        "vg autoindex --workflow giraffe -r {input.genome} -v {input.variants} -p results/vg/autoindex/idx -t {threads}"

#new mapping strategy:
#first, align against with bwa against the primary assembly
#second, extract all reads that belong to HLA loci, use samtools view and then fastq
#last, map only to those extracted reads

#Step 1: map reads by bwa
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
    params:
        index="results/bwa-index/hs_genome",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",             
        sort_order="coordinate", 
    threads: 30
    wrapper:
        "v2.0.0/bio/bwa/mem"

rule samtools_index_bwa:
    input:
        "results/bwa_alignment/{sample}_mapped.bam"
    output:
        "results/bwa_alignment/{sample}_mapped.bai"
    log:
        "logs/samtools_index_bwa/{sample}.log"
    threads: 4
    wrapper:
        "v2.0.0/bio/samtools/index"

#Step 2: extract reads that map to HLA loci
rule samtools_view:
    input:
        "results/bwa_alignment/{sample}_mapped.bam",
        regions="resources/HLA_regions/regions.bed" #only contains the regions for HLA-A,HLA-B, HLA-C and HLA-DQB1
    output:
        bam="results/bwa_alignment/{sample}_mapped_extracted.bam",
        idx="results/bwa_alignment/{sample}_mapped_extracted.bai"
    log:
        "logs/samtools-view/{sample}.log",
    params:
        extra=lambda wc, input: "-L {}".format(input.regions)
    threads: 10
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
    params:
        sort="-m 4G",
        fastq="-n",
    threads: 10
    wrapper:
        "v2.0.0/bio/samtools/fastq/separate"

rule vg_giraffe:
    input:
        reads_1 = "results/extracted_reads/{sample}.1.fq",
        reads_2 = "results/extracted_reads/{sample}.2.fq",
        idx = "resources/vg-pangenome/hprc-v1.0-mc-grch38.xg"
    output:
        "results/vg/alignment/{sample}_vg.bam"
    log:
        "logs/vg/alignment/{sample}.log"
    conda:
        "../envs/vg.yaml"
    threads: 30
    shell:
        "vg giraffe -x {input.idx} -f {input.reads_1} -f {input.reads_2} --max-multimaps 3 --output-format BAM -t {threads}  > {output} 2> {log}"

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

rule reheader:
    input:
        "results/vg/alignment/{sample}_vg.sorted.bam"
    output:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam"
    log:
        "logs/samtools_reheader/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools view -H {input} | sed 's/GRCh38.chr//g' | samtools reheader - {input} > {output} 2> {log}"
    
rule samtools_index_after_reheader:
    input:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam"
    output:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam.bai"
    log:
        "logs/samtools_index_after_reheader/{sample}.log"
    threads: 4
    wrapper:
        "v1.22.0/bio/samtools/index"


rule keep_only_primary_chr:
    input:
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam",
        "results/vg/alignment/{sample}_vg.sorted.reheadered.bam.bai",
    output:
        bam="results/vg/alignment/{sample}_vg.sorted.reheadered.extracted.bam",
        idx="results/vg/alignment/{sample}_vg.sorted.reheadered.extracted.bam.bai"
    log:
        "logs/samtools-view-primary-chr/{sample}.log",
    params:
        region="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
    threads: 10
    wrapper:
        "v2.0.0/bio/samtools/view"
