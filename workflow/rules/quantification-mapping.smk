# ruleorder: orthanq_candidates > varlociraptor_preprocess

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

rule vg_autoindex:
    input:
        genome="results/refs/hs_genome.fasta",
        variants= "results/orthanq-candidates/merged.vcf"
    output:
        "results/vg/autoindex/idx.giraffe.gbz"
    log:
        "logs/vg/autoindex.log"
    conda:
        "../envs/vg.yaml"
    threads: 20
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
    threads: 18
    shell:
        "vg giraffe -Z results/vg/autoindex/idx.giraffe.gbz -f {input.reads[0]} -f {input.reads[1]} --output-format BAM -t {threads} > {output} 2> {log}"

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

rule kallisto_index:
    input:
        fasta = "resources/IMGTHLA-3.32.0-alpha/fasta/{hla}_gen.fasta"
    output:
        index = "results/kallisto-index/{hla}.idx"
    params:
        extra = "--kmer-size=31"
    log:
        "logs/kallisto/index/{hla}.log"
    threads: 10
    wrapper:
        "v1.25.0/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq = get_fastq_input,
        index = "results/kallisto-index/{hla}.idx"
    output:
        directory('results/kallisto/quant_results_{sample}_{hla}')
    params:
        extra = "-b 5 --seed=42 --pseudobam"
    log:
        "logs/kallisto/kallisto_quant_{sample}_{hla}.log"
    threads: 20
    wrapper:
        "v1.25.0/bio/kallisto/quant"