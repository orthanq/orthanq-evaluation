ruleorder: orthanq_candidates > varlociraptor_preprocess

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
        "v1.3.2/bio/reference/ensembl-sequence"

rule genome_index:
    input:
        "results/refs/hs_genome.fasta"
    output:
        "results/refs/hs_genome.fasta.fai"
    log:
        "logs/genome_index/hs_genome.log"
    threads: 4
    wrapper:
        "v1.3.2/bio/samtools/faidx" #update wrapper version

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
    log:
        "logs/orthanq-candidates/candidates.log"
    shell:
        "~/orthanq/target/release/orthanq candidates --alleles {input.alleles} --genome {input.genome} "
        "--xml {input.xml} --allele-freq {input.allele_freqs} --wes --output results/orthanq-candidates 2> {log}" # --wes option for protein level hla type variant generation, --wgs for individual types 

# rule kallisto_index:
#     input:
#         fasta = "results/orthanq-candidates/{hla}.fasta"
#     output:
#         index = "results/kallisto-index/{hla}.idx"
#     params:
#         extra = "--kmer-size=31"
#     log:
#         "logs/kallisto/index/{hla}.log"
#     threads: 2
#     wrapper:
#         "v1.3.2/bio/kallisto/index"

# rule kallisto_quant:
#     input:
#         fastq = get_fastq_input,
#         index = "results/kallisto-index/{hla}.idx"
#     output:
#         directory('results/kallisto/quant_results_{sample}_{hla}')
#     params:
#         extra = "-b 20 --seed=42"
#     log:
#         "logs/kallisto/kallisto_quant_{sample}_{hla}.log"
#     threads: 20
#     wrapper:
#         "v1.3.2/bio/kallisto/quant"

# rule salmon_index:
#     input:
#         sequences="results/orthanq-candidates/{hla}.fasta",
#     output:
#         directory("results/salmon-index/{hla}")
#     log:
#         "logs/salmon-index/{hla}.log",
#     threads: 10
#     params:
#         # optional parameters
#         extra="",
#     wrapper:
#         "v1.19.1-1-g03532856/bio/salmon/index"

# rule salmon_quant_reads:
#     input:
#         index="results/salmon-index/{hla}",
#         r1=get_fastq_input_1,
#         r2=get_fastq_input_2,
#     output:
#         quant="results/salmon-quant/{sample}_{hla}/quant.sf",
#         lib="results/salmon-quant/{sample}_{hla}/lib_format_counts.json",
#     log:
#         "logs/salmon-quant/{sample}_{hla}.log",
#     params:
#         # optional parameters
#         libtype="A",
#         numBootstraps=20
#     threads: 10
#     wrapper:
#         "v1.19.1-1-g03532856/bio/salmon/quant"

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
        "v1.4.0/bio/bwa/index" #update wrapper version

rule bwa_mem:
    input:
        reads = get_fastq_input,
        #reads = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"],
        idx = multiext("results/bwa-index/hs_genome", ".amb", ".ann", ".bwt", ".pac", ".sa")
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
        "v1.4.0/bio/bwa/mem" #update wrapper version

rule samtools_index:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    threads: 4
    wrapper:
        "v1.4.0/bio/samtools/index"
