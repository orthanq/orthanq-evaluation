ruleorder: orthanq_candidates > varlociraptor_preprocess

rule bwakit_index:
    output:
        "results/bwakit-genome/hs38DH.fa",
        "results/bwakit-genome/hs38DH.fa.alt"
    log:
        "logs/bwakit/bwakit-index.log"
    conda:
        "../envs/bwakit.yaml"
    threads: 20
    shell:
        "mkdir -p results/bwakit-genome && cd results/bwakit-genome && run-gen-ref hs38DH && bwa index hs38DH.fa"

rule bwakit_alignment:
    input:
        reads = get_fastq_input,
        genome = "results/bwakit-genome/hs38DH.fa"
    output:
        "results/mapped/{sample}.aln.bam"
    log:
        "logs/bwakit/bwakit-align/{sample}.log"
    conda:
        "../envs/bwakit.yaml"
    threads: 20
    shell:
        "run-bwamem -o results/mapped/{wildcards.sample} -s -t {threads} -H {input.genome} {input.reads} |sh 2> {log}"

rule samtools_view:
    input:
        "results/mapped/{sample}.aln.bam"
    output:
        "results/mapped/{sample}.aln.sam"
    log:
        "logs/samtools_view/{sample}.log"
    threads: 6
    params:
        extra="-h",  #include header
    wrapper:
        "v1.22.0/bio/samtools/view"

rule bwakit_post_process:
    input:
        sample="results/mapped/{sample}.aln.sam",
        genome_alt="results/bwakit-genome/hs38DH.fa.alt",
        bwakit="resources/bwa/bwakit"
    output:
        "results/processed_mapped/{sample}.aln.sam"
    log:
        "logs/bwakit/bwakit-postalt/{sample}.log"
    conda:
        "../envs/bwakit.yaml"
    threads:10
    shell:
        #should be fixed.
        "k8 {input.bwakit}/bwa-postalt.js {input.genome_alt} {input.sample} > {output} 2> {log}"

rule realignment:
    input:
        sam="results/processed_mapped/{sample}.aln.sam",
        genome_alt="results/bwakit-genome/hs38DH.fa",
        realigner="resources/realignment-after-postalt/target/debug/realignment"
    output:
        "results/realignment/{sample}.realigned.bam"
    log:
        "logs/realignment/{sample}.log"
    shell:
        "{input.realigner} --alignment {input.sam} --reference {input.genome_alt} --output {output} 2> {log}"

rule genome_index:
    input:
        "results/bwakit-genome/hs38DH.fa"
    output:
        "results/bwakit-genome/hs38DH.fa.fai"
    log:
        "logs/genome_index/bwa_hs_genome.log"
    threads: 4
    wrapper:
        "v1.22.0/bio/samtools/faidx" #update wrapper version

rule samtools_sort:
    input:
        "results/realignment/{sample}.realigned.bam"
    output:
        "results/realignment/{sample}.realigned.sorted.bam"
    log:
        "logs/samtools_sort/{sample}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v1.22.0/bio/samtools/sort"

rule samtools_index:
    input:
        "results/realignment/{sample}.realigned.sorted.bam"
    output:
        "results/realignment/{sample}.realigned.sorted.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    threads: 4
    wrapper:
        "v1.22.0/bio/samtools/index"

rule prepare_orthanq_genome:
    input:
        genome="results/bwakit-genome/hs38DH.fa",
        chromosomes="resources/genome-prep/chromosomes_to_extract.txt"
    output:
        "results/orthanq-genome/hg38_genome_extracted.fa"
    log:
        "logs/seqtk/extract-chromosomes.log"
    threads: 2
    shell:
        "seqtk subseq {input.genome} {input.chromosomes} > {output} 2> {log}"

##orthanq_candidates generates both candidate variants and individual genome sequences for each locus(to be used for quantification indices in the downstream processing)
rule orthanq_candidates:
    input:
        alleles = config["allele_db"],
        genome = "results/orthanq-genome/hg38_genome_extracted.fa",
        xml = config["allele_db_xml"],
        allele_freqs = "resources/allele_freqs/allele_frequencies.csv"
    output:
        candidate_variants = expand("results/orthanq-candidates/{hla}.vcf",hla=loci),
        #fasta = expand("results/orthanq-candidates/{hla}.fasta",hla=loci)
    log:
        "logs/orthanq-candidates/candidates.log"
    shell:
        "../orthanq/target/release/orthanq candidates --alleles {input.alleles} --genome {input.genome} "
        "--xml {input.xml} --allele-freq {input.allele_freqs} --wes --output results/orthanq-candidates 2> {log}" # --wes option for protein level hla type variant generation, --wgs for individual types 
