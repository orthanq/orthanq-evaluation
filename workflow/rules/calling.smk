rule varlociraptor_preprocess:
    input:
        ref="results/refs/hs_genome.fasta",
        fai = "results/refs/hs_genome.fasta.fai",
        candidates = "results/orthanq-candidates/{hla}.vcf",
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai"
    output:
        "results/observations/{sample}_{hla}.bcf"
    log:
        "logs/varlociraptor/preprocess/{sample}_{hla}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --report-fragment-ids --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"

rule varlociraptor_call:
    input:
        obs="results/observations/{sample}_{hla}.bcf",
        scenario="resources/scenarios/scenario.yaml",
    output:
        "results/calls/{sample}_{hla}.bcf"
    log:
        "logs/varlociraptor/call/{sample}_{hla}.log",
    conda:
        "../envs/varlociraptor.yaml" 
    shell:
        "varlociraptor call variants --omit-strand-bias --omit-read-position-bias --omit-read-orientation-bias --omit-softclip-bias --omit-homopolymer-artifact-detection --omit-alt-locus-bias generic --obs sample={input.obs} " ##varlociraptor v5.3.0
        "--scenario {input.scenario} > {output} 2> {log}"

rule orthanq_call:
    input:
        calls = "results/calls/{sample}_{hla}.bcf",
        counts = "results/kallisto/quant_results_{sample}_{hla}",
        candidate_variants = "results/orthanq-candidates/{hla}.vcf",
    output:
        report(
            "results/orthanq/{sample}_{hla}/{sample}_{hla}.tsv",
            caption="../report/haplotype_abundances.rst",
        ) #this directory also contains json files for solutions
    log:
        "logs/orthanq-call/{sample}_{hla}.log"
    shell:
        "~/orthanq/target/release/orthanq call --haplotype-counts {input.counts}/abundance.h5 --haplotype-calls {input.calls} "
        "--haplotype-variants {input.candidate_variants} --min-norm-counts 0.0 --max-haplotypes 5 --output {output} 2> {log}"

rule arcasHLA_reference:
    output:
        "foo.txt"
    log:
        "logs/arcashla/reference.log"
    conda:
        "../envs/arcasHLA.yaml"
    params:
        version = config["arcasHLA_version"]
    shell:
        "touch foo.txt; "
        "arcasHLA reference --version {params}"

rule arcasHLA_extract:
    input:
        "foo.txt",
        bam="results/mapped/{sample}.bam",
    output:
        extracted_read1="results/arcasHLA/{sample}/{sample}.extracted.1.fq.gz",
        extracted_read2="results/arcasHLA/{sample}/{sample}.extracted.2.fq.gz",
    log:
        "logs/arcashla/extract/{sample}.log"
    conda:
        "../envs/arcasHLA.yaml"
    threads: 8
    shell:
        "arcasHLA extract {input.bam} -o results/arcasHLA/{wildcards.sample} -t {threads} -v 2> {log}"

rule arcasHLA_genotype:
    input:
        extracted_read1="results/arcasHLA/{sample}/{sample}.extracted.1.fq.gz",
        extracted_read2="results/arcasHLA/{sample}/{sample}.extracted.2.fq.gz",
    output:
        "results/arcasHLA/{sample}_{hla}/{sample}.genotype.json"
    log:
        "logs/arcashla/genotype/{sample}_{hla}.log",
    conda:
        "../envs/arcasHLA.yaml"
    threads: 4
    params:
        locus = "{hla}"
    shell:
        "arcasHLA genotype {input} -g {params} -o results/arcasHLA/{wildcards.sample}_{wildcards.hla} -t {threads} 2> {log}"


