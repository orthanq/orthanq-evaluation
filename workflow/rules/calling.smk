ruleorder: orthanq_candidates > varlociraptor_preprocess

rule orthanq_candidates:
    input:
        alleles = config["allele_db"],
        genome = "results/refs/hs_genome.fasta",
    output:
        candidate_variants = expand("results/orthanq-candidates/{hla}.vcf",hla=loci) #to be changed later
    log:
        "logs/orthanq-candidates/candidates.log"
    shell:
        "~/orthanq/target/release/orthanq candidates --alleles {input.alleles} "
        "--genome {input.genome} --wes --output results/orthanq-candidates 2> {log}" # --wes option for protein level hla type variant generation, --wgs for individual types 

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
        "varlociraptor/target/release/varlociraptor preprocess variants --report-fragment-ids --candidates {input.candidates} "
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
        "varlociraptor/target/release/varlociraptor call variants --omit-alt-locus-bias generic --obs sample={input.obs} " ##varlociraptor v5.3.0
        "--scenario {input.scenario} > {output} 2> {log}"

rule orthanq_call:
    input:
        calls = "results/calls/{sample}_{hla}.bcf",
        candidate_variants = "results/orthanq-candidates/{hla}.vcf",
    output:
        report(
            "results/orthanq/{sample}_{hla}.tsv",
            caption="../report/haplotype_abundances.rst",
        )
    log:
        "logs/orthanq-call/{sample}_{hla}.log"
    shell:
        "~/orthanq/target/release/orthanq call --haplotype-calls {input.calls} "
        "--haplotype-variants {input.candidate_variants} --max-haplotypes 5 --output {output} 2> {log}"

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


