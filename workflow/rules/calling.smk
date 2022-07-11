ruleorder: orthanq_candidates > varlociraptor_preprocess

rule orthanq_candidates:
    input:
        alleles = "/vol/compute/hamdiyes_project/orthanq-evaluation/resources/HLA-alleles/hla_all/hla_gen.fasta",
        genome = "results/refs/hs_genome.fasta",
    output:
        candidate_variants = "results/orthanq-candidates/{hla}.vcf" #to be changed later
    log:
        "logs/orthanq-candidates/{hla}.log"
    shell:
        "~/orthanq/target/release/orthanq candidates --alleles {input.alleles} "
        "--genome {input.genome} --wes --output results/orthanq-candidates" # --wes option for protein level hla type variant generation, --wgs for individual types 

rule varlociraptor_preprocess:
    input:
        ref="results/refs/hs_genome.fasta",
        fai = "results/refs/hs_genome.fasta.fai",
        candidates = expand("results/orthanq-candidates/{hla}.vcf", hla=loci),
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai"
    output:
        "results/observations/{sample}.bcf"
    log:
        "logs/varlociraptor/preprocess/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor/target/release/varlociraptor preprocess variants --report-fragment-ids --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"

rule varlociraptor_call:
    input:
        obs="results/observations/{sample}.bcf",
        scenario="resources/scenarios/scenario.yaml",
    output:
        "results/calls/{sample}.bcf"
    log:
        "logs/varlociraptor/call/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor/target/release/varlociraptor call variants generic --obs sample={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"

rule orthanq_call:
    input:
        calls = "results/calls/{sample}.bcf",
        obs = "results/observations/{sample}.bcf",
        candidate_variants = "results/orthanq-candidates/{hla}.vcf",
    output:
        report(
            "results/orthanq/{sample}_{hla}.tsv",
            caption="../report/haplotype_abundances.rst",
        )
    log:
        "logs/orthanq-call/{sample}_{hla}.log"
    shell:
        "~/orthanq/target/release/orthanq call --haplotype-calls {input.calls} --observations {input.obs} "
        "--haplotype-variants {input.candidate_variants} --max-haplotypes 5 --output {output} 2> {log}"

rule arcasHLA_extract:
    input:
        bam="results/mapped/{sample}.bam",
    output:
        extracted_read1="results/arcasHLA/{sample}/{sample}.extracted.1.fq.gz",
        extracted_read2="results/arcasHLA/{sample}/{sample}.extracted.2.fq.gz",
    log:
        "logs/arcashla/extract/{sample}.log"
    conda:
        "../envs/arcasHLA.yaml"
    threads: 8
    params:
        version = "3.24.0"
    shell:
        "arcasHLA reference --version {params}; "
        "arcasHLA extract {input} -o results/arcasHLA/{wildcards.sample} -t {threads} -v 2> {log}"

rule arcasHLA_genotype:
    input:
        extracted_read1="results/arcasHLA/{sample}/{sample}.extracted.1.fq.gz",
        extracted_read2="results/arcasHLA/{sample}/{sample}.extracted.2.fq.gz",
    output:
        "results/arcasHLA/{sample}/{sample}.genotype.json"
    log:
        "logs/arcashla/genotype/{sample}.log",
    conda:
        "../envs/arcasHLA.yaml"
    threads: 4
    params:
        locus = "DQA1"
    shell:
        "arcasHLA genotype {input} -g {params} -o results/arcasHLA/{wildcards.sample} -t {threads} 2> {log}"


