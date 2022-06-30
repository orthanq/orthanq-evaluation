rule varlociraptor_preprocess:
    input:
        ref="results/refs/hs_genome.fasta",
        fai = "results/refs/hs_genome.fasta.fai",
        candidates="resources/hla-allele-variants_v4.vcf.gz",
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

rule orthanq_candidates:
    input:
        alleles = "resources/HLA-alleles/hla_all/hla_gen.fasta",
        genome = "results/refs/hs_genome.fasta",
    output:
        candidate_variants = "results/orthanq-candidates/{hla}.vcf" #to be changed later
    log:
        "logs/orthanq-candidates/{hla}.log"
    shell:
        "~/orthanq/target/release/orthanq candidates --alleles {input.alleles} "
        "--genome {input.genome} --wes --output results/orthanq-candidates" # --wes option for protein level hla type variant generation, --wgs for individual types 

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
        "--haplotype-variants {input.candidate_variants} "
        "--max-haplotypes 5 --output {output} 2> {log}"
