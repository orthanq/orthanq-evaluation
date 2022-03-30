rule varlociraptor_preprocess:
    input:
        ref="resources/hs_genome/hs_genome.fasta",
        candidates="resources/hla-allele-variants_v3.vcf.gz",
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai"
    output:
        "results/observations/{sample}.bcf"
    log:
        "logs/varlociraptor/preprocess/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
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
        "varlociraptor call variants generic --obs sample={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"
