rule varlociraptor_preprocess:
    input:
        ref = "resources/reference/NC_045512.2.fasta",
        fai = "resources/reference/NC_045512.2.fasta.fai",
        candidates = "results/orthanq-candidates/candidates.vcf",
        bam="results/vg/alignment/{sample}_vg.sorted.bam",
        bai="results/vg/alignment/{sample}_vg.sorted.bam.bai"
    output:
        "results/observations/{sample}.bcf"
    log:
        "logs/varlociraptor/preprocess/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants "
        "--report-fragment-ids --omit-mapq-adjustment --atomic-candidate-variants --candidates {input.candidates} "
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
        "varlociraptor call variants --omit-strand-bias --omit-read-position-bias --omit-read-orientation-bias --omit-softclip-bias --omit-homopolymer-artifact-detection --omit-alt-locus-bias generic --obs sample={input.obs} " ##varlociraptor v5.3.0
        "--scenario {input.scenario} > {output} 2> {log}"

rule orthanq_call:
    input:
        calls = "results/calls/{sample}.bcf",
        candidate_variants = "results/orthanq-candidates/candidates.vcf",
    output:
        "results/orthanq/{sample}/{sample}.tsv"
    conda:
        "../envs/orthanq.yaml"
    log:
        "logs/orthanq-call/{sample}.log"
    params:
        prior = config["orthanq_prior"]
    shell:
        "../../orthanq/target/debug/orthanq call virus --haplotype-calls {input.calls} --haplotype-variants {input.candidate_variants} "
        "--prior {params} --output {output} 2> {log}"