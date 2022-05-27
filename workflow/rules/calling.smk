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

rule orthanq_call:
    input:
        calls = "results/calls/{sample}.bcf",
        obs = "results/observations/{sample}.bcf",
        candidate_variants = "resources/hla-allele-variants_v4.vcf.gz",
        counts = "results/kallisto/quant_results_{sample}_{hla}"
    output:
        report(
            "results/orthanq/{sample}_{hla}.tsv",
            caption="../report/haplotype_abundances.rst",
        )
    log:
        "logs/orthanq/{sample}_{hla}.log"
    shell:
        "~/orthanq/target/release/orthanq call --haplotype-calls {input.calls} --observations {input.obs} "
        "--haplotype-variants {input.candidate_variants} --haplotype-counts {input.counts}/abundance.h5 "
        "--min-norm-counts 0.01 --max-haplotypes 2 --use-evidence both --output {output} 2> {log}" #--use-evidence, for easier debugging (available options: varlociraptor, kallisto or both.)
