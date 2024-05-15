# This smk calls variants by varlociraptor and quantifies haplotypes by orthanq.

rule varlociraptor_preprocess:
    input:
        ref="results/refs/hs_genome.fasta",
        fai = "results/refs/hs_genome.fasta.fai",
        candidates = "results/orthanq-candidates-intersected/{hla}.vcf",
        bam="results/vg/alignment/{sample}_vg.sorted.reheadered.extracted.bam",
        bai="results/vg/alignment/{sample}_vg.sorted.reheadered.extracted.bam.bai"
    output:
        "results/observations/{sample}_{hla}.bcf"
    log:
        "logs/varlociraptor/preprocess/{sample}_{hla}.log",
    benchmark:    
        "benchmarks/varlociraptor_preprocess/{sample}_{hla}.tsv"  
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants " #for testing the new mnv improvement in varlociraptor
        "--report-fragment-ids --omit-mapq-adjustment --atomic-candidate-variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"

rule varlociraptor_call:
    input:
        obs="results/observations/{sample}_{hla}.bcf",
        scenario="resources/scenarios/scenario.yaml",
    output:
        "results/calls/{sample}_{hla}.bcf"
    log:
        "logs/varlociraptor/call/{sample}_{hla}.log",
    benchmark:    
        "benchmarks/varlociraptor_call/{sample}_{hla}.tsv"  
    conda:
        "../envs/varlociraptor.yaml" 
    shell:
        "varlociraptor call variants --omit-strand-bias --omit-read-position-bias --omit-read-orientation-bias --omit-softclip-bias --omit-homopolymer-artifact-detection --omit-alt-locus-bias generic --obs sample={input.obs} " ##varlociraptor v5.3.0
        "--scenario {input.scenario} > {output} 2> {log}"

rule orthanq_call:
    input:
        calls = "results/calls/{sample}_{hla}.bcf",
        candidate_variants = "results/orthanq-candidates-intersected/{hla}.vcf",
        xml = config["allele_db_xml"]
    output:
        table="results/orthanq/{sample}_{hla}/{sample}_{hla}.csv",
        three_field_solutions="results/orthanq/{sample}_{hla}/3_field_solutions.json",
        two_field_solutions="results/orthanq/{sample}_{hla}/2_field_solutions.json",
        final_solutions="results/orthanq/{sample}_{hla}/final_solution.json",
        lp_solution="results/orthanq/{sample}_{hla}/lp_solution.json",
        two_field_table="results/orthanq/{sample}_{hla}/2-field.csv",
        g_groups="results/orthanq/{sample}_{hla}/G_groups.csv"
    conda:
        "../envs/orthanq.yaml"
    log:
        "logs/orthanq-call/{sample}_{hla}.log"
    benchmark:    
        "benchmarks/orthanq_call/{sample}_{hla}.tsv"
    params:
        prior = config["orthanq_prior"]
    shell:
        "../orthanq/target/release/orthanq call hla --haplotype-calls {input.calls} --haplotype-variants "
        " {input.candidate_variants} --xml {input.xml} --prior {params} --output {output.table} 2> {log}"
