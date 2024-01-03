rule art_simulation:
    input:
        "resources/lineages/{lineage}.fasta"
    output:
        "results/art/{lineage}_1.fq",
        "results/art/{lineage}_2.fq"
    log:
        "logs/art/{lineage}.log",
    #threads: config["threads"]["art"]
    conda:
        "../envs/art.yaml"
    params: config["f_coverage"]
    shell:
        "art_illumina -ss HS25 -i {input} -p -l 150 -s 10 -m 200 -f {params} --noALN --rndSeed 31303889 -o"
        " results/art/{wildcards.lineage}_ 2> {log}"

rule get_fractions:
    input:
        fq1="results/art/{lineage}_1.fq",
        fq2="results/art/{lineage}_2.fq" 
    output:
        out_fq1="results/fractions/{sample}-{lineage}-{num}_1.fq",
        out_fq2="results/fractions/{sample}-{lineage}-{num}_2.fq"
    log:
        "logs/seqtk/{sample}-{lineage}-{num}.log",
    params:
        "{num}"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample -s100 {input.fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.fq2} {params} > {output.out_fq2} 2> {log}"
 
rule concat_fractions: 
    input:
        fq1=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_1.fq",
            zip,
            lineage=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['lineage'],
            num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
            ),
        fq2=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_2.fq",
            zip,
            lineage=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['lineage'],
            num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
            ),
    output:
        out_fq1="results/mixed/{sample}_1.fq",
        out_fq2="results/mixed/{sample}_2.fq"
    log:
        "logs/mixed/{sample}.log",
    shell:
        "cat {input.fq1} > {output.out_fq1}; "
        "cat {input.fq2} > {output.out_fq2}"


