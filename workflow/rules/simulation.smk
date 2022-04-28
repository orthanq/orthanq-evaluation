rule art_simulation:
    input:
        "resources/HLA-alleles/{allele}_nuc.fasta"
    output:
        "results/art/{allele}_1.fq",
        "results/art/{allele}_2.fq"
    log:
        "logs/art/{allele}.log",
    #threads: config["threads"]["art"]
    conda:
        "../envs/art.yaml"
    params: config["f_coverage"]
    shell:
        "art_illumina -ss HS25 -i {input} -p -l 150 -s 10 -m 200 -f {params} --noALN --rndSeed 31303889 -o"
        " results/art/{wildcards.allele}_ 2> {log}"

rule get_fractions:
    input:
        fq1="results/art/{allele}_1.fq",
        fq2="results/art/{allele}_2.fq" 
    output:
        out_fq1="results/fractions/{sample}-{allele}-{num}_1.fq",
        out_fq2="results/fractions/{sample}-{allele}-{num}_2.fq"
    log:
        "logs/seqtk/{sample}-{allele}-{num}.log",
    params:
        "{num}"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample -s100 {input.fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.fq2} {params} > {output.out_fq2} 2> {log}"

rule simulate_nonhla:
    input:
        "resources/chr21.fa"
    output:
        "results/art_nonhla/chr21_1.fq",
        "results/art_nonhla/chr21_2.fq"
    log:
        "logs/art_nonhla/nonhla.log",
    conda:
        "../envs/art.yaml"
    params: config["f_coverage"]
    shell:
        "art_illumina -ss HS25 -i {input} -p -l 150 -s 10 -m 200 -f {params} --noALN --rndSeed 31303889 -o"
        " results/art_nonhla/chr21_ 2> {log}"

rule get_fractions_nonhla:
    input:
        fq1="results/art_nonhla/chr21_1.fq",
        fq2="results/art_nonhla/chr21_2.fq" 
    output:
        out_fq1="results/fraction_nonhla/chr21_1.fq",
        out_fq2="results/fraction_nonhla/chr21_2.fq"
    log:
        "logs/nonhla_seqtk/nonhla.log",
    params:
        config["n_nonhla_reads"]
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample -s100 {input.fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.fq2} {params} > {output.out_fq2} 2> {log}"
  
rule concat_fractions: 
    input:
        fq1=lambda wc: expand("results/fractions/{{sample}}-{allele}-{num}_1.fq",
            zip,
            allele=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['hla'],
            num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
            ),
        fq2=lambda wc: expand("results/fractions/{{sample}}-{allele}-{num}_2.fq",
            zip,
            allele=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['hla'],
            num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
            ),
        nonhla_fq1 = "results/fraction_nonhla/chr21_1.fq",
        nonhla_fq2 = "results/fraction_nonhla/chr21_2.fq"
    output:
        out_fq1="results/mixed/{sample}_1.fq",
        out_fq2="results/mixed/{sample}_2.fq"
    log:
        "logs/mixed/{sample}.log",
    shell:
        "cat {input.fq1} {input.nonhla_fq1} > {output.out_fq1}; "
        "cat {input.fq2} {input.nonhla_fq2} > {output.out_fq2}"


