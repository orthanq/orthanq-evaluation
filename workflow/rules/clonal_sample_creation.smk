rule concat_hapmap: 
    input:
        fq1= samples["fq1"],
        fq2= samples["fq2"]
    output:
        out_fq1="results/subclonal/subclonal_1.fq",
        out_fq2="results/subclonal/subclonal_2.fq"
    log:
        "logs/clonal_samples/subclonal_sample.log",
    shell:
        "cat {input.fq1} > {output.out_fq1}; "
        "cat {input.fq2} > {output.out_fq2}"

