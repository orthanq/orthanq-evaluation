# This script evaluates orthanq together with two other HLA typers: arcasHLA and HLA-LA.
# Therefore, it outputs individual tables for three tools as well as a performance table for comparison.

rule arcasHLA_reference:
    output:
        "foo.txt"
    log:
        "logs/arcashla/reference.log"
    benchmark:    
        "benchmarks/arcasHLA/reference/reference.tsv"  
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
        bam="results/bwa_alignment/{sample}_mapped.bam" #note: arcasHLA reads sample name as 'SRR01234_mapped', as a result it creates extracted fastqs with this name, e.g. 'SRR01234_mapped.extracted.1.fq.gz'
    output:
        extracted_read1="results/arcasHLA/{sample}/{sample}_mapped.extracted.1.fq.gz",
        extracted_read2="results/arcasHLA/{sample}/{sample}_mapped.extracted.2.fq.gz",
    log:
        "logs/arcashla/extract/{sample}.log"
    benchmark:    
        "benchmarks/arcasHLA/extract/{sample}.tsv"  
    conda:
        "../envs/arcasHLA.yaml"
    threads: 20
    shell:
        "arcasHLA extract {input.bam} -o results/arcasHLA/{wildcards.sample} -t {threads} -v 2> {log}"

rule arcasHLA_genotype:
    input:
        extracted_read1="results/arcasHLA/{sample}/{sample}_mapped.extracted.1.fq.gz",
        extracted_read2="results/arcasHLA/{sample}/{sample}_mapped.extracted.2.fq.gz",
    output:
        "results/arcasHLA/{sample}_{hla}/{sample}_mapped.genotype.json"
    log:
        "logs/arcashla/genotype/{sample}_{hla}.log",
    benchmark:    
        "benchmarks/arcasHLA/genotype/{sample}_{hla}.tsv"  
    conda:
        "../envs/arcasHLA.yaml"
    threads: 20
    params:
        locus = "{hla}"
    shell:
        "arcasHLA genotype {input} -g {params} -o results/arcasHLA/{wildcards.sample}_{wildcards.hla} -t {threads} 2> {log}"

rule HLA_LA:
    input:
        bam="results/bwa_alignment/{sample}_mapped.bam",
        bai="results/bwa_alignment/{sample}_mapped.bai",
        index="resources/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH", #V3.32
    output:
        "results/HLA-LA/{sample}/hla/R1_bestguess_G.txt",
    threads: 40
    log:
        "logs/HLA-LA/{sample}.log",
    benchmark:    
        "benchmarks/hla-la/{sample}.tsv"  
    params:
        graph=lambda w, input: os.path.basename(os.path.dirname(input.index)),
        graphdir=lambda w, input: os.path.dirname(os.path.dirname(input.index)),
        workdir=lambda w, output: os.path.dirname(
            os.path.dirname(os.path.dirname(output[0]))
        )
    conda:
        "../envs/hla-la.yaml"
    shell:
        "HLA-LA.pl --bam {input.bam} --sampleID {wildcards.sample} --graph {params.graph} --customGraphDir {params.graphdir} --workingDir {params.workdir} --maxThreads {threads} > {log} 2>&1"

#use razers3 before optiype

#razers3 uses a lot of memory and it will not run for wgs samples with the memory that we have in our clusters.
#for this reason, we will subset fastqs in 30 batches and align them to the reference.
split_numbers = list(range(1,31))
rule fastq_split:
    input:
        get_fastq_input
    output:
        fq1=expand("results/fastq_split/{{sample}}_piece{split_no}_1.fastq.gz",split_no=split_numbers),
        fq2=expand("results/fastq_split/{{sample}}_piece{split_no}_2.fastq.gz",split_no=split_numbers)
    log:
        "logs/fastqsplitter/{sample}.log"
    benchmark:    
        "benchmarks/fastq_split/{sample}.tsv"  
    conda:
        "../envs/fastqsplitter.yaml"
    params:
        fq1=lambda w, output: [ "-o " + fq for fq in output.fq1], #the output of fastqsplitter requires "-o and output name" as many as the number of pieces that is desired
        fq2=lambda w, output: [ "-o " + fq for fq in output.fq2],
    threads: 10
    resources:
        mem_mb=50000
    shell:
        "fastqsplitter -i {input[0]} -t {threads} {params.fq1} 2> {log} && "
        "fastqsplitter -i {input[1]} -t {threads} {params.fq2} 2>> {log}"

#reads are recommended to be aligned separately. If both are supplied, core dumps or the process is killed at some point.
rule razers3:
    input:
        genome="results/refs/hs_genome.fasta",
        fq1="results/fastq_split/{sample}_piece{split_no}_1.fastq.gz",
        fq2="results/fastq_split/{sample}_piece{split_no}_2.fastq.gz",
    output:
        b1="results/razers3/mapped/{sample}_piece{split_no}_1.bam",
        b2="results/razers3/mapped/{sample}_piece{split_no}_2.bam"
    log:
        "logs/razers3/map/{sample}_{split_no}.log"
    benchmark:    
        "benchmarks/razers3/{sample}_{split_no}.tsv"  
    threads: 30
    # resources:
    #     mem_mb=lambda w, input: input.size
    conda:
        "../envs/razers3.yaml"
    shell:
        "razers3 -i 95 -m 1 -dr 0 -tc {threads} -o {output.b1} {input.genome} {input.fq1} 2> {log} && "
        "razers3 -i 95 -m 1 -dr 0 -tc {threads} -o {output.b2} {input.genome} {input.fq2} 2>> {log}"

#sort the resulting bam files to be able to use samtools merge in the next step
rule samtools_sort_razers3:
    input:
        "results/razers3/mapped/{sample}_piece{split_no}_{pair}.bam",
    output:
        bam="results/razers3/mapped/{sample}_piece{split_no}_{pair}_sorted.bam",
        idx="results/razers3/mapped/{sample}_piece{split_no}_{pair}_sorted.bam.bai",
    log:
        "logs/samtools_sort_razers3/{sample}_{split_no}_{pair}.log",
    benchmark:    
        "benchmarks/samtools_sort_razers3/{sample}_{split_no}_{pair}.tsv" 
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v1.22.0/bio/samtools/sort"

#merge the pieces into a single bam file
rule merge_bam:
    input:
        expand("results/razers3/mapped/{{sample}}_piece{split_no}_{{pair}}_sorted.bam", split_no=split_numbers)
    output:
        "results/razers3/merged/{sample}_{pair}.bam",
        idx="results/razers3/merged/{sample}_{pair}.bam.bai"
    log:
        "logs/samtools_merge/{sample}_{pair}.log",
    benchmark:    
        "benchmarks/samtools_merge/{sample}_{pair}.tsv" 
    threads: 10
    wrapper:
        "v2.3.2/bio/samtools/merge"


rule razers3_bam_to_fastq:
    input:
        "results/razers3/merged/{sample}_{pair}.bam",
    output:
        "results/razers3/reads/{sample}.{pair}.fq",
    log:
        "logs/razers3/extract/{sample}.{pair}.separate.log",
    benchmark:    
        "benchmarks/razers3_bam_to_fastq/{sample}.{pair}.tsv"  
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        "samtools fastq {input} > {output}"

rule optitype:
    input:
        # list of input reads
        # reads=get_fastq_input
        reads=["results/razers3/reads/{sample}.1.fq", "results/razers3/reads/{sample}.2.fq"]
    output:
        pdf="results/optitype/{sample}_coverage_plot.pdf",
        tsv="results/optitype/{sample}_result.tsv",
    log:
        "logs/optitype/{sample}.log"
    benchmark:    
        "benchmarks/optitype/{sample}.tsv"  
    params:
        # Type of sequencing data. Can be 'dna' or 'rna'. Default is 'dna'.
        sequencing_type="dna",
        # optiype config file, optional
        config="",
        # additional parameters
        extra=""
    threads: 20
    conda:
        "../envs/optitype.yaml"
    shell: #in case user configs have both uppercase and lowercase no_proxy values (optitype throws errors in this case)
        "unset http_proxy ftp_proxy https_proxy no_proxy; "
        "OptiTypePipeline.py -i {input.reads[0]} {input.reads[1]} --dna --outdir results/optitype --prefix {wildcards.sample}"
        
rule parse_HLAs:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/{sample}_{hla}.tsv", 
        sample=samples.sample_name,
        hla=loci
        ),
        hla_la=expand("results/HLA-LA/{sample}/hla/R1_bestguess.txt",
        sample=samples.sample_name
        ),
        arcasHLA=expand("results/arcasHLA/{sample}_{hla}/{sample}_mapped.genotype.json",
        sample=samples.sample_name,
        hla=loci
        )
    output:
        orthanq=
            "results/orthanq/final_report.csv",
        hla_la=
            "results/HLA-LA/final_report.csv",
        arcasHLA=
            "results/arcasHLA/final_report.csv",
    log:
        "logs/parse_HLAs/parse_HLA_alleles.log"
    script:
        "../scripts/parse_HLA_alleles.py"

rule compare_tools:
    input:
        orthanq="results/orthanq/orthanq_validation.tsv",
        hla_la="results/HLA-LA/final_report.csv",
        arcasHLA="results/arcasHLA/final_report.csv",
        ground_truth="resources/ground_truth/HLA-ground-truth-CEU-for-paper.tsv"
    output:
        validation="results/comparison/validation.tsv"
    log:
        "logs/comparison/compare_tools.log"
    script:
        "../scripts/validation.py"

rule datavzrd_config:
    input:
        template="resources/datavzrd.yaml",
        validation="results/comparison/validation.tsv",
        orthanq="results/orthanq/final_report.csv",
        hla_la="results/HLA-LA/final_report.csv",
        arcasHLA="results/arcasHLA/final_report.csv",
        ground_truth="resources/ground_truth/HLA-ground-truth-CEU-for-paper.tsv"
    output:
        "results/datavzrd/validation_datavzrd.yaml"
    log:
        "logs/datavzrd-config/template.log"
    template_engine:
        "yte"

rule datavzrd:
    input:
        config="results/datavzrd/validation_datavzrd.yaml",
        validation="results/comparison/validation.tsv",
        orthanq="results/orthanq/final_report.csv",
        hla_la="results/HLA-LA/final_report.csv",
        arcasHLA="results/arcasHLA/final_report.csv",
        ground_truth="resources/ground_truth/HLA-ground-truth-CEU-for-paper.tsv"
    output:
        report(
            directory("results/datavzrd-report/validation"),
            htmlindex="index.html",
            category="Validation_results",
        ),
    log:
        "logs/datavzrd/validation.log",
    wrapper:
        "v1.21.4/utils/datavzrd"
