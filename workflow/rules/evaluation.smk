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
    resources:
        mem_mb=120000
    threads: 40
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
    threads: 40
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
    threads: 40
    resources:
        mem_mb=60000
    shell:
        "fastqsplitter -i {input[0]} -t {threads} {params.fq1} 2> {log} && "
        "fastqsplitter -i {input[1]} -t {threads} {params.fq2} 2>> {log}"

#reads are recommended to be aligned separately. If both are supplied, core dumps or the process is killed at some point.
rule razers3:
    input:
        genome="resources/optitype/hla_reference_dna.fasta",
        fq1="results/fastq_split/{sample}_piece{split_no}_1.fastq.gz",
        fq2="results/fastq_split/{sample}_piece{split_no}_2.fastq.gz",
    output:
        b1="results/razers3/mapped/{sample}_piece{split_no}_1.bam",
        b2="results/razers3/mapped/{sample}_piece{split_no}_2.bam"
    log:
        "logs/razers3/map/{sample}_{split_no}.log"
    benchmark:    
        "benchmarks/razers3/{sample}_{split_no}.tsv"  
    threads: 40
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
    threads: 40
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
    threads: 40
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
    threads: 40
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
    threads: 40
    conda:
        "../envs/optitype.yaml"
    shell: #in case user configs have both uppercase and lowercase no_proxy values (optitype throws errors in this case)
        "unset http_proxy ftp_proxy https_proxy no_proxy; "
        "OptiTypePipeline.py -i {input.reads[0]} {input.reads[1]} --dna --outdir results/optitype --prefix {wildcards.sample}"

rule merge_sample_sheets:
    input:
        samples_evaluated="resources/ground_truth/1K_CEU_evaluated.tsv",
        samples_incomplete_truth="resources/ground_truth/1K_CEU_incomplete_HLA.tsv",
        samples_low_coverage="resources/ground_truth/1K_CEU_low_coverage.tsv",
    output:
        merged="resources/ground_truth/merged_sample_sheet.csv"
    conda:
        "../envs/altair.yaml"
    log:
        "logs/merge_sample_sheets/merge_sample_sheets.log",
    script:
        "../scripts/merge_sample_sheets.py"

rule find_read_length_and_count:
    input:
        sample_sheet="resources/ground_truth/merged_sample_sheet.csv",
        fastq_dir=config["all_fastq_dir"],
        fastq_dir_2=config["all_fastq_dir_additional"]
    output:
        sample_sheet="resources/ground_truth/merged_sample_sheet_w_read_info.csv"
    conda:
        "../envs/altair.yaml"
    threads: 20
    log:
        "logs/insert_read_info/insert_read_info.log",
    script:
        "../scripts/insert_read_info.py"

rule parse_HLAs:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/{sample}_{hla}.csv", 
        sample=samples.sample_name,
        hla=loci
        ),
        hla_la=expand("results/HLA-LA/{sample}/hla/R1_bestguess.txt",
        sample=samples.sample_name
        ),
        arcasHLA=expand("results/arcasHLA/{sample}_{hla}/{sample}_mapped.genotype.json",
        sample=samples.sample_name,
        hla=loci
        ),
        optitype=expand("results/optitype/{sample}_result.tsv",
        sample=samples.sample_name
        ),
    output:
        orthanq=
            "results/orthanq/final_report.csv",
        hla_la=
            "results/HLA-LA/final_report.csv",
        arcasHLA=
            "results/arcasHLA/final_report.csv",
        optitype=
            "results/optitype/final_report.csv"
    log:
        "logs/parse_HLAs/parse_HLA_alleles.log"
    script:
        "../scripts/parse_HLA_alleles.py"

rule validate_orthanq:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/{sample}_{hla}.csv", 
        sample=samples.sample_name,
        hla=loci
        ),
        ground_truth_evaluated="resources/ground_truth/1K_CEU_evaluated.tsv",
        ground_truth_low_coverage="resources/ground_truth/1K_CEU_low_coverage.tsv",
        allele_freqs="resources/allele_freqs/allele_frequencies.csv",
        orthanq_final_table="results/orthanq/final_report.csv"
    output:
        validation_low="results/validation/orthanq_validation_low.tsv",
        validation_high="results/validation/orthanq_validation_high.tsv",
        validation_all="results/validation/orthanq_validation_all.tsv",
        tp_fp_table="results/validation/tp_fp_table_all.tsv",
        orthanq_A_tp_fp="results/orthanq/A_tp_fp.tsv",
        orthanq_B_tp_fp="results/orthanq/B_tp_fp.tsv",
        orthanq_C_tp_fp="results/orthanq/C_tp_fp.tsv",
        orthanq_DQB1_tp_fp="results/orthanq/DQB1_tp_fp.tsv",
    log:
        "logs/validate_orthanq/validate_orthanq.log"
    script:
        "../scripts/parse_and_validate_diploid.py"

rule plot_tp_fp:
    input:
        tp_fp_table="results/validation/tp_fp_table_all.tsv",
    output:
        tp_fp_plot="results/validation/tp_fp_table_all.json",
    conda:
        "../envs/altair.yaml"
    log:
        "logs/tp_fp_plot/tp_fp_plot.log"
    script:
        "../scripts/tp_fp_plot.py"

rule validate_tools:
    input:
        orthanq_low="results/validation/orthanq_validation_low.tsv",
        orthanq_high="results/validation/orthanq_validation_high.tsv",
        orthanq_all="results/validation/orthanq_validation_all.tsv",
        orthanq_A_tp_fp="results/orthanq/A_tp_fp.tsv",
        orthanq_B_tp_fp="results/orthanq/B_tp_fp.tsv",
        orthanq_C_tp_fp="results/orthanq/C_tp_fp.tsv",
        orthanq_DQB1_tp_fp="results/orthanq/DQB1_tp_fp.tsv",
        hla_la="results/HLA-LA/final_report.csv",
        arcasHLA="results/arcasHLA/final_report.csv",
        ground_truth="resources/ground_truth/1K_CEU_all.tsv",
        optitype="results/optitype/final_report.csv",
        samples_evaluated="resources/ground_truth/1K_CEU_evaluated.tsv",
        samples_low_coverage="resources/ground_truth/1K_CEU_low_coverage.tsv",
        allele_freqs="resources/allele_freqs/allele_frequencies.csv"
    output:
        validation_low="results/validation/validation_low.tsv",
        validation_high="results/validation/validation_high.tsv",
        validation_all="results/validation/validation_all.tsv",
        A_tp_fp="results/validation/A_tp_fp.tsv",
        B_tp_fp="results/validation/B_tp_fp.tsv",
        C_tp_fp="results/validation/C_tp_fp.tsv",
        DQB1_tp_fp="results/validation/DQB1_tp_fp.tsv"
    log:
        "logs/validation/validation.log"
    script:
        "../scripts/validation.py"

rule evaluation_plot:
    input:
        template="resources/templates/evaluation_plot.json",
        validation_low="results/validation/validation_low.tsv",
        validation_high="results/validation/validation_high.tsv",
    output:
        plot_low="results/evaluation/evaluation_plot_low.json",
        plot_high="results/evaluation/evaluation_plot_high.json"
    log:
        "logs/evaluation/evaluation_plot.log"
    script:
        "../scripts/evaluation_plot.py"

rule gather_benchmark:
    input: 
        benchmarks = "benchmarks" #to avoid too many inputs
    output:
        runtimes_table = "results/runtimes/runtimes.csv",
        runtimes_plot = "results/runtimes/runtimes.json",
    conda:
        "../envs/gather_benchmarks.yaml"
    log:
        "logs/runtimes/runtimes.log"
    script:
        "../scripts/runtimes.py"

rule vg2svg_orthanq:
    input:
        three_field="results/orthanq/{sample}_{hla}/3_field_solutions.json",
        two_field="results/orthanq/{sample}_{hla}/2_field_solutions.json",
        lp_solution="results/orthanq/{sample}_{hla}/lp_solution.json",
        final_solution="results/orthanq/{sample}_{hla}/final_solution.json"
    output:
        three_field=report("results/orthanq/{sample}_{hla}/3_field_solutions.html",category="Orthanq detailed solutions", subcategory="{sample}_{hla}",labels={
            "sample": "{sample}",
            "locus": "{hla}",
            "figure": "3-field solutions"
        }),
        two_field=report("results/orthanq/{sample}_{hla}/2_field_solutions.html",category="Orthanq detailed solutions", subcategory="{sample}_{hla}", labels={
            "sample": "{sample}",
            "locus": "{hla}",
            "figure": "2-field solutions"
        }),
        lp_solution=report("results/orthanq/{sample}_{hla}/lp_solution.html",category="Orthanq detailed solutions", subcategory="{sample}_{hla}", labels={
            "sample": "{sample}",
            "locus": "{hla}",
            "figure": "lp solution"
        }),
        final_solution=report("results/orthanq/{sample}_{hla}/final_solution.html",category="Orthanq detailed solutions", subcategory="{sample}_{hla}", labels={
            "sample": "{sample}",
            "locus": "{hla}",
            "figure": "final solution"
        })        
    log:
        "logs/vg2svg/orthanq/{sample}_{hla}.log",
    conda:
        "../envs/vega.yaml"
    shell:
        "vl-convert vl2html --input {input.three_field} --output {output.three_field} 2> {log} && "
        "vl-convert vl2html --input {input.two_field} --output {output.two_field} 2>> {log} && "
        "vl-convert vl2html --input {input.lp_solution} --output {output.lp_solution} 2>> {log} && "
        "vl-convert vl2html --input {input.final_solution} --output {output.final_solution} 2>> {log}"

rule vg2svg_evaluation:
    input:
        runtimes_plot = "results/runtimes/runtimes.json",
        evaluation_low = "results/evaluation/evaluation_plot_low.json",
        evaluation_high = "results/evaluation/evaluation_plot_high.json",
        tp_fp_plot = "results/validation/tp_fp_table_all.json"
    output:
        runtimes_svg="results/vega/runtimes.svg",
        evaluation_low_svg="results/vega/evaluation_low.svg",
        evaluation_high_svg="results/vega/evaluation_high.svg",
        tp_fp_plot_svg = "results/vega/tp_fp_table_all.svg",
        runtimes_html=report("results/vega/runtimes.html", category="Runtime performance", labels={
            "type": "figure"
        }),
        evaluation_low_html=report("results/vega/evaluation_low.html", category="Accuracy", labels={
            "name": "accuracy comparison (low coverage)",
            "type": "figure"
        }),
        evaluation_high_html=report("results/vega/evaluation_high.html", category="Accuracy", labels={
            "name": "accuracy comparison",
            "type": "figure"
        }),
        tp_fp_plot_html=report("results/vega/tp_fp_table_all.html", category="Orthanq density accuracy", labels={
            "figure": "log-scaled density plot"
        })
    log:
        "logs/vg2svg/evaluation_plots.log",
    conda:
        "../envs/vega.yaml"
    shell:
        "vl2svg {input.runtimes_plot} {output.runtimes_svg} 2> {log} && "
        "vl-convert vl2html --input {input.runtimes_plot} --output {output.runtimes_html} 2>> {log} && "
        "vl2svg {input.evaluation_low} {output.evaluation_low_svg} 2>> {log} && "
        "vl-convert vl2html --input {input.evaluation_low} --output {output.evaluation_low_html} 2>> {log} && "
        "vl2svg {input.evaluation_high} {output.evaluation_high_svg} 2>> {log} && "
        "vl-convert vl2html --input {input.evaluation_high} --output {output.evaluation_high_html} 2>> {log} && "
        "vl2svg {input.tp_fp_plot} {output.tp_fp_plot_svg} 2>> {log} && "
        "vl-convert vl2html --input {input.tp_fp_plot} --output {output.tp_fp_plot_html} 2>> {log} "

rule datavzrd_config_runtimes:
    input:
        template="resources/datavzrd/runtimes.yaml",
        runtimes_table = "results/runtimes/runtimes.csv",
    output:
        "results/datavzrd/runtimes.yaml"
    log:
        "logs/datavzrd-config/runtimes.log"
    template_engine:
        "yte"

rule datavzrd_runtimes:
    input:
        config="results/datavzrd/runtimes.yaml",
        runtimes_table = "results/runtimes/runtimes.csv",
    output:
        report(
            directory("results/datavzrd-report/runtimes"),
            htmlindex="index.html",
            category="Runtime performance", labels={
            "type": "table"
        }
        ),
    log:
        "logs/datavzrd/runtimes.log",
    wrapper:
        "v2.13.0/utils/datavzrd"

rule datavzrd_config_orthanq:
    input:
        template="resources/datavzrd/orthanq.yaml",
        orthanq="results/orthanq/final_report.csv",
    output:
        "results/datavzrd/orthanq.yaml"
    log:
        "logs/datavzrd-config/orthanq.log"
    template_engine:
        "yte"

rule datavzrd_orthanq:
    input:
        config="results/datavzrd/orthanq.yaml",
        orthanq="results/orthanq/final_report.csv",
    output:
        report(
            directory("results/datavzrd-report/orthanq"),
            htmlindex="index.html",
            category="Accuracy", labels={
            "type": "figure",
            "name": "orthanq predictions"
        }
        ),
    log:
        "logs/datavzrd/orthanq.log",
    wrapper:
        "v2.13.0/utils/datavzrd"

rule datavzrd_config_optitype:
    input:
        template="resources/datavzrd/optitype.yaml",
        optitype="results/optitype/final_report.csv",
    output:
        "results/datavzrd/optitype.yaml"
    log:
        "logs/datavzrd-config/optitype.log"
    template_engine:
        "yte"

rule datavzrd_optitype:
    input:
        config="results/datavzrd/optitype.yaml",
        optitype="results/optitype/final_report.csv",
    output:
        report(
            directory("results/datavzrd-report/optitype"),
            htmlindex="index.html",
            category="Accuracy", labels={
            "type": "table",
            "name": "optitype predictions"
        }
        ),
    log:
        "logs/datavzrd/optitype.log",
    wrapper:
        "v2.13.0/utils/datavzrd"

rule datavzrd_config_arcashla:
    input:
        template="resources/datavzrd/arcashla.yaml",
        arcasHLA="results/arcasHLA/final_report.csv",
    output:
        "results/datavzrd/arcashla.yaml"
    log:
        "logs/datavzrd-config/arcashla.log"
    template_engine:
        "yte"

rule datavzrd_arcashla:
    input:
        config="results/datavzrd/arcashla.yaml",
        arcasHLA="results/arcasHLA/final_report.csv",
    output:
        report(
            directory("results/datavzrd-report/arcashla"),
            htmlindex="index.html",
            category="Accuracy", labels={
            "type": "table",
            "name": "arcasHLA predictions"
        }
        ),
    log:
        "logs/datavzrd/arcashla.log",
    wrapper:
        "v2.13.0/utils/datavzrd"


rule datavzrd_config_hla_la:
    input:
        template="resources/datavzrd/hla_la.yaml",
        hla_la="results/HLA-LA/final_report.csv",
    output:
        "results/datavzrd/hla_la.yaml"
    log:
        "logs/datavzrd-config/hla_la.log"
    template_engine:
        "yte"

rule datavzrd_hla_la:
    input:
        config="results/datavzrd/hla_la.yaml",
        hla_la="results/HLA-LA/final_report.csv",
    output:
        report(
            directory("results/datavzrd-report/hla_la"),
            htmlindex="index.html",
            category="Accuracy", labels={
            "type": "table",
            "name": "HLA-LA predictions"
        }
        ),
    log:
        "logs/datavzrd/hla_la.log",
    wrapper:
        "v2.13.0/utils/datavzrd"

rule datavzrd_config_comparison_low:
    input:
        template="resources/datavzrd/accuracy_comparison_low.yaml",
        validation="results/validation/validation_low.tsv",
    output:
        "results/datavzrd/accuracy_comparison_low.yaml"
    log:
        "logs/datavzrd-config/accuracy_comparison_low.log"
    template_engine:
        "yte"

rule datavzrd_comparison_low:
    input:
        config="results/datavzrd/accuracy_comparison_low.yaml",
        validation="results/validation/validation_low.tsv",
    output:
        report(
            directory("results/datavzrd-report/accuracy_comparison_low"),
            htmlindex="index.html",
            category="Accuracy", labels={
            "name": "accuracy comparison (low coverage)",
            "type": "table"
        }
        ),
    log:
        "logs/datavzrd/accuracy_comparison_low.log",
    wrapper:
        "v2.13.0/utils/datavzrd"

rule datavzrd_config_comparison_high:
    input:
        template="resources/datavzrd/accuracy_comparison_high.yaml",
        validation="results/validation/validation_high.tsv",
    output:
        "results/datavzrd/accuracy_comparison_high.yaml"
    log:
        "logs/datavzrd-config/accuracy_comparison_high.log"
    template_engine:
        "yte"

rule datavzrd_comparison_high:
    input:
        config="results/datavzrd/accuracy_comparison_high.yaml",
        validation="results/validation/validation_high.tsv",
    output:
        report(
            directory("results/datavzrd-report/accuracy_comparison_high"),
            htmlindex="index.html",
            category="Accuracy", labels={
            "name": "accuracy comparison",
            "type": "table"
        }
        ),
    log:
        "logs/datavzrd/accuracy_comparison_high.log",
    wrapper:
        "v2.13.0/utils/datavzrd"

rule datavzrd_config_evaluation_tables:
    input:
        template="resources/datavzrd/evaluation_tables.yaml",
        A_tp_fp="results/validation/A_tp_fp.tsv",
        B_tp_fp="results/validation/B_tp_fp.tsv",
        C_tp_fp="results/validation/C_tp_fp.tsv",
        DQB1_tp_fp="results/validation/DQB1_tp_fp.tsv"
    output:
        "results/datavzrd/evaluation_tables.yaml"
    log:
        "logs/datavzrd-config/evaluation_tables.log"
    template_engine:
        "yte"

rule datavzrd_evaluation_tables:
    input:
        config="results/datavzrd/evaluation_tables.yaml",
        A_tp_fp="results/validation/A_tp_fp.tsv",
        B_tp_fp="results/validation/B_tp_fp.tsv",
        C_tp_fp="results/validation/C_tp_fp.tsv",
        DQB1_tp_fp="results/validation/DQB1_tp_fp.tsv"
    output:
        report(
            directory("results/datavzrd-report/evaluation_tables"),
            htmlindex="index.html",
            category="Evaluation (all tools)", labels={
            "name": "evaluation table",
            "type": "table"
        }
        ),
    log:
        "logs/datavzrd/evaluation_tables.log",
    wrapper:
        "v2.13.0/utils/datavzrd"

rule datavzrd_config_sample_sheet:
    input:
        template="resources/datavzrd/sample_sheet.yaml",
        sample_sheet="resources/ground_truth/merged_sample_sheet_w_read_info.csv"
    output:
        "results/datavzrd/sample_sheet.yaml"
    log:
        "logs/datavzrd-config/template.log"
    template_engine:
        "yte"

rule datavzrd_sample_sheet:
    input:
        config="results/datavzrd/sample_sheet.yaml",
        sample_sheet="resources/ground_truth/merged_sample_sheet_w_read_info.csv"
    output:
        report(
            directory("results/datavzrd-report/sample_sheet"),
            htmlindex="index.html",
            category="Sample sheet",
            labels={
            "name": "sample sheet",
            "type": "table"
        }
        ),
    log:
        "logs/datavzrd/sample_sheet.log",
    wrapper:
        "v2.13.0/utils/datavzrd"
        