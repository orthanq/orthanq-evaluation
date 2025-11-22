rule parse_HLAs_orthanq:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/predictions.csv", 
        sample=samples.sample_name,
        hla=loci
        )
    output:
        orthanq=
            "results/orthanq/final_report.csv",
    log:
        "logs/parse_HLAs/parse_HLA_alleles.log"
    script:
        "../scripts/parse_HLA_alleles_orthanq.py"

rule validate_orthanq:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/predictions.csv", 
        sample=samples.sample_name,
        hla=loci
        ),
        ground_truth_evaluated="resources/ground_truth/1K_CEU_evaluated.tsv",
        allele_freqs="resources/allele_freqs/allele_frequencies.csv",
        orthanq_final_table="results/orthanq/final_report.csv"
    output:
        validation_all="results/validation/orthanq_validation_all.tsv",
        orthanq_A_tp_fp="results/orthanq/A_tp_fp.tsv",
        orthanq_B_tp_fp="results/orthanq/B_tp_fp.tsv",
        orthanq_C_tp_fp="results/orthanq/C_tp_fp.tsv",
        orthanq_DQB1_tp_fp="results/orthanq/DQB1_tp_fp.tsv",
        threshold_results="results/threshold_results/threshold_results.tsv"
    log:
        "logs/validate_orthanq/validate_orthanq.log"
    script:
        "../scripts/parse_and_validate_diploid.py"

rule gather_benchmark:
    input: 
        benchmarks = "benchmarks", #to avoid too many inputs
        sample_sheet = config["samples"]
    output:
        runtimes_table = "results/runtimes/runtimes.csv",
        runtimes_plot = "results/runtimes/runtimes.json",
    conda:
        "../envs/altair.yaml"
    log:
        "logs/runtimes/runtimes.log"
    script:
        "../scripts/runtimes.py"

rule calculate_lp_pruned:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/predictions.csv", 
        sample=samples.sample_name,
        hla=loci
        ),
    output:
        pruned_haplotypes="results/calculate_lp_pruned/pruned_table.tsv"
    log:
        "logs/calculate_lp_pruned/calculate_lp_pruned.log"
    script:
        "../scripts/calculate_lp_pruned.py"    

rule plot_lp_pruned:
    input:
        pruned_table="results/calculate_lp_pruned/pruned_table.tsv"
    output:
        pruned_json="results/calculate_lp_pruned/pruned.json"
    conda:
        "../envs/altair.yaml"
    log:
        "logs/plot_lp_pruned/plot_lp_pruned.log"
    script:
        "../scripts/plot_lp_pruned.py"   

rule plot_validation_threshold:
    input:
        threshold_results="results/threshold_results/threshold_results.tsv"
    output:
        threshold_results_line_plot="results/threshold_results/threshold_results_plot.json"
    conda:
        "../envs/altair.yaml"
    log:
        "logs/plot_threshold_accuracy/plot_threshold_accuracy.log"
    script:
        "../scripts/plot_threshold_accuracy.py"   

# find read counts per locus for mapq > 20 (just for demonstration purposes, to be put in the rebuttal)
rule extract_read_count_hla_genes:
    input:
        bam="results/vg/alignment/{sample}_vg.sorted.reheadered.extracted.bam",
        region="resources/HLA_regions/genes/{gene}.bed"
    output:
        bed="results/count_coverage_per_gene/{sample}_{gene}.bed"
    conda:
        "../envs/bedops.yaml"
    log:
        "logs/count_coverage_per_gene/{sample}_{gene}.log"
    shell:
        "set +o pipefail; "
        " samtools view -h -q 20 {input.bam} | "
        "bam2bed --reduced - | bedmap --echo --count --delim '\t' {input.region} - > {output} 2> {log}"

rule export_read_count_to_table:
    input:
        read_counts=expand("results/count_coverage_per_gene/{sample}_{gene}.bed", sample=samples.sample_name, gene=["A", "B", "C", "DQB1"])
    output:
        table="results/count_coverage_per_gene/read_counts.csv"
    log:
        "logs/count_coverage_per_gene/export_read_count_to_table.log"
    script:
        "../scripts/export_read_count_to_table.py"

rule vg2svg_orthanq:
    input:
        three_field="results/orthanq/{sample}_{hla}/3_field_solutions.json",
        two_field="results/orthanq/{sample}_{hla}/2_field_solutions.json",
        lp_solution="results/orthanq/{sample}_{hla}/lp_solution.json",
        final_solution="results/orthanq/{sample}_{hla}/best_solution.json",
        arrow_plot="results/orthanq/{sample}_{hla}/arrow_plot.json"
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
        final_solution=report("results/orthanq/{sample}_{hla}/best_solution.html",category="Orthanq detailed solutions", subcategory="{sample}_{hla}", labels={
            "sample": "{sample}",
            "locus": "{hla}",
            "figure": "final solution"
        }),
        arrow_plot= report("results/orthanq/{sample}_{hla}/arrow_plot.html",category="Orthanq detailed solutions", subcategory="{sample}_{hla}", labels={
            "sample": "{sample}",
            "locus": "{hla}",
            "figure": "arrow plot"
        })      
    log:
        "logs/vg2svg/orthanq/{sample}_{hla}.log",
    conda:
        "../envs/vega.yaml"
    shell:
        "vl-convert vl2html --input {input.three_field} --output {output.three_field} 2> {log} && "
        "vl-convert vl2html --input {input.two_field} --output {output.two_field} 2>> {log} && "
        "vl-convert vl2html --input {input.lp_solution} --output {output.lp_solution} 2>> {log} && "
        "vl-convert vl2html --input {input.final_solution} --output {output.final_solution} 2>> {log} && "
        "vl-convert vl2html --input {input.arrow_plot} --output {output.arrow_plot} 2>> {log} "


rule datavzrd_runtimes:
    input:
        config="resources/datavzrd/runtimes.yaml",
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
    group: "runtimes"
    wrapper:
        "v3.10.2/utils/datavzrd"

rule datavzrd_orthanq:
    input:
        config="resources/datavzrd/orthanq.yaml",
        orthanq="results/orthanq/final_report.csv",
    output:
        report(
            directory("results/datavzrd-report/orthanq"),
            htmlindex="index.html",
            category="Accuracy", labels={
            "type": "table",
            "name": "orthanq predictions"
        }
        ),
    log:
        "logs/datavzrd/orthanq.log",
    group: "orthanq"
    wrapper:
        "v3.10.2/utils/datavzrd"
