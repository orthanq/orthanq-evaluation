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
        "--report-fragment-ids --omit-mapq-adjustment --candidates {input.candidates} "
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
        "results/orthanq/{sample}_{hla}/{sample}_{hla}.tsv"
    conda:
        "../envs/orthanq.yaml"
    log:
        "logs/orthanq-call/{sample}_{hla}.log"
    benchmark:    
        "benchmarks/orthanq_call/{sample}_{hla}.tsv"
    shell:
        "../orthanq/target/release/orthanq call --haplotype-calls {input.calls} --haplotype-variants {input.candidate_variants} "
        "--xml {input.xml} --max-haplotypes 5 --prior diploid --output {output} 2> {log}"

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

rule parse_and_validate:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/{sample}_{hla}.tsv", 
        sample=samples.sample_name,
        hla=loci
        ),
        ground_truth="resources/ground_truth/HLA-ground-truth-CEU-for-paper.tsv"
    output:
        validation="results/comparison/comparison.tsv"
    log:
        "logs/comparison/compare_and_validate.log"
    script:
        "../scripts/parse_and_validate_diploid_subclonal.py"

rule parse_HLAs:
    input:
        orthanq=expand("results/orthanq/{sample}_{hla}/{sample}_{hla}.tsv", 
        sample=samples.sample_name,
        hla=loci
        ),
        # hla_la=expand("results/HLA-LA/{sample}/hla/R1_bestguess.txt",
        # sample=samples.sample_name
        # ),
        # arcasHLA=expand("results/arcasHLA/{sample}_{hla}/{sample}_mapped.genotype.json",
        # sample=samples.sample_name,
        # hla=loci
        # )
    output:
        orthanq=
            "results/orthanq/final_report.csv",
        # hla_la=
        #     "results/HLA-LA/final_report.csv",
        # arcasHLA=
        #     "results/arcasHLA/final_report.csv",
    log:
        "logs/parse_HLAs/parse_HLA_alleles.log"
    script:
        "../scripts/parse_HLA_alleles.py"

# rule compare_tools:
#     input:
#         orthanq="results/orthanq/final_report.csv",
#         # hla_la="results/HLA-LA/final_report.csv",
#         # arcasHLA="results/arcasHLA/final_report.csv",
#         ground_truth="resources/ground_truth/HLA-ground-truth-CEU-for-paper.tsv"
#     output:
#         validation="results/comparison/comparison.tsv"
#     log:
#         "logs/comparison/compare_tools.log"
#     script:
#         "../scripts/validation.py"

rule datavzrd_config:
    input:
        template="resources/datavzrd.yaml",
        validation="results/comparison/comparison.tsv",
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
        validation="results/comparison/comparison.tsv",
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
