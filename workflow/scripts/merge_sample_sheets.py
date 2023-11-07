import pandas as pd

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    #read the tables into a pd df
    # samples_evaluated = "ground_truth_for_paper/samples_evaluated.tsv"
    # indiv_no_fastq = "ground_truth_for_paper/individuals_wo_fastq_data.tsv"
    # samples_low_coverage = "ground_truth_for_paper/samples_with_low_coverage.tsv"

    samples_evaluated = snakemake.input.samples_evaluated
    indiv_no_fastq = snakemake.input.indiv_no_fastq
    samples_low_coverage = snakemake.input.samples_low_coverage

    samples_evaluated = pd.read_csv(samples_evaluated, sep = "\t", keep_default_na=False)
    indiv_no_fastq = pd.read_csv(indiv_no_fastq, sep = "\t")
    samples_low_coverage = pd.read_csv(samples_low_coverage, sep = "\t")

    #add NAs for the column "Run Accession" tp the table with no FASTQ data
    indiv_no_fastq["Run_Accession"] = "NA"

    #add a column for inclusion to the evaluation
    samples_evaluated["Inclusion"] = "yes"
    indiv_no_fastq["Inclusion"] = "no (no FASTQ data)"
    samples_low_coverage["Inclusion"] = "no (vg issues - low coverage)"

    #write SRR071132 the exclusion reason with "Inclusion" to "no"
    samples_evaluated.loc[samples_evaluated["Run_Accession"] == "SRR071132", "Inclusion"] = "no (discrepancy in ground truth)"

    #concetanate all tables together
    samples_evaluated = pd.concat([samples_evaluated,indiv_no_fastq,samples_low_coverage])
    print(samples_evaluated)

    #write concatenated table to file
    # samples_evaluated.to_csv("whole_sample_sheet.csv", index=False)
    samples_evaluated.to_csv(snakemake.output.merged, index=False)
