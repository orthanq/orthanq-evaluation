import pandas as pd

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #read the tables into a pd df
    samples_evaluated = snakemake.input.samples_evaluated
    samples_incomplete = snakemake.input.samples_incomplete_truth

    samples_evaluated = pd.read_csv(samples_evaluated, sep = "\t", keep_default_na=False)
    samples_incomplete = pd.read_csv(samples_incomplete, sep = "\t")

    #add NAs for the column "Run Accession" tp the table with no FASTQ data
    samples_incomplete["Run_Accession"] = "NA"

    #add a column for inclusion to the evaluation
    samples_evaluated["Inclusion"] = "yes"
    samples_incomplete["Inclusion"] = "no (incompletely typed)"
    #add the second reason for NA11894
    samples_incomplete.loc[samples_incomplete["Sample_ID"] == "NA11894", "Inclusion"] = "no (incompletely typed & no WES or WGS data)"

    #write SRR071132 the exclusion reason with "Inclusion" to "no"
    samples_evaluated.loc[samples_evaluated["Run_Accession"] == "SRR071132", "Inclusion"] = "no (discrepancy in ground truth)"
    samples_evaluated.loc[samples_evaluated["Run_Accession"] == "SRR1601876", "Inclusion"] = "no (discrepancy in ground truth)"
    
    #concetanate all tables together
    # samples_evaluated = pd.concat([samples_evaluated,indiv_no_fastq,samples_low_coverage])
    samples_evaluated = pd.concat([samples_evaluated,samples_incomplete])

    # add a column for sequencing method
    samples_evaluated["Sequencing method"] = "WES"

    # correct the sequencing method for the following WGS samples, SRR1601854, ERR194147, ERR194160, ERR194161
    samples_evaluated.loc[samples_evaluated["Run_Accession"] == "SRR1601854", "Sequencing method"] = "WGS"
    samples_evaluated.loc[samples_evaluated["Run_Accession"] == "ERR194147", "Sequencing method"] = "WGS"
    samples_evaluated.loc[samples_evaluated["Run_Accession"] == "ERR194160", "Sequencing method"] = "WGS"
    samples_evaluated.loc[samples_evaluated["Run_Accession"] == "ERR194161", "Sequencing method"] = "WGS"

    print(samples_evaluated)

    #write concatenated table to file
    # samples_evaluated.to_csv("whole_sample_sheet.csv", index=False)
    samples_evaluated.to_csv(snakemake.output.merged, index=False)
