import pandas as pd
import os

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #read sample sheet
    sample_sheet = pd.read_csv(snakemake.input.sample_sheet)

    #initialize columns for read length and count 
    sample_sheet["Read Length"] = ""
    sample_sheet["Read Count"] = ""

    for file in snakemake.input.sample_calculation:
        #read calculation file
        calc = pd.read_csv(file, names=["calc"])

        #parse count and length
        read_count = calc.loc[0, "calc"] 
        read_length = calc.loc[1, "calc"]

        #parse sample name from the file
        print(file)
        basename=os.path.basename(file)
        splitted = basename.split("_")
        if "D1_S1_L001" in basename:
            sample_name = "SRR2962669"
        else:
            sample_name = splitted[0].split(".")[0]
        print(sample_name)

        #insert the read length and read count to the sample sheet
        sample_sheet.loc[sample_sheet["Run Accession"] == sample_name, "Read Length"] = read_length
        sample_sheet.loc[sample_sheet["Run Accession"] == sample_name, "Read Count"] = read_count

        print(sample_sheet)

    sample_sheet.to_csv(snakemake.output.updated_sample_sheet, index=False)



