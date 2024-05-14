import pandas as pd
import os

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #
    new_dict = {'sample': [], 'locus': [], 'read_count': []}

    #read sample sheet
    read_count_paths = snakemake.input.read_counts

    #loop through each report
    for (index, _) in enumerate(read_count_paths):

        #retrieve sample and locus names
        splitted = os.path.basename(read_count_paths[index]).split("_")
        sample_name = splitted[0]
        locus_name = splitted[1].split(".")[0]

        #get read count
        read_count_pd = pd.read_csv(read_count_paths[index], header=None, sep="\t")
        read_count = read_count_pd.iloc[0, 3]
        print(read_count)

        #then add the row to the dataframe
        new_dict['sample'].append(sample_name)
        new_dict['locus'].append(locus_name)
        new_dict['read_count'].append(read_count)

    new_df = pd.DataFrame(new_dict)

    #write table to output
    new_df.to_csv(snakemake.output.table, index=False)