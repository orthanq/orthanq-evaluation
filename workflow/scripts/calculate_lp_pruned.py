import pandas as pd
import os
import json
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #orthanq results
    orthanq_input = snakemake.input.orthanq

    #number of haplotypes per locus
    n_haplotypes = {'A': 142, 'B': 193, 'C': 121, 'DQB1': 42}

    #initialize a dictionary to store percentages
    percentage_dict = {'record': [], 'percentage': []}

    #iterate
    for (index, _) in enumerate(orthanq_input):
        #retrieve sample and locus names
        splitted = os.path.basename(orthanq_input[index]).split("_")
        sample_name = splitted[0]
        locus_name = splitted[-1].split(".")[0]

        #read orthanq result
        results = pd.read_csv(orthanq_input[index])

        #read predictions
        haplotypes = [col for col in results if col.startswith(locus_name)]

        #find percentage of haplotypes per locus
        print(sample_name)
        print(locus_name)
        if locus_name != 'DQA1':
            percentage = (n_haplotypes[locus_name] - len(haplotypes))/n_haplotypes[locus_name]*100
            print(percentage)

            #insert percentage dict
            percentage_dict['record'].append(sample_name + "_" + locus_name)
            percentage_dict['percentage'].append(percentage)
    
    print(percentage_dict)
    percentage_df = pd.DataFrame(percentage_dict)

    percentage_df.to_csv(
        snakemake.output.pruned_haplotypes, sep="\t", index=False, header=True
    )