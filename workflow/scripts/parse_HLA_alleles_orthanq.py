import pandas as pd
import os
import json

#orthanq
orthanq_input = snakemake.input.orthanq

#initialize the dataframe
orthanq_final_table = pd.DataFrame(columns=('sample', 'A', 'B', 'C', 'DQA1','DQB1'))

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    for index in range(len(orthanq_input)):
        if "D1_S1_L001" in orthanq_input[index]:
            splitted = os.path.dirname(orthanq_input[index]).split("_")
            sample_name = splitted[0] + splitted[1] + splitted[2] 
            locus_name = splitted[3].split(".")[0]
            #rename sample name for giab to the accesssion id
            sample_name = "SRR2962669"
        else:
            splitted = os.path.dirname(orthanq_input[index]).split("_")
            sample_name = splitted[0]
            locus_name = splitted[1].split(".")[0]
        if not sample_name in orthanq_final_table['sample'].tolist():
            # new_row = {'sample': sample_name, 'A': [], 'B': [], 'C': [], 'DQA1': [], 'DQB1': []}
            new_row = pd.DataFrame([[sample_name, '', '', '', '', '']],
                    columns=['sample', 'A', 'B', 'C', 'DQA1', 'DQB1'])
            # orthanq_final_table = orthanq_final_table.append(new_row, ignore_index=True)
            orthanq_final_table = pd.concat([orthanq_final_table, new_row], ignore_index=True)

        ##more than one best record
        results = pd.read_csv(orthanq_input[index])

        ##find the records that have the same density
        best_odds = [1, 1.0, 1.00]
        best_results = results[results.odds.isin(best_odds)]
        print(sample_name)
        print(locus_name)
        print("best results: ", best_results)

        all_combinations=[]
        if not best_results.empty: #orthanq doesn't have predictions for some samples
            #retrieve the predicted haplotypes
            #collect locus names
            filtered_cols = []
            filtered_cols = [col for col in results if col.startswith(locus_name)]
            # print(filtered_cols)
            #loop over best results
            for (i,result_row) in best_results.iterrows():
                value_to_add = []     
                #collect first two fields and cumulative fractions of haplotypes
                for col in filtered_cols:
                    if best_results[col][i] == 0.5:
                        value_to_add.append(col)
                    elif best_results[col][i] == 1.0:
                        value_to_add.append(col) #two times the value, homozygous sample handling
                        value_to_add.append(col)
                value_to_add = '/'.join(value_to_add)
                all_combinations.append(value_to_add)
        print(all_combinations) 
        orthanq_final_table.loc[orthanq_final_table['sample'] == sample_name, locus_name] = ','.join(all_combinations)

    orthanq_final_table.to_csv(
        snakemake.output.orthanq, sep="\t", index=False, header=True
    )
