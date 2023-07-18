# This script parses the results of HLA typing that is performed by HLA-LA, arcasHLA and orthanq.
# It then creates separate tables for each in a unified manner.abs
# Sidenote for Orthanq: orthanq may output several combinations of HLA genotypes per locus. We output all of them separated by commas.

import pandas as pd
import os
import json

#orthanq
orthanq_input = snakemake.input.orthanq

#initialize the dataframe
orthanq_final_table = pd.DataFrame(columns=('sample', 'A', 'B', 'C', 'DQA1','DQB1'))

for index in range(len(orthanq_input)):
    splitted = os.path.basename(orthanq_input[index]).split("_")
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
    print("best results: ", best_results)
    
    #retrieve the predicted haplotypes
    #collect locus names
    filtered_cols = []
    filtered_cols = [col for col in results if col.startswith(locus_name)]
    print(filtered_cols)
    #loop over best results
    all_combinations=[]
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
    
    orthanq_final_table.loc[orthanq_final_table['sample'] == sample_name, locus_name] = ','.join(all_combinations)

orthanq_final_table.to_csv(
    snakemake.output.orthanq, sep="\t", index=False, header=True
)

# #arcasHLA
# arcasHLA_input = snakemake.input.arcasHLA

# #initialize the dataframe
# arcasHLA_final_table = pd.DataFrame(columns=('sample', 'A', 'B', 'C', 'DQA1','DQB1'))
# for index in range(len(arcasHLA_input)):
#     splitted = os.path.basename(os.path.dirname(arcasHLA_input[index])).split("_")
#     sample_name = splitted[0]
#     locus_name = splitted[1].split(".")[0]
#     if not sample_name in arcasHLA_final_table['sample'].tolist():
#         new_row = {'sample': sample_name, 'A': [], 'B': [], 'C': [], 'DQA1': [], 'DQB1': []}
#         arcasHLA_final_table = arcasHLA_final_table.append(new_row, ignore_index=True)
#     # Opening JSON file
#     f = open(arcasHLA_input[index])
    
#     # returns JSON object as 
#     # a dictionary
#     data = json.load(f)

#     filtered_cols = data[locus_name]
#     value_to_add = []
#     for col in filtered_cols:
#         value_to_add.append(col)
#     value_to_add = '/'.join(value_to_add)
#     arcasHLA_final_table.loc[arcasHLA_final_table['sample'] == sample_name, locus_name] = value_to_add

# arcasHLA_final_table.to_csv(
#     snakemake.output.arcasHLA, sep="\t", index=False, header=True
# )

# #hla-la
# hla_la_input = snakemake.input.hla_la

# #initialize the dataframe
# hla_la_final_table = pd.DataFrame(columns=('sample', 'A', 'B', 'C', 'DQA1','DQB1'))
# for index in range(len(hla_la_input)):
#     sample_name = os.path.basename(os.path.dirname(os.path.dirname(hla_la_input[index]))).split("_")[0]
#     # reading
#     data = pd.read_csv(hla_la_input[index], sep = "\t")
#     # print(selected)
#     A = []
#     B = []
#     C = []
#     DQA1 = []
#     DQB1 = []

#     selected1 = data.loc[(data['Locus'] == 'A') & (data['Chromosome'] == 1)]
#     A.append(selected1.loc[0,'Allele'])
#     selected2 = data.loc[(data['Locus'] == 'A') & (data['Chromosome'] == 2)]
#     selected2 = selected2.reset_index()
#     A.append(selected2.loc[0,'Allele'])
#     A = '/'.join(A)

#     selected1 = data.loc[(data['Locus'] == 'B') & (data['Chromosome'] == 1)]
#     selected1 = selected1.reset_index()
#     B.append(selected1.loc[0,'Allele'])
#     selected2 = data.loc[(data['Locus'] == 'B') & (data['Chromosome'] == 2)]
#     selected2 = selected2.reset_index()
#     B.append(selected2.loc[0,'Allele'])
#     B = '/'.join(B)

#     selected1 = data.loc[(data['Locus'] == 'C') & (data['Chromosome'] == 1)]
#     selected1 = selected1.reset_index()
#     C.append(selected1.loc[0,'Allele'])
#     selected2 = data.loc[(data['Locus'] == 'C') & (data['Chromosome'] == 2)]
#     selected2 = selected2.reset_index()
#     C.append(selected2.loc[0,'Allele'])
#     C = '/'.join(C)

#     selected1 = data.loc[(data['Locus'] == 'DQA1') & (data['Chromosome'] == 1)]
#     selected1 = selected1.reset_index()
#     DQA1.append(selected1.loc[0,'Allele'])
#     selected2 = data.loc[(data['Locus'] == 'DQA1') & (data['Chromosome'] == 2)]
#     selected2 = selected2.reset_index()
#     DQA1.append(selected2.loc[0,'Allele'])
#     DQA1 = '/'.join(DQA1)

#     selected1 = data.loc[(data['Locus'] == 'DQB1') & (data['Chromosome'] == 1)]
#     selected1 = selected1.reset_index()
#     DQB1.append(selected1.loc[0,'Allele'])
#     selected2 = data.loc[(data['Locus'] == 'DQB1') & (data['Chromosome'] == 2)]
#     selected2 = selected2.reset_index()
#     DQB1.append(selected2.loc[0,'Allele'])
#     DQB1 = '/'.join(DQB1)

#     new_row = {'sample': sample_name, 'A': A, 'B': B, 'C': C, 'DQA1': DQA1, 'DQB1': DQB1}
#     hla_la_final_table = hla_la_final_table.append(new_row, ignore_index=True)

# hla_la_final_table.to_csv(
#     snakemake.output.hla_la, sep="\t", index=False, header=True
# )