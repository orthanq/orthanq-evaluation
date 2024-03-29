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

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    for index in range(len(orthanq_input)):
        if "D1_S1_L001" in orthanq_input[index]:
            splitted = os.path.basename(orthanq_input[index]).split("_")
            sample_name = splitted[0] + splitted[1] + splitted[2] 
            locus_name = splitted[3].split(".")[0]
            #rename sample name for giab to the accesssion id
            sample_name = "SRR2962669"
        else:
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

    #arcasHLA
    arcasHLA_input = snakemake.input.arcasHLA

    #initialize the dataframe
    arcasHLA_final_table = pd.DataFrame(columns=('sample', 'A', 'B', 'C', 'DQA1','DQB1'))
    for index in range(len(arcasHLA_input)):
        if "D1_S1_L001" in arcasHLA_input[index]:
            splitted = os.path.basename(os.path.dirname(arcasHLA_input[index])).split("_")
            sample_name = splitted[0] + splitted[1] + splitted[2] 
            locus_name = splitted[3].split(".")[0]
            #rename sample name for giab to the accesssion id
            sample_name = "SRR2962669"
        else:
            splitted = os.path.basename(os.path.dirname(arcasHLA_input[index])).split("_")
            sample_name = splitted[0]
            locus_name = splitted[1].split(".")[0]
        if not sample_name in arcasHLA_final_table['sample'].tolist():
            new_row = pd.DataFrame([[sample_name, '', '', '', '', '']],
                    columns=['sample', 'A', 'B', 'C', 'DQA1', 'DQB1'])
            arcasHLA_final_table = pd.concat([arcasHLA_final_table, new_row], ignore_index=True)

        # Opening JSON file
        f = open(arcasHLA_input[index])
        
        # returns JSON object as 
        # a dictionary
        data = json.load(f)
        
        if data != {}: #arcasHLA may not produce predictions for some loci
            filtered_cols = data[locus_name]
            value_to_add = []
            for col in filtered_cols:
                value_to_add.append(col)
            value_to_add = '/'.join(value_to_add)
        else:
            value_to_add = None
        arcasHLA_final_table.loc[arcasHLA_final_table['sample'] == sample_name, locus_name] = value_to_add

    arcasHLA_final_table.to_csv(
        snakemake.output.arcasHLA, sep="\t", index=False, header=True
    )

    #hla-la
    hla_la_input = snakemake.input.hla_la

    #initialize the dataframe
    hla_la_final_table = pd.DataFrame(columns=('sample', 'A', 'B', 'C', 'DQA1','DQB1'))
    for index in range(len(hla_la_input)):
        if "D1_S1_L001" in hla_la_input[index]:
            splitted = os.path.basename(os.path.dirname(os.path.dirname(hla_la_input[index])))
            sample_name = splitted[0] + splitted[1] + splitted[2]
            #rename sample name for giab to the accesssion id
            sample_name = "SRR2962669"
        else:
            sample_name = os.path.basename(os.path.dirname(os.path.dirname(hla_la_input[index]))).split("_")[0]
        # reading
        data = pd.read_csv(hla_la_input[index], sep = "\t")
        A = []
        B = []
        C = []
        DQA1 = []
        DQB1 = []

        selected1 = data.loc[(data['Locus'] == 'A') & (data['Chromosome'] == 1)]
        A.append(selected1.loc[0,'Allele'])
        selected2 = data.loc[(data['Locus'] == 'A') & (data['Chromosome'] == 2)]
        selected2 = selected2.reset_index()
        A.append(selected2.loc[0,'Allele'])
        A = '/'.join(A)

        selected1 = data.loc[(data['Locus'] == 'B') & (data['Chromosome'] == 1)]
        selected1 = selected1.reset_index()
        B.append(selected1.loc[0,'Allele'])
        selected2 = data.loc[(data['Locus'] == 'B') & (data['Chromosome'] == 2)]
        selected2 = selected2.reset_index()
        B.append(selected2.loc[0,'Allele'])
        B = '/'.join(B)

        selected1 = data.loc[(data['Locus'] == 'C') & (data['Chromosome'] == 1)]
        selected1 = selected1.reset_index()
        C.append(selected1.loc[0,'Allele'])
        selected2 = data.loc[(data['Locus'] == 'C') & (data['Chromosome'] == 2)]
        selected2 = selected2.reset_index()
        C.append(selected2.loc[0,'Allele'])
        C = '/'.join(C)

        selected1 = data.loc[(data['Locus'] == 'DQA1') & (data['Chromosome'] == 1)]
        selected1 = selected1.reset_index()
        DQA1.append(selected1.loc[0,'Allele'])
        selected2 = data.loc[(data['Locus'] == 'DQA1') & (data['Chromosome'] == 2)]
        selected2 = selected2.reset_index()
        DQA1.append(selected2.loc[0,'Allele'])
        DQA1 = '/'.join(DQA1)

        selected1 = data.loc[(data['Locus'] == 'DQB1') & (data['Chromosome'] == 1)]
        selected1 = selected1.reset_index()
        DQB1.append(selected1.loc[0,'Allele'])
        selected2 = data.loc[(data['Locus'] == 'DQB1') & (data['Chromosome'] == 2)]
        selected2 = selected2.reset_index()
        DQB1.append(selected2.loc[0,'Allele'])
        DQB1 = '/'.join(DQB1)

        new_row = pd.DataFrame([[sample_name, A, B, C, DQA1, DQB1]],
                    columns=['sample', 'A', 'B', 'C', 'DQA1', 'DQB1'])
        hla_la_final_table = pd.concat([hla_la_final_table, new_row], ignore_index=True)

    hla_la_final_table.to_csv(
        snakemake.output.hla_la, sep="\t", index=False, header=True
    )

    #optitype
    optitype_input = snakemake.input.optitype   

    #initialize the dataframe
    optitype_final_table = pd.DataFrame(columns=('sample', 'A', 'B', 'C'))
    for index in range(len(optitype_input)):
        if "D1_S1_L001" in optitype_input[index]:
            splitted = os.path.basename(optitype_input[index])
            sample_name = splitted[0] + splitted[1] + splitted[2]
            #rename sample name for giab to the accesssion id
            sample_name = "SRR2962669"
        else:
            sample_name = os.path.basename(optitype_input[index]).split("_")[0]
        #read
        data = pd.read_csv(optitype_input[index], sep = "\t")

        A = []
        B = []
        C = []

        if not data['A1'].isnull().values.any(): #optitype does not may not output predictions for some samples and loci
            A.append(data.loc[0,'A1'])  
        if not data['A2'].isnull().values.any():
            A.append(data.loc[0,'A2'])
        A = '/'.join(A)

        if not data['B1'].isnull().values.any():
            B.append(data.loc[0,'B1'])
        if not data['B2'].isnull().values.any():
            B.append(data.loc[0,'B2'])
        B = '/'.join(B)

        if not data['C1'].isnull().values.any():
            C.append(data.loc[0,'C1'])
        if not data['C2'].isnull().values.any():
            C.append(data.loc[0,'C2'])
        C = '/'.join(C)

        new_row = pd.DataFrame([[sample_name, A, B, C]],
                    columns=['sample', 'A', 'B', 'C'])
        optitype_final_table = pd.concat([optitype_final_table, new_row], ignore_index=True)

    optitype_final_table.to_csv(
        snakemake.output.optitype, sep="\t", index=False, header=True
    )
