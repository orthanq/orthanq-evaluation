##this script does the following:
# It parses and validates orthanq results compared to the ground truth (applicable diploid priors)
# Finally, create a table for accuracy values
# In addition, this script also writes a table for TP and FP predictions of Orthanq later to be used for plotting purposes.
# Side note: validation from final report of orthanq inside validation.py cannot be performed because we also need fraction info for TP-FP table, we do both that and validation at once here.
# For accuracy, we count, for each tool: TP / (TP+FP).
# For call rate, we count, for each tool: (TP + FP) / (TP + FP + 'no call')

import pandas as pd
import os
import json
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #orthanq results
    orthanq_input = snakemake.input.orthanq

    # read sample sheets with truth values for evaluated and low coverage samples
    ground_truth_evaluated = pd.read_csv(snakemake.input.ground_truth_evaluated, sep = "\t")

    def truth_for_sample(locus, values_in_truth, ground_truth, sample_name):
        for (chr_index,chr_number) in enumerate([1, 2]):
            locus_in_truth= "HLA-" + locus + " " + str(chr_number)
            value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_in_truth].array[0]

            ##split in case of alles that have e.g. 23:01/02/04 in the truth
            tmp_alleles = []
            if "/" in value_in_truth:
                first_split = value_in_truth.split("/")
                first_field = first_split[0].split(":")[0]
                tmp_alleles.append(locus + "*" + first_split[0])
                for splitted in first_split[1:]:
                    tmp_alleles.append(locus + "*" + first_field + ":" +splitted)
                values_in_truth["{0}".format(chr_index)] = tmp_alleles
            elif value_in_truth.endswith("*"):
                value_in_truth=value_in_truth.rstrip("*")
                values_in_truth["{0}".format(chr_index)] = [locus + "*" + value_in_truth]
            else:
                values_in_truth["{0}".format(chr_index)] = [locus + "*" + value_in_truth]
        return values_in_truth

    def validate_orthanq(orthanq_input, threshold_check, threshold_density, threshold_number_of_haplotypes, validation_table, orthanq_tp_fp_table, ground_truth, sample_list, orthanq_tp_fp_A, orthanq_tp_fp_B, orthanq_tp_fp_C, orthanq_tp_fp_DQB1):
        #loop over loci and orthanq results
        loci = ['A', 'B', 'C', 'DQB1']
        for locus in loci:   
            samples_collected = [] 
            #initialize collected hits
            collected = 0.0
            samples_called = 0 #some samples in orthanq have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, _) in enumerate(orthanq_input):
                #collect ground truth in values_in_truth, handle special cases e.g. 23:01/02/04

                #retrieve sample and locus names
                splitted = os.path.basename(orthanq_input[index]).split("_")
                sample_name = splitted[0]

                if sample_name in sample_list and sample_name != "D1":
                    values_in_truth = {}
                    values_in_truth = truth_for_sample(locus, values_in_truth, ground_truth, sample_name)

                    samples_collected.append(sample_name) #keep track of samples that are evaluated
                    locus_in_orthanq = splitted[1].split(".")[0]
                    print("sample name: ", sample_name)
                                     
                    ##more than one best record
                    results = pd.read_csv(orthanq_input[index])
                    # results['odds'] = results['odds'].map('{:,.3f}'.format).astype(float)

                    ##find the records that have the same density
                    best_odds = [1, 1.0, 1.00]
                    best_results = results[results.odds.isin(best_odds)]

                    #fill in the tp fp table

                    #collecting values_in_truth_w_separators is necessary as some truth values have multiple values (required for the tp_fp table)
                    values_in_truth_w_separators = []
                    for key_chr,values in values_in_truth.items():
                        values_with_sep = []
                        for value in values:
                            values_with_sep.append(value)
                        values_with_sep = ';'.join(values_with_sep)
                        values_in_truth_w_separators.append(values_with_sep)
                    print(values_in_truth_w_separators)

                    if locus == "A" and sample_name not in orthanq_tp_fp_A['sample'].tolist():
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_w_separators), '']],
                            columns=['sample', 'ground', 'orthanq evaluation'])
                        orthanq_tp_fp_A = pd.concat([orthanq_tp_fp_A, new_row_tp], ignore_index=True)
                    if locus == "B" and sample_name not in orthanq_tp_fp_B['sample'].tolist():
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_w_separators), '']],
                            columns=['sample', 'ground', 'orthanq evaluation'])
                        orthanq_tp_fp_B = pd.concat([orthanq_tp_fp_B, new_row_tp], ignore_index=True)
                    if locus == "C" and sample_name not in orthanq_tp_fp_C['sample'].tolist():
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_w_separators), '']],
                            columns=['sample', 'ground', 'orthanq evaluation'])
                        orthanq_tp_fp_C = pd.concat([orthanq_tp_fp_C, new_row_tp], ignore_index=True)
                    if locus == "DQB1" and sample_name not in orthanq_tp_fp_DQB1['sample'].tolist():
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_w_separators), '']],
                            columns=['sample', 'ground', 'orthanq evaluation'])
                        orthanq_tp_fp_DQB1 = pd.concat([orthanq_tp_fp_DQB1, new_row_tp], ignore_index=True)

                    if not best_results.empty: #orthanq has no predictions for some samples

                        #retrieve the predicted haplotypes
                        haplotypes = []
                        haplotypes = [col for col in results if col.startswith(locus)]
            
                        #enter the best density to the orthanq_tp_fp_table
                        best_density = results['density'][0]

                        #initialize the orthanq tp fp table 
                        if not (((orthanq_tp_fp_table["Sample"] == sample_name) & (orthanq_tp_fp_table["Locus"] == locus_in_orthanq)).any()):
                            row_to_add = pd.DataFrame([[sample_name, locus_in_orthanq, "FP", best_density]],
                                columns=['Sample', 'Locus', 'Prediction', 'Best_Density'])
                            orthanq_tp_fp_table = pd.concat([orthanq_tp_fp_table, row_to_add], ignore_index=True)

                        #only include predictions having density above determined threshold
                        sum_of_densities = 0
                        sum_of_densities = best_results["density"].sum()
                        print("sum_of_densities: ", sum_of_densities)
                        print(best_results)

                        #loop over best results (rows that share the same best odds score, i.e. same density and solutions but different haplotypes)
                        counter_fp = 0
                        for (i,row) in best_results.iterrows():

                            #collect first two fields
                            predictions = []
                            for h in haplotypes:
                                h_splt = h.split(":")
                                first_field = h_splt[0].split("*")[1]
                                two_field = locus_in_orthanq + "*" + first_field + ":" + h_splt[1]
                                if row[h] == 0.5:
                                    predictions.append(two_field)
                                if row[h] == 1.0:
                                    predictions.append(two_field)
                                    predictions.append(two_field)
 
                            print("predictions: ", predictions)
                            print("values_in_truth: ", values_in_truth)

                            allele_match_check = 0 #should be 2 to pass
                            values_in_truth_clone = values_in_truth.copy() #an explicit copy is necessary as we need to keep values in truth for the next iteration
                            for p in predictions:
                                for c_key, c_values in values_in_truth_clone.copy().items():
                                    print(p)
                                    print(c_values)
                                    if p in c_values:
                                        allele_match_check += 1
                                        values_in_truth_clone.pop(c_key)
                                        break #necessary to avoid checking for the other chr just after the match here

                            print("allele_match_check: ", allele_match_check)

                            #threshold check
                            if threshold_check:
                                threshold_haplotypes_bool = best_results.shape[0] <= threshold_number_of_haplotypes
                                threshold_density_bool = sum_of_densities > threshold_density
                            else:
                                threshold_haplotypes_bool = True
                                threshold_density_bool = True
                            print("threshold_number_of_haplotypes: ", threshold_number_of_haplotypes)
                            print("threshold_density: ", threshold_density)
                            print("threshold_haplotypes: ",threshold_haplotypes_bool)
                            print("threshold_density_bool: ",threshold_density_bool)

                            if locus == locus_in_orthanq:
                                if allele_match_check == 2 and (threshold_density_bool and threshold_haplotypes_bool): #this number is configurable:: 
                                    print("first")
                                    #count the number of samples that are called
                                    if counter_fp == 0:
                                        samples_called += 1

                                    #count TP for accuracy calculation
                                    collected += 1.0
                                    print("inner loop collected: ", collected)

                                    #fill up the orthanq_tp_fp_table with TP
                                    orthanq_tp_fp_table.loc[(orthanq_tp_fp_table.Sample == sample_name) & (orthanq_tp_fp_table.Locus == locus_in_orthanq) & (orthanq_tp_fp_table.Best_Density == best_density),'Prediction'] = "TP"
                                    
                                    #fill in the table for tp fp info
                                    if locus_in_orthanq == "A":
                                        orthanq_tp_fp_A.loc[orthanq_tp_fp_A['sample'] == sample_name, 'orthanq evaluation'] = 'TP'
                                    if locus_in_orthanq == "B":
                                        orthanq_tp_fp_B.loc[orthanq_tp_fp_B['sample'] == sample_name, 'orthanq evaluation'] = 'TP'
                                    if locus_in_orthanq == "C":
                                        orthanq_tp_fp_C.loc[orthanq_tp_fp_C['sample'] == sample_name, 'orthanq evaluation'] = 'TP'
                                    if locus_in_orthanq == "DQB1":
                                        orthanq_tp_fp_DQB1.loc[orthanq_tp_fp_DQB1['sample'] == sample_name, 'orthanq evaluation'] = 'TP'
                                    break ##break as soon as the collected gets +1 
                                
                                elif allele_match_check == 2  and (not (threshold_density_bool and threshold_haplotypes_bool)): #this number is configurable:: 
                                    print("second")
                                    #fractions exceed the threshold but but it's an considered uncalled, locus==locus_in_orthanq check avoids multiple assignments
                                    for index,row in orthanq_tp_fp_A.iterrows():
                                        if locus==locus_in_orthanq and locus=="A" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_A.loc[index, "orthanq evaluation"] = "true & uncalled"
                                    for index,row in orthanq_tp_fp_B.iterrows():
                                        if locus==locus_in_orthanq and locus=="B" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_B.loc[index, "orthanq evaluation"] = "true & uncalled"
                                    for index,row in orthanq_tp_fp_C.iterrows():
                                        if locus==locus_in_orthanq and locus=="C" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_C.loc[index, "orthanq evaluation"] = "true & uncalled"
                                    for index,row in orthanq_tp_fp_DQB1.iterrows():
                                        if locus==locus_in_orthanq and locus=="DQB1" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_DQB1.loc[index, "orthanq evaluation"] = "true & uncalled"
                                    #also, make the general tp_fp table true & uncalled
                                    orthanq_tp_fp_table.loc[(orthanq_tp_fp_table.Sample == sample_name) & (orthanq_tp_fp_table.Locus == locus_in_orthanq) & (orthanq_tp_fp_table.Best_Density == best_density),'Prediction'] = "true & uncalled"
                                    break

                                elif (allele_match_check != 2 ) and (threshold_density_bool and threshold_haplotypes_bool): #this number is configurable:: 
                                    print("third")
                                    if counter_fp == 0: #this counter is to avoid same density solutions to add to samples_called all the time
                                        #count the number of samples that are called
                                        samples_called += 1
                                        counter_fp += 1
                                    if locus_in_orthanq == "A":
                                        orthanq_tp_fp_A.loc[orthanq_tp_fp_A['sample'] == sample_name, 'orthanq evaluation'] = 'FP'
                                    if locus_in_orthanq == "B":
                                        orthanq_tp_fp_B.loc[orthanq_tp_fp_B['sample'] == sample_name, 'orthanq evaluation'] = 'FP'
                                    if locus_in_orthanq == "C":
                                        orthanq_tp_fp_C.loc[orthanq_tp_fp_C['sample'] == sample_name, 'orthanq evaluation'] = 'FP'
                                    if locus_in_orthanq == "DQB1":
                                        orthanq_tp_fp_DQB1.loc[orthanq_tp_fp_DQB1['sample'] == sample_name, 'orthanq evaluation'] = 'FP'
                                    #fill up the orthanq_tp_fp_table with FP
                                    orthanq_tp_fp_table.loc[(orthanq_tp_fp_table.Sample == sample_name) & (orthanq_tp_fp_table.Locus == locus_in_orthanq) & (orthanq_tp_fp_table.Best_Density == best_density),'Prediction'] = "FP"
                            
                                elif (allele_match_check != 2 ) and (not (threshold_density_bool and threshold_haplotypes_bool)): #this number is configurable:: 
                                    print("fourth")
                                    #if sum of densities remain under the threshold, then make it a no call, locus==locus_in_orthanq check avoids multiple assignments
                                    for index,row in orthanq_tp_fp_A.iterrows():
                                        if locus==locus_in_orthanq and locus=="A" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_A.loc[index, "orthanq evaluation"] = "false & uncalled"
                                    for index,row in orthanq_tp_fp_B.iterrows():
                                        if locus==locus_in_orthanq and locus=="B" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_B.loc[index, "orthanq evaluation"] = "false & uncalled"
                                    for index,row in orthanq_tp_fp_C.iterrows():
                                        if locus==locus_in_orthanq and locus=="C" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_C.loc[index, "orthanq evaluation"] = "false & uncalled"
                                    for index,row in orthanq_tp_fp_DQB1.iterrows():
                                        if locus==locus_in_orthanq and locus=="DQB1" and row["orthanq evaluation"] == "":
                                            orthanq_tp_fp_DQB1.loc[index, "orthanq evaluation"] = "false & uncalled"
                                    #ill up the orthanq_tp_fp_table with false & uncalled
                                    orthanq_tp_fp_table.loc[(orthanq_tp_fp_table.Sample == sample_name) & (orthanq_tp_fp_table.Locus == locus_in_orthanq) & (orthanq_tp_fp_table.Best_Density == best_density),'Prediction'] = "false & uncalled"
                    else:
                        print("empty df sample name", sample_name)
                        #if the df is empty for a result then, then make it a no call, locus==locus_in_orthanq check avoids multiple assignments
                        for index,row in orthanq_tp_fp_A.iterrows():
                            if locus==locus_in_orthanq and locus=="A" and row["orthanq evaluation"] == "":
                                orthanq_tp_fp_A.loc[index, "orthanq evaluation"] = "no call"
                        for index,row in orthanq_tp_fp_B.iterrows():
                            if locus==locus_in_orthanq and locus=="B" and row["orthanq evaluation"] == "":
                                orthanq_tp_fp_B.loc[index, "orthanq evaluation"] = "no call"
                        for index,row in orthanq_tp_fp_C.iterrows():
                            if locus==locus_in_orthanq and locus=="C" and row["orthanq evaluation"] == "":
                                orthanq_tp_fp_C.loc[index, "orthanq evaluation"] = "no call"
                        for index,row in orthanq_tp_fp_DQB1.iterrows():
                            if locus==locus_in_orthanq and locus=="DQB1" and row["orthanq evaluation"] == "":
                                orthanq_tp_fp_DQB1.loc[index, "orthanq evaluation"] = "no call"
                        #also, make the general tp_fp table no call
                        orthanq_tp_fp_table.loc[(orthanq_tp_fp_table.Sample == sample_name) & (orthanq_tp_fp_table.Locus == locus_in_orthanq) & (orthanq_tp_fp_table.Best_Density == best_density),'Prediction'] = "no call"
            print("collected: ",collected)

            #calculate accuracy and call rate

            n_samples_collected = len(list(set(samples_collected)))
            print("n_samples_collected: ",n_samples_collected)
            # to be fair to all tools, the denominator should be the number of samples that the tool results in a prediction.
            if samples_called == 0:
                accuracy = "NA"
            else:
                accuracy = 100*collected/samples_called
            print("samples_called: ",samples_called)
            call_rate = samples_called/n_samples_collected

            # add the numbers in parenthesis for call rate and accuracy
            if accuracy != "NA":
                accuracy = "{:.2f}".format(accuracy) + " (" + str(int(collected)) + "/" + str(samples_called) + ")"
            call_rate = "{:.2f}".format(call_rate) + " (" + str(samples_called) + "/" + str(n_samples_collected) + ")"

            #concatenate the locus-wise statistics to the validation table
            new_row = pd.DataFrame([[locus, n_samples_collected, call_rate, accuracy]],
                            columns=['Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'])

            validation_table = pd.concat([validation_table, new_row], ignore_index=True)
        
        #if the evaluation field is empty, make all empty to FP
        for index,row in orthanq_tp_fp_A.iterrows():
            if row["orthanq evaluation"] == "":
                orthanq_tp_fp_A.loc[index, "orthanq evaluation"] = "FP"
        for index,row in orthanq_tp_fp_B.iterrows():
            if row["orthanq evaluation"] == "":
                orthanq_tp_fp_B.loc[index, "orthanq evaluation"] = "FP"
        for index,row in orthanq_tp_fp_C.iterrows():
            if row["orthanq evaluation"] == "":
                orthanq_tp_fp_C.loc[index, "orthanq evaluation"] = "FP"
        for index,row in orthanq_tp_fp_DQB1.iterrows():
            if row["orthanq evaluation"] == "":
                orthanq_tp_fp_DQB1.loc[index, "orthanq evaluation"] = "FP"
        print("orthanq_tp_fp_B")
        print(orthanq_tp_fp_B)
        return validation_table, orthanq_tp_fp_table, orthanq_tp_fp_A, orthanq_tp_fp_B, orthanq_tp_fp_C, orthanq_tp_fp_DQB1
    
    #check if the predicted alleles are in the allele freq table and below or above the freq filter
    def check_alleles_in_database(predictions_table, allele_freq_table, locus_name):
        predictions = predictions_table[locus_name].to_list()
        for i in range(len(predictions)):
            if predictions[i] != '':
            #in case there are multiple solutions, split by ,
                splitted = predictions[i].split(",")
                for allele_pair in splitted:
                    #check for two field resolution of orthanq prediction, three field doesn't make sense
                    first_allele_three_field = allele_pair.split("/")[0]
                    first_allele_splitted = first_allele_three_field.split(":")
                    first_allele_two_field = first_allele_splitted[0] + ":" + first_allele_splitted[1]

                    sec_allele_three_field = allele_pair.split("/")[1]
                    sec_allele_splitted = sec_allele_three_field.split(":")
                    sec_allele_two_field = sec_allele_splitted[0] + ":" + sec_allele_splitted[1]
                
                #one of the alleles in allele pair are labeled if one of them meets the criteria.
                for allele, freq in zip(allele_freq_table["var"], allele_freq_table["frequency"]):
                    if freq > 0.05:
                        if (first_allele_three_field == allele or first_allele_two_field == allele) or (sec_allele_three_field == allele or sec_allele_two_field == allele):
                            break

                # the check for the first two fields of a three field or four field containing record, only if first two fields don't match is not necessary.
                # because the info in the db belongs to a subset of the two field allele, so that wouldn't be a valid comparison.
                else: # if the break is not entered above, then else is executed in this for/else construct.
                    predictions_table.loc[i, "orthanq evaluation"] = "allele not considered in the truth set" #it's enough if one the alleles don't break the if above
        return predictions_table

    #initialize the table to output TP and FP samples
    orthanq_validation_table_all= pd.DataFrame(columns=('Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'))
    orthanq_tp_fp_table_all = pd.DataFrame(columns=('Sample', 'Locus', 'Prediction', 'Best_Density'))

    #initialize locus-wise tp fp table for orthanq
    orthanq_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
    orthanq_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
    orthanq_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
    orthanq_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))

    #validate all samples
    sample_list_all = ground_truth_evaluated["Run Accession"].to_list()
    threshold_density_in_paper = 0.7
    threshold_n_haplotypes_in_paper = 5
    threshold_haplotype_check = True

    orthanq_validation_results_all = validate_orthanq(orthanq_input, threshold_haplotype_check, threshold_density_in_paper, threshold_n_haplotypes_in_paper, orthanq_validation_table_all, orthanq_tp_fp_table_all, ground_truth_evaluated, sample_list_all, orthanq_tp_fp_A, orthanq_tp_fp_B, orthanq_tp_fp_C, orthanq_tp_fp_DQB1)

    final_orthanq_all = orthanq_validation_results_all[0]
    final_tp_fp_table_all = orthanq_validation_results_all[1]

    #merge orthanq predictions from the final table of orthanq, here it wasn't possible because we break the loop during validation so we don't get to see the other solutions.
    orthanq_parsed_table = pd.read_csv(snakemake.input.orthanq_final_table, sep="\t", keep_default_na=False)
    merged_with_orthanq_predictions_A = pd.merge(orthanq_validation_results_all[2], orthanq_parsed_table[["sample", "A"]], how='left', on=['sample'])
    merged_with_orthanq_predictions_B = pd.merge(orthanq_validation_results_all[3], orthanq_parsed_table[["sample", "B"]], how='left', on=['sample'])
    merged_with_orthanq_predictions_C = pd.merge(orthanq_validation_results_all[4], orthanq_parsed_table[["sample", "C"]], how='left', on=['sample'])
    merged_with_orthanq_predictions_DQB1 = pd.merge(orthanq_validation_results_all[5], orthanq_parsed_table[["sample", "DQB1"]], how='left', on=['sample'])

    #finally, for tp fp tables for each locus, check if they contain alleles that are below the criteria according to the Abi-Rached 2018 paper.
    #if that's the case label them as "allele not considered in the truth set".
    allele_freqs = pd.read_csv(snakemake.input.allele_freqs)
    final_orthanq_tp_fp_A = check_alleles_in_database(merged_with_orthanq_predictions_A, allele_freqs, "A") #2 for A, 3 for B, 4 for C, 5 for DQB1orthanq_validation_results_all[5]
    final_orthanq_tp_fp_B = check_alleles_in_database(merged_with_orthanq_predictions_B, allele_freqs, "B")
    final_orthanq_tp_fp_C = check_alleles_in_database(merged_with_orthanq_predictions_C, allele_freqs, "C")
    final_orthanq_tp_fp_DQB1 = check_alleles_in_database(merged_with_orthanq_predictions_DQB1, allele_freqs, "DQB1")

    #rename the locus columns to "orthanq"
    final_orthanq_tp_fp_A = final_orthanq_tp_fp_A.rename(columns={"A": "orthanq"})
    final_orthanq_tp_fp_B = final_orthanq_tp_fp_B.rename(columns={"B": "orthanq"})
    final_orthanq_tp_fp_C = final_orthanq_tp_fp_C.rename(columns={"C": "orthanq"})
    final_orthanq_tp_fp_DQB1 = final_orthanq_tp_fp_DQB1.rename(columns={"DQB1": "orthanq"})

    #if the "orthanq" field is empty, make "orthanq evaluation" to "no call"
    for index,row in final_orthanq_tp_fp_A.iterrows():
        if row["orthanq"] == "":
            final_orthanq_tp_fp_A.loc[index, "orthanq evaluation"] = "no call"
    for index,row in final_orthanq_tp_fp_B.iterrows():
        if row["orthanq"] == "":
            final_orthanq_tp_fp_B.loc[index, "orthanq evaluation"] = "no call"
    for index,row in final_orthanq_tp_fp_C.iterrows():
        if row["orthanq"] == "":
            final_orthanq_tp_fp_C.loc[index, "orthanq evaluation"] = "no call"
    for index,row in final_orthanq_tp_fp_DQB1.iterrows():
        if row["orthanq"] == "":
            final_orthanq_tp_fp_DQB1.loc[index, "orthanq evaluation"] = "no call"

    #reorder columns
    final_orthanq_tp_fp_A = final_orthanq_tp_fp_A[['sample', 'ground', 'orthanq', 'orthanq evaluation']]
    final_orthanq_tp_fp_B = final_orthanq_tp_fp_B[['sample', 'ground', 'orthanq', 'orthanq evaluation']]
    final_orthanq_tp_fp_C = final_orthanq_tp_fp_C[['sample', 'ground', 'orthanq', 'orthanq evaluation']]
    final_orthanq_tp_fp_DQB1 = final_orthanq_tp_fp_DQB1[['sample', 'ground', 'orthanq', 'orthanq evaluation']]

    #write the validation table

    final_orthanq_all.to_csv(
        snakemake.output.validation_all, sep="\t", index=False, header=True
    )

    #write the orthanq_tp_fp_table table for all
    final_tp_fp_table_all.to_csv(
        snakemake.output.tp_fp_table, sep="\t", index=False, header=True
    )

    #write locus wise tables
    final_orthanq_tp_fp_A.to_csv(
        snakemake.output.orthanq_A_tp_fp, sep="\t", index=False, header=True
    )

    final_orthanq_tp_fp_B.to_csv(
        snakemake.output.orthanq_B_tp_fp, sep="\t", index=False, header=True
    )

    final_orthanq_tp_fp_C.to_csv(
        snakemake.output.orthanq_C_tp_fp, sep="\t", index=False, header=True
    )

    final_orthanq_tp_fp_DQB1.to_csv(
        snakemake.output.orthanq_DQB1_tp_fp, sep="\t", index=False, header=True
    )

    #### Addition: write accuracy-callrate-threshold table to file (request from first review) ###
    
    #initialize dataframe
    threshold_results = {'threshold_density': [], 'threshold_haplotypes': [] ,'locus': [], 'call_rate': [], 'accuracy': []}

    #try varying thresholds from 0.0 to 1.0
    thresholds_density = (x * 0.1 for x in range(0, 11))
    # thresholds_density = (x * 0.1 for x in range(0, 2)) #for testing

    thresholds_n_haplotypes = [5, 6, 7, 8, 9, 10]
    # thresholds_n_haplotypes = [5, 6] #for testing

    print("thresholds:", thresholds_density)
    for t_density in thresholds_density:
        for t_haplotype in thresholds_n_haplotypes:
            print("t_density: ", t_density)
            print("t_haplotype: ", t_haplotype)

            #first, initialize required dfs
            orthanq_validation_table = pd.DataFrame(columns=('Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'))
            orthanq_tp_fp_table_all = pd.DataFrame(columns=('Sample', 'Locus', 'Prediction', 'Best_Density'))
            orthanq_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
            orthanq_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
            orthanq_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
            orthanq_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))

            #get validation results for the given threshold
            threshold_haplotype_check = True
            validation_all_tables = validate_orthanq(orthanq_input, threshold_haplotype_check, t_density, t_haplotype, orthanq_validation_table_all, orthanq_tp_fp_table_all, ground_truth_evaluated, sample_list_all, orthanq_tp_fp_A, orthanq_tp_fp_B, orthanq_tp_fp_C, orthanq_tp_fp_DQB1)

            #get validation table with accuracy and call rate
            validation_results = validation_all_tables[0]

            #loop over loci and append values to the results table
            for row in validation_results.itertuples():
                threshold_results['threshold_density'].append(t_density)
                threshold_results['threshold_haplotypes'].append(t_haplotype)
                threshold_results['locus'].append(row.Locus)
                threshold_results['call_rate'].append(row.Orthanq_Call_Rate)
                threshold_results['accuracy'].append(row.Orthanq_Accuracy)

    print("threshold_results: ", threshold_results)

    ##lastly, add the case with no threshold at all, values for threshold_density adn threshold_haplotypes do not matter as long as they're boyh False.
    orthanq_validation_table = pd.DataFrame(columns=('Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'))
    orthanq_tp_fp_table_all = pd.DataFrame(columns=('Sample', 'Locus', 'Prediction', 'Best_Density'))
    orthanq_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
    orthanq_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
    orthanq_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
    orthanq_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'orthanq evaluation'))
    threshold_haplotype_check = False
    validation_all_tables = validate_orthanq(orthanq_input, threshold_haplotype_check, 0.0, 0, orthanq_validation_table_all, orthanq_tp_fp_table_all, ground_truth_evaluated, sample_list_all, orthanq_tp_fp_A, orthanq_tp_fp_B, orthanq_tp_fp_C, orthanq_tp_fp_DQB1)
    validation_results = validation_all_tables[0]
    for row in validation_results.itertuples():
        threshold_results['threshold_density'].append("no threshold")
        threshold_results['threshold_haplotypes'].append("no threshold")
        threshold_results['locus'].append(row.Locus)
        threshold_results['call_rate'].append(row.Orthanq_Call_Rate)
        threshold_results['accuracy'].append(row.Orthanq_Accuracy)

    threshold_results_df = pd.DataFrame(threshold_results)

    threshold_results_df.to_csv(
        snakemake.output.threshold_results, sep="\t", index=False, header=True
    )