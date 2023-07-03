import pandas as pd
import os
import json
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #orthanq results
    orthanq_input = snakemake.input.orthanq

    #ground truth
    ground_truth = pd.read_csv(snakemake.input.ground_truth, sep = "\t")

    #initialize orthanq validation table
    orthanq_validation_table = pd.DataFrame(columns=('Locus', 'N', 'Orthanq - Call Rate', 'Orthanq - Accuracy'))

    #loop over loci and orthanq results
    loci = ['A', 'B', 'C', 'DQB1', 'DRB1']
    samples_collected = []
    for locus in loci:    
        #initialize collected hits
        collected = 0.0
        for (index, _) in enumerate(orthanq_input):
            #retrieve sample and locus names
            splitted = os.path.basename(orthanq_input[index]).split("_")
            sample_name = splitted[0]
            samples_collected.append(sample_name)
            locus_name = splitted[1].split(".")[0]
            print("sample name: ", sample_name)
            
            #take the best record
            # best_result = pd.read_csv(orthanq_input[index]).iloc[[0]] #best result, first row
            
            ##more than one best record
            results = pd.read_csv(orthanq_input[index])
            # results['odds'] = results['odds'].map('{:,.3f}'.format).astype(float)

            ##find the records that have the same density
            best_odds = [1, 1.0, 1.00]
            best_results = results[results.odds.isin(best_odds)]
            print("best results: ", best_results)
            
            #retrieve the predicted haplotypes
            filtered_cols = []
            filtered_cols = [col for col in results if col.startswith(locus)]
   
            #loop over best results
            for (i,result_row) in best_results.iterrows():

                #collect first two fields and cumulative fractions of haplotypes
                two_field_fractions = {}
                for col in filtered_cols:
                    col_splt = col.split(":")
                    two_field = col_splt[0] + ":" + col_splt[1]
                    if two_field in two_field_fractions:
                        two_field_fractions[two_field] += best_results[col][i]
                    else:
                        two_field_fractions[two_field] = best_results[col][i]
                print("two_field_fractions: ",two_field_fractions)

                #collect ground truth in values_in_truth, handle special cases e.g. 23:01/02/04
                values_in_truth = {}
                for (chr_index,chr_number) in enumerate([1, 2]):

                    #get ground truth values for the sample and locus
                    locus_in_truth = "HLA-" + locus + " " + str(chr_number)
                    value_in_truth = locus+"*"+ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_in_truth].array[0] ##array[0] to reach what's inside the Series
                    # print("values_in_truth: ",values_in_truth)

                    ##split in case of alleles that have 
                    #e.g. 23:01/02/04 in the truth
                    ##the sample that has multiple possible alleles, the maximum allele is being added to the cumulative fractions (1st condition below)
                    tmp_alleles = []
                    if "/" in value_in_truth:
                        first_split = value_in_truth.split("/")
                        first_field = first_split[0].split(":")[0]
                        tmp_alleles.append(first_split[0])
                        for splitted in first_split[1:]:
                            tmp_alleles.append(first_field + ":" +splitted)
                        values_in_truth["{0}".format(chr_index)] = tmp_alleles
                    elif value_in_truth.endswith("*"):  #also handle the cases where samples have '*' character in the end; remove the character.
                        value_in_truth=value_in_truth.rstrip("*")
                        values_in_truth["{0}".format(chr_index)] = [value_in_truth]
                    else: #the normal case
                        values_in_truth["{0}".format(chr_index)] = [value_in_truth]
                print("values_in_truth: ",values_in_truth)

                #collect the fractions of haplotypes matching with ground truth
                #loop over values in truth and stop when one matches an allele from the output of orthanq
                #prev_value is used to stop collecting the fraction from the identical allele match
                ##1-)the ground truth might have multiple possible alleles for single haplotypes (chromosomes)
                ##  in that case, collect all allele-fraction combinations to collected_fractions
                ##  and take account the allele that has the highest fraction
                
                prev_value = ""
                fractions = 0.0
                for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
                    for chr_value in chr_values:
                        collected_fractions = {}
                        if chr_value in two_field_fractions.keys() and chr_value != prev_value:
                            collected_fractions[chr_value]=two_field_fractions[chr_value]
                        if collected_fractions:
                            fractions += collected_fractions[max(collected_fractions)] ##get the maximum one
                            prev_value = max(collected_fractions)
                print("fractions: ",fractions)

                #score only haplotypes that make up more than 50%
                if fractions > 0.5:    
                    collected += 1.0
                    print("inner loop collected: ", collected)
                    break ##break as soon as the collected gets +1 
        print("collected: ",collected)

        #calculate accuracy
        accuracy = 100*collected/len(list(set(samples_collected)))

        #concatenate the locus-wise statistics to the validation table
        new_row = pd.DataFrame([[locus, len(list(set(samples_collected))), 1.00, accuracy]],
                        columns=['Locus', 'N', 'Orthanq - Call Rate', 'Orthanq - Accuracy'])
        orthanq_validation_table = pd.concat([orthanq_validation_table, new_row], ignore_index=True)

    #write the validation table
    orthanq_validation_table.to_csv(
        snakemake.output.validation, sep="\t", index=False, header=True
    )