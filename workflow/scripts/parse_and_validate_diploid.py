##this script does the following:
# It parses and validates orthanq results compared to the ground truth (applicable diploid priors)
# Finally, create a table for accuracy values
# In addition, this script also writes a table for TP and FP predictions of Orthanq later to be used for plotting purposes.
# Side note: validation from final report of orthanq inside validation.py cannot be performed because we also need fraction info for TP-FP table, we do both that and validation at once here.

import pandas as pd
import os
import json
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #orthanq results
    orthanq_input = snakemake.input.orthanq
    print(orthanq_input)
    #ground truth
    # ground_truth = pd.read_csv(snakemake.input.ground_truth, sep = "\t")

    # read sample sheets with truth values for evaluated and low coverage samples
    samples_evaluated = pd.read_csv(snakemake.input.samples_evaluated, sep = "\t")
    samples_low_coverage = pd.read_csv(snakemake.input.samples_low_coverage, sep = "\t")
    samples_all = pd.concat([samples_evaluated, samples_low_coverage], ignore_index=True)

    def validate_orthanq(orthanq_input, validation_table, orthanq_tp_fp_table, ground_truth, sample_list):
        #loop over loci and orthanq results
        loci = ['A', 'B', 'C', 'DQB1']
        for locus in loci:   
            samples_collected = [] 
            #initialize collected hits
            collected = 0.0
            samples_called = 0 #some samples in orthanq have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, _) in enumerate(orthanq_input):

                #retrieve sample and locus names
                splitted = os.path.basename(orthanq_input[index]).split("_")
                sample_name = splitted[0]
                if sample_name in sample_list:
                    samples_collected.append(sample_name) #keep track of samples that are evaluated
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

                    if not best_results.empty: #orthanq has no predictions for some samples

                        #retrieve the predicted haplotypes
                        haplotypes = []
                        haplotypes = [col for col in results if col.startswith(locus)]
            
                        #enter the best density to the orthanq_tp_fp_table
                        best_density = results['density'][0]

                        #fill up the orthanq_tp_fp_table with 0 for FP as default
                        if not (((orthanq_tp_fp_table["Sample"] == sample_name) & (orthanq_tp_fp_table["Locus"] == locus_name)).any()):
                            row_to_add = pd.DataFrame([[sample_name, locus_name, 0, best_density]],
                                columns=['Sample', 'Locus', 'TP', 'Best_Density'])
                            orthanq_tp_fp_table = pd.concat([orthanq_tp_fp_table, row_to_add], ignore_index=True)

                        #only include predictions having density above determined threshold
                        threshold = 0.7
                        sum_of_densities = 0
                        sum_of_densities = best_results["density"].sum()
                        print(sample_name)
                        print(best_results)
                        print(sum_of_densities)
                        if sum_of_densities > threshold:
                            #count the number of samples that are called
                            if locus == locus_name:
                                samples_called += 1
                            #loop over best results (rows that share the same best odds score, i.e. same density and solutions but different haplotypes)
                            for (i,_) in best_results.iterrows():
                                #collect first two fields and cumulative fractions of haplotypes
                                #required to group haplotypes with identical two field values (ground truth contains two field information)
                                two_field_fractions = {}
                                for haplo in haplotypes:
                                    haplo_splt = haplo.split(":")
                                    two_field = haplo_splt[0] + ":" + haplo_splt[1]
                                    if two_field in two_field_fractions:
                                        two_field_fractions[two_field] += best_results[haplo][i]
                                    else:
                                        two_field_fractions[two_field] = best_results[haplo][i]
                                print("two_field_fractions: ",two_field_fractions)

                                #collect ground truth in values_in_truth, handle special cases e.g. 23:01/02/04
                                values_in_truth = {}
                                for (chr_index,chr_number) in enumerate([1, 2]):

                                    #get ground truth values for the sample and locus
                                    locus_in_truth = "HLA-" + locus + " " + str(chr_number)
                                    value_in_truth = locus+"*"+ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_in_truth].array[0] ##array[0] to reach what's inside the Series
                                    print("values_in_truth: ",values_in_truth)

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
                                            fractions += collected_fractions[max(collected_fractions)] ##get a random one, max() can be used (orders allele names alphabetically), it doesn't matter anyways in case of diploid priors
                                            prev_value = max(collected_fractions)
                                print("fractions: ",fractions)

                                #score only haplotypes that make up more than 50%
                                #also enter the entry to the tp_fp table to use in the plots in further rules
                                if fractions > 0.5: 
                                    collected += 1.0

                                    #fill up the orthanq_tp_fp_table with 1 for TP
                                    orthanq_tp_fp_table.loc[(orthanq_tp_fp_table.Sample == sample_name) & (orthanq_tp_fp_table.Locus == locus_name) & (orthanq_tp_fp_table.Best_Density == best_density),'TP'] = 1
                                    print("inner loop collected: ", collected)
                                    break ##break as soon as the collected gets +1 
            print("collected: ",collected)

            #calculate accuracy and call rate
            n_samples_collected = len(list(set(samples_collected)))
            print("n_samples_collected: ",n_samples_collected)
            # to be fair to all tools, the denominator should be the number of samples that the tool results in a prediction.
            accuracy = 100*collected/samples_called
            print("samples_called: ",samples_called)
            call_rate = samples_called/n_samples_collected

            #concatenate the locus-wise statistics to the validation table
            new_row = pd.DataFrame([[locus, n_samples_collected, call_rate, accuracy]],
                            columns=['Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'])
            validation_table = pd.concat([validation_table, new_row], ignore_index=True)
        return validation_table, orthanq_tp_fp_table

    ##low
    #initialize the table to output TP and FP samples
    orthanq_validation_table_low= pd.DataFrame(columns=('Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'))
    orthanq_tp_fp_table_low = pd.DataFrame(columns=('Sample', 'Locus', 'TP', 'Best_Density'))

    #validate low samples
    sample_list_low = samples_low_coverage["Run Accession"].to_list()
    orthanq_validation_table_low = validate_orthanq(orthanq_input, orthanq_validation_table_low, orthanq_tp_fp_table_low, samples_low_coverage, sample_list_low)
    print("orthanq_validation_table_low")
    print(orthanq_validation_table_low)
    final_orthanq_low = orthanq_validation_table_low[0]
    print("final_orthanq_low")
    print(final_orthanq_low)
    
    ##high 
    #initialize the table to output TP and FP samples
    orthanq_validation_table_high= pd.DataFrame(columns=('Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'))
    orthanq_tp_fp_table_high = pd.DataFrame(columns=('Sample', 'Locus', 'TP', 'Best_Density'))

    #validate high samples
    sample_list_evaluated = samples_evaluated["Run Accession"].to_list()
    orthanq_validation_table_high = validate_orthanq(orthanq_input, orthanq_validation_table_high, orthanq_tp_fp_table_high, samples_evaluated, sample_list_evaluated)
    
    final_orthanq_high = orthanq_validation_table_high[0]

    ##all
    #initialize the table to output TP and FP samples
    orthanq_validation_table_all= pd.DataFrame(columns=('Locus', 'N', 'Orthanq_Call_Rate', 'Orthanq_Accuracy'))
    orthanq_tp_fp_table_all = pd.DataFrame(columns=('Sample', 'Locus', 'TP', 'Best_Density'))

    #validate all samples
    sample_list_all = samples_all["Run Accession"].to_list()
    orthanq_validation_results_all = validate_orthanq(orthanq_input, orthanq_validation_table_all, orthanq_tp_fp_table_all, samples_all, sample_list_all)

    final_orthanq_all = orthanq_validation_results_all[0]
    final_tp_fp_table_all = orthanq_validation_results_all[1]

    #write the validation table


    final_orthanq_high.to_csv(
        snakemake.output.validation_high, sep="\t", index=False, header=True
    )

    final_orthanq_low.to_csv(
        snakemake.output.validation_low, sep="\t", index=False, header=True
    )

    final_orthanq_all.to_csv(
        snakemake.output.validation_all, sep="\t", index=False, header=True
    )

    #write the orthanq_tp_fp_table table for all
    final_tp_fp_table_all.to_csv(
        snakemake.output.tp_fp_table, sep="\t", index=False, header=True
    )