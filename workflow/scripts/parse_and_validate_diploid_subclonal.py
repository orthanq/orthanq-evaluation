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
    loci = ['A', 'B', 'C', 'DQB1']
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
            best_result = pd.read_csv(orthanq_input[index]).iloc[[0]] #best result, first row

            #retrieve the predicted haplotypes
            filtered_cols = []
            filtered_cols = [col for col in best_result if col.startswith(locus)]

            #collect first two fields and cumulative fractions of haplotypes
            two_field_fractions = {}
            for col in filtered_cols:
                col_splt = col.split(":")
                two_field = col_splt[0] + ":" + col_splt[1]
                if two_field in two_field_fractions:
                    two_field_fractions[two_field] += best_result[col][0]
                else:
                    two_field_fractions[two_field] = best_result[col][0]
            print("two_field_fractions: ",two_field_fractions)

            #get ground truth values for the sample and locus
            values_in_truth = []
            locus_in_truth_1 = "HLA-" + locus + " " + str(1)
            locus_in_truth_2 = "HLA-" + locus + " " + str(2)
            values_in_truth.append(locus + "*" + ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_in_truth_1].array[0])
            values_in_truth.append(locus + "*" + ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_in_truth_2].array[0])
            print("values_in_truth: ",values_in_truth)

            #collect the fractions of haplotypes matching with ground truth
            prev_value = ""
            fractions = 0.0
            for value in values_in_truth:
                if value in two_field_fractions.keys() and value != prev_value:
                    prev_value = value
                    fractions += two_field_fractions[value]
            print("fractions: ",fractions)

            #score only haplotypes that make up more than 50%
            if fractions > 0.5:    
                collected += 1.0
            print("inner loop collected: ", collected)
        print(collected)
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