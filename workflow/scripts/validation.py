# Stratify validation for high coverage and low coverage samples. Generate separate tables for each tool.
# this script also creates a table for combining all results locuswise and containing information about TP, FP, no call and there is a call but the allele is not supposed to be found in the ground truth.
# For accuracy, we count, for each tool: TP / (TP+FP).
# For call rate, we count, for each tool: (TP + FP) / (TP + FP + 'no call')

import pandas as pd
import os 

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #read tables
    #orthanq
    ground_truth = pd.read_csv(snakemake.input.ground_truth, sep = "\t")
    orthanq_input_low = pd.read_csv(snakemake.input.orthanq_low, sep = "\t")
    orthanq_input_high = pd.read_csv(snakemake.input.orthanq_high, sep = "\t")
    orthanq_input_all = pd.read_csv(snakemake.input.orthanq_all, sep = "\t")
    orthanq_A_tp_fp = pd.read_csv(snakemake.input.orthanq_A_tp_fp, sep = "\t")
    orthanq_B_tp_fp = pd.read_csv(snakemake.input.orthanq_B_tp_fp, sep = "\t")
    orthanq_C_tp_fp = pd.read_csv(snakemake.input.orthanq_C_tp_fp, sep = "\t")
    orthanq_DQB1_tp_fp = pd.read_csv(snakemake.input.orthanq_DQB1_tp_fp, sep = "\t")
    allele_freqs = pd.read_csv(snakemake.input.allele_freqs)

    #other tools
    arcasHLA_input = pd.read_csv(snakemake.input.arcasHLA, sep = "\t")
    hla_la_input = pd.read_csv(snakemake.input.hla_la, sep = "\t")
    optitype_input = pd.read_csv(snakemake.input.optitype, sep = "\t")

    # read sample sheets for evaluated and low coverage samples
    samples_evaluated = pd.read_csv(snakemake.input.samples_evaluated, sep = "\t")
    samples_low_coverage = pd.read_csv(snakemake.input.samples_low_coverage, sep = "\t")
    
    #rename sample column
    samples_evaluated = samples_evaluated.rename(columns={"Run Accession": "sample"})
    samples_low_coverage = samples_low_coverage.rename(columns={"Run Accession": "sample"})

    # split tables for high and low coverage samples
    #low
    # orthanq_input_low = pd.merge(orthanq_input_low, samples_low_coverage[["sample"]], on = "sample", how = "left")
    arcasHLA_input_low = pd.merge(arcasHLA_input, samples_low_coverage[["sample"]], on = "sample")
    hla_la_input_low = pd.merge(hla_la_input, samples_low_coverage[["sample"]], on = "sample")
    optitype_input_low = pd.merge(optitype_input, samples_low_coverage[["sample"]], on = "sample")
    
    #high
    # orthanq_input_high = pd.merge(orthanq_input_high, samples_evaluated[["sample"]], on = "sample", how = "left")
    arcasHLA_input_high = pd.merge(arcasHLA_input, samples_evaluated[["sample"]], on = "sample")
    hla_la_input_high = pd.merge(hla_la_input, samples_evaluated[["sample"]], on = "sample")
    optitype_input_high = pd.merge(optitype_input, samples_evaluated[["sample"]], on = "sample")    

    def truth_for_sample(locus, values_in_truth, ground_truth, sample_name):
        for (chr_index,chr_number) in enumerate([1, 2]):
            locus_name_in_truth = "HLA-" + locus + " " + str(chr_number)
            value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_name_in_truth].array[0] ##array[0] to reach what's inside the Series
            
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

    def arcashla_validation(arcasHLA_results, ground_truth, arcasHLA_validation_table, arcasHLA_tp_fp_A, arcasHLA_tp_fp_B, arcasHLA_tp_fp_C, arcasHLA_tp_fp_DQB1):
        #arcasHLA
        #loop over all loci
        loci = ['A', 'B', 'C', 'DQB1']
        for locus in loci:
            values = arcasHLA_results[locus]
            collected = 0
            samples_called = 0 #some samples in arcashla have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, value_in_arcasHLA) in enumerate(values):
                sample_name = arcasHLA_results.loc[index,'sample']
                values_in_truth = {}
                values_in_truth = truth_for_sample(locus, values_in_truth, ground_truth, sample_name)
                
                values_in_truth_list = []
                # for chr_number in [1, 2]:
                #     locus_name_in_truth = "HLA-" + locus + " " + str(chr_number)
                #     value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_name_in_truth].array[0]
                #     value_in_truth=value_in_truth.rstrip("*") #remove the alleles ending with "*"
                #     values_in_truth_list.append(locus + "*" + value_in_truth)
                for key,values in values_in_truth.items():
                    values_in_truth_list.append(";".join(values))
                print(values_in_truth_list)

                if pd.isnull(value_in_arcasHLA):
                    #fill in the tp fp table
                    if locus == "A":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_A = pd.concat([arcasHLA_tp_fp_A, new_row_tp], ignore_index=True)
                    if locus == "B":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_B = pd.concat([arcasHLA_tp_fp_B, new_row_tp], ignore_index=True)
                    if locus == "C":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_C = pd.concat([arcasHLA_tp_fp_C, new_row_tp], ignore_index=True)
                    if locus == "DQB1":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_DQB1 = pd.concat([arcasHLA_tp_fp_DQB1, new_row_tp], ignore_index=True)
                    continue
                else:
                    sample_name = arcasHLA_results.loc[index,'sample']
                    alleles = []
                    for (chr_index,chr_number) in enumerate([1, 2]):
                        #value in arcasHLA 
                        first_allele = value_in_arcasHLA.split("/")[chr_index]
                        # first_allele_no_locus = first_allele.split("*")[1]
                        allele_in_arcasHLA = first_allele.split(":")[0] + ":" +first_allele.split(":")[1]
                        alleles.append(allele_in_arcasHLA)

                    values_in_truth_clone = values_in_truth
                    print("sample name: " ,sample_name)
                    print("arcasHLA prediction: " + ''.join(alleles))
                    print("truth_values: " + str(values_in_truth))

                    #fill in the tp fp table
                    if locus == 'A':
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_arcasHLA, '']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_A = pd.concat([arcasHLA_tp_fp_A, new_row_tp], ignore_index=True)
                    if locus == 'B':
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_arcasHLA, '']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_B = pd.concat([arcasHLA_tp_fp_B, new_row_tp], ignore_index=True)
                    if locus == 'C':
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_arcasHLA, '']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_C = pd.concat([arcasHLA_tp_fp_C, new_row_tp], ignore_index=True)
                    if locus == 'DQB1':
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_arcasHLA, '']],
                            columns=['sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'])
                        arcasHLA_tp_fp_DQB1 = pd.concat([arcasHLA_tp_fp_DQB1, new_row_tp], ignore_index=True)

                    allele_present=0
                    for allele in alleles: ##???
                        values_in_truth = values_in_truth_clone
                        for chr_index, chr_values in values_in_truth.items():
                            if allele in chr_values:
                                allele_present+= 1 #both alleles should be correct
                                #to stop homozygous alleles inaccurately match, we should remove the one that matches
                                if len(values_in_truth) > 1:
                                    values_in_truth_clone.pop(chr_index)
                                break
                    samples_called+=1
                    if allele_present==2: #both alleles should be correct
                        collected+=1

                        #fill in the table for tp fp info
                        if locus == "A":
                            arcasHLA_tp_fp_A.loc[arcasHLA_tp_fp_A['sample'] == sample_name, 'arcasHLA evaluation'] = 'TP'
                        if locus == "B":
                            arcasHLA_tp_fp_B.loc[arcasHLA_tp_fp_B['sample'] == sample_name, 'arcasHLA evaluation'] = 'TP'
                        if locus == "C":
                            arcasHLA_tp_fp_C.loc[arcasHLA_tp_fp_C['sample'] == sample_name, 'arcasHLA evaluation'] = 'TP'
                        if locus == "DQB1":
                            arcasHLA_tp_fp_DQB1.loc[arcasHLA_tp_fp_DQB1['sample'] == sample_name, 'arcasHLA evaluation'] = 'TP'
                    print("collected: " + str(collected))
            call_rate = samples_called/len(arcasHLA_results.index)
            # to be fair to all tools, the denominator should be the number of samples that the tool results in a prediction.
            if samples_called != 0:
                accuracy = 100*collected/samples_called
            else:
                 accuracy = 0
            print("len(arcasHLA_results.index):"+str(len(arcasHLA_results.index)))
            
            # add the numbers in parenthesis for call rate and accuracy
            accuracy = "{:.2f}".format(accuracy) + " (" + str(int(collected)) + "/" + str(samples_called) + ")"
            call_rate = "{:.2f}".format(call_rate) + " (" + str(samples_called) + "/" + str(len(arcasHLA_results.index)) + ")"

            # concatenate the row the validation table 
            new_row = pd.DataFrame([[locus, len(arcasHLA_results.index), call_rate, accuracy]],
                        columns=['Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'])
            arcasHLA_validation_table = pd.concat([arcasHLA_validation_table, new_row], ignore_index=True)
                
        #if the evaluation field is empty, make all empty to FP
        for index,row in arcasHLA_tp_fp_A.iterrows():
            if row["arcasHLA evaluation"] == "":
                arcasHLA_tp_fp_A.loc[index, "arcasHLA evaluation"] = "FP"
        for index,row in arcasHLA_tp_fp_B.iterrows():
            if row["arcasHLA evaluation"] == "":
                arcasHLA_tp_fp_B.loc[index, "arcasHLA evaluation"] = "FP"
        for index,row in arcasHLA_tp_fp_C.iterrows():
            if row["arcasHLA evaluation"] == "":
                arcasHLA_tp_fp_C.loc[index, "arcasHLA evaluation"] = "FP"
        for index,row in arcasHLA_tp_fp_DQB1.iterrows():
            if row["arcasHLA evaluation"] == "":
                arcasHLA_tp_fp_DQB1.loc[index, "arcasHLA evaluation"] = "FP"

        return arcasHLA_validation_table, arcasHLA_tp_fp_A, arcasHLA_tp_fp_B, arcasHLA_tp_fp_C, arcasHLA_tp_fp_DQB1

    def hla_la_validation(hla_la_results, ground_truth, hla_la_validation_table, hla_la_tp_fp_A, hla_la_tp_fp_B, hla_la_tp_fp_C, hla_la_tp_fp_DQB1):
        # #HLA-LA
        # #loop over all loci
        loci = ['A', 'B', 'C', 'DQB1']
        for locus in loci:
            values = hla_la_results[locus]
            collected = 0
            samples_called = 0 #some samples in hla-la have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, value_in_hla_la) in enumerate(values):
                sample_name = hla_la_results.loc[index,'sample']
                values_in_truth = {}
                values_in_truth = truth_for_sample(locus, values_in_truth, ground_truth, sample_name)
                
                values_in_truth_list = []
                # for chr_number in [1, 2]:
                #     locus_name_in_truth = "HLA-" + locus + " " + str(chr_number)
                #     value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_name_in_truth].array[0]
                #     value_in_truth=value_in_truth.rstrip("*") #remove the alleles ending with "*"                    
                #     values_in_truth_list.append(locus + "*" + value_in_truth)
                for key,values in values_in_truth.items():
                    values_in_truth_list.append(";".join(values))
                print(values_in_truth_list)

                if pd.isnull(value_in_hla_la):
                    #fill in the tp fp table
                    if locus == "A":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_A = pd.concat([hla_la_tp_fp_A, new_row_tp], ignore_index=True)
                    if locus == "B":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_B = pd.concat([hla_la_tp_fp_B, new_row_tp], ignore_index=True)
                    if locus == "C":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_C = pd.concat([hla_la_tp_fp_C, new_row_tp], ignore_index=True)
                    if locus == "DQB1":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_DQB1 = pd.concat([hla_la_tp_fp_DQB1, new_row_tp], ignore_index=True)
                    continue
                else:
                    sample_name = hla_la_results.loc[index,'sample']
                    alleles = {}
                    for (chr_index,chr_number) in enumerate([1, 2]):
                        #values in hla-la are as in the following: A*02:01:02;A*03:01:04/A*02:01:02;A*03:01:04
                        all_alleles=[]
                        first_allele = value_in_hla_la.split("/")[chr_index]
                        splitted_alleles = first_allele.split(";")
                        for splitted in splitted_alleles:
                            # allele_no_locus = splitted.split("*")[1]
                            allele_in_hla_la = splitted.split(":")[0] + ":" +splitted.split(":")[1]
                            all_alleles.append(allele_in_hla_la)
                        alleles["{0}".format(chr_index)] = all_alleles

                    values_in_truth_clone = values_in_truth
                    print("sample name: " ,sample_name)
                    print("hla-la prediction: " + ''.join(sum(alleles.values(), [])))

                    #collecting alleles_w_separators is necessary as some truth values have multiple values
                    alleles_w_separators = []
                    for key_chr,value_alleles in alleles.items():
                        alleles_with_sep = []
                        for value in value_alleles:
                            alleles_with_sep.append(value)
                        alleles_with_sep = ';'.join(alleles_with_sep)
                        alleles_w_separators.append(alleles_with_sep)
                    print(alleles_w_separators)

                    #fill in the tp fp table
                    if locus == "A":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_hla_la, '']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_A = pd.concat([hla_la_tp_fp_A, new_row_tp], ignore_index=True)
                    if locus == "B":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_hla_la, '']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_B = pd.concat([hla_la_tp_fp_B, new_row_tp], ignore_index=True)
                    if locus == "C":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_hla_la, '']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_C = pd.concat([hla_la_tp_fp_C, new_row_tp], ignore_index=True)
                    if locus == "DQB1":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_hla_la, '']],
                            columns=['sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'])
                        hla_la_tp_fp_DQB1 = pd.concat([hla_la_tp_fp_DQB1, new_row_tp], ignore_index=True)

                    allele_present=0 #both alleles should be correct
                    for _, hla_la_alleles in alleles.items():
                        for allele in hla_la_alleles:
                            values_in_truth = values_in_truth_clone
                            for chr_index, chr_values in values_in_truth.items():
                                if allele in chr_values:
                                    allele_present+= 1
                                    #to stop homozygous alleles inaccurately match, we should remove the one that matches
                                    if len(values_in_truth) > 1:
                                        values_in_truth_clone.pop(chr_index)
                                    break
                            else:
                                continue
                            break
                    if allele_present==2: #both alleles should be correct
                        collected+=1

                        #fill in the table for tp fp info
                        if locus == "A":
                            hla_la_tp_fp_A.loc[hla_la_tp_fp_A['sample'] == sample_name, 'HLA-LA evaluation'] = 'TP'
                        if locus == "B":
                            hla_la_tp_fp_B.loc[hla_la_tp_fp_B['sample'] == sample_name, 'HLA-LA evaluation'] = 'TP'
                        if locus == "C":
                            hla_la_tp_fp_C.loc[hla_la_tp_fp_C['sample'] == sample_name, 'HLA-LA evaluation'] = 'TP'
                        if locus == "DQB1":
                            hla_la_tp_fp_DQB1.loc[hla_la_tp_fp_DQB1['sample'] == sample_name, 'HLA-LA evaluation'] = 'TP'
                    print("collected: " + str(collected))
                    samples_called+=1
            call_rate = samples_called/len(hla_la_results.index)
            # to be fair to all tools, the denominator should be the number of samples that the tool results in a prediction.
            accuracy = 100*collected/samples_called
            print("len(hla_la_results.index):"+str(len(hla_la_results.index)))      

            # add the numbers in parenthesis for call rate and accuracy
            accuracy = "{:.2f}".format(accuracy) + " (" + str(int(collected)) + "/" + str(samples_called) + ")"
            call_rate = "{:.2f}".format(call_rate) + " (" + str(samples_called) + "/" + str(len(hla_la_results.index)) + ")"

            # concatenate the row the validation table  
            new_row = pd.DataFrame([[locus, len(hla_la_results.index), call_rate, accuracy]],
                        columns=['Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'])
            hla_la_validation_table = pd.concat([hla_la_validation_table, new_row], ignore_index=True)
                        
        #if the evaluation field is empty, make all empty to FP
        for index,row in hla_la_tp_fp_A.iterrows():
            if row["HLA-LA evaluation"] == "":
                hla_la_tp_fp_A.loc[index, "HLA-LA evaluation"] = "FP"
        for index,row in hla_la_tp_fp_B.iterrows():
            if row["HLA-LA evaluation"] == "":
                hla_la_tp_fp_B.loc[index, "HLA-LA evaluation"] = "FP"
        for index,row in hla_la_tp_fp_C.iterrows():
            if row["HLA-LA evaluation"] == "":
                hla_la_tp_fp_C.loc[index, "HLA-LA evaluation"] = "FP"
        for index,row in hla_la_tp_fp_DQB1.iterrows():
            if row["HLA-LA evaluation"] == "":
                hla_la_tp_fp_DQB1.loc[index, "HLA-LA evaluation"] = "FP"
        return hla_la_validation_table, hla_la_tp_fp_A, hla_la_tp_fp_B, hla_la_tp_fp_C, hla_la_tp_fp_DQB1
    
    def optitype_validation(optitype_results, ground_truth, optitype_validation_table, optitype_tp_fp_A, optitype_tp_fp_B, optitype_tp_fp_C, optitype_tp_fp_DQB1):
        #optitype
        #loop over all loci
        loci = ['A', 'B', 'C']
        print("optitype")
        for locus in loci:
            values = optitype_results[locus]
            collected = 0
            samples_called = 0 #some samples in optitype have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, value_in_optitype) in enumerate(values):
                sample_name = optitype_results.loc[index,'sample']
                values_in_truth = {}
                values_in_truth = truth_for_sample(locus, values_in_truth, ground_truth, sample_name)
                
                values_in_truth_list = []
                # for chr_number in [1, 2]:
                #     locus_name_in_truth = "HLA-" + locus + " " + str(chr_number)
                #     value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_name_in_truth].array[0]
                #     value_in_truth=value_in_truth.rstrip("*") #remove the alleles ending with "*"
                #     values_in_truth_list.append(locus + "*" + value_in_truth)
                for key,values in values_in_truth.items():
                    values_in_truth_list.append(";".join(values))
                print(values_in_truth_list)

                if pd.isnull(value_in_optitype):
                    #fill in the tp fp table
                    if locus == "A":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_A = pd.concat([optitype_tp_fp_A, new_row_tp], ignore_index=True)
                    if locus == "B":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_B = pd.concat([optitype_tp_fp_B, new_row_tp], ignore_index=True)
                    if locus == "C":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_C = pd.concat([optitype_tp_fp_C, new_row_tp], ignore_index=True)
                    if locus == "DQB1":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), '', 'no call']],
                            columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_DQB1 = pd.concat([optitype_tp_fp_DQB1, new_row_tp], ignore_index=True)
                    continue
                else:
                    alleles = []
                    for (chr_index,chr_number) in enumerate([1, 2]):
                        #value in optitype 
                        first_allele = value_in_optitype.split("/")[chr_index]
                        # allele_in_arcasHLA = first_allele.split("*")[1]
                        alleles.append(first_allele)

                    values_in_truth_clone = values_in_truth
                    print("sample name: " ,sample_name)
                    print("optitype prediction: " + ''.join(alleles))
                    print("truth_values: " + str(values_in_truth))

                    #fill in the tp fp table
                    if locus == "A":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_optitype, '']],
                        columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_A = pd.concat([optitype_tp_fp_A, new_row_tp], ignore_index=True)
                    if locus == "B":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_optitype, '']],
                        columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_B = pd.concat([optitype_tp_fp_B, new_row_tp], ignore_index=True)
                    if locus == "C":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_optitype, '']],
                        columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_C = pd.concat([optitype_tp_fp_C, new_row_tp], ignore_index=True)
                    if locus == "DQB1":
                        new_row_tp = pd.DataFrame([[sample_name, '/'.join(values_in_truth_list), value_in_optitype, '']],
                        columns=['sample', 'ground', 'optitype', 'optitype evaluation'])
                        optitype_tp_fp_DQB1 = pd.concat([optitype_tp_fp_DQB1, new_row_tp], ignore_index=True)

                    allele_present = 0
                    for allele in alleles: ##???
                        values_in_truth = values_in_truth_clone
                        for chr_index,chr_values in values_in_truth.items():
                            if allele in chr_values:
                                allele_present+= 1
                                #to stop homozygous alleles inaccurately match, we should remove the one that matches
                                if len(values_in_truth) > 1:
                                    values_in_truth_clone.pop(chr_index)
                                break
                    if allele_present==2:
                        collected+=1

                        #fill in the table for tp fp info
                        if locus == "A":
                            optitype_tp_fp_A.loc[optitype_tp_fp_A['sample'] == sample_name, 'optitype evaluation'] = 'TP'
                        if locus == "B":
                            optitype_tp_fp_B.loc[optitype_tp_fp_B['sample'] == sample_name, 'optitype evaluation'] = 'TP'
                        if locus == "C":
                            optitype_tp_fp_C.loc[optitype_tp_fp_C['sample'] == sample_name, 'optitype evaluation'] = 'TP'
                        if locus == "DQB1":
                            optitype_tp_fp_DQB1.loc[optitype_tp_fp_DQB1['sample'] == sample_name, 'optitype evaluation'] = 'TP'
                    print("collected: " + str(collected))
                    samples_called+=1
            call_rate = samples_called/len(optitype_results.index)
            # to be fair to all tools, the denominator should be the number of samples that the tool results in a prediction.
            accuracy = 100*collected/samples_called
            print("len(optitype_results.index):"+str(len(optitype_results.index)))      

            # add the numbers in parenthesis for call rate and accuracy
            accuracy = "{:.2f}".format(accuracy)  + " (" + str(int(collected)) + "/" + str(samples_called) + ")"
            call_rate = "{:.2f}".format(call_rate)  + " (" + str(samples_called) + "/" + str(len(optitype_results.index)) + ")"

            # concatenate the row the validation table
            new_row = pd.DataFrame([[locus, len(optitype_results.index), call_rate, accuracy]],
                        columns=['Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'])
            optitype_validation_table = pd.concat([optitype_validation_table, new_row], ignore_index=True)
                                   
        #if the evaluation field is empty, make all empty to FP
        for index,row in optitype_tp_fp_A.iterrows():
            if row["optitype evaluation"] == "":
                optitype_tp_fp_A.loc[index, "optitype evaluation"] = "FP"
        for index,row in optitype_tp_fp_B.iterrows():
            if row["optitype evaluation"] == "":
                optitype_tp_fp_B.loc[index, "optitype evaluation"] = "FP"
        for index,row in optitype_tp_fp_C.iterrows():
            if row["optitype evaluation"] == "":
                optitype_tp_fp_C.loc[index, "optitype evaluation"] = "FP"
        for index,row in optitype_tp_fp_DQB1.iterrows():
            if row["optitype evaluation"] == "":
                optitype_tp_fp_DQB1.loc[index, "optitype evaluation"] = "FP"
        return optitype_validation_table, optitype_tp_fp_A, optitype_tp_fp_B, optitype_tp_fp_C, optitype_tp_fp_DQB1

    #check if the predicted alleles are in the allele freq table and below or above the freq filter
    def check_alleles_in_database(predictions_table, allele_freq_table):
        tool_name = predictions_table.columns[2]
        predictions = predictions_table[tool_name].to_list()
        for i in range(len(predictions)):
            if predictions[i] != '':
                if tool_name == "HLA-LA":
                    splitted = predictions[i].split("/")
                    for full_field_allele in splitted:
                        #check for two field resolution of hla-la prediction, three field doesn't make sense
                        count_of_alleles = 0
                        for single_allele in full_field_allele.split(";"):
                            splitted = single_allele.split(":")
                            first_allele_three_field = ""
                            first_allele_two_field = ""
                            if len(splitted) == 4:
                                first_allele_three_field = splitted[0] + ":" + splitted[1] + ":" + splitted[2] 
                            if len(splitted) == 3:
                                first_allele_three_field = splitted[0] + ":" + splitted[1]
                            #all alleles should be below 0.05 to be labeled as below.
                            for allele, freq in zip(allele_freq_table["var"], allele_freq_table["frequency"]):
                                if freq > 0.05:
                                    if single_allele == allele:
                                        count_of_alleles += 1
                                        break
                                    elif first_allele_three_field != "" and first_allele_three_field == allele:
                                        count_of_alleles += 1
                                        break
                                    elif first_allele_two_field != "" and first_allele_two_field == allele:
                                        count_of_alleles += 1
                                        break
                            # the check for the first two fields of a three field or four field containing record, only if first two fields don't match is not necessary.
                            # because the info in the db belongs to a subset of the two field allele, so that wouldn't be a valid comparison.
                            # else: # if the break is not entered above, then else is executed in this for/else construct.
                        if count_of_alleles == 0: # all alleles should be below the freq
                            predictions_table.loc[i, "HLA-LA evaluation"] = "allele not considered in the truth set"
                else: # the case for arcas (3 field) and optitype (two field)
                    splitted = predictions[i].split("/")
                    for allele_in_tool in splitted: 
                        #all alleles should be below 0.05 to be labeled as below.
                        two_field = ""
                        count_of_alleles = 0
                        if len(allele_in_tool.split(":")) == 3:
                            splitted_three_field = allele_in_tool.split(":")
                            two_field = splitted_three_field[0] + ":" + splitted_three_field[1]
                        for allele, freq in zip(allele_freq_table["var"], allele_freq_table["frequency"]):
                            if freq > 0.05:
                                if allele_in_tool == allele:
                                    count_of_alleles += 1
                                    break
                                elif two_field != "" and two_field == allele:
                                    count_of_alleles += 1
                                    break
                        # the check for the first two fields of a three field or four field containing record, only if first two fields don't match is not necessary.
                        # because the info in the db belongs to a subset of the two field allele, so that wouldn't be a valid comparison.
                        else: # if the break is not entered above, then else is executed in this for/else construct.
                            predictions_table.loc[i,  tool_name + " evaluation"] = "allele not considered in the truth set"
        print("predictions_table")
        print(predictions_table)
        return predictions_table

    #low coverage samples

    ##initialize performance tables for tools 
    arcasHLA_validation_table_low = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'))
    hla_la_validation_table_low = pd.DataFrame(columns=('Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'))
    optitype_validation_table_low = pd.DataFrame(columns=('Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'))

    #initialize locus-wise tp fp tables locus wise
    arcasHLA_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))

    optitype_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))

    hla_la_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))

    #merge tables to get the final table, the zeroth index gives the validation table for all tools
    arcasHLA_validation_table_low = arcashla_validation(arcasHLA_input_low, ground_truth, arcasHLA_validation_table_low,  arcasHLA_tp_fp_A, arcasHLA_tp_fp_B, arcasHLA_tp_fp_C, arcasHLA_tp_fp_DQB1)
    first_merge = pd.merge(orthanq_input_low, arcasHLA_validation_table_low[0], how='left', on=['Locus', 'N'])
    hla_la_validation_table_low = hla_la_validation(hla_la_input_low, ground_truth, hla_la_validation_table_low, hla_la_tp_fp_A, hla_la_tp_fp_B, hla_la_tp_fp_C, hla_la_tp_fp_DQB1)
    second_merge = pd.merge(first_merge, hla_la_validation_table_low[0], how='left', on=['Locus', 'N'])
    optitype_validation_table_low = optitype_validation(optitype_input_low, ground_truth, optitype_validation_table_low, optitype_tp_fp_A, optitype_tp_fp_B, optitype_tp_fp_C, optitype_tp_fp_DQB1)
    low_final_table = pd.merge(second_merge, optitype_validation_table_low[0], how='left', on=['Locus', 'N'])
    print(low_final_table)
    #high coverage samples

    ##initialize performance tables for tools 
    arcasHLA_validation_table_high = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'))
    hla_la_validation_table_high = pd.DataFrame(columns=('Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'))
    optitype_validation_table_high = pd.DataFrame(columns=('Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'))

    #initialize locus-wise tp fp tables locus wise
    arcasHLA_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))

    optitype_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))

    hla_la_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    
    #merge tables to get the final table
    arcasHLA_validation_table_high = arcashla_validation(arcasHLA_input_high, ground_truth, arcasHLA_validation_table_high, arcasHLA_tp_fp_A, arcasHLA_tp_fp_B, arcasHLA_tp_fp_C, arcasHLA_tp_fp_DQB1)
    first_merge = pd.merge(orthanq_input_high, arcasHLA_validation_table_high[0], how='left', on=['Locus', 'N'])
    hla_la_validation_table_high = hla_la_validation(hla_la_input_high, ground_truth, hla_la_validation_table_high, hla_la_tp_fp_A, hla_la_tp_fp_B, hla_la_tp_fp_C, hla_la_tp_fp_DQB1)
    second_merge = pd.merge(first_merge, hla_la_validation_table_high[0], how='left', on=['Locus', 'N'])
    optitype_validation_table_high = optitype_validation(optitype_input_high, ground_truth, optitype_validation_table_high, optitype_tp_fp_A, optitype_tp_fp_B, optitype_tp_fp_C, optitype_tp_fp_DQB1)
    high_final_table = pd.merge(second_merge, optitype_validation_table_high[0], how='left', on=['Locus', 'N'])
    #all samples

    ##initialize performance tables for tools 
    arcasHLA_validation_table_all = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'))
    hla_la_validation_table_all = pd.DataFrame(columns=('Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'))
    optitype_validation_table_all = pd.DataFrame(columns=('Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'))

    #initialize locus-wise tp fp tables locus wise
    arcasHLA_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))
    arcasHLA_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'arcasHLA', 'arcasHLA evaluation'))

    optitype_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))
    optitype_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'optitype', 'optitype evaluation'))

    hla_la_tp_fp_A = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_B = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_C = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))
    hla_la_tp_fp_DQB1 = pd.DataFrame(columns=('sample', 'ground', 'HLA-LA', 'HLA-LA evaluation'))

    #merge validation tables to get the final table
    arcasHLA_validation_output = arcashla_validation(arcasHLA_input, ground_truth, arcasHLA_validation_table_all, arcasHLA_tp_fp_A, arcasHLA_tp_fp_B, arcasHLA_tp_fp_C, arcasHLA_tp_fp_DQB1)
    first_merge = pd.merge(orthanq_input_all, arcasHLA_validation_output[0], how='left', on=['Locus', 'N'])
    hla_la_validation_output = hla_la_validation(hla_la_input, ground_truth, hla_la_validation_table_all, hla_la_tp_fp_A, hla_la_tp_fp_B, hla_la_tp_fp_C, hla_la_tp_fp_DQB1)
    second_merge = pd.merge(first_merge, hla_la_validation_output[0], how='left', on=['Locus', 'N'])
    optitype_validation_output = optitype_validation(optitype_input, ground_truth, optitype_validation_table_all, optitype_tp_fp_A, optitype_tp_fp_B, optitype_tp_fp_C, optitype_tp_fp_DQB1)
    all_final_table = pd.merge(second_merge, optitype_validation_output[0], how='left', on=['Locus', 'N'])

    #then, for tp fp tables for each locus, check if they contain alleles that are below the criteria according to the Abi-Rached 2018 paper.
    #if that's the case label them as "allele not considered in the truth set".
    #arcas
    arcasHLA_freq_checked_A = check_alleles_in_database(arcasHLA_validation_output[1], allele_freqs) 
    arcasHLA_freq_checked_B = check_alleles_in_database(arcasHLA_validation_output[2], allele_freqs) 
    arcasHLA_freq_checked_C = check_alleles_in_database(arcasHLA_validation_output[3], allele_freqs) 
    arcasHLA_freq_checked_DQB1 = check_alleles_in_database(arcasHLA_validation_output[4], allele_freqs) 

    #optitype
    optitype_freq_checked_A = check_alleles_in_database(optitype_validation_output[1], allele_freqs) 
    optitype_freq_checked_B = check_alleles_in_database(optitype_validation_output[2], allele_freqs) 
    optitype_freq_checked_C = check_alleles_in_database(optitype_validation_output[3], allele_freqs) 

    #HLA-LA
    hla_la_freq_checked_A = check_alleles_in_database(hla_la_validation_output[1], allele_freqs) 
    hla_la_freq_checked_B = check_alleles_in_database(hla_la_validation_output[2], allele_freqs) 
    hla_la_freq_checked_C = check_alleles_in_database(hla_la_validation_output[3], allele_freqs) 
    hla_la_freq_checked_DQB1 = check_alleles_in_database(hla_la_validation_output[4], allele_freqs) 

    #merge locus wise tp fp tables 1:A, 2:B, 3:C, 4:DQB1

    #A
    A_tp_fp = pd.merge(arcasHLA_freq_checked_A, optitype_freq_checked_A, on=['sample', 'ground'])
    A_tp_fp = pd.merge(A_tp_fp, hla_la_freq_checked_A, on=['sample', 'ground'])
    A_tp_fp = pd.merge(A_tp_fp, orthanq_A_tp_fp, on=['sample', 'ground'])

    #B
    B_tp_fp = pd.merge(arcasHLA_freq_checked_B, optitype_freq_checked_B, on=['sample', 'ground'])
    B_tp_fp = pd.merge(B_tp_fp, hla_la_freq_checked_B, on=['sample', 'ground'])
    B_tp_fp = pd.merge(B_tp_fp, orthanq_B_tp_fp, on=['sample', 'ground'])

    #C
    C_tp_fp = pd.merge(arcasHLA_freq_checked_C, optitype_freq_checked_C, on=['sample', 'ground'])
    C_tp_fp = pd.merge(C_tp_fp, hla_la_freq_checked_C, on=['sample', 'ground'])
    C_tp_fp = pd.merge(C_tp_fp, orthanq_C_tp_fp, on=['sample', 'ground'])

    #DQB1 
    DQB1_tp_fp = pd.merge(arcasHLA_freq_checked_DQB1, hla_la_freq_checked_DQB1, on=['sample', 'ground'])
    DQB1_tp_fp = pd.merge(DQB1_tp_fp, orthanq_DQB1_tp_fp, on=['sample', 'ground'])

    #finally,
    A_tp_fp["coverage"] = ""
    B_tp_fp["coverage"] = ""
    C_tp_fp["coverage"] = ""
    DQB1_tp_fp["coverage"] = ""

    def label_coverage(table, low_sample_sheet):
        for i, row in table[["sample"]].iterrows():
            if row["sample"] in low_sample_sheet["sample"].to_list():
                table.loc[i, "coverage"] = "low"
            else:
                table.loc[i, "coverage"] = "high"
        return table

    A_tp_fp_final = label_coverage(A_tp_fp, samples_low_coverage)
    B_tp_fp_final = label_coverage(B_tp_fp, samples_low_coverage)
    C_tp_fp_final = label_coverage(C_tp_fp, samples_low_coverage)
    DQB1_tp_fp_final = label_coverage(DQB1_tp_fp, samples_low_coverage)

    #write all tables to csv 

    #validation tables
    low_final_table.to_csv(
        snakemake.output.validation_low, sep="\t", index=False, header=True
    )

    high_final_table.to_csv(
        snakemake.output.validation_high, sep="\t", index=False, header=True
    )

    all_final_table.to_csv(
        snakemake.output.validation_all, sep="\t", index=False, header=True
    )

    #tp fp tables
    A_tp_fp_final.to_csv(
        snakemake.output.A_tp_fp, sep="\t", index=False, header=True
    )

    B_tp_fp_final.to_csv(
        snakemake.output.B_tp_fp, sep="\t", index=False, header=True
    )

    C_tp_fp_final.to_csv(
        snakemake.output.C_tp_fp, sep="\t", index=False, header=True
    )
    
    DQB1_tp_fp_final.to_csv(
        snakemake.output.DQB1_tp_fp, sep="\t", index=False, header=True
    )