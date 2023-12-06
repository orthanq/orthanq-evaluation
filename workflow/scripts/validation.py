# Stratify validation for high coverage and low coverage samples. Generate separate tables for each tool.

import pandas as pd
import os 

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #read tables
    ground_truth = pd.read_csv(snakemake.input.ground_truth, sep = "\t")
    orthanq_input_low = pd.read_csv(snakemake.input.orthanq_low, sep = "\t")
    orthanq_input_high = pd.read_csv(snakemake.input.orthanq_high, sep = "\t")
    orthanq_input_all = pd.read_csv(snakemake.input.orthanq_all, sep = "\t")

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

    def arcashla_validation(arcasHLA_results, ground_truth, arcasHLA_validation_table):
        #arcasHLA
        #loop over all loci
        loci = ['A', 'B', 'C', 'DQB1']
        for locus in loci:
            values = arcasHLA_results[locus]
            collected = 0
            samples_called = 0 #some samples in arcashla have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, value_in_arcasHLA) in enumerate(values):
                if pd.isnull(value_in_arcasHLA):
                    continue
                else:
                    sample_name = arcasHLA_results.loc[index,'sample']
                    values_in_truth = {}
                    alleles = []
                    for (chr_index,chr_number) in enumerate([1, 2]):
                        #value in arcasHLA 
                        first_allele = value_in_arcasHLA.split("/")[chr_index]
                        first_allele_no_locus = first_allele.split("*")[1]
                        allele_in_arcasHLA = first_allele_no_locus.split(":")[0] + ":" +first_allele_no_locus.split(":")[1]
                        alleles.append(allele_in_arcasHLA)
                        locus_name_in_truth_1 = "HLA-" + locus + " " + str(chr_number)
                        value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_name_in_truth_1].array[0] ##array[0] to reach what's inside the Series
                        ##split in case of alles that have e.g. 23:01/02/04 in the truth
                        tmp_alleles = []
                        if "/" in value_in_truth:
                            first_split = value_in_truth.split("/")
                            first_field = first_split[0].split(":")[0]
                            tmp_alleles.append(first_split[0])
                            for splitted in first_split[1:]:
                                tmp_alleles.append(first_field + ":" +splitted)
                            values_in_truth["{0}".format(chr_index)] = tmp_alleles
                            
                        elif value_in_truth.endswith("*"):
                            value_in_truth=value_in_truth.rstrip("*")
                            values_in_truth["{0}".format(chr_index)] = [value_in_truth]
                        else:
                            values_in_truth["{0}".format(chr_index)] = [value_in_truth]

                    values_in_truth_clone = values_in_truth
                    print("sample name: " ,sample_name)
                    print("arcasHLA prediction: " + ''.join(alleles))
                    print("truth_values: " + str(values_in_truth))
                    allele_present=0
                    for allele in alleles: ##???
                        values_in_truth = values_in_truth_clone
                        for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
                            if allele in chr_values:
                                allele_present+= 1 #both alleles should be correct
                                #to stop homozygous alleles inaccurately match, we should remove the one that matches
                                if len(values_in_truth) > 1:
                                    values_in_truth_clone.pop(chr_index, None)
                                break
                    samples_called+=1
                    if allele_present==2: #both alleles should be correct
                        collected+=1
                    print("collected: " + str(collected))
            call_rate = samples_called/len(arcasHLA_results.index)
            accuracy = 100*collected/len(arcasHLA_results.index)
            print("len(arcasHLA_results.index):"+str(len(arcasHLA_results.index)))
            # new_row = {'Locus': locus, 'N': len(arcasHLA_input.index), 'arcasHLA - Call Rate': 1.00, 'arcasHLA - Accuracy': accuracy}
            new_row = pd.DataFrame([[locus, len(arcasHLA_results.index), call_rate, accuracy]],
                        columns=['Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'])
            arcasHLA_validation_table = pd.concat([arcasHLA_validation_table, new_row], ignore_index=True)
        return arcasHLA_validation_table

    def hla_la_validation(hla_la_results, ground_truth, hla_la_validation_table):
        # #HLA-LA
        # #loop over all loci
        loci = ['A', 'B', 'C', 'DQB1']
        for locus in loci:
            values = hla_la_results[locus]
            collected = 0
            samples_called = 0 #some samples in hla-la have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, value_in_hla_la) in enumerate(values):
                if pd.isnull(value_in_hla_la):
                    continue
                else:
                    sample_name = hla_la_results.loc[index,'sample']
                    values_in_truth = {}
                    alleles = {}
                    for (chr_index,chr_number) in enumerate([1, 2]):
                        #values in hla-la are as in the following: A*02:01:02;A*03:01:04/A*02:01:02;A*03:01:04
                        all_alleles=[]
                        first_allele = value_in_hla_la.split("/")[chr_index]
                        splitted_alleles = first_allele.split(";")
                        for splitted in splitted_alleles:
                            allele_no_locus = splitted.split("*")[1]
                            allele_in_hla_la = allele_no_locus.split(":")[0] + ":" +allele_no_locus.split(":")[1]
                            all_alleles.append(allele_in_hla_la)
                        alleles["{0}".format(chr_index)] = all_alleles

                        locus_name_in_truth_1 = "HLA-" + locus + " " + str(chr_number)
                        value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_name_in_truth_1].array[0] ##array[0] to reach what's inside the Series
                        ##split in case of alles that have e.g. 23:01/02/04 in the truth
                        tmp_alleles = []
                        if "/" in value_in_truth:
                            first_split = value_in_truth.split("/")
                            first_field = first_split[0].split(":")[0]
                            tmp_alleles.append(first_split[0])
                            for splitted in first_split[1:]:
                                tmp_alleles.append(first_field + ":" +splitted)
                            values_in_truth["{0}".format(chr_index)] = tmp_alleles
                            
                        elif value_in_truth.endswith("*"):
                            value_in_truth=value_in_truth.rstrip("*")
                            values_in_truth["{0}".format(chr_index)] = [value_in_truth]
                        else:
                            values_in_truth["{0}".format(chr_index)] = [value_in_truth]

                    values_in_truth_clone = values_in_truth
                    print("sample name: " ,sample_name)
                    print("hla-la prediction: " + ''.join(alleles))
                    print("truth_values: " + str(values_in_truth))
                    allele_present=0 #both alleles should be correct
                    for _, hla_la_alleles in alleles.items():
                        for allele in hla_la_alleles:
                            values_in_truth = values_in_truth_clone
                            for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
                                if allele in chr_values:
                                    allele_present+= 1
                                    #to stop homozygous alleles inaccurately match, we should remove the one that matches
                                    if len(values_in_truth) > 1:
                                        values_in_truth_clone.pop(chr_index, None)
                                    break
                            else:
                                continue
                            break
                    if allele_present==2: #both alleles should be correct
                        collected+=1
                    print("collected: " + str(collected))
                    samples_called+=1
            call_rate = samples_called/len(hla_la_results.index)
            accuracy = 100*collected/len(hla_la_results.index)
            print("len(hla_la_results.index):"+str(len(hla_la_results.index)))        
            # accuracy = 100*collected/(2*len(hla_la_results.index))
            # new_row = {'Locus': locus, 'N': len(hla_la_results.index), 'HLA-LA - Call Rate': 1.00, 'HLA-LA - Accuracy': accuracy}
            # hla_la_validation_table = hla_la_validation_table.append(new_row, ignore_index=True)
            new_row = pd.DataFrame([[locus, len(hla_la_results.index), call_rate, accuracy]],
                        columns=['Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'])
            hla_la_validation_table = pd.concat([hla_la_validation_table, new_row], ignore_index=True)
        return hla_la_validation_table
 
    def optitype_validation(optitype_results, ground_truth, optitype_validation_table):
        #optitype
        #loop over all loci
        loci = ['A', 'B', 'C']
        print("optitype")
        for locus in loci:
            values = optitype_results[locus]
            collected = 0
            samples_called = 0 #some samples in optitype have some loci that is uncalled, so we need to count the number of samples that are called here.
            for (index, value_in_optitype) in enumerate(values):
                if pd.isnull(value_in_optitype):
                    continue
                else:
                    sample_name = optitype_results.loc[index,'sample']
                    values_in_truth = {}
                    alleles = []
                    for (chr_index,chr_number) in enumerate([1, 2]):
                        #value in optitype 
                        first_allele = value_in_optitype.split("/")[chr_index]
                        allele_in_arcasHLA = first_allele.split("*")[1]
                        alleles.append(allele_in_arcasHLA)
                        
                        locus_name_in_truth_1 = "HLA-" + locus + " " + str(chr_number)
                        value_in_truth = ground_truth.loc[ground_truth['Run Accession'] == sample_name, locus_name_in_truth_1].array[0] ##array[0] to reach what's inside the Series
                        ##split in case of alles that have e.g. 23:01/02/04 in the truth
                        tmp_alleles = []
                        if "/" in value_in_truth:
                            first_split = value_in_truth.split("/")
                            first_field = first_split[0].split(":")[0]
                            tmp_alleles.append(first_split[0])
                            for splitted in first_split[1:]:
                                tmp_alleles.append(first_field + ":" +splitted)
                            values_in_truth["{0}".format(chr_index)] = tmp_alleles
                            
                        elif value_in_truth.endswith("*"):
                            value_in_truth=value_in_truth.rstrip("*")
                            values_in_truth["{0}".format(chr_index)] = [value_in_truth]
                        else:
                            values_in_truth["{0}".format(chr_index)] = [value_in_truth]

                    values_in_truth_clone = values_in_truth
                    print("sample name: " ,sample_name)
                    print("optitype prediction: " + ''.join(alleles))
                    print("truth_values: " + str(values_in_truth))
                    allele_present = 0
                    for allele in alleles: ##???
                        values_in_truth = values_in_truth_clone
                        for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
                            if allele in chr_values:
                                allele_present+= 1
                                #to stop homozygous alleles inaccurately match, we should remove the one that matches
                                if len(values_in_truth) > 1:
                                    values_in_truth_clone.pop(chr_index, None)
                                break
                    if allele_present==2:
                        collected+=1
                    print("collected: " + str(collected))
                    samples_called+=1
            call_rate = samples_called/len(optitype_results.index)
            accuracy = 100*collected/len(optitype_results.index)
            print("len(optitype_results.index):"+str(len(optitype_results.index)))        
            new_row = pd.DataFrame([[locus, len(optitype_results.index), call_rate, accuracy]],
                        columns=['Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'])
            optitype_validation_table = pd.concat([optitype_validation_table, new_row], ignore_index=True)
        return optitype_validation_table

    #low coverage samples

    ##initialize performance tables for tools 
    arcasHLA_validation_table_low = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'))
    hla_la_validation_table_low = pd.DataFrame(columns=('Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'))
    optitype_validation_table_low = pd.DataFrame(columns=('Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'))

    #merge tables to get the final table
    arcasHLA_validation_table_low = arcashla_validation(arcasHLA_input_low, ground_truth, arcasHLA_validation_table_low)
    first_merge = pd.merge(orthanq_input_low, arcasHLA_validation_table_low, how='left', on=['Locus', 'N'])
    hla_la_validation_table_low = hla_la_validation(hla_la_input_low, ground_truth, hla_la_validation_table_low)
    second_merge = pd.merge(first_merge, hla_la_validation_table_low, how='left', on=['Locus', 'N'])
    optitype_validation_table_low = optitype_validation(optitype_input_low, ground_truth, optitype_validation_table_low)
    low_final_table = pd.merge(second_merge, optitype_validation_table_low, how='left', on=['Locus', 'N'])

    #high coverage samples

    ##initialize performance tables for tools 
    arcasHLA_validation_table_high = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'))
    hla_la_validation_table_high = pd.DataFrame(columns=('Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'))
    optitype_validation_table_high = pd.DataFrame(columns=('Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'))

    #merge tables to get the final table
    arcasHLA_validation_table_high = arcashla_validation(arcasHLA_input_high, ground_truth, arcasHLA_validation_table_high)
    first_merge = pd.merge(orthanq_input_high, arcasHLA_validation_table_high, how='left', on=['Locus', 'N'])
    hla_la_validation_table_high = hla_la_validation(hla_la_input_high, ground_truth, hla_la_validation_table_high)
    second_merge = pd.merge(first_merge, hla_la_validation_table_high, how='left', on=['Locus', 'N'])
    optitype_validation_table_high = optitype_validation(optitype_input_high, ground_truth, optitype_validation_table_high)
    high_final_table = pd.merge(second_merge, optitype_validation_table_high, how='left', on=['Locus', 'N'])

    #all samples

    ##initialize performance tables for tools 
    arcasHLA_validation_table_all = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'))
    hla_la_validation_table_all = pd.DataFrame(columns=('Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'))
    optitype_validation_table_all = pd.DataFrame(columns=('Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'))

    #merge tables to get the final table
    arcasHLA_validation_table_all = arcashla_validation(arcasHLA_input, ground_truth, arcasHLA_validation_table_all)
    first_merge = pd.merge(orthanq_input_all, arcasHLA_validation_table_all, how='left', on=['Locus', 'N'])
    hla_la_validation_table_all = hla_la_validation(hla_la_input, ground_truth, hla_la_validation_table_all)
    second_merge = pd.merge(first_merge, hla_la_validation_table_all, how='left', on=['Locus', 'N'])
    optitype_validation_table_all = optitype_validation(optitype_input, ground_truth, optitype_validation_table_all)
    all_final_table = pd.merge(second_merge, optitype_validation_table_all, how='left', on=['Locus', 'N'])

    #write all tables to csv 

    low_final_table.to_csv(
        snakemake.output.validation_low, sep="\t", index=False, header=True
    )

    high_final_table.to_csv(
        snakemake.output.validation_high, sep="\t", index=False, header=True
    )

    all_final_table.to_csv(
        snakemake.output.validation_all, sep="\t", index=False, header=True
    )