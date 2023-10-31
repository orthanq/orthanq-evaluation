import pandas as pd
import os 

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #read tables
    ground_truth = pd.read_csv(snakemake.input.ground_truth, sep = "\t")
    orthanq_input = pd.read_csv(snakemake.input.orthanq, sep = "\t")
    arcasHLA_input = pd.read_csv(snakemake.input.arcasHLA, sep = "\t")
    hla_la_input = pd.read_csv(snakemake.input.hla_la, sep = "\t")
    optitype_input = pd.read_csv(snakemake.input.optitype, sep = "\t")

    ##initialize performance tables for tools
    arcasHLA_validation_table = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'))
    hla_la_validation_table = pd.DataFrame(columns=('Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'))
    optitype_validation_table = pd.DataFrame(columns=('Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'))

    #arcasHLA
    #loop over all loci
    loci = ['A', 'B', 'C', 'DQB1']
    for locus in loci:
        values = arcasHLA_input[locus]
        collected = 0
        samples_called = 0 #some samples in arcashla have some loci that is uncalled, so we need to count the number of samples that are called here.
        for (index, value_in_arcasHLA) in enumerate(values):
            if pd.isnull(value_in_arcasHLA):
                continue
            else:
                sample_name = arcasHLA_input.loc[index,'sample']
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
        call_rate = samples_called/len(arcasHLA_input.index)
        accuracy = 100*collected/samples_called
        # new_row = {'Locus': locus, 'N': len(arcasHLA_input.index), 'arcasHLA - Call Rate': 1.00, 'arcasHLA - Accuracy': accuracy}
        new_row = pd.DataFrame([[locus, len(arcasHLA_input.index), call_rate, accuracy]],
                    columns=['Locus', 'N', 'arcasHLA_Call_Rate', 'arcasHLA_Accuracy'])
        arcasHLA_validation_table = pd.concat([arcasHLA_validation_table, new_row], ignore_index=True)

    first_merge = pd.merge(orthanq_input, arcasHLA_validation_table, how='left', on=['Locus', 'N'])

    # #HLA-LA
    # #loop over all loci
    loci = ['A', 'B', 'C', 'DQB1']
    for locus in loci:
        values = hla_la_input[locus]
        collected = 0
        samples_called = 0 #some samples in arcashla have some loci that is uncalled, so we need to count the number of samples that are called here.
        for (index, value_in_hla_la) in enumerate(values):
            if pd.isnull(value_in_arcasHLA):
                continue
            else:
                sample_name = hla_la_input.loc[index,'sample']
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
        call_rate = samples_called/len(arcasHLA_input.index)
        accuracy = 100*collected/samples_called
        # accuracy = 100*collected/(2*len(hla_la_input.index))
        # new_row = {'Locus': locus, 'N': len(hla_la_input.index), 'HLA-LA - Call Rate': 1.00, 'HLA-LA - Accuracy': accuracy}
        # hla_la_validation_table = hla_la_validation_table.append(new_row, ignore_index=True)
        new_row = pd.DataFrame([[locus, len(hla_la_input.index), call_rate, accuracy]],
                    columns=['Locus', 'N', 'HLA_LA_Call_Rate', 'HLA_LA_Accuracy'])
        hla_la_validation_table = pd.concat([hla_la_validation_table, new_row], ignore_index=True)

    second_merge = pd.merge(first_merge, hla_la_validation_table, how='left', on=['Locus', 'N'])

    #optitype
    #loop over all loci
    loci = ['A', 'B', 'C']
    print("optitype")
    for locus in loci:
        values = optitype_input[locus]
        collected = 0
        samples_called = 0 #some samples in arcashla have some loci that is uncalled, so we need to count the number of samples that are called here.
        for (index, value_in_optitype) in enumerate(values):
            if pd.isnull(value_in_optitype):
                continue
            else:
                sample_name = optitype_input.loc[index,'sample']
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
        call_rate = samples_called/len(optitype_input.index)
        accuracy = 100*collected/samples_called
        new_row = pd.DataFrame([[locus, len(optitype_input.index), call_rate, accuracy]],
                    columns=['Locus', 'N', 'Optitype_Call_Rate', 'Optitype_Accuracy'])
        optitype_validation_table = pd.concat([optitype_validation_table, new_row], ignore_index=True)

    final_table = pd.merge(second_merge, optitype_validation_table, how='left', on=['Locus', 'N'])

    final_table.to_csv(
        snakemake.output.validation, sep="\t", index=False, header=True
    )