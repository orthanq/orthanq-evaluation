import pandas as pd

#read tables
ground_truth = pd.read_csv(snakemake.input.ground_truth, sep = "\t")
orthanq_input = pd.read_csv(snakemake.input.orthanq, sep = "\t")
arcasHLA_input = pd.read_csv(snakemake.input.arcasHLA, sep = "\t")
hla_la_input = pd.read_csv(snakemake.input.hla_la, sep = "\t")

##initialize performance tables for tools
orthanq_validation_table = pd.DataFrame(columns=('Locus', 'N', 'Orthanq - Call Rate', 'Orthanq - Accuracy - 1', 'Orthanq - Accuracy - 2'))
arcasHLA_validation_table = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA - Call Rate', 'arcasHLA - Accuracy - 1', 'arcasHLA - Accuracy - 2'))
hla_la_validation_table = pd.DataFrame(columns=('Locus', 'N', 'HLA-LA - Call Rate', 'HLA-LA - Accuracy - 1', 'HLA-LA - Accuracy - 2'))

#ORTHANQ
#loop over all loci
loci = ['A', 'B', 'C', 'DQB1']
for locus in loci:
    values = orthanq_input[locus]
    collected_1 = 0
    collected_2 = 0
    for (index, value_in_orthanq) in enumerate(values):
        sample_name = orthanq_input.loc[index,'sample']
        values_in_truth = {}
        alleles = []
        for (chr_index,chr_number) in enumerate([1, 2]):
            #value in orthanq 
            first_allele = value_in_orthanq.split("/")[chr_index]
            first_allele_no_locus = first_allele.split("*")[1]
            allele_in_orthanq = first_allele_no_locus.split(":")[0] + ":" +first_allele_no_locus.split(":")[1]
            alleles.append(allele_in_orthanq)
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
        alleles_clone = alleles
        for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
            alleles = alleles_clone
            for allele in alleles:
                if allele in chr_values:
                    if chr_index == 0:
                        collected_1+= 1
                    if chr_index == 1:
                        collected_2+= 1
                    #to stop homozygous alleles inaccurately match, we should remove the one that matches
                    if len(values_in_truth) > 1:
                        alleles_clone.remove(allele)
                    break
    accuracy_1 = 100*collected_1/(len(orthanq_input.index))
    accuracy_2 = 100*collected_2/(len(orthanq_input.index))
    new_row = {'Locus': locus, 'N': len(orthanq_input.index), 'Orthanq - Call Rate': 1.00, 'Orthanq - Accuracy - 1': accuracy_1, 'Orthanq - Accuracy - 2': accuracy_2}
    orthanq_validation_table = orthanq_validation_table.append(new_row, ignore_index=True)

#arcasHLA
#loop over all loci
loci = ['A', 'B', 'C', 'DQB1']
for locus in loci:
    values = arcasHLA_input[locus]
    collected_1 = 0
    collected_2 = 0
    for (index, value_in_arcasHLA) in enumerate(values):
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
        print(alleles)
        print(values_in_truth)
        alleles_clone = alleles
        for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
            alleles = alleles_clone
            for allele in alleles:
                if allele in chr_values:
                    if chr_index == 0:
                        collected_1+= 1
                    if chr_index == 1:
                        collected_2+= 1
                    #to stop homozygous alleles inaccurately match, we should remove the one that matches
                    if len(values_in_truth) > 1:
                        alleles_clone.remove(allele)
                    break
        print("collected_1:" + str(collected_1))
        print("collected_2:" + str(collected_2))
    accuracy_1 = 100*collected_1/(len(arcasHLA_input.index))
    accuracy_2 = 100*collected_2/(len(arcasHLA_input.index))
    new_row = {'Locus': locus, 'N': len(arcasHLA_input.index), 'arcasHLA - Call Rate': 1.00, 'arcasHLA - Accuracy - 1': accuracy_1, 'arcasHLA - Accuracy - 2': accuracy_2}
    arcasHLA_validation_table = arcasHLA_validation_table.append(new_row, ignore_index=True)

first_merge = pd.merge(orthanq_validation_table, arcasHLA_validation_table, how='left', on=['Locus', 'N'])

# #HLA-LA
# #loop over all loci
loci = ['A', 'B', 'C', 'DQB1']
for locus in loci:
    values = hla_la_input[locus]
    collected_1 = 0
    collected_2 = 0
    for (index, value_in_hla_la) in enumerate(values):
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
        alleles_clone = alleles
        for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
            alleles = alleles_clone
            for (allele_index,allele_list) in alleles.items():
                for allele in allele_list:
                    if allele in chr_values:
                        if chr_index == 0:
                            collected_1+= 1
                        if chr_index == 1:
                            collected_2+= 1
                        #to stop homozygous alleles inaccurately match, we should remove the one that matches
                        if len(values_in_truth) > 1:
                            alleles_clone.pop(allele_index, None)
                        break
                else:
                    continue  # only executed if the inner loop did NOT break
                break  # only executed if the inner loop DID break
    accuracy_1 = 100*collected_1/(len(hla_la_input.index))
    accuracy_2 = 100*collected_2/(len(hla_la_input.index))
    new_row = {'Locus': locus, 'N': len(hla_la_input.index), 'HLA-LA - Call Rate': 1.00, 'HLA-LA - Accuracy - 1': accuracy_1, 'HLA-LA - Accuracy - 2': accuracy_2}
    hla_la_validation_table = hla_la_validation_table.append(new_row, ignore_index=True)

final_table = pd.merge(first_merge, hla_la_validation_table, how='left', on=['Locus', 'N'])

final_table.to_csv(
    snakemake.output.validation, sep="\t", index=False, header=True
)