import pandas as pd

orthanq_input = snakemake.input.orthanq
hla_la_input = snakemake.input.hla_la
arcasHLA_input = snakemake.input.arcasHLA
ground_truth = snakemake.input.ground_truth


#initialize the dataframe
validation_table = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA - Call Rate', 'arcasHLA - Accuracy', 'HLA-LA - Call Rate','HLA-LA - Accuracy', 'Orthanq - Call Rate', 'Orthanq - Accuracy'))

#read tables
ground_truth = pd.read_csv(ground_truth, sep = "\t")
orthanq_input = pd.read_csv(orthanq_input, sep = "\t")
arcasHLA_input = pd.read_csv(arcasHLA_input, sep = "\t")
hla_la_input = pd.read_csv(hla_la_input, sep = "\t")

##initialize performance tables for tools
orthanq_validation_table = pd.DataFrame(columns=('Locus', 'N', 'Orthanq - Call Rate', 'Orthanq - Accuracy'))
arcasHLA_validation_table = pd.DataFrame(columns=('Locus', 'N', 'arcasHLA - Call Rate', 'arcasHLA - Accuracy'))
hla_la_validation_table = pd.DataFrame(columns=('Locus', 'N', 'HLA-LA - Call Rate', 'HLA-LA - Accuracy'))

#ORTHANQ
#loop over all loci
loci = ['A', 'B', 'C', 'DQB1']
for locus in loci:
    values = orthanq_input[locus]
    collected = 0
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

        values_in_truth_clone = values_in_truth
        for allele in alleles: ##???
            values_in_truth = values_in_truth_clone
            for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
                if allele in chr_values:
                    collected+= 1
                    #to stop homozygous alleles inaccurately match, we should remove the one that matches
                    if len(values_in_truth) > 1:
                        values_in_truth_clone.pop(chr_index, None)
                    break
    accuracy = 100*collected/(2*len(orthanq_input.index))
    new_row = {'Locus': locus, 'N': len(orthanq_input.index), 'Orthanq - Call Rate': 1.00, 'Orthanq - Accuracy': accuracy}
    orthanq_validation_table = orthanq_validation_table.append(new_row, ignore_index=True)

#arcasHLA
#loop over all loci
loci = ['A', 'B', 'C', 'DQB1']
for locus in loci:
    values = arcasHLA_input[locus]
    collected = 0
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

        values_in_truth_clone = values_in_truth
        for allele in alleles: ##???
            values_in_truth = values_in_truth_clone
            for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
                if allele in chr_values:
                    collected+= 1
                    #to stop homozygous alleles inaccurately match, we should remove the one that matches
                    if len(values_in_truth) > 1:
                        values_in_truth_clone.pop(chr_index, None)
                    break
    accuracy = 100*collected/(2*len(arcasHLA_input.index))
    new_row = {'Locus': locus, 'N': len(arcasHLA_input.index), 'arcasHLA - Call Rate': 1.00, 'arcasHLA - Accuracy': accuracy}
    arcasHLA_validation_table = arcasHLA_validation_table.append(new_row, ignore_index=True)

first_merge = pd.merge(orthanq_validation_table, arcasHLA_validation_table, how='left', on=['Locus', 'N'])

# #HLA-LA
# #loop over all loci
loci = ['A', 'B', 'C', 'DQB1']
for locus in loci:
    values = hla_la_input[locus]
    collected = 0
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

        values_in_truth_clone = values_in_truth
        for _, hla_la_alleles in alleles.items():
            for allele in hla_la_alleles:
                values_in_truth = values_in_truth_clone
                for chr_index, (_,chr_values) in enumerate(values_in_truth.items()):
                    if allele in chr_values:
                        collected+= 1
                        #to stop homozygous alleles inaccurately match, we should remove the one that matches
                        if len(values_in_truth) > 1:
                            values_in_truth_clone.pop(chr_index, None)
                        break
                else:
                    continue
                break
    accuracy = 100*collected/(2*len(hla_la_input.index))
    new_row = {'Locus': locus, 'N': len(hla_la_input.index), 'HLA-LA - Call Rate': 1.00, 'HLA-LA - Accuracy': accuracy}
    hla_la_validation_table = hla_la_validation_table.append(new_row, ignore_index=True)

final_table = pd.merge(first_merge, hla_la_validation_table, how='left', on=['Locus', 'N'])

final_table.to_csv(
    snakemake.output.validation, sep="\t", index=False, header=True
)