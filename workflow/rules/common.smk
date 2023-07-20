# This script contains functions for creating simulated samples and collecting fastq inputs.

import pandas as pd

configfile: "config/config.yaml"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# Simulation of samples from HLA alleles and haplotype quantification for samples that are simulated #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #                                                                                                   #

# the HLA loci that have our focus in this workflow.
loci = ["A","B","C","DQA1","DQB1"]

# reading alleles to be simulated
n_reads = config["n_reads"]
alleles = pd.read_csv(config["hla"], sep ="\t")

# input function to create a simulated sample
def create_sample(): #n: number of samples in the end, k: number of fractions
    alleles['sample_name'] = "Sample" + "_" + '_'.join(str(row['hla']) + "-" + str(row['fraction']) for _, row in alleles.iterrows())
    # print(alleles)
    alleles['num_reads'] = alleles['fraction']*n_reads #assume the mixture will have 2000 reads
    alleles['num_reads'] = alleles['num_reads'].astype(int).astype('str')
    del alleles['fraction'] 
    return alleles

# the simulated sample is created below
simulated_sample = create_sample()
# print(simulated_sample)


# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Direct haplotype quantification for given samples #
# # # # # # # # # # # # # # # # # # # # # # # # # # #

# reading samples
samples = pd.read_csv(config["samples"], sep="\t").set_index(
    ["sample_name"], drop=False
)

# input function to retrieve fastq samples
def get_fastq_input(wildcards):
    if config["simulation"] == False:
        sample = samples.loc[wildcards.sample]
        return [sample["fq1"], sample["fq2"]]
    else:
        sample = wildcards.sample
        simulated = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"]
        return simulated
