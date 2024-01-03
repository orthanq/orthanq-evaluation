import pandas as pd

configfile: "config/config.yaml"
n_reads = config["n_reads"]

#in case of simulation:
lineages = pd.read_csv(config["viruses"], sep ="\t")

#the following line needs to be implemented for the analysis part of the workflow
# samples = pd.read_csv(config["samples"], sep="\t").set_index(
#     ["sample_name"], drop=False
# )

# loci = ["A","B","C","DQB1"]

# input function for simulation sample
def create_sample(): #n: number of samples in the end, k: number of fractions
    lineages['sample_name'] = "Sample" + "_" + '_'.join(str(row['lineage']) + "-" + str(row['fraction']) for _, row in lineages.iterrows())
    # print(lineages)
    lineages['num_reads'] = lineages['fraction']*n_reads #assume the mixture will have 2000 reads
    lineages['num_reads'] = lineages['num_reads'].astype(int).astype('str')
    del lineages['fraction'] 
    return lineages

simulated_sample = create_sample()
print(simulated_sample)

# input function to retrieve fastq samples
def get_fastq_input(wildcards):
    if config["simulation"] == False:
        sample = samples.loc[wildcards.sample]
        return [sample["fq1"], sample["fq2"]]
    else:
        sample = wildcards.sample
        simulated = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"]
        return simulated

def get_fastq_input_1(wildcards):
    if config["simulation"] == False:
        sample = samples.loc[wildcards.sample]
        return sample["fq1"]
    else:
        sample = wildcards.sample
        simulated = "results/mixed/{sample}_1.fq"
        return simulated

def get_fastq_input_2(wildcards):
    if config["simulation"] == False:
        sample = samples.loc[wildcards.sample]
        return sample["fq2"]
    else:
        sample = wildcards.sample
        simulated = "results/mixed/{sample}_2.fq"
        return simulated

