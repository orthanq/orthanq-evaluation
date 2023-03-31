import pandas as pd

configfile: "config/config.yaml"
n_reads = config["n_reads"]
alleles = pd.read_csv(config["hla"], sep ="\t")
samples = pd.read_csv(config["samples"], sep="\t").set_index(
    ["sample_name"], drop=False
)
#loci = ["A","DPA1","DRB4","V","B","DPB1","DRB5","W","C","DQA1", "E", "DQA2", "F", "S", "DMA", "DQB1", "G", "TAP1", "DMB", "DRA", "HFE", "TAP2", "DOA", "DRB1", "T", "DOB", "DRB3", "MICA", "U"]
loci = ["A","B","C","DQA1","DQB1"]
loci = ["A","B","C"]

# input function for simulation sample
def create_sample(): #n: number of samples in the end, k: number of fractions
    alleles['sample_name'] = "Sample" + "_" + '_'.join(str(row['hla']) + "-" + str(row['fraction']) for _, row in alleles.iterrows())
    # print(alleles)
    alleles['num_reads'] = alleles['fraction']*n_reads #assume the mixture will have 2000 reads
    alleles['num_reads'] = alleles['num_reads'].astype(int).astype('str')
    del alleles['fraction'] 
    return alleles

simulated_sample = create_sample()
# print(simulated_sample)

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

