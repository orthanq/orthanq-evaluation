import pandas as pd

samples = pd.read_csv("../config/samples_final.tsv", sep="\t")
sample_infos = pd.read_csv("../resources/ground_truth/merged_sample_sheet_w_read_info.csv")

below_threshold = sample_infos[sample_infos["Read Count 1"] <= 10000000]
below_threshold = below_threshold.rename(columns={"Run Accession": "sample_name"})

final_below_threshold=pd.merge(samples, below_threshold, on="sample_name")
final_below_threshold=final_below_threshold.iloc[:,0:3]

above_threshold = sample_infos[sample_infos["Read Count 1"] > 10000000]
above_threshold = above_threshold.rename(columns={"Run Accession": "sample_name"})

final_above_threshold=pd.merge(samples, above_threshold, on="sample_name")
final_above_threshold=final_above_threshold.iloc[:,0:3]

final_below_threshold.to_csv("../config/samples_below_threshold.tsv", index=False, sep="\t")
final_above_threshold.to_csv("../config/samples_above_threshold.tsv", index=False, sep="\t")
