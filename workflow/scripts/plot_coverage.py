#coverage plot excluding giab sample

import os
import glob 
import pandas as pd
import altair as alt

# pd.set_option('display.max_rows', 500)


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #import read cound and tp fp tables first
    # tp_fp_table = pd.read_csv("tp_fp_table_all.tsv", sep="\t")
    tp_fp_table = pd.read_csv(snakemake.input.tp_fp_table, sep="\t")

    #read counts global
    # read_counts_global = pd.read_csv("merged_sample_sheet_w_read_info_2.csv", sep=",")
    read_counts_global = pd.read_csv(snakemake.input.read_counts_global, sep=",")

    #select only the required columns from read count table
    read_counts_global = read_counts_global[["Run Accession", "Read Length", "Sequencing method"]]

    #read counts per locus
    # read_counts_per_locus = pd.read_csv("read_counts.csv", sep=",")
    read_counts_per_locus = pd.read_csv(snakemake.input.read_counts_per_locus, sep=",")

    #renme locus column
    read_counts_per_locus = read_counts_per_locus.rename(columns={'locus': 'Locus'}) 

    #read per locus gene lengths into one to calculate locus sizes
    # A_region = pd.read_csv("genes/A.bed", sep="\t", header=None)
    # B_region = pd.read_csv("genes/B.bed", sep="\t", header=None)
    # C_region = pd.read_csv("genes/C.bed", sep="\t", header=None)
    # DQB1_region = pd.read_csv("genes/DQB1.bed", sep="\t", header=None)

    A_region = pd.read_csv(snakemake.input.A_coord, sep="\t", header=None)
    B_region = pd.read_csv(snakemake.input.B_coord, sep="\t", header=None)
    C_region = pd.read_csv(snakemake.input.C_coord, sep="\t", header=None)
    DQB1_region = pd.read_csv(snakemake.input.DQB1_coord, sep="\t", header=None)

    A_length = A_region.iloc[0, 2]-A_region.iloc[0, 1]
    B_length = B_region.iloc[0, 2]-B_region.iloc[0, 1]
    C_length = C_region.iloc[0, 2]-C_region.iloc[0, 1]
    DQB1_length = DQB1_region.iloc[0, 2]-DQB1_region.iloc[0, 1]

    #get cds lengths of HLA genes (retrieved from Biomart GRCh38.p14)
    A_cds_length = 1098
    B_cds_length = 1089
    C_cds_length = 1020
    DQB1_cds_length = 675

    #merge global read count and per locus read counts file to add read length and sequencing method information 
    read_counts_all_information = read_counts_per_locus.merge(read_counts_global, how = "inner", left_on='sample', right_on='Run Accession')
    print(read_counts_all_information)

    #combine true & uncalled and false & uncalled into called: no
    tp_fp_table["called"] = "no input"
    for index, row in tp_fp_table.iterrows():
        if (row["Prediction"] == "true & uncalled") or (row["Prediction"] == "false & uncalled"):
            tp_fp_table.loc[index, "called"] = "no"
        if (row["Prediction"] == "TP") or (row["Prediction"] == "FP"):
            tp_fp_table.loc[index, "called"] = "yes"

    #remove DQA1
    tp_fp_table.drop(tp_fp_table[tp_fp_table.Locus == "DQA1"].index, inplace=True)

    #then merge tp_fp with read count table
    merged = tp_fp_table.merge(read_counts_all_information, how = "inner", left_on=['Sample', 'Locus'], right_on=['Run Accession', 'Locus'])
    print("mergeddd")
    print(merged)

    #calculate coverage
    merged["coverage"] = 0

    for index, row in merged.iterrows():
        print("Sample", row["Sample"])
        print("Locus", row["Locus"])
        print("Sequencing method,", row["Sequencing method"])
        print("Read Length,", row["Read Length"])
        print("read_count,", row["read_count"])
        if row["Sequencing method"] == "WGS":
            if row["Locus"] == "A":
                coverage = row["read_count"]*(2*row["Read Length"])/A_length
            elif row["Locus"] == "B":
                coverage = row["read_count"]*(2*row["Read Length"])/B_length
            elif row["Locus"] == "C":
                coverage = row["read_count"]*(2*row["Read Length"])/C_length
            elif row["Locus"] == "DQB1":
                coverage = row["read_count"]*(2*row["Read Length"])/DQB1_length        
        elif row["Sequencing method"] == "WES":
            if row["Locus"] == "A":
                coverage = row["read_count"]*(2*row["Read Length"])/A_cds_length
            if row["Locus"] == "B":
                coverage = row["read_count"]*(2*row["Read Length"])/B_cds_length
            if row["Locus"] == "C":
                coverage = row["read_count"]*(2*row["Read Length"])/C_cds_length
            if row["Locus"] == "DQB1":
                coverage = row["read_count"]*(2*row["Read Length"])/DQB1_cds_length

        print("coverage: ", coverage)
        merged.loc[index, "coverage"] = coverage

    print("merged updated")
    print(merged)
    def plot_layered(df):
        coverage_boxplot = alt.Chart(df).mark_boxplot(ticks=True).encode(
            alt.X("Locus:N"),
            alt.Y("coverage:Q"),
            alt.Color("called:N"),
            alt.XOffset("called:N"),
            )
        coverage_scatterplot = alt.Chart(df).mark_circle(size=9,opacity=0.2).transform_calculate(jitter="random()").encode(
                x='Locus:N',
                y='coverage:Q',
                xOffset='jitter:Q',
            color=alt.value("black"),
            )
        coverage_plot = alt.layer(coverage_boxplot,coverage_scatterplot)     
        return coverage_plot

    A_locus = merged[merged["Locus"] == "A"]
    B_locus = merged[merged["Locus"] == "B"]
    C_locus = merged[merged["Locus"] == "C"]
    DQB1_locus = merged[merged["Locus"] == "DQB1"]

    plot = plot_layered(A_locus) | plot_layered(B_locus) | plot_layered(C_locus) | plot_layered(DQB1_locus)

    # plot.save("boxplot.json")
    plot.save(snakemake.output.boxplot_json)

