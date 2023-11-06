# import altair with an abbreviated alias
# %%

import altair as alt
import pandas as pd

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    #read the table and drop DQA1
    # orthanq_table = pd.read_csv("orthanq_tp.tsv", sep = "\t")
    orthanq_table = pd.read_csv(snakemake.input.tp_fp_table, sep = "\t")

    orthanq_table.drop(orthanq_table[orthanq_table.Locus == "DQA1"].index, inplace=True)

    print(orthanq_table[orthanq_table.TP == 0])

    #create True Positive and False Positive dataframes
    fp=orthanq_table[orthanq_table.TP == 0]
    tp=orthanq_table[orthanq_table.TP == 1]

    #histogram for FP
    hist_fp = alt.Chart(fp).transform_calculate(FP="datum.TP == 0").transform_aggregate(fp_count="sum(FP)",
                                        groupby=["Best_Density"]).mark_bar(color='red', 
                                    orient='horizontal',
                                    height=1,
                                    baseline='middle',
                                    binSpacing=0,
                                    tooltip=True).encode(
                                            alt.X('fp_count:Q',
                                                    axis=alt.Axis(values=[0, 1, 2, 3],
                                                                tickMinStep=1,
                                                                title="false predictions"),
                                                    scale=alt.Scale(reverse=True)),
                                            alt.Y('Best_Density:Q',
                                                axis=alt.Axis(orient="right",
                                                                tickCount=10,
                                                                tickOffset=0),
                                                scale=alt.Scale(domain=[0,1]),
                                                bin="binned",
                                                title=None)
                                                )
    hist_fp

    # #histogram for TP
    hist_tp = alt.Chart(tp).mark_bar(color='black', 
                                    tooltip=True).encode(
                                        alt.X('count():Q', 
                                            title="correct predictions"),
                                        alt.Y("Best_Density:Q",
                                            bin=alt.Bin(maxbins=20), 
                                                title=None, 
                                                scale=alt.Scale(domain=[0,1]), 
                                                axis=alt.Axis(labels=False))
    )
    hist_tp

    chart = alt.hconcat(hist_fp, 
                        hist_tp).resolve_axis(y='shared').configure(
        concat=alt.CompositionConfig(spacing=0)).configure_title().configure_view(continuousHeight=300, 
                                                continuousWidth=300).configure_axis(grid=False)

    # chart.save('tp_fp_densities.json')
    chart.save(snakemake.output.tp_fp_plot)
    # chart
# %%
