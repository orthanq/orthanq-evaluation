# import altair with an abbreviated alias
# %%

import altair as alt
import pandas as pd

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    #read the table
    # orthanq_table = pd.read_csv("tp_fp_table_all.tsv", sep = "\t")
    orthanq_table = pd.read_csv(snakemake.input.tp_fp_table, sep = "\t")

    #remove DQA1
    orthanq_table.drop(orthanq_table[orthanq_table.Locus == "DQA1"].index, inplace=True)

    #stratify samples
    false_predictions=orthanq_table[(orthanq_table.Prediction == "FP") | (orthanq_table.Prediction == "false & uncalled")]
    true_predictions=orthanq_table[(orthanq_table.Prediction == "TP") | (orthanq_table.Prediction == "true & uncalled")]

    #assign numerical values for TP (1), FP(0), no call(2) (conditions do not work with the strings (don't know why))
    true_predictions.loc[true_predictions.Prediction == "TP", "Prediction"] = 1
    true_predictions.loc[true_predictions.Prediction == "true & uncalled", "Prediction"] = 2
    false_predictions.loc[false_predictions.Prediction == "FP", "Prediction"] = 3
    false_predictions.loc[false_predictions.Prediction == "false & uncalled", "Prediction"] = 4

    print(false_predictions.to_string())
    print(true_predictions.to_string())

    #histogram for FP
    domain = [3, 4]
    range_ = ["red", "gray"]
    hist_fp = alt.Chart(false_predictions).transform_calculate(false_pred=("datum.Prediction == 3 | datum.Prediction == 4")).mark_bar(
                                    orient='horizontal',
                                    height=1,
                                    baseline='middle',
                                    binSpacing=0,
                                    tooltip=True).encode(
                                            x=alt.X('count():Q',
                                                    axis=alt.Axis(values=[0, 1, 2, 3, 4],
                                                                tickMinStep=1,
                                                                title="false predictions"),
                                                    scale=alt.Scale(reverse=True)),
                                            y= alt.Y('Best_Density:Q',
                                                axis=alt.Axis(orient="right",
                                                                tickCount=10,
                                                                tickOffset=0),
                                                scale=alt.Scale(domain=[0,1]),
                                                bin="binned",
                                                title=None),
                                            color=alt.Color('Prediction:N',scale=alt.Scale(domain=domain, range=range_), legend=None)
                                        )

                                        


    # #histogram for TP
    domain = [1, 2]
    range_ = ["black", "gray"]
    hist_tp = alt.Chart(true_predictions).transform_calculate(true_pred=("datum.Prediction == 1 | datum.Prediction == 2")).mark_bar(
                                    tooltip=True).encode(
                                        x=alt.X('count():Q', 
                                            title="true predictions"),
                                        y=alt.Y("Best_Density:Q",
                                            bin=alt.Bin(maxbins=20), 
                                                title=None, 
                                                scale=alt.Scale(domain=[0,1]), 
                                                axis=alt.Axis(labels=False)),
                                        color=alt.Color('Prediction:N',scale=alt.Scale(domain=domain, range=range_), legend=None)
                                        )
    hist_tp

    chart = alt.hconcat(hist_fp, 
                        hist_tp).resolve_axis(y='shared').resolve_scale(color='independent').configure(
        concat=alt.CompositionConfig(spacing=0)).configure_title().configure_view(continuousHeight=300, 
                                                continuousWidth=300).configure_axis(grid=False)

    chart.save('tp_fp_local.json')
    chart.save(snakemake.output.tp_fp_plot)
# %%
