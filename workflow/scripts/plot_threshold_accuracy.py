import pandas as pd
import altair as alt

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    def take_first(value):
        splitted=value.split(' ')
        return splitted[0]

    def convert_decimal(value):
        value=value/100
        return value

    # table = pd.read_csv("threshold_results.tsv", sep = "\t")
    table = pd.read_csv(snakemake.input.threshold_results, sep = "\t")

    print(table)
    print(table.dtypes)

    #first convert to string for string operation in the following step
    table["accuracy"] = table["accuracy"].astype(str)

    #take only the call rate and accuracies as float 
    table["call_rate"] = table["call_rate"].apply(take_first).astype(float)
    table["accuracy"] = table["accuracy"].apply(take_first).astype(float)

    #convert percentages to decimals
    table["accuracy"] = table["accuracy"].apply(convert_decimal)

    #separate no rows with no threshold on the threshold_density
    no_threshold = table[table["threshold_density"] == "no threshold"]
    rest = table[table["threshold_density"] != "no threshold"]

    #convert threshold densities for rest to float
    rest["threshold_density"] = rest["threshold_density"].astype(float)
    print(rest.dtypes)

    #remove 1.0 rows
    rest.drop(rest[rest.threshold_density == 1.0].index, inplace=True)

    def line_plot(table, locus, column_name):
        chart = alt.Chart(table).mark_line(point=True, tooltip=True).encode(
                                    x=alt.X('threshold_density',
                                            axis=alt.Axis(
                                                        title="min density", labels=True),
                                            scale=alt.Scale(domain=[0.0,0.9])),
                                    y= alt.Y(column_name, scale=alt.Scale(domain=[0.6,1.0])),
                                    color=alt.Color("threshold_haplotypes:N",legend=alt.Legend(title="max ambiguous solutions") )).transform_filter(
                                        alt.datum.locus == locus).properties(title=locus)
        return chart

    #add density values to the no_threshold table to make a lineplot
    thresholds_density = (x * 0.1 for x in range(0, 10))
    new_dict = {"threshold_density": [], "threshold_haplotypes": [], "locus": [], "call_rate": [], "accuracy": []}
    for density_i in thresholds_density:
        for i, row in enumerate(no_threshold.itertuples(), 1):
            new_dict['threshold_density'].append(density_i)
            new_dict['threshold_haplotypes'].append("no threshold")
            new_dict['locus'].append(row.locus)
            new_dict['call_rate'].append(row.call_rate)
            new_dict['accuracy'].append(row.accuracy)

    new_df = pd.DataFrame(new_dict)
    print("new_df", new_df)

    #plots 
    #A - call rate
    rest_call_rate_A = line_plot(rest, "A", "call_rate")
    no_threshold_call_rate_A = line_plot(new_df, "A", "call_rate")
    A_call_rate = alt.layer(rest_call_rate_A, no_threshold_call_rate_A)

    #A - accuracy
    rest_accuracy_A = line_plot(rest, "A", "accuracy")
    no_threshold_accuracy_A = line_plot(new_df, "A", "accuracy")
    A_accuracy = alt.layer(rest_accuracy_A, no_threshold_accuracy_A)

    #B - call rate
    rest_call_rate_B = line_plot(rest, "B", "call_rate")
    no_threshold_call_rate_B = line_plot(new_df, "B", "call_rate")
    B_call_rate = alt.layer(rest_call_rate_B, no_threshold_call_rate_B)

    #B - accuracy
    rest_accuracy_B = line_plot(rest, "B", "accuracy")
    no_threshold_accuracy_B = line_plot(new_df, "B", "accuracy")
    B_accuracy = alt.layer(rest_accuracy_B, no_threshold_accuracy_B)

    #C - call rate
    rest_call_rate_C = line_plot(rest, "C", "call_rate")
    no_threshold_call_rate_C = line_plot(new_df, "C", "call_rate")
    C_call_rate = alt.layer(rest_call_rate_C, no_threshold_call_rate_C)

    #C - accuracy
    rest_accuracy_C = line_plot(rest, "C", "accuracy")
    no_threshold_accuracy_C = line_plot(new_df, "C", "accuracy")
    C_accuracy = alt.layer(rest_accuracy_C, no_threshold_accuracy_C)

    #DQB1 - call rate
    rest_call_rate_DQB1 = line_plot(rest, "DQB1", "call_rate")
    no_threshold_call_rate_DQB1 = line_plot(new_df, "DQB1", "call_rate")
    DQB1_call_rate = alt.layer(rest_call_rate_DQB1, no_threshold_call_rate_DQB1)

    #C - accuracy
    rest_accuracy_DQB1 = line_plot(rest, "DQB1", "accuracy")
    no_threshold_accuracy_DQB1 = line_plot(new_df, "DQB1", "accuracy")
    DQB1_accuracy = alt.layer(rest_accuracy_DQB1, no_threshold_accuracy_DQB1)

    final_chart = (A_call_rate & A_accuracy) | (B_call_rate & B_accuracy) | (C_call_rate & C_accuracy) | (DQB1_call_rate & DQB1_accuracy)

    # final_chart.save('plot.json')
    final_chart.save(snakemake.output.threshold_results_line_plot)

    final_chart
#%%