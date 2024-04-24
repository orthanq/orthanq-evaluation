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

    #remove 1.0 rows
    table.drop(table[table.threshold == 1.0].index, inplace=True)

    #take only the rate
    table["call_rate"] = table["call_rate"].apply(take_first).astype(float)
    table["accuracy"] = table["accuracy"].apply(take_first).astype(float)

    print(table)

    #convert percentages to decimals
    table["accuracy"] = table["accuracy"].apply(convert_decimal)

    print(table)

    def plot_line_plot(table, locus, column_name, color):
        chart = alt.Chart(table).mark_line(point=True, tooltip=True, color=color).encode(
                                        x=alt.X('threshold',
                                                axis=alt.Axis(
                                                            title="threshold", labels=True),
                                                scale=alt.Scale(domain=[0.0,0.9])),
                                        y= alt.Y(column_name, scale=alt.Scale(domain=[0.0,1.0]))).transform_filter(
                                            alt.datum.locus == locus).properties(title=locus)
        return chart

    chart_call_rate_A = plot_line_plot(table, "A", "call_rate", "orange")
    chart_accuracy_A = plot_line_plot(table, "A", "accuracy", "blue")

    chart_call_rate_B = plot_line_plot(table, "B", "call_rate", "orange")
    chart_accuracy_B = plot_line_plot(table, "B", "accuracy", "blue")

    chart_call_rate_C = plot_line_plot(table, "C", "call_rate", "orange")
    chart_accuracy_C = plot_line_plot(table, "C", "accuracy", "blue")

    chart_call_rate_DQB1 = plot_line_plot(table, "DQB1", "call_rate", "orange")
    chart_accuracy_DQB1 = plot_line_plot(table, "DQB1", "accuracy", "blue")

    final_chart = alt.hconcat((chart_call_rate_A + chart_accuracy_A), (chart_call_rate_B + chart_accuracy_B), (chart_call_rate_C + chart_accuracy_C), (chart_call_rate_DQB1 + chart_accuracy_DQB1)).configure_axis(grid=False).resolve_scale(y='shared')

    # final_chart.save('plot.json')
    final_chart.save(snakemake.output.threshold_results_line_plot)
#%%