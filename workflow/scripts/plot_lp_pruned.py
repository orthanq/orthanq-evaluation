import altair as alt
import pandas as pd
import sys

# %
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    # orthanq_table = pd.read_csv("pruned.tsv", sep = "\t")
    orthanq_table = pd.read_csv(snakemake.input.pruned_table, sep = "\t")

    chart = alt.Chart(orthanq_table).mark_bar(tooltip=True).encode(
                                        x=alt.X('percentage:Q', title="LP reduction [%]"),
                                        y=alt.Y('count()', title="count"))
    chart
    # chart.save('pruned_haplotypes.json')
    chart.save(snakemake.output.pruned_json)
#%%