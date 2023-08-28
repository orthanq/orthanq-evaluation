import json
import pandas as pd

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    # Opening JSON file
    f = open(snakemake.input.template)
    # f = open("resources/templates/evaluation_plot.json")

    # open the validation table
    # validation = pd.read_table()
    validation = pd.read_table(snakemake.input.validation, sep="\t", keep_default_na=False)
    # validation = pd.read_table("results/comparison/validation.tsv", sep="\t")
    
    # returns JSON object as a dictionary
    template = json.load(f)

    #initialize a list to collect tool_accuracies
    tool_accuracies = []

    #loop through the validation table and fill in the json template
    for row in validation.itertuples():
        arcashla = {
            "category": row.Locus,
            "group": "arcasHLA",
            "value": row.arcasHLA_Accuracy
        }
        hla_la = {
            "category": row.Locus,
            "group": "HLA-LA",
            "value": row.HLA_LA_Accuracy
        }
        optitype = {
            "category": row.Locus,
            "group": "Optitype",
            "value": row.Optitype_Accuracy
        }
        orthanq = {
            "category": row.Locus,
            "group": "Orthanq",
            "value": row.Orthanq_Accuracy
        }
        tool_accuracies.append(arcashla)
        tool_accuracies.append(hla_la)
        tool_accuracies.append(optitype)
        tool_accuracies.append(orthanq)

    template['datasets']['tool_accuracies'] = tool_accuracies

    jsonData = json.dumps(template)
    print(jsonData)

    #write json to file
    # output_path = "results/comparison/evaluation_plot.json"
    output_path = snakemake.output.plot

    # with open(output_path, 'w', encoding='utf-8') as f:
    #     json.dump(jsonData, f, ensure_ascii=False, indent=4)

    with open(output_path, "w") as text_file:
        text_file.write(jsonData)
