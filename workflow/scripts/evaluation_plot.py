import json
import pandas as pd

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    # Opening JSON file
    f = open(snakemake.input.template)

    # open the validation table
    # validation_low = pd.read_table(snakemake.input.validation_low, sep="\t", keep_default_na=False)
    # validation_high = pd.read_table(snakemake.input.validation_high, sep="\t", keep_default_na=False)
    validation_all = pd.read_table(snakemake.input.validation_all, sep="\t", keep_default_na=False)

    # validation = pd.read_table("results/comparison/validation.tsv", sep="\t")
    
    def create_evaluation_plot(validation_table, template_file):
        # returns JSON object as a dictionary
        template = json.load(f)

        #initialize a list to collect tool_accuracies
        tool_accuracies = []

        #loop through the validation table and fill in the json template
        for row in validation_table.itertuples():

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
        return jsonData

    # #evaluation plot for low coverage
    # json_data_low = create_evaluation_plot(validation_low, f)

    # #write json to file
    # output_path_low = snakemake.output.plot_low

    # with open(output_path_low, "w") as text_file:
    #     text_file.write(json_data_low)

    #evaluation plot
    json_data = create_evaluation_plot(validation_all, f)

    #write json to file
    output_path_all = snakemake.output.plot_all

    with open(output_path_all, "w") as text_file:
        text_file.write(json_data)
