import json
import pandas as pd

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    # Opening JSON file
    f = open(snakemake.input.template)

    # open the validation table
    validation_low = pd.read_table(snakemake.input.validation_low, sep="\t", keep_default_na=False)
    validation_high = pd.read_table(snakemake.input.validation_high, sep="\t", keep_default_na=False)

    # validation = pd.read_table("results/comparison/validation.tsv", sep="\t")
    
    def create_evaluation_plot(validation_table, template_file):
        # returns JSON object as a dictionary
        template = json.load(f)

        #initialize a list to collect tool_accuracies
        tool_accuracies = []

        #loop through the validation table and fill in the json template
        for row in validation_table.itertuples():
            #arcasHLA
            arcasHLA_accuracy_with_fraction = row.arcasHLA_Accuracy.split("(")[0] #parse the accuracy value
            arcasHLA_accuracy = arcasHLA_accuracy_with_fraction.rstrip() #remove the trailing whitespace

            #HLA-LA
            hla_la_accuracy_with_fraction = row.HLA_LA_Accuracy.split("(")[0] #parse the accuracy value
            hla_la_accuracy = hla_la_accuracy_with_fraction.rstrip() #remove the trailing whitespace

            #optitype
            optitype_accuracy_with_fraction = row.Optitype_Accuracy.split("(")[0] #parse the accuracy value
            optitype_accuracy = optitype_accuracy_with_fraction.rstrip() #remove the trailing whitespace

            #orthanq
            orthanq_accuracy_with_fraction = row.Orthanq_Accuracy.split("(")[0] #parse the accuracy value
            orthanq_accuracy = orthanq_accuracy_with_fraction.rstrip() #remove the trailing whitespace

            arcashla = {
                "category": row.Locus,
                "group": "arcasHLA",
                "value": arcasHLA_accuracy
            }
            hla_la = {
                "category": row.Locus,
                "group": "HLA-LA",
                "value": hla_la_accuracy
            }
            optitype = {
                "category": row.Locus,
                "group": "Optitype",
                "value": optitype_accuracy
            }
            orthanq = {
                "category": row.Locus,
                "group": "Orthanq",
                "value": orthanq_accuracy
            }
            tool_accuracies.append(arcashla)
            tool_accuracies.append(hla_la)
            tool_accuracies.append(optitype)
            tool_accuracies.append(orthanq)

        template['datasets']['tool_accuracies'] = tool_accuracies

        jsonData = json.dumps(template)
        print(jsonData)
        return jsonData

    #evaluation plot for low coverage
    json_data_low = create_evaluation_plot(validation_low, f)

    #write json to file
    output_path_low = snakemake.output.plot_low

    with open(output_path_low, "w") as text_file:
        text_file.write(json_data_low)

    #evaluation plot for high coverage
    f = open(snakemake.input.template)
    json_data_high = create_evaluation_plot(validation_high, f)

    #write json to file
    output_path_high = snakemake.output.plot_high

    with open(output_path_high, "w") as text_file:
        text_file.write(json_data_high)
