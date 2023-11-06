import os
import glob 
import pandas as pd
import altair as alt


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #create the df
    d = {'sample': [], 'tool': [], "max_rss": [], "s": []}
    tool_runtimes = pd.DataFrame(data=d)
    # path = os.getcwd()
    path = snakemake.input.benchmarks
    print("path is: ")
    print(path)
    #the following function inserts values given in df to tool_runtimes.
    #It does this by firstly checking the sample and tool names and if they exist, the values are summed up and 
    #max_rss and s are updated in GB and minutes, if not new entries are added.
    def insert_values(tool_runtimes, df, sample_name, tool_name, filter1_sample_list):
        if (tool_name in tool_runtimes['tool']) or (sample_name in filter1_sample_list):
            idx = tool_runtimes.index[(tool_runtimes['tool'] == tool_name) & (tool_runtimes['sample'] == sample_name)]
            max_rss = tool_runtimes.iloc[idx[0]]['max_rss'] #we know that there will only be one index

            # new_max_rs = max_rss
            # new_max_rs += (df.iloc[0]['max_rss'])/1024
            
            max_max_rss = max(max_rss, (df.iloc[0]['max_rss'])/1024) # get the max of each processes
            
            minutes = tool_runtimes.iloc[idx]['s'] 

            new_minutes = minutes
            new_minutes += (df.iloc[0]['s'])/60
            
            tool_runtimes.loc[idx, ['max_rss']] = max_max_rss
            tool_runtimes.loc[idx, ['s']] = new_minutes
        else:
            new_record = pd.DataFrame([{'sample': sample_name, 'tool': tool_name, 'max_rss': (df.iloc[0]["max_rss"])/1024 , 's': (df.iloc[0]["s"])/60}])
            tool_runtimes = pd.concat([tool_runtimes, new_record], ignore_index=True)
        return tool_runtimes

    #collect_measurements loops in given list of tsv_files, collects sample and locus names (if present),
    #checks for multiple entries for each sample (only for orthanq (both preprocessing and genotyping) and arcasHLA)
    #checks for only the relevant loci (A, B, C, DQB1)
    #then updates or adds new values to too_runtimes
    def collect_measurements(tool_name,tsv_files, tool_runtimes, separator):
        # loop over the list of csv files 
        for f in tsv_files: 
            # read the csv file 
            df = pd.read_csv(f, sep = "\t")
            file_name = os.path.basename(f)
            splitted = file_name.split(separator) #there is no "_" in file name in case of arcasHLA and orthanq, but then whole name is returned here.
            sample_name = splitted[0]
            sample_name = sample_name.split("_")[0] #split by "_" again because of file names in razers3_bam_to_fastq, e.g. D1_S1_L001.1.tsv, we want the D1 which is the name comes from other tools. This op. will have no effect on the other sample names
            sample_name = sample_name.rsplit('.', 1)[0] #in case of file names with no locus. e.g. SRR702070.tsv and the "separator" in input function is "_".

            locus_name = splitted[-1]
            locus_name = locus_name.rsplit('.', 1)[0]
            filter1 = tool_runtimes[tool_runtimes["tool"] == tool_name]
            filter1_sample_list = filter1["sample"].tolist()
            
            #check only loci that we evaluate in the paper
            if (tool_name == "orthanq (excl. vg)") or (tool_name == "arcasHLA"):
                if locus_name == 'A' or locus_name == 'B' or locus_name == 'C' or locus_name == 'DQB1': #varlociraptor from preprocessing, orthanq_call and arcasHLA
                    # update the value of max_rss and s in the final table, for rows having the same tool and sample name
                    tool_runtimes = insert_values(tool_runtimes,df, sample_name, tool_name, filter1_sample_list)
                elif locus_name == "DQA1" or locus_name == "DRB1":
                    continue
                else:
                    tool_runtimes = insert_values(tool_runtimes, df, sample_name, tool_name, filter1_sample_list)
            else:
                tool_runtimes = insert_values(tool_runtimes, df, sample_name, tool_name, filter1_sample_list)
        return tool_runtimes

    #create a plot for final table
    def create_plot(tool_runtimes):
        seconds = alt.Chart(tool_runtimes).mark_line(point=True, tooltip=True).encode(
            x=alt.X('sample:N', axis=alt.Axis(labels=False, title=None)).sort(field="order", order= "descending"),
            y=alt.Y('s:Q', title="runtime [m]").scale(type="log"),
            color='tool:N',
        ).transform_joinaggregate(
            order='min(s)',
            groupby=["sample"]
        )
        max_rss = alt.Chart(tool_runtimes).mark_line(point=True, tooltip=True).encode(
            x=alt.X('sample:N', axis=alt.Axis(labels=False)).sort(field="order", order= "descending"),
            y=alt.Y('max_rss:Q', title="RSS memory consumption [GB]").scale(type="log"),
            color='tool:N',
        ).transform_joinaggregate(
            order='min(s)',
            groupby=["sample"]
        )
        final = seconds & max_rss
        return final

    ## genotyping/calling with preprocessing combined ##
    for subdir, dirs, files in os.walk(path):
            for dir in dirs:
                tsv_files = glob.glob(os.path.join(subdir + "/" + dir, "*.tsv"))  
                if dir == "hla-la":
                    tool_runtimes = collect_measurements(dir,tsv_files,tool_runtimes, "_")
                if dir == "optitype" or dir == "fastq_split" or dir == "razers3" or dir == "samtools_merge" or dir == "samtools_sort_razers3": #optitype
                    tool_runtimes = collect_measurements("optitype",tsv_files,tool_runtimes, "_")
                if dir == "razers3_bam_to_fastq": #also optitype
                    tool_runtimes = collect_measurements("optitype",tsv_files,tool_runtimes, ".")
                if dir == "vg_giraffe" or dir == "samtools_view_extract_HLA" or dir ==  "bam_to_fastq": #orthanq + vg
                    tool_runtimes = collect_measurements("vg",tsv_files,tool_runtimes, "_")
                if dir == "orthanq_call" or dir == "samtools_index_after_reheader" or dir == "samtools_reheader" or dir == "samtools_sort" or dir == "samtools_view_primary_chr" or dir == "varlociraptor_preprocess" or dir == "varlociraptor_call":
                    tool_runtimes = collect_measurements("orthanq (excl. vg)",tsv_files,tool_runtimes, "_")
                if dir == "genotype" or dir == "extract":
                    print(dir)
                    tool_runtimes = collect_measurements("arcasHLA",tsv_files,tool_runtimes, "_")
    print(tool_runtimes)

    #separate vg from (orthanq excl. vg)
    vg = tool_runtimes[tool_runtimes["tool"] == "vg"]
    orthanq = tool_runtimes[tool_runtimes["tool"] == "orthanq (excl. vg)"]

    #take the sum of vg runtimes from orthanq and vg
    orthanq_vg_runtime = pd.concat([vg, orthanq]).groupby(["sample"], as_index=False).sum()
    orthanq_vg_runtime.to_csv("orthanq_plus_vg_runtime.csv", index=False)

    #take the max of orthanq and vg
    orthanq_vg_max_memory = pd.concat([vg, orthanq]).groupby(["sample"], as_index=False)['max_rss'].max()
    orthanq_vg_max_memory.to_csv("orthanq_vg_max_memory.csv", index=False)

    #create a new dataframe for vg included orthanq results with new s and max_rss
    d = {'sample': orthanq_vg_runtime["sample"], 
        'tool': "orthanq (incl. vg)", 
        "max_rss": orthanq_vg_max_memory["max_rss"], 
        "s": orthanq_vg_runtime["s"]}
    orthanq_plus_vg = pd.DataFrame(data=d)

    #append the vg included df to the tool runtimes df
    tool_runtimes = pd.concat([tool_runtimes,orthanq_plus_vg])

    #remove rows with "vg"
    tool_runtimes = tool_runtimes[tool_runtimes["tool"] != "vg"]
    # tool_runtimes.to_csv("concatenated.csv", index=False)

    #write the final table to csv
    tool_runtimes.to_csv(snakemake.output.runtimes_table, index=False)
    # tool_runtimes.to_csv("runtimes.csv", index=False)

    #create plot
    runtimes_plot = create_plot(tool_runtimes)

    #write to json
    runtimes_plot.save(snakemake.output.runtimes_plot)
    # runtimes_plot.save("runtimes.json")
