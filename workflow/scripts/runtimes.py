#some considerations to be done manually for files:
#remove the same that was removed rom sample list: SRR071132
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
    #the following function updates tool_runtimes by inserting new rows with new records from tsv files
    #if the record is already in the tool runtimes, it adds up the max_rss and s and updates the tool_runtimes in place
    def collect_measurements_preprocessing(tool_name,tsv_files, tool_runtimes, separator):
        for f in tsv_files: 
            # read the csv file 
            df = pd.read_csv(f, sep = "\t")
            file_name = os.path.basename(f)
            # print(file_name)
            splitted = file_name.split(separator) #there is no "_" in file name in case of arcasHLA and orthanq, but then whole name is returned here.
            sample_name = splitted[0]
            sample_name = sample_name.split("_")[0] #split by "_" again because of file names in razers3_bam_to_fastq, e.g. D1_S1_L001.1.tsv, we want the D1 which is the expected name for this sample, it will have no effect on the other samples and tools
            sample_name = sample_name.rsplit('.', 1)[0] #only in case of D1_S1_L001 this is not applied
            # update the value of max_rss and s in the final table, for rows having the same tool and sample name
            filter1 = tool_runtimes[tool_runtimes["tool"] == tool_name]
            filter1_sample_list = filter1["sample"].tolist()
            if (tool_name in tool_runtimes['tool']) or (sample_name in filter1_sample_list):
                idx = tool_runtimes.index[(tool_runtimes['tool'] == tool_name) & (tool_runtimes['sample'] == sample_name)]
                max_rs = tool_runtimes.iloc[idx]['max_rss']
                new_max_rs = max_rs

                new_max_rs += (df.iloc[0]['max_rss'])/1024 # in gb

                seconds = tool_runtimes.iloc[idx]['s'] 
                new_seconds = seconds

                new_seconds += (df.iloc[0]['s'])/60 # in minutes

                tool_runtimes.loc[idx, ['max_rss']] = new_max_rs
                tool_runtimes.loc[idx, ['s']] = new_seconds
            else:
                new_record = pd.DataFrame([{'sample': sample_name, 'tool': tool_name, 'max_rss': (df.iloc[0]["max_rss"])/1024 , 's': (df.iloc[0]["s"])/60}])
                tool_runtimes = pd.concat([tool_runtimes, new_record], ignore_index=True)
        return tool_runtimes

    def collect_measurements_genotyping(tool_name,tsv_files, tool_runtimes):
        # loop over the list of csv files 
        # check initially to update the value of max_rss and s in the final table for arcas and orthanq as they output by perlocus 
        if (tool_name == "orthanq") or (tool_name == "arcasHLA"):
            for f in tsv_files: 
                # read the csv file 
                df = pd.read_csv(f, sep = "\t")
                file_name = os.path.basename(f)
                splitted = file_name.split("_")
                sample_name = splitted[0]
                locus_name = splitted[-1]
                locus_name = locus_name.rsplit('.', 1)[0]
                #check only loci that we evaluate in the paper
                if locus_name == 'A' or locus_name == 'B' or locus_name == 'C' or locus_name == 'DQB1':
                    filter1 = tool_runtimes[tool_runtimes["tool"] == tool_name]
                    filter1_sample_list = filter1["sample"].tolist()
                    # update the value of max_rss and s in the final table, for rows having the same tool and sample name
                    if (tool_name in tool_runtimes['tool']) or (sample_name in filter1_sample_list):
                        idx = tool_runtimes.index[(tool_runtimes['tool'] == tool_name) & (tool_runtimes['sample'] == sample_name)]
                        max_rs = tool_runtimes.iloc[idx]['max_rss']
                        new_max_rs = max_rs

                        new_max_rs += (df.iloc[0]['max_rss'])/1024 # in gb

                        seconds = tool_runtimes.iloc[idx]['s'] 
                        new_seconds = seconds

                        new_seconds += (df.iloc[0]['s'])/60 # in minutes

                        tool_runtimes.loc[idx, ['max_rss']] = new_max_rs
                        tool_runtimes.loc[idx, ['s']] = new_seconds
                    else:
                        new_record = pd.DataFrame([{'sample': sample_name, 'tool': tool_name, 'max_rss': (df.iloc[0]["max_rss"])/1024 , 's': (df.iloc[0]["s"])/60}])
                        tool_runtimes = pd.concat([tool_runtimes, new_record], ignore_index=True)
        else:
            for f in tsv_files: 
                df = pd.read_csv(f, sep = "\t")
                file_name = os.path.basename(f)
                splitted = file_name.split("_")
                sample_name = splitted[0]
                sample_name = sample_name.rsplit('.', 1)[0]
                new_record = pd.DataFrame([{'sample': sample_name, 'tool': tool_name, 'max_rss': (df.iloc[0]["max_rss"])/1024 , 's': (df.iloc[0]["s"])/60}])
                tool_runtimes = pd.concat([tool_runtimes, new_record], ignore_index=True)
        return tool_runtimes
    #create a plot for final table
    def create_plot(tool_runtimes):
        seconds = alt.Chart(tool_runtimes).mark_line(point=True).encode(
            x=alt.X('sample:N', axis=alt.Axis(labels=True, title=None)).sort(field="order", order= "descending"),
            y=alt.Y('s:Q', title="runtime [m]").scale(type="log"),
            color='tool:N',
        ).transform_joinaggregate(
            order='min(s)',
            groupby=["sample"]
        )
        max_rss = alt.Chart(tool_runtimes).mark_line(point=True).encode(
            x=alt.X('sample:N', axis=alt.Axis(labels=True)).sort(field="order", order= "descending"),
            y=alt.Y('max_rss:Q', title="RSS memory consumption [GB]").scale(type="log"),
            color='tool:N',
        ).transform_joinaggregate(
            order='min(s)',
            groupby=["sample"]
        )
        final = seconds & max_rss
        return final

    ## only preprocessing ##
    for subdir, dirs, files in os.walk(path):
            # tool_name = os.path.basename(subdir) #subdir is /home/uzuner/Documents/orthanq_paper/benchmarks/arcasHLA
            for dir in dirs:
                tsv_files = glob.glob(os.path.join(subdir + "/" + dir, "*.tsv"))  
                if dir == "extract": #arcasHLA
                    tool_runtimes = collect_measurements_preprocessing("arcasHLA",tsv_files,tool_runtimes, "_")
                if dir == "samtools_index_after_reheader" or dir == "samtools_reheader" or dir == "samtools_sort" or dir == "samtools_view_primary_chr": #orthanq
                    tool_runtimes = collect_measurements_preprocessing("orthanq",tsv_files,tool_runtimes, "_")
                if dir == "vg_giraffe" or dir == "samtools_view_extract_HLA" or dir ==  "bam_to_fastq": #orthanq + vg
                    tool_runtimes = collect_measurements_preprocessing("orthanq+vg",tsv_files,tool_runtimes, "_")
                if dir == "fastq_split" or dir == "razers3" or dir == "samtools_merge" or dir == "samtools_sort_razers3": #optitype
                    tool_runtimes = collect_measurements_preprocessing("optitype",tsv_files,tool_runtimes, "_")
                if dir == "razers3_bam_to_fastq":
                    tool_runtimes = collect_measurements_preprocessing("optitype",tsv_files,tool_runtimes, ".")

    #write the final table to csv
    tool_runtimes.to_csv(snakemake.output.preprocessing_table, index=False)

    # #sort by max_rs
    # tool_runtimes = tool_runtimes.sort_values('s', ascending=False) 

    preprocessing = create_plot(tool_runtimes)
    preprocessing.save(snakemake.output.preprocessing_plot)

    ## only calling or genotyping ##

    d = {'sample': [], 'tool': [], "max_rss": [], "s": []}
    tool_runtimes = pd.DataFrame(data=d)

    for subdir, dirs, files in os.walk(path):
            # tool_name = os.path.basename(subdir) #subdir is /home/uzuner/Documents/orthanq_paper/benchmarks/arcasHLA
            for dir in dirs:
                tsv_files = glob.glob(os.path.join(subdir + "/" + dir, "*.tsv"))  
                if dir == "hla-la" or dir == "optitype":
                    tool_runtimes = collect_measurements_genotyping(dir,tsv_files,tool_runtimes)
                if dir == "orthanq_call":
                    tool_runtimes = collect_measurements_genotyping("orthanq",tsv_files,tool_runtimes)
                if dir == "genotype":
                    tool_runtimes = collect_measurements_genotyping("arcasHLA",tsv_files,tool_runtimes)
    print(tool_runtimes)
    #write the final table to csv
    tool_runtimes.to_csv(snakemake.output.genotyping_table, index=False)

    # #sort by max_rs
    # tool_runtimes = tool_runtimes.sort_values('s', ascending=False) 

    genotyping = create_plot(tool_runtimes)
    genotyping.save(snakemake.output.genotyping_plot)

