import pandas as pd
import glob
import os
import subprocess
import gzip
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #list input fastq - part 1
    # input_fastq = os.getcwd()
    input_fastq = snakemake.input.fastq_dir
    print(input_fastq)
    input_fastq = os.path.join(input_fastq, '*.fastq.gz')
    print(input_fastq)
    input_fastq_list = glob.glob(input_fastq)
    print(input_fastq_list)

    #list input fastq - part 2 (additional) - files downloaded later
    input_fastq_2 = snakemake.input.fastq_dir_2

    input_fastq_2 = os.path.join(input_fastq_2, '*.fastq.gz')
    input_fastq_list_2 = glob.glob(input_fastq_2)
    print(input_fastq_list_2)

    #read sample sheet
    # sample_sheet = pd.read_csv("merge_sample_sheet.csv")
    sample_sheet = pd.read_csv(snakemake.input.sample_sheet)
    print(sample_sheet)

    #initialize columns for read length and count 
    sample_sheet["Read Length 1"] = ""
    sample_sheet["Read Length 2"] = ""
    sample_sheet["Read Count"] = ""

    def insert_read_info(sample_sheet, fq):       
        #find read length and read count of each fastq file
        zcat = subprocess.Popen(('zcat', fq), stdout=subprocess.PIPE)
        output = subprocess.check_output(('awk', '-b', 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'), stdin=zcat.stdout)
        zcat.wait()
        # read_length_dict = {}
        # fastq_file = gzip.open(fq, 'rb')
        # seq = ''
        # count=1
        # for line in fastq_file:
        #     line=line.decode("utf-8")
        #     line = line.rstrip('\n')
        #     if not line.startswith('+'): #avoid counting the separator (+) sign
        #         if line.startswith('@'):
        #             if seq:
        #                 read_length_dict[str(len(seq))]=count
        #                 count+=1
        #                 seq = ''
        #         else:
        #             seq = line 
        # fastq_file.close()
        #
        # for key,value in read_length_dict.items():
        #     read_length = key
        #     read_count = value
        #     print("read_length: "+read_length)
        #     print("read_count: "+str(read_count))

        read_length=output.decode('utf-8').split()[0]
        read_count=output.decode('utf-8').split()[1]

        #find sample name and make sure there is only one read length
        splitted = os.path.basename(fq).split("_")
        sample_name = splitted[0]
        read_pair = splitted[1].split(".")[0]
        print(sample_name)
        print(read_pair)
        read_length_pair_name = "Read Length " + read_pair
        print(read_length_pair_name)
        # print(len(read_length_dict))
        # len(read_length_dict) == 1 # there should only be one read length
        # print(read_length_dict)

        #append the read length and read count to the sample sheet
        sample_sheet.loc[sample_sheet["Run_Accession"] == sample_name, read_length_pair_name] = read_length
        sample_sheet.loc[sample_sheet["Run_Accession"] == sample_name, "Read Count"] = read_count

        return sample_sheet

    for fq in input_fastq_list:
        sample_sheet = insert_read_info(sample_sheet, fq)

    for fq in input_fastq_list_2:
        sample_sheet = insert_read_info(sample_sheet, fq)

    # sample_sheet.to_csv("sample_sheet_w_read_info.csv", index=False)
    sample_sheet.to_csv(snakemake.output.sample_sheet, index=False)
