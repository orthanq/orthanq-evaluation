#download and write the snapshot of allele frequencies database with a commit id of "11bde30", ?raw=true added later to the link to make the download possible
# download.file("https://github.com/Genentech/midasHLA/blob/11bde30cbbf11b34f2dea29a6284371a9c1e9440/data/allele_frequencies.rda?raw=true", snakemake@output[["rda_file"]])
# load(snakemake@output[["rda_file"]])
# write.csv(allele_frequencies,file=snakemake@output[["csv"]])

download.file("https://github.com/Genentech/midasHLA/blob/11bde30cbbf11b34f2dea29a6284371a9c1e9440/data/allele_frequencies.rda?raw=true", snakemake@output[["rda_file"]])
load(snakemake@output[["rda_file"]])
write.csv(allele_frequencies, snakemake@output[["csv"]])
