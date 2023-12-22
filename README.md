# Orthanq-evaluation [![DOI](https://zenodo.org/badge/475809624.svg)](https://zenodo.org/badge/latestdoi/475809624)

An evaluation workflow that performs haplotype quantification (i.e. more specifically, HLA typing, at the moment. Haplotype quantification is achieved for given samples or simulated samples based on the user configuration on the [config](config/config.yaml). It is carried out for 3 scenarios that are possible in the workflow:
1) simulation of samples given HLA alleles in [alleles sheet](config/alleles.tsv). 
2) generation of subclonal samples with samples given in [sample sheet](config/samples.tsv).
3) quantification based on given samples in [sample sheet](config/samples.tsv).
