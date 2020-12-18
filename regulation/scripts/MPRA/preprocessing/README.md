
This README contains information on how to pre-process and analyse fq.gz files from Illumina containing MPRA barcode reads from an assay involving stable viral transduction of the plasmid library (reporter: GFP, but could be any other).

One needs to:

1) download the MPRA libraries with the accession number E-MTAB-9952. Within all the following scripts, adjust the working directory.
2) `I_deML.sh` - demultiplex libraries
3) `II_filtering.sh` - this script needs two helper files, also present here - III_quality.filtering_counting_pipeline_GFP.sh and IV_reads.through.filtering_all.samples.sh. The first one filters the demultiplexed read files based on base quality (Q>10 for all 10 barcode bases --> means 90% probability of calling the right base). In addition, it checks whether the constant region (here, the beginning of the GFP gene) is present and extracts the 10 bases before, thereby making sure to avoid wrong barcodes due to reads that did not have the correct starting point during amplification. The second script summarises filtering statistics into one, R-friendly file.
4) `V_Count_table_matching_metadata_MPRA.R` - generation of count tables.

The output of these scripts can be found under 


