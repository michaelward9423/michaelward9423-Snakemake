# Using Snakemake with Kallisto

The snakefile included here has 4 rules. One counts the total number of reads. The trim rule removes the low quality reads. Then a quant rule quantifies the abundances of transcripts. Finally, an index rule generates an index in a FASTA file or files.
