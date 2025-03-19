This folder contains subsampled and edited plasmid QC data from the 300125_novaseq_qc_50048 run:

* sample_sheet/tsv: saample sheet that should be used when testing whether > 40%, 30% - 40% and < 40% library reads lead to PASS or FAIL in the html report. This is a downsampled and edited
version of the original sample sheet.

* Unedited chromosome files.

* Unedited query_counts and lib_counts files.

* 7135STDY14473712.lib_counts_ds.tsv.gz: downsampled file in order to get < 30% library reads. It was generated the following way:
zcat data/plasmid_300125_novaseq_qc_50048_TD778/SGE_QC15395666.lib_counts.tsv.gz | head -n 1200 | gzip > data/plasmid_300125_novaseq_qc_50048_TD778/SGE_QC15395666.lib_counts_ds.tsv.gz

NB, no samples contain exact 40% library reads in this test dataset.