## Data ASOptimizer
This data was retrieved from the ASOptimizer repository on GitHub - https://github.com/Spidercores/ASOptimizer

Hwang G, Kwon M, Seo D, Kim DH, Lee D, Lee K, Kim E, Kang M, Ryu JH. ASOptimizer: Optimizing antisense 
oligonucleotides through deep learning for IDO1 gene regulation. Mol Ther Nucleic Acids. 
2024 Apr 6;35(2):102186. doi: 10.1016/j.omtn.2024.102186. PMID: 38706632; PMCID: PMC11066473.

## Latest version of data
data_asoptimizer_05_08.csv

## NEW UPDATE FOR ALL CSV's
A new column was added to each version of the data - 'index'. Each row has a unique index to keep order
of data and further features. Please make sure to run your scripts again with the indices to ensure
all data is appropriately saved and organized.

## Added this version
1. Columns for flags for future invalidation of rows (not relevant for features and other information, please ignore 
columns if working on those subjects). 
- aso_volume_high_check: all rows with ASO volume higher than 20,000 nM
- cell_line_check: all rows with unknown cell line
- density_zero_check: all rows with density zero
- invalidating: rows suspected in invalidating
2. New column - cell_line_after_check, added cell line names for unknown cell lines (as found in patents)
3. Updated cell line organism for unknown cell lines
4. Correction for Transcription, Location_in_sequence and Location_div_by_length - if not found, Transcription is Nan,
Location and div by length are -1, in order to differentiate between sequences not found and sequences appear in the
first nucleotide of the gene.
5. SNP gene was found for several sequences and location was updated accordingly.



## Previous versions (latest to oldest)
1. data_updated_inhibition.csv

2. data_asoptimizer_updated.csv (mos_scan, true_length_of_seq and removing no inhibition)

3. data_updated_18.5.csv (columns about the location of the sequence in the gene mRNA or transcript were added)

4. data_from_article_fixed.csv (first version after analyzing)


