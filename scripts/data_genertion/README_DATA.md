## Data ASOptimizer
This data was retrieved from the ASOptimizer repository on GitHub - https://github.com/Spidercores/ASOptimizer

Hwang G, Kwon M, Seo D, Kim DH, Lee D, Lee K, Kim E, Kang M, Ryu JH. ASOptimizer: Optimizing antisense 
oligonucleotides through deep learning for IDO1 gene regulation. Mol Ther Nucleic Acids. 
2024 Apr 6;35(2):102186. doi: 10.1016/j.omtn.2024.102186. PMID: 38706632; PMCID: PMC11066473.

## Latest version of data
data_updated_inhibition.csv

## NEW UPDATE FOR ALL CSV's
1. A new column was added to each version of the data - 'index'. Each row has a unique index to keep order
of data and further features. Please make sure to run your scripts again with the indices to ensure
all data is appropriately saved and organized.
2. A new column was added to combine cell line A431 and A-431 - 'cell_line_uniform'. For the rest of the
cell lines, it is the same in the new column as in 'Cell_Line'.

## Added this version
1. Column 'corrected_inhibition' - average of inhibition for duplicated rows that has the same values 
for all column except inhibition and index.

2. Column 'mutual_index' - for duplicated rows, a mutual index was made, that is the index of the
first row that appear with this data. 



## Previous versions (latest to oldest)
1. data_asoptimizer_updated.csv (mos_scan, true_length_of_seq and removing no inhibition)

2. data_updated_18.5.csv (columns about the location of the sequence in the gene mRNA or transcript were added)

3. data_from_article_fixed.csv (first version after analyzing)

