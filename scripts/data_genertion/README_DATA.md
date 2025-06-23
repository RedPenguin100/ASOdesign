## Data ASOptimizer
This data was retrieved from the ASOptimizer repository on GitHub - https://github.com/Spidercores/ASOptimizer

Hwang G, Kwon M, Seo D, Kim DH, Lee D, Lee K, Kim E, Kang M, Ryu JH. ASOptimizer: Optimizing antisense 
oligonucleotides through deep learning for IDO1 gene regulation. Mol Ther Nucleic Acids. 
2024 Apr 6;35(2):102186. doi: 10.1016/j.omtn.2024.102186. PMID: 38706632; PMCID: PMC11066473.

## Latest version of data
data_asoptimizer_updated.csv

## Added this version
1. Column 'mod_scan' - binary column states whether the experiment is part of modification scan (1) 
or not (0). An experiment is considered as mod_scan if there is a subset of rows (more than 1) 
that has identical values in the following columns - ['Sequence', 'Density(cells/well)', 
'Transfection', 'Treatment_Period(hours)', 'ASO_volume(nM)']. If there is more than one unique
value in the following columns, the whole subset is considered a mod_scan - ['Modification', 
'Location', 'Chemical_Pattern', 'Linkage', 'Linkage_Location', 'Smiles'].
This is used to filtered out rows that has identical experiment features and sequence, but different inhibition
due to changes in the modifications type, location or chemical pattern.

2. Column 'true_length_of_seq' - several sequences have errors in the 'seq_length' column,
so the true length of the sequence itself was added.

3. Rows with no inhibitions where removed.


## Previous versions (latest to oldest)
1. data_updated_18.5.csv (columns about the location of the sequence in the gene mRNA or transcript were added)

2. data_from_article_fixed.csv (first version after analyzing)

