import pandas as pd
import numpy as np
from pandas.core.interchange.dataframe_protocol import DataFrame

path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/ACH-001328_transcriptome_top500.csv'
df = pd.read_csv(path)
# print(df.columns.tolist())

yes = df[df['Gene'].str.contains('LSM5', na=False)]
# print(yes['expression_TPM'])


blah = pd.DataFrame({
    'Original Transcript sequence': ['ATCGTTA', 'TTATGG', np.nan, 'GGCAT'],
    'Mutated Transcript sequence': ['ATTGCT', np.nan, 'TTTT', 'CGT']
})

cols = ['Original Transcript sequence', 'Mutated Transcript sequence']
for col in cols:
    blah[col] = blah[col].apply(lambda x: x.replace('T', 'U') if isinstance(x, str) else x)

print(blah)