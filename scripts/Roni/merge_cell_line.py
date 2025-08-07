import pandas as pd
import os
from glob import glob
from get_premRNA_sequences_II import cell_line2dataII



path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/'
cell_line_list = ['ACH-000232', 'ACH-000463', 'ACH-000739', 'ACH-001086', 'ACH-001188', 'ACH-001328']
addition = '_transcriptome_premRNA'

n=5
print(n)
df1 = pd.read_csv(path + cell_line_list[n] + addition + '.csv')
print(df1.shape)
#df2 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_5_8.premRNA.completed.csv')
#df3 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_9_12.premRNA.completed.csv')
df4 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_13_16.premRNA.completed.csv')
df5 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_17_20.premRNA.completed.csv')
df6 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_21_22_XY.premRNA.completed.csv')
df7 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_M.premRNA.completed.csv')
df8 = pd.read_csv(path + cell_line_list[n] + addition + '.chr_rest.premRNA.completed.csv')

df_list = [df4, df5, df6, df7, df8]  # exclude df1 itself

df1.set_index("Transcript_ID", inplace=True)
for df in df_list:
    print(df.shape)
    df.set_index("Transcript_ID", inplace=True)
    for col in df.columns:
        if col in df1.columns:
            df1[col] = df1[col].combine_first(df[col])
        else:
            df1[col] = df[col]
    print('merged')

df1.reset_index(inplace=True)
df1_sorted = df1.sort_values(by=['expression_TPM'], ascending=False)
print(df1_sorted.shape)
df1_sorted.to_csv(path + cell_line_list[n] + addition + '.merged.csv', index=False)