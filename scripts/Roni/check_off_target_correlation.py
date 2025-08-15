import pandas as pd
import matplotlib.pyplot as plt

may_df = pd.read_csv("/home/oni/ASOdesign/scripts/data_genertion/data_updated_inhibition.csv")
off_target_df = pd.read_csv("/home/oni/ASOdesign/scripts/data_genertion/off_target.top500.cutoff1400.premRNA.csv")

merged = pd.merge(may_df[['index', 'Inhibition(%)']], off_target_df[['index', 'off_target_score_log']], on='index', how='inner')

correlation = merged['Inhibition(%)'].corr(merged['off_target_score_log'])
print(correlation)

plt.scatter(merged['Inhibition(%)'], merged['off_target_score_log'])
plt.xlabel('Inhibition(%)')
plt.ylabel('off_target_score_log')
plt.title('Inhibition(%) vs off_target_score_log')
plt.grid(True)
plt.show()