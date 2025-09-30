from utils import *
from feature_creation import *

def get_n_best_res(genes_lst,n,tp=24,col_to_plot = ['sense','score']):
    res ={}
    gene_to_data =  create_gene_to_data(genes_lst)
    dfs = fill_dfs(genes_lst,gene_to_data,tp = tp)
    model = create_and_load_model()
    for gene, df in dfs.items():
        df = dfs[gene].copy()
        df["score"] = model.predict(df[selected_features].values)
        df["sense"] = df[SEQUENCE].astype(str).str.translate(tbl).str[::-1]
        top_n = df.nlargest(n, 'score')
        res[gene] = top_n[col_to_plot].copy()
    return res    
        
    