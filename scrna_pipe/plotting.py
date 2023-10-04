import requests
import numpy as np
import pandas as pd
import scanpy as sc


def plot_string_network(gene_list):
    """
    """
    string_api_url = "https://version-11-5.string-db.org/api"

    output_format = "highres_image"
    method = "network"

    params = {

        "identifiers" : "\r".join(gene_list),
        "species" : 9606, 
        "limit" : 1, 
        "echo_query" : 1, 
        "caller_identity" : "www.awesome_app.org" 

    }

    request_url = "/".join([string_api_url, output_format, method])

    response = requests.post(request_url, data=params)

    return response.content


def extract_cluster_map(de_df, adata, cell_type, contrast):
    """
    """
    de_df = de_df[abs(de_df['logFC']) > 1].copy()
    de_df.sort_values(by=['FDR', 'logFC'], inplace=True)
    
    try:
        de_genes = list(de_df[0:50]['gene'])
    except:
        de_genes = list(de_df['gene'])

    adata.obs['sample'] = adata.obs.apply(lambda x: x['sample'].replace('-','_'), axis=1)

    adata.obs["cell_type"] = [ct.replace(" ", "_") for ct in adata.obs["cell_type"]]
    adata.obs["cell_type"] = [ct.replace("+", "") for ct in adata.obs["cell_type"]]
    adata.obs["replicate"] = adata.obs.apply(lambda x: x['sample'].split('_')[-1], axis=1)
    adata.obs["label"] = adata.obs.apply(lambda x: '_'.join(x['sample'].split('_')[:-1]), axis=1)

    adata.obs["replicate"] = adata.obs["replicate"].astype("category")
    adata.obs["label"] = adata.obs["label"].astype("category")
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print(contrast, cell_type)
    print(set(adata.obs['label']))
    print(set(adata.obs['cell_type']))

    plot_adata = adata[(adata.obs['cell_type'] == cell_type) &\
                       (adata.obs['label'].isin(list(contrast)))]
    
    plot_adata = plot_adata[:, plot_adata.var.index.isin(de_genes)]
    
    cg = sc.pl.clustermap(plot_adata, use_raw=False,
                      obs_keys='label', cmap='mako', show=False)
    
    row_order = cg.dendrogram_row.reordered_ind
    col_order = cg.dendrogram_col.reordered_ind

    plot_adata.obs['cell_id'] = plot_adata.obs.apply(lambda x: x.label + x.name, axis=1)
    df = pd.DataFrame(plot_adata.X.toarray(), columns=plot_adata.var.index, index=plot_adata.obs['cell_id'])
    plot_df = df.iloc[row_order, col_order]
    
    plot_df['label'] = plot_df.apply(lambda x: np.max(plot_df) if contrast[0] in x.name else 0, axis=1)
    
    return plot_df