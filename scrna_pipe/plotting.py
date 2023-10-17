import requests
import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px
import plotly.graph_objects as go

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram


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

    plot_adata = adata[(adata.obs['cell_type'] == cell_type) &\
                       (adata.obs['label'].isin(list(contrast)))]
    
    plot_adata = plot_adata[:, plot_adata.var.index.isin(de_genes)]

    row_linkage = linkage(pdist(plot_adata.X.todense()), method='average')
    col_linkage = linkage(pdist(plot_adata.X.T.todense()), method='average')

    dendro_row = dendrogram(row_linkage, no_plot=True)
    dendro_col = dendrogram(col_linkage, no_plot=True)

    plot_adata.obs['cell_id'] = plot_adata.obs.apply(lambda x: x.label + x.name, axis=1)
    df = pd.DataFrame(plot_adata.X.toarray(), columns=plot_adata.var.index, index=plot_adata.obs['cell_id'])
    plot_df = df.iloc[dendro_row['leaves'], dendro_col['leaves']]
    
    plot_df['label'] = plot_df.apply(lambda x: np.max(plot_df) if contrast[0] in x.name else 0, axis=1)
    
    return plot_df


def assign_opacity(row, treatment, control):
    """
    """
    if (row[f'Significant_{treatment}'] and not row[f'Significant_{control}']):
        return .8
    elif (row[f'Significant_{treatment}'] and row[f'Significant_{control}']):
        return .4
    elif (row[f'Significant_{control}']):
        return .2
    else:
        return .1


def assign_size(row, treatment, control):
    """
    """
    if (row[f'Significant_{treatment}'] and not row[f'Significant_{control}']):
        return 10
    elif (row[f'Significant_{treatment}'] and row[f'Significant_{control}']):
        return 5
    elif (row[f'Significant_{control}']):
        return 3
    else:
        return 1


def plot_static_umap(adata_dim, samples, gene):
    """
    """
    try:
        gene_idx = adata_dim.var.index.get_loc(gene)
        max_val = np.percentile(adata_dim.X[:, gene_idx], q=98)
    except:
        max_val = 5

    ax = sc.pl.umap(
        adata_dim[adata_dim.obs['sample'].isin(samples)],
        color=[gene],
        frameon=False,
        sort_order=False,
        wspace=1,
        return_fig=True,
        vmin=0.1,
        vmax=max_val,
        cmap='viridis'
    )

    return ax


def plot_pathway_volcano(plot_gsva_df):
    """
    """
    fig = px.scatter(plot_gsva_df, x='logFC', y='-Log10(FDR)', opacity=0.5,
                    color="Significant", color_discrete_map={True: "blue", False: "red"},
                    size='point_size', size_max=5, template='plotly_white',
                    labels={"logFC": "Log2(FoldChange)"},
                    hover_data=['Pathway'])

    fig.update_layout(legend_font=dict(size=18))

    return fig


def plot_gene_volcano(plot_gene_df):
    """
    """
    fig = px.scatter(plot_gene_df, x='logFC', y='-Log10(FDR)', opacity=0.5,
                    color="Significant", color_discrete_map={True: "blue", False: "red"},
                    size='point_size', size_max=5, template='plotly_white',
                    labels={"logFC": "Log2(FoldChange)"},
                    hover_data=['gene'])

    fig.update_layout(legend_font=dict(size=18))

    return fig


def plot_diff_scatter(top_dfs, contrast_map, cell_type, dataset_abbr, contrast1, contrast2, 
                      diff_type='gene', algorithm='edgeR'):
    """
    """
    cell_type = cell_type.replace(' ', '').replace("_", "")
    
    key1 = f'{dataset_abbr}__{algorithm}__{contrast_map[contrast1]}__{cell_type}'
    key2 = f'{dataset_abbr}__{algorithm}__{contrast_map[contrast2]}__{cell_type}'

    df1 = top_dfs[key1]
    df2 = top_dfs[key2]

    treatment = contrast1.split('-')[0].strip()
    control = contrast2.split('-')[0].strip()

    if diff_type == 'gene':
        df1[f'Significant'] = df1.apply(lambda x: x.FDR < .05 and np.abs(x.logFC) > 1, axis=1)
        df2[f'Significant'] = df2.apply(lambda x: x.FDR < .05 and np.abs(x.logFC) > 1, axis=1)
    else:
        df1[f'Significant'] = df1.apply(lambda x: x.FDR < .05, axis=1)
        df2[f'Significant'] = df2.apply(lambda x: x.FDR < .05, axis=1)
        
    scatter_df = df1.merge(df2, left_on=diff_type, right_on=diff_type, suffixes=(f'_{treatment}', f'_{control}'))

    scatter_df['size'] = scatter_df.apply(lambda x:  assign_size(x, treatment, control), axis=1)
    scatter_df['opacity'] = scatter_df.apply(lambda x:  assign_opacity(x, treatment, control), axis=1)
    
    scatter = px.scatter(scatter_df, x=f'logFC_{treatment}', y=f'logFC_{control}', color=f'Significant_{treatment}', 
                        color_discrete_map={True: "blue", False: "red"},
                        size_max=20, template='plotly_white', opacity=scatter_df['opacity'], size=scatter_df['size'],
                        hover_data=[diff_type, 'opacity', f'Significant_{treatment}', f'Significant_{control}'])

    scatter.update_layout(title=f"Contrast {diff_type} scatter")

    scatter_df = scatter_df[scatter_df[f'Significant_{treatment}'] & ~scatter_df[f'Significant_{control}']]
    
    
    scatter_cols = [diff_type, f'logFC_{treatment}', f'FDR_{treatment}', f'Significant_{treatment}',
                        f'logFC_{control}', f'FDR_{control}', f'Significant_{control}']
    
    if diff_type == 'Pathway':
        scatter_cols += [f'leading_edge_{treatment}',
                         f'DE Genes In Pathway_{treatment}',
                         f'Total Genes In Pathway_{treatment}']
                         
    print(scatter_cols, scatter_df.shape)
    return scatter, scatter_df[scatter_cols]


def plot_umap(plot_df, color_map, title):
    """
    """
    fig = px.scatter(plot_df, x=0, y=1, opacity=0.5,
                    color="cell_type", color_discrete_map=color_map,
                    labels={
                        "0": "UMAP 2",
                        "1": "UMAP 1",
                    },)

    fig.update_layout(legend_font=dict(size=24), title=title)

    return fig


def overview_count_table(contrast_tuples, gene_de_dfs, dataset_name):
    """
    """
    all_contrasts = [f'{x[0]}-{x[1]}' for x in contrast_tuples]
    gene_contrasts = [x for x in list(gene_de_dfs.keys()) if any([c in x for c in all_contrasts]) and 'subset' in x]
    cell_types = [x.split('_')[-1] for x in gene_contrasts]

    count_dict = {}

    for contrast in all_contrasts:
        if 'None' in contrast: continue
        count_dict[contrast] = {}
        for cell_type in cell_types:
            key = f'{dataset_name}__edgeR__{contrast}__{cell_type}'
            if key in gene_de_dfs.keys():
                sig_df = gene_de_dfs[key]
                x = sig_df[(sig_df['FDR'] < .05) & (np.abs(sig_df['logFC']) > 1)]
                up_cnt = x[x['logFC'] > 0].shape[0]
                dwn_cnt = x[x['logFC'] < 0].shape[0]
                count_dict[contrast][(cell_type, 'Up')] = up_cnt
                count_dict[contrast][(cell_type, 'Down')] = dwn_cnt

    df = pd.DataFrame(count_dict).T
    df = df.astype(str)

    return df

