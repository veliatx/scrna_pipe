import requests
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
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
    elif (not row[f'Significant_{treatment}'] and row[f'Significant_{control}']):
        return 3
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
                    size='point_size', size_max=10, template='plotly_white',
                    labels={"logFC": "Log2(FoldChange)"},
                    hover_data=['Pathway'])

    fig.update_layout(legend_font=dict(size=18))

    return fig


def plot_gene_volcano(plot_gene_df):
    """
    """
    fig = px.scatter(plot_gene_df, x='logFC', y='-Log10(FDR)', opacity=0.5,
                    color="Significant", color_discrete_map={True: "blue", False: "red"},
                    size='point_size', size_max=10, template='plotly_white',
                    labels={"logFC": "Log2(FoldChange)"},
                    hover_data=['gene'], render_mode='webgl')

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
        df1 = df1[df1['logCPM'] > 5]
        df2 = df2[df2['logCPM'] > 5]
    else:
        df1[f'Significant'] = df1.apply(lambda x: x.FDR < .05, axis=1)
        df2[f'Significant'] = df2.apply(lambda x: x.FDR < .05, axis=1)
        
    scatter_df = df1.merge(df2, left_on=diff_type, right_on=diff_type, suffixes=(f'_{treatment}', f'_{control}'))

    scatter_df['size'] = scatter_df.apply(lambda x: assign_size(x, treatment, control), axis=1)
    #scatter_df['opacity'] = scatter_df.apply(lambda x:  assign_opacity(x, treatment, control), axis=1)
    
    scatter = px.scatter(scatter_df, x=f'logFC_{treatment}', y=f'logFC_{control}', color=f'Significant_{treatment}', 
                        color_discrete_map={True: "blue", False: "red"},
                        template='plotly_white', opacity=.6, size=scatter_df['size'],
                        hover_data=[diff_type, f'Significant_{treatment}', f'Significant_{control}'])

    scatter.update_layout(title=f"Contrast {diff_type} scatter")

    scatter_df = scatter_df[scatter_df[f'Significant_{treatment}'] & ~scatter_df[f'Significant_{control}']]
    
    scatter_cols = [diff_type, f'logFC_{treatment}', f'FDR_{treatment}', f'Significant_{treatment}',
                        f'logFC_{control}', f'FDR_{control}', f'Significant_{control}']
    
    if diff_type == 'Pathway':
        scatter_cols += [f'leading_edge_{treatment}',
                         f'DE Genes In Pathway_{treatment}',
                         f'Total Genes In Pathway_{treatment}']
                         
    return scatter, scatter_df[scatter_cols]


def plot_umap(plot_df, color_map, title):
    """
    """
    fig = px.scatter(plot_df, x=0, y=1, opacity=0.5,
                    color="cell_type_display", color_discrete_map=color_map,
                    labels={
                        "0": "UMAP 2",
                        "1": "UMAP 1",
                    }, height=800)

    fig.update_layout(legend_font=dict(size=24), title=title)

    return fig


def plot_gene_boxplot(cpm_df, gene, cell_type, contrast=None):
    """
    """
    plot_df = pd.DataFrame(cpm_df.loc[gene])

    plot_df = np.log2(plot_df)
    plot_df.replace([np.inf, -np.inf, np.nan], .5, inplace=True)

    plot_df['label'] = plot_df.apply(lambda x: '_'.join(x.name.split('_')[1:-3]).rstrip('_'), axis=1)
    
    if contrast:
        plot_df = plot_df[plot_df.apply(lambda x: any(c in x.name for c in contrast), axis=1)]

    ax = sns.boxplot(data=plot_df, y=gene, x='label', palette='vlag')
    sns.stripplot(data=plot_df, y=gene, x='label', ax=ax)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=14)
    ax.set_ylabel(f'{gene} Expression - Log2(CPM)', fontsize=14)
    ax.set_title(f'{gene} expression - {cell_type} 5hr', fontsize=16)
    ax.figure.set_size_inches((8,6))

    return ax.figure

