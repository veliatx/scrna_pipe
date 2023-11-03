import numpy as np
import pandas as pd


def overview_de_count_table(contrast_tuples, gene_de_dfs, dataset_name):
    """
    """
    all_contrasts = [f'{x[0]}-{x[1]}' for x in contrast_tuples]
    gene_contrasts = [x for x in list(gene_de_dfs.keys()) if any([c in x for c in all_contrasts])]
    cell_types = [x.split('_')[-1] for x in gene_contrasts]

    count_dict = {}

    for contrast in all_contrasts:
        if 'None' in contrast: continue
        count_dict[contrast] = {}
        for cell_type in cell_types:
            key = f'{dataset_name}__edgeR__{contrast}__{cell_type}'

            if key in gene_de_dfs.keys():

                sig_df = gene_de_dfs[key]
                x = sig_df[(sig_df['FDR'] < .05) &
                           (np.abs(sig_df['logFC']) > 1) &
                           (sig_df['logCPM'] > 5)]
                up_cnt = x[x['logFC'] > 0].shape[0]
                dwn_cnt = x[x['logFC'] < 0].shape[0]
                count_dict[contrast][(cell_type, 'Up')] = up_cnt
                count_dict[contrast][(cell_type, 'Down')] = dwn_cnt

    df = pd.DataFrame(count_dict).T
    df = df.astype(float)

    return df


def overview_cell_type_table(adata_dim):
    """
    """
    cell_cnt_df = pd.DataFrame(adata_dim.obs.groupby(['label', 'cell_type']).size())
    cell_cnt_df.reset_index(inplace=True)
    cell_cnt_pivot_df = cell_cnt_df.pivot_table(values=0, columns='cell_type', index='label')
    cell_cnt_pivot_df = cell_cnt_pivot_df[cell_cnt_pivot_df.describe().loc['mean'].sort_values(ascending=False).index]
    cell_cnt_pivot_df['Total Cells'] = cell_cnt_pivot_df.sum(axis=1)
    cell_cnt_pivot_df.sort_values(by='Total Cells', ascending=False, inplace=True)

    return cell_cnt_pivot_df