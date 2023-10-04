import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad

import matplotlib.pyplot as plt

from itertools import product
from celltypist import models
from pathlib import Path


def load_count_data(output_dir, sample_sheet_path, keys=['PBMC', '5hr']):
    """
    """
    sample_df = pd.read_csv(sample_sheet_path)

    mask = sample_df['sample'].apply(lambda x: all(k in x for k in keys))

    sample_df = sample_df[mask]
    samples = list(sample_df['sample'])

    adatas = {}

    for sample in samples:
        sample_dir = output_dir.joinpath('samples', f'{sample}.filtered')
        adat = sc.read_mtx(sample_dir.joinpath('matrix.mtx')).transpose()
        adatas[sample] = adat

    sample_name = '_'.join(keys)
    adata = ad.concat(adatas, label="sample")

    feat_df = pd.read_csv(sample_dir.joinpath('features.tsv'), sep='\t', 
                      names=['ensembl_gene_id', 'hgnc_name', 'tag'])

    adata.var = feat_df

    adata.var.set_index('hgnc_name', inplace=True)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def load_count_data_2plate(output_dir, sample_sheet_path, 
    in_keys=['PBMC', '5hr'], out_keys=[]):
    """
    This function will merge the samples from the second distribution 
    plate. It assumes sample names from the first plate end with 'a'
    and sample names from the second plate end with 'b'.
    """
    sample_df = pd.read_csv(sample_sheet_path)

    sample_df = sample_df[sample_df['sample'].str.endswith('a')].copy()

    mask = sample_df['sample'].apply(lambda x: all(k in x for k in in_keys) \
                                       and not any(k in x for k in out_keys))

    sample_df = sample_df[mask]
    samples = list(sample_df['sample'])

    adatas = {}

    for sample in samples:
        sample_a_dir = output_dir.joinpath('samples', f'{sample}.filtered')
        sample_b_dir = output_dir.joinpath('samples', f'{sample[:-1]}b.filtered')
            
        adat_a = sc.read_mtx(sample_a_dir.joinpath('matrix.mtx')).transpose()
        adat_b = sc.read_mtx(sample_b_dir.joinpath('matrix.mtx')).transpose()
        adat = ad.concat([adat_a, adat_b])
        adatas[sample[:-3]] = adat

    sample_name = '_'.join(in_keys)
    adata = ad.concat(adatas, label="sample")

    feat_df = pd.read_csv(sample_a_dir.joinpath('features.tsv'), sep='\t', 
                      names=['ensembl_gene_id', 'hgnc_name', 'tag'])

    adata.var = feat_df

    adata.var.set_index('hgnc_name', inplace=True)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def filter_cells(adata, min_genes=200, max_counts=2e4, min_cells=3):
    """
    """
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_counts=max_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    adata.var['mt'] = adata.var.apply(lambda x: x.name.startswith('MT-'), axis=1)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt < 2, :]

    adata.raw = adata

    return adata


def run_celltypist(adata, model_name="Immune_All_High.pkl"):
    """
    """
    adata_celltypist = adata.copy()  # make a copy of our adata
    adata_celltypist.X = adata.raw.X  # set adata.X to raw counts
    sc.pp.normalize_per_cell(
        adata_celltypist, counts_per_cell_after=10**4
    )  # normalize to 10,000 counts per cell
    sc.pp.log1p(adata_celltypist)  # log-transform
    # make .X dense instead of sparse, for compatibility with celltypist:
    adata_celltypist.X = adata_celltypist.X.toarray()

    models.download_models(
        force_update=False, model=[model_name]
    )

    model = models.Model.load(model=model_name)

    predictions = celltypist.annotate(
        adata_celltypist, model=model, majority_voting=True
    )

    predictions_adata = predictions.to_adata()

    adata.obs["cell_type"] = predictions_adata.obs.loc[
        adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score"] = predictions_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]

    return adata


def dim_reduction(adata, min_dist=.4, spread=.9):
    """
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    ax_pca = sc.tl.pca(adata, svd_solver='arpack')

    sc.pl.pca_variance_ratio(adata, log=True)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

    sc.tl.leiden(adata, resolution=.4)

    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False)

    sc.tl.umap(adata, init_pos='paga', min_dist=min_dist, spread=spread)

    return adata, ax_pca