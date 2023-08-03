import warnings

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
import random
import sc_toolbox

import rpy2.rinterface_lib.callbacks
import anndata2ri
import logging

from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, numpy2ri
from scipy.sparse import csc_matrix

sc.settings.verbosity = 0
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)


def aggregate_and_filter(adata, cell_identity, obs_to_keep=[]):
    """
    """
    
    sample_key = 'sample'
    condition_key = 'label'
    cell_identity_key = 'cell_type'
    num_cell_per_sample = 30
    replicates_per_sample = 1
    
    # subset adata to the given cell identity
    print(cell_identity_key)
    print(adata.obs.columns)
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # check which samples to keep according to the number of cells specified with num_cell_per_sample
    size_by_sample = adata_cell_pop.obs.groupby([sample_key]).size()
    samples_to_drop = [
        sample
        for sample in size_by_sample.index
        if size_by_sample[sample] <= num_cell_per_sample
    ]
    if len(samples_to_drop) > 0:
        print("Dropping the following samples:")
        print(samples_to_drop)
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])

    adata_cell_pop.obs[sample_key] = adata_cell_pop.obs[sample_key].astype("category")
    for i, sample in enumerate(samples := adata_cell_pop.obs[sample_key].cat.categories):
        print(f"\tProcessing sample {i+1} out of {len(samples)}...", end="\r")
        if sample not in samples_to_drop:
            adata_sample = adata_cell_pop[adata_cell_pop.obs[sample_key] == sample]
            # create replicates for each sample
            indices = list(adata_sample.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_sample)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_sample[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each sample and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"
                # create a df with all genes, sample and condition info
                df_sample = pd.DataFrame(adata_replicate.X.A)
                df_sample.index = adata_replicate.obs_names
                df_sample.columns = adata_replicate.var_names
                df_sample = df_sample.join(adata_replicate.obs[obs_to_keep])
                # aggregate
                df_sample = df_sample.groupby(sample_key).agg(agg_dict)
                df_sample[sample_key] = sample
                df.loc[f"{sample}_{i}"] = df_sample.loc[sample]
    print("\n")
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names)
    )
    
    adata_cell_pop.X = csc_matrix(adata_cell_pop.X.astype(float))
    adata_cell_pop.obs['label'] = adata_cell_pop.obs['label'].astype('category')
    adata_cell_pop.obs['cell_type'] = adata_cell_pop.obs['cell_type'].astype('category')
    adata_cell_pop.obs['replicate'] = adata_cell_pop.obs['replicate'].astype('category')
    adata_cell_pop.obs['sample'] = adata_cell_pop.obs['sample'].astype('category')
    
    return adata_cell_pop


def load_anndata_pseudobulk(adata_path, overwrite=False):
    """
    """

    adata_pb_path = adata_path.parent.joinpath(f'{adata_path.stem}_pseudobulk.h5ad')
    
    if adata_pb_path.exists() and not overwrite:
        adata_pb = sc.read(adata_pb_path)
    
    else:
        adata = sc.read(adata_path)
        adata.obs['sample'] = adata.obs.apply(lambda x: x['sample'].replace('-','_'), axis=1)

        adata.layers["counts"] = adata.raw.X.copy()
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)

        adata.obs["cell_type"] = [ct.replace(" ", "_") for ct in adata.obs["cell_type"]]
        adata.obs["cell_type"] = [ct.replace("+", "") for ct in adata.obs["cell_type"]]
        adata.obs["replicate"] = adata.obs.apply(lambda x: x['sample'].split('_')[-1], axis=1)
        adata.obs["label"] = adata.obs.apply(lambda x: '_'.join(x['sample'].split('_')[:-1]), axis=1)

        adata.obs["replicate"] = adata.obs["replicate"].astype("category")
        adata.obs["label"] = adata.obs["label"].astype("category")
        adata.obs["sample"] = adata.obs["sample"].astype("category")
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

        obs_to_keep = ["label", "cell_type", "replicate", "sample"]
        adata.X = adata.layers["counts"].copy()

        # process first cell type separately...
        cell_type = adata.obs["cell_type"].cat.categories[0]
        print(
            f'Processing {cell_type} (1 out of {len(adata.obs["cell_type"].cat.categories)})...'
        )
        adata_pb = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
        for i, cell_type in enumerate(adata.obs["cell_type"].cat.categories[1:]):

            print(
                f'Processing {cell_type} ({i+2} out of {len(adata.obs["cell_type"].cat.categories)})...'
            )
            adata_cell_type = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
            adata_pb = adata_pb.concatenate(adata_cell_type)

        adata_pb.write_h5ad(adata_pb_path)

    return adata_pb


def main():
    importr('edgeR')


if __name__ == "__main__":
    main()



