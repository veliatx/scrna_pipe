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

from rpy2 import robjects
from rpy2.robjects import r, pandas2ri, numpy2ri
from rpy2.robjects.vectors import StrVector, ListVector
from rpy2.robjects.packages import importr
from scipy.sparse import csc_matrix

from pathlib import Path

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


def process_adata_pbmc(adata_pb):
    """
    """
    adata_pb.obs['cell_type'] = adata_pb.obs.apply(lambda x: x.cell_type.replace('_', ''), axis=1)
    adata_pb.obs['new_index'] = adata_pb.obs.apply(lambda x: x.cell_type + '_' + '_'.join(x.name.split('_')[:-1]), axis=1)
    adata_pb.obs.set_index('new_index', inplace=True)

    adata_pb.layers['counts'] = adata_pb.X.copy()

    sc.pp.normalize_total(adata_pb, target_sum=1e6)
    sc.pp.log1p(adata_pb)
    sc.pp.pca(adata_pb)

    adata_pb.obs["lib_size"] = np.sum(adata_pb.layers["counts"], axis=1)
    adata_pb.obs["lib_size"] = adata_pb.obs["lib_size"].astype(float)
    adata_pb.obs["log_lib_size"] = np.log(adata_pb.obs["lib_size"])
    adata_pb.X = adata_pb.layers['counts'].copy()

    return adata_pb


def process_adata_cancer(adata_pb):
    """
    """
    adata_pb.obs['cell_type'] = adata_pb.obs.apply(lambda x: x.cell_type.replace('_', ''), axis=1)

    adata_pb.obs.index.name = 'new_index'
    adata_pb.layers['counts'] = adata_pb.X.copy()

    sc.pp.normalize_total(adata_pb, target_sum=1e6)
    sc.pp.log1p(adata_pb)
    sc.pp.pca(adata_pb)

    adata_pb.obs["lib_size"] = np.sum(adata_pb.layers["counts"], axis=1)
    adata_pb.obs["lib_size"] = adata_pb.obs["lib_size"].astype(float)
    adata_pb.obs["log_lib_size"] = np.log(adata_pb.obs["lib_size"])
    adata_pb.X = adata_pb.layers['counts'].copy()

    return adata_pb


def r_fit_model(adata_pb):
    """
    """

    importr('edgeR')
    importr('base')
    importr('stats')
    importr('limma')
    importr('GSVA')
    
    fit_model_code = """
    fit_model <- function(adata_){
        # create an edgeR object with counts and grouping factor
        y <- DGEList(assay(adata_, "X"), group = colData(adata_)$label)
        # filter out genes with low counts
        print("Dimensions before subsetting:")
        print(dim(y))
        print("")
        keep <- filterByExpr(y)
        y <- y[keep, , keep.lib.sizes=FALSE]
        print("Dimensions after subsetting:")
        print(dim(y))
        print("")
        # normalize
        y <- calcNormFactors(y)
        # create a vector that is concatentation of condition and cell type that we will later use with contrasts
        group <- paste0(colData(adata_)$label, ".", colData(adata_)$cell_type)
        replicate <- colData(adata_)$replicate
        # create a design matrix: here we have multiple donors so also consider that in the design matrix
        design <- model.matrix(~ 0 + group + replicate)
        # estimate dispersion
        y <- estimateDisp(y, design = design)
        # fit the model
        fit <- glmQLFit(y, design)
        return(list("fit"=fit, "design"=design, "y"=y))
    }
    """

    run_model_code = """
    outs <- fit_model(adata_pb)

    adata_pb_fit <- outs$fit
    adata_pb_y <- outs$y
    """

    cpm_model_code = """
    cpm_df <- edgeR::cpm(adata_pb_y)
    rownames(cpm_df) <- rownames(adata_pb_y)
    cpm_py_df <- as.data.frame(cpm_df)
    """

    r_adata_pb = anndata2ri.py2rpy(adata_pb)
    robjects.r.assign("adata_pb", r_adata_pb)

    _ = robjects.r(fit_model_code)

    _ = robjects.r(run_model_code)

    _ = robjects.r(cpm_model_code)

    cpm_df = pandas2ri.rpy2py_dataframe(robjects.r['cpm_py_df'])

    return cpm_df


def load_gene_set():
    """
    """
    msig_dir = Path('/home/ec2-user/velia-data-dev/VDC_004_annotation/msigdb/v2023.1_json/')

    pathway_df = pd.read_json(msig_dir.joinpath('c2.cp.v2023.1.Hs.json')).T

    gene_sets = {}
    for i, row in pathway_df.iterrows():
        gene_sets[row.name] = StrVector(row.geneSymbols)

    return gene_sets


def r_run_gsva(adata_pb, cpm_df, focus_contrasts, adata_path, gene_sets):
    """
    """
            
    r_list = ListVector(gene_sets)

    r_list = robjects.r.assign("gene_sets", r_list)

    gsva_code = """
    gsva_df <- gsva(cpm_df, gene_sets, method='gsva')
    """

    _ = robjects.r(gsva_code)

    contrast_str = ''
    contrast_list = []

    all_contrasts = set(['_'.join(x.split('_')[0:-1]) for x in cpm_df.columns])
    cell_types = set(adata_pb.obs['cell_type'])

    for c1, c2 in focus_contrasts:
        for cell_type in cell_types:
            c1_str = f'{cell_type}_{c1}'
            c2_str = f'{cell_type}_{c2}'
            if c1_str in all_contrasts and c2_str in all_contrasts:
                contrast_str += f'"{c1_str}-{c2_str}" = conditions{c1_str} - conditions{c2_str},\n '
                contrast_list.append(f'"{c1_str}-{c2_str}"')
                
    conditions = StrVector(['_'.join(x.split('_')[0:-1]) for x in cpm_df.columns])
    robjects.r.assign("conditions", conditions)

    # Construct design matrix
    _ = robjects.r('''
    design <- model.matrix(~ 0 + conditions)
    ''')

    # Create contrasts for pairwise comparisons
    _ = robjects.r(f'''
    contrast_matrix <- makeContrasts(
        {contrast_str}
        levels = design
    )
    ''')

    _ = robjects.r(f'''
    limma_fit <- lmFit(gsva_df, design)
    limma_fit2 <- contrasts.fit(limma_fit, contrast_matrix)
    limma_fit2 <- eBayes(limma_fit2)
    ''')

    con_names = robjects.r('colnames(contrast_matrix)')

    gsva_dfs = {}

    gsva_data_path = adata_path.parent

    for cell_type in cell_types:
        gsva_dfs[cell_type] = {}
        for contrast in focus_contrasts:
            
            limma_contrast = f'{cell_type}_{contrast[0]}-{cell_type}_{contrast[1]}'
            
            if limma_contrast not in con_names: 
                continue

            _ = robjects.r(f'tt <- topTable(limma_fit2, coef="{limma_contrast}" , n=Inf)')

            df = pandas2ri.rpy2py_dataframe(robjects.r['tt'])
            gsva_dfs[cell_type][limma_contrast] = df
            
            gsva_file_name = f'{adata_path.stem}_gsva_{contrast[0]}-{contrast[1]}_{cell_type}.csv'

            df.to_csv(gsva_data_path.joinpath('gsva', gsva_file_name))

    return gsva_dfs


def r_run_edger(adata_pb, focus_contrasts, adata_path):
    """
    """
            
    fit_cols = robjects.r('fit_cols <- colnames(adata_pb_fit)')

    contrast_str = ''
    contrast_map = {}
    cell_types = set(adata_pb.obs['cell_type'])
    cons_fine = [x[5:] for x in fit_cols]

    for c1, c2 in focus_contrasts:

        for cell_type in cell_types:
            c1_str = f'{c1}.{cell_type}'
            c2_str = f'{c2}.{cell_type}'
            if c1_str in cons_fine and c2_str in cons_fine:
                contrast_str += f'"{c1_str}-{c2_str}" = fit_colsgroup{c1_str} - fit_colsgroup{c2_str},\n '
                contrast_map[f"{c1_str}-{c2_str}"] = f"group{c1_str} - group{c2_str}"
            

    # Construct design matrix
    _ = robjects.r('''
    design_gene <- model.matrix(~ 0 + fit_cols)
    ''')

    # Create contrasts for pairwise comparisons
    _ = robjects.r(f'''
    contrast_matrix_gene <- makeContrasts(
        {contrast_str}
        levels = design_gene
    )
    ''')

    edgeR_dfs = {}

    for contrast_name, contrast in contrast_map.items():

        edgeR_code = f"""

        contrast2 <- makeContrasts("{contrast}", levels=colnames(adata_pb_fit$design))

        qlf <- glmQLFTest(adata_pb_fit, contrast=contrast2)
        tt_gene <- topTags(qlf, n = Inf)
        tt_table = tt_gene$table
        """
        _ = robjects.r(edgeR_code)

        df = pandas2ri.rpy2py_dataframe(robjects.r['tt_table'])
        
        edgeR_dfs[contrast_name] = df
        
        edgeR_file_name = f'{adata_path.stem}_edgeR_{contrast_name}.csv'

        outpath = adata_path.parent.joinpath('edgeR')
        if not outpath.exists():
            outpath.mkdir()

        df.to_csv(outpath.joinpath('edgeR_file_name'))
    
    return edgeR_dfs


def load_gsva_dfs(adata_path):
    """
    """
    pass


def main():
    """
    """

    adata_pb = load_anndata_pseudobulk(adata_path, overwrite=False)

    adata_pb = process_adata_pbmc(adata_pb)

    r_adata_pb = anndata2ri.py2rpy(adata_pb)
    robjects.r.assign("adata_pb", r_adata_pb)




if __name__ == "__main__":
    main()



