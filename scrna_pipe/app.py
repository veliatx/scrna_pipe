import altair as alt
import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc

from pathlib import Path
from scrna_pipe import differential
from streamlit_echarts import st_echarts
from streamlit_plotly_events import plotly_events
import plotly.express as px

st.set_page_config(layout="wide")


@st.cache_data
def load_anndata(adata_path):
    """
    """
    adata = sc.read(adata_path)
    adata_dim = sc.read(adata_path.parent.joinpath(f'{adata_path.stem}_umap.h5ad'))
    
    adata_pb = differential.load_anndata_pseudobulk(adata_path, overwrite=False)
    adata_pb = differential.process_adata_pbmc(adata_pb)    

    return adata, adata_dim, adata_pb


@st.cache_data
def load_differential_dfs(adata_path):
    gsva_dfs = {}
    for df_path in adata_path.parent.joinpath('gsva').glob('*.csv'):
        gsva_dfs[df_path.stem] = pd.read_csv(df_path) 

    return gsva_dfs


output_dir = Path('/home/ec2-user/velia-analyses-dev/VAP_20230711_single_cell_moa')
adata_path = output_dir.joinpath('outputs', 'run_5', 'analysis', 'pbmc_5hr.h5ad')

focus_contrasts = [
    ('None', 'None'),
    ('PBMC_5hr_LPS', 'PBMC_5hr_Mock'),
    ('PBMC_5hr_sORF2184_0', 'PBMC_5hr_Mock'), 
    ('PBMC_5hr_sORF2341_0', 'PBMC_5hr_Mock'),
    ('PBMC_5hr_LPS_sORF2184_0', 'PBMC_5hr_LPS'),
    ('PBMC_5hr_LPS_sORF2341_0', 'PBMC_5hr_LPS'),
    ('PBMC_5hr_sORF2341_0', 'PBMC_5hr_sORF2184_0'),
    ('PBMC_5hr_LPS_sORF2341_0', 'PBMC_5hr_LPS_sORF2184_0'),
]

adata, adata_dim, adata_pb = load_anndata(adata_path)
gsva_dfs = load_differential_dfs(adata_path)

cell_types = list(set(adata.obs['cell_type']))
cell_types.sort()
cell_types.insert(0, 'All')

with st.sidebar:

    cell_type = st.selectbox(
        'Choose a cell type',
        cell_types, index=0
    )

    st.divider()

    contrast = st.selectbox(
        'Choose a contrast',
        [f'{x[0]} - {x[1]}' for x in focus_contrasts]
    )

color_map =  {
    "Macrophages": "#1f77b4",  
    "B cells": "#ff7f0e", 
    "T cells": "#2ca02c", 
    "Monocytes":  "#d62728",  
    "pDC": "#9467bd", 
    "ILC": "#8c564b", 
}

if cell_type and cell_type != 'All':
    set_color = color_map[cell_type]
    color_map = {c: "#D3D3D3" for c in color_map.keys()}
    color_map[cell_type] = set_color

plot_df = pd.DataFrame(adata_dim.obsm['X_umap'])
cell_df = pd.DataFrame(adata_dim.obs['cell_type'])
cell_df.reset_index(inplace=True)
plot_df = pd.concat([plot_df, cell_df], axis=1)


with st.expander(label='Cell Annotation', expanded=True):
    fig = px.scatter(plot_df, x=0, y=1, opacity=0.5,
                    color="cell_type", color_discrete_map=color_map,
                    labels={
                        "0": "UMAP 2",
                        "1": "UMAP 1",
                    },)

    st.plotly_chart(fig, theme="streamlit")


with st.expander(label='Differential Expression', expanded=True):
    if cell_type != 'All' and contrast != 'None - None':
        cell_type = cell_type.replace(' ', '')
        contrast_map = {f'{x[0]} - {x[1]}': f'{x[0]}-{x[1]}' for x in focus_contrasts}
        key = f'{adata_path.stem}_gsva_{contrast_map[contrast]}_{cell_type}'


        plot_df = gsva_dfs[key]
        plot_df.rename(columns={'Unnamed: 0': 'Pathway'}, inplace=True)
        plot_df['-Log10(FDR)'] = -1*np.log10(plot_df['adj.P.Val'])
        plot_df['significant'] = plot_df.apply(lambda x: x['adj.P.Val'] < .05, axis=1)

        plot_df['FDR'] = plot_df['adj.P.Val'].astype(float)

        alt.data_transformers.disable_max_rows()

        col1, col2 = st.columns(2)

        with col1:
            st.altair_chart(alt.Chart(plot_df).mark_circle(size=60).encode(
                x='logFC',
                y='-Log10(FDR)',
                color='significant',
                tooltip=['Pathway', 'FDR', 'logFC']
            ))

        with col2:
            st.dataframe(plot_df[plot_df['significant']].sort_values(by='FDR'))
