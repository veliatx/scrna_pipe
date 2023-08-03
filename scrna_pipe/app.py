import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc

from pathlib import Path
from scrna_pipe import differential
from streamlit_echarts import st_echarts
from streamlit_plotly_events import plotly_events
import plotly.express as px


@st.cache_data
def load_anndata(adata_path):
    """
    """
    adata = sc.read(adata_path)
    adata_dim = sc.read(adata_path.parent.joinpath(f'{adata_path.stem}_umap.h5ad'))
    adata_pb = differential.load_anndata_pseudobulk(adata_path, overwrite=False)

    return adata, adata_dim, adata_pb


output_dir = Path('/home/ec2-user/velia-analyses-dev/VAP_20230711_single_cell_moa')
adata_path = output_dir.joinpath('outputs', 'run_4', 'analysis', 'pbmc_5hr.h5ad')

adata, adata_dim, adata_pb = load_anndata(adata_path)

cell_types = list(set(adata.obs['cell_type']))
cell_types.insert(0, 'All')

with st.sidebar:
    cell_type = st.selectbox(
        'Choose a cell type',
        cell_types, index=0)


st.write('You selected:', cell_type)

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

fig = px.scatter(plot_df, x=0, y=1, opacity=0.5,
                 color="cell_type", color_discrete_map=color_map)

st.plotly_chart(fig, theme="streamlit")


