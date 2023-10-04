import altair as alt
import streamlit as st
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc

from pathlib import Path
from scrna_pipe import differential
from streamlit_echarts import st_echarts
from streamlit_plotly_events import plotly_events

import plotly.express as px
import plotly.graph_objects as go

from scrna_pipe import plotting

st.set_page_config(layout="wide")


@st.cache_data
def load_anndata(adata_path):
    """
    """
    adata = sc.read(adata_path)
    adata_dim = sc.read(adata_path.parent.joinpath(f'{adata_path.stem}_umap.h5ad'))
    
    #adata_pb = differential.load_anndata_pseudobulk(adata_path, overwrite=False)
    #adata_pb = differential.process_adata_pbmc(adata_pb)    

    return adata, adata_dim#, adata_pb


@st.cache_data
def load_differential_dfs(adata_path):
    gene_de_dfs = {}
    gsva_de_dfs = {}

    for df_path in adata_path.parent.joinpath('edgeR').glob('*.csv'):
        gene_de_dfs[df_path.stem] = pd.read_csv(df_path, index_col=0) 

    for df_path in adata_path.parent.joinpath('gsva').glob('*.csv'):
        gsva_de_dfs[df_path.stem] = pd.read_csv(df_path, index_col=0) 

    return gene_de_dfs, gsva_de_dfs


def convert_df(df):
   return df.to_csv(index=False).encode('utf-8')


#output_dir = Path('/home/ubuntu/scrna_pipe/data')
#output_dir = Path('/home/ec2-user/scrna_pipe/data')
#output_dir = Path('/home/ec2-user/velia-analyses-dev/VAP_20230711_single_cell_moa/outputs/run_6')
output_dir = Path('/home/ec2-user/velia-analyses-dev/VAP_20230919_single_cell_pbmc_hits/outputs/run_2')


adata_paths = {
    #'PBMC 5hr': output_dir.joinpath('analysis', '5hr_coarse.h5ad'),
    #'PBMC 24hr': output_dir.joinpath('analysis', '24hr_coarse.h5ad'),
    '5hr': output_dir.joinpath('analysis', '5hr_coarse.h5ad'),
    '24hr': output_dir.joinpath('analysis', '24hr_coarse.h5ad'),
    #'PBMC': output_dir.joinpath('analysis', 'PBMC_coarse.h5ad'),
    #'HCT116': output_dir.joinpath('outputs', 'run_5', 'analysis', 'HCT116.h5ad'),
    #'A549': output_dir.joinpath('outputs', 'run_5', 'analysis', 'A549.h5ad'),
}

focus_contrasts = {
    'PBMC 5hr': [
        ('None', 'None'),
        ('PBMC_5hr_LPS', 'PBMC_5hr_Mock'),
        ('PBMC_5hr_sORF2184_0', 'PBMC_5hr_Mock'), 
        ('PBMC_5hr_sORF2341_0', 'PBMC_5hr_Mock'),
        ('PBMC_5hr_LPS_sORF2184_0', 'PBMC_5hr_LPS'),
        ('PBMC_5hr_LPS_sORF2341_0', 'PBMC_5hr_LPS'),
        ('PBMC_5hr_sORF2341_0', 'PBMC_5hr_sORF2184_0'),
        ('PBMC_5hr_LPS_sORF2341_0', 'PBMC_5hr_LPS_sORF2184_0'),
    ],
    'PBMC 24hr': [
        ('PBMC_24hr_LPS', 'PBMC_24hr_Mock'),
        ('PBMC_24hr_sORF2184_0', 'PBMC_24hr_Mock'), 
        ('PBMC_24hr_sORF2341_0', 'PBMC_24hr_Mock'),
        ('PBMC_24hr_LPS_sORF2184_0', 'PBMC_24hr_LPS'),
        ('PBMC_24hr_LPS_sORF2341_0', 'PBMC_24hr_LPS'),
        ('PBMC_24hr_sORF2341_0', 'PBMC_24hr_sORF2184_0'),
        ('PBMC_24hr_LPS_sORF2341_0', 'PBMC_24hr_LPS_sORF2184_0'),
        ('PBMC_24hr_LPS_sORF2406', 'PBMC_24hr_LPS'),
    ],
    'PBMC': [
        ('None', 'None'),
        ('PBMC_5hr_LPS', 'PBMC_5hr_Mock'),
        ('PBMC_5hr_sORF2184_0', 'PBMC_5hr_Mock'), 
        ('PBMC_5hr_sORF2341_0', 'PBMC_5hr_Mock'),
        ('PBMC_5hr_LPS_sORF2184_0', 'PBMC_5hr_LPS'),
        ('PBMC_5hr_LPS_sORF2341_0', 'PBMC_5hr_LPS'),
        ('PBMC_5hr_sORF2341_0', 'PBMC_5hr_sORF2184_0'),
        ('PBMC_5hr_LPS_sORF2341_0', 'PBMC_5hr_LPS_sORF2184_0'),
        ('PBMC_24hr_LPS', 'PBMC_24hr_Mock'),
        ('PBMC_24hr_sORF2184_0', 'PBMC_24hr_Mock'), 
        ('PBMC_24hr_sORF2341_0', 'PBMC_24hr_Mock'),
        ('PBMC_24hr_LPS_sORF2184_0', 'PBMC_24hr_LPS'),
        ('PBMC_24hr_LPS_sORF2341_0', 'PBMC_24hr_LPS'),
        ('PBMC_24hr_sORF2341_0', 'PBMC_24hr_sORF2184_0'),
        ('PBMC_24hr_LPS_sORF2341_0', 'PBMC_24hr_LPS_sORF2184_0'),
        ('PBMC_24hr_LPS_sORF2406', 'PBMC_24hr_LPS'),
        ('PBMC_5hr_Mock', 'PBMC_24hr_Mock'),
        ('PBMC_5hr_LPS', 'PBMC_24hr_LPS'),
        ('PBMC_5hr_sORF2184_0', 'PBMC_24hr_sORF2184_0'),
        ('PBMC_5hr_LPS_sORF2184_0', 'PBMC_24hr_LPS_sORF2184_0'),
    ],
    '5hr': [
        ('None', 'None'),
        ('Mock__R848_5hr_', 'Mock__Fc_1uM_5hr_'),
        ('VTX0851359__5hr_', 'Mock__Fc_1uM_5hr_'),
        ('VTX0851359__R848_5hr_', 'Mock__R848_5hr_'), 
        ('IL10__5hr_', 'Mock__Fc_1uM_5hr_'),
        ('IL10__R848_5hr_', 'Mock__R848_5hr_'),
        ('VTX0852555__5hr_', 'Mock__Fc_1uM_5hr_'),
        ('VTX0852555__R848_5hr_', 'Mock__R848_5hr_'),
        ('VTX0851359__R848_5hr_', 'IL10__R848_5hr_'), 
        ('VTX0851359__5hr_', 'IL10__5hr_'),

        ('VTX0852555__5hr_', 'IL10__5hr_'),
        ('VTX0852555__R848_5hr_', 'IL10__R848_5hr_'),
        
        ('VTX0852488__5hr_', 'IL10__5hr_'),
        ('VTX0852488__R848_5hr_', 'IL10__R848_5hr_'),
        
        ('VTX0851359__5hr_', 'VTX0852555__5hr_'),
        ('VTX0851359__R848_5hr_', 'VTX0852555__R848_5hr_'),
    ],
    '24hr': [
        ('None', 'None'),
        ('Mock__R848_24hr_', 'Mock__Fc_1uM_24hr_'),
        ('VTX0851359__24hr_', 'Mock__Fc_1uM_24hr_'),
        ('VTX0851359__R848_24hr_', 'Mock__R848_24hr_'),
        ('IL10__24hr_', 'Mock__Fc_1uM_24hr_'),
        ('IL10__R848_24hr_', 'Mock__R848_24hr_'),
        ('VTX0852555__24hr_', 'Mock__Fc_1uM_24hr_'),
        ('VTX0852555__R848_24hr_', 'Mock__R848_24hr_'),        
    ],
}


datasets = ['5hr']#'PBMC', 'PBMC 5hr', 'PBMC 24hr']#, 'A549', 'HCT116']

with st.sidebar:

    dataset = st.selectbox(
        'Choose a dataset',
        datasets, index=0,
    )
    adata, adata_dim = load_anndata(adata_paths[dataset])
    
    st.divider()

    cell_types = ['Macrophages', 'T cells', 'B cells', 'Monocytes', 'ILC', 'DC']
    #cell_types = list(set(adata.obs['cell_type']))
    cell_types.sort()
    cell_types.insert(0, 'All')

    cell_type = st.selectbox(
        'Choose a cell type',
        cell_types, index=0
    )

    st.divider()

    contrast = st.selectbox(
        'Choose a contrast',
        [f'{x[0]} - {x[1]}' for x in focus_contrasts[dataset]]
    )


gene_de_dfs, gsva_de_dfs = load_differential_dfs(adata_paths[dataset])

color_map =  {
    "Macrophages": "#1f77b4",  
    "B cells": "#ff7f0e", 
    "T cells": "#2ca02c", 
    "Monocytes":  "#d62728",  
    "pDC": "#9467bd", 
    "ILC": "#8c564b",
    "DC": "#e377c2",
    "Endothelial cells": "#7f7f7f",
    "Plasma cells": "#bcbd22",
    "ETP": "#1f77b4",
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

    fig.update_layout(legend_font=dict(size=24))

    st.plotly_chart(fig, theme="streamlit")


with st.expander(label='Differential Expression', expanded=True):
    if cell_type != 'All' and contrast != 'None - None':
        cell_type = cell_type.replace(' ', '')
        contrast_map = {f'{x[0]} - {x[1]}': f'{x[0]}-{x[1]}' for x in focus_contrasts[dataset]}
        
        gsva_key = f'{adata_paths[dataset].stem}__gsva__{contrast_map[contrast]}__{cell_type.replace("_", "")}'
        edger_key = f'{adata_paths[dataset].stem}__edgeR__{contrast_map[contrast]}__{cell_type.replace("_", "")}'
        
        st.write(edger_key)
        if gsva_key in gsva_de_dfs.keys():

            plot_gene_df = gene_de_dfs[edger_key]
            plot_gene_df['point_size'] = 5
            plot_gene_df['Significant'] = plot_gene_df.apply(lambda x: abs(x.logFC) > 1 and x.FDR < .05, axis=1)

            plot_gsva_df = gsva_de_dfs[gsva_key]
            plot_gsva_df['point_size'] = 5
                        
            col1, col2 = st.columns(2)

            with col1:
                st.subheader('Gene Volcano Plot')
                select_gene_df = plot_gene_df.copy()

                fig = px.scatter(select_gene_df, x='logFC', y='-Log10(FDR)', opacity=0.5,
                                 color="Significant", size='point_size', size_max=5, template='plotly_white',
                                 labels={"logFC": "Log2(FoldChange)"},
                                 hover_data=['gene'])


                fig.update_layout(legend_font=dict(size=18))
                
                selected_points = plotly_events(fig)


            with col2:
                st.subheader('Gene DE Table')
                gene_cols = ['gene', 'logFC', 'FDR', 'logCPM', 
                             'PValue', '-Log10(FDR)', 'F']


                sig_df = plot_gene_df[plot_gene_df['Significant']].sort_values(by='FDR')
                st.dataframe(sig_df[gene_cols].style.format({"FDR": "{:.2E}", "PValue": "{:.2E}"}))
                csv = convert_df(sig_df[gene_cols])

                st.download_button(
                    "Download Table",
                    csv,
                    f"gene_{edger_key}.csv",
                    "text/csv",
                    key='download-csv-gene'
                )

            x = contrast.split(' - ')
            try:
                plot_df = plotting.extract_cluster_map(plot_gene_df.copy(), adata, cell_type, x)

            
                heatmap = go.Figure(data=go.Heatmap(z=plot_df.values, 
                                                    x=plot_df.columns, 
                                                    y=plot_df.index,
                                                    colorscale='Viridis'))
                heatmap.update_layout(xaxis=dict(nticks=len(plot_df.columns)))

                st.plotly_chart(heatmap, theme="streamlit", use_container_width=True)
            except:
                None


            if selected_points:
                st.subheader('UMAP projection of selected gene')
                
                select_gene_df = select_gene_df[select_gene_df['Significant']]

                gene = select_gene_df.iloc[selected_points[0]["pointIndex"]]['gene']

                c1, c2 = contrast.split('-')
                
                c1 = c1.replace('_', '-').strip()
                c2 = c2.replace('_', '-').strip()

                samples1 = [f'{c1}-{i}' for i in range(1, 4)]
                samples2 = [f'{c2}-{i}' for i in range(1, 4)]
                col3, col4 = st.columns(2)

                try:
                    gene_idx = adata_dim.var.index.get_loc(gene)
                    max_val = np.percentile(adata_dim.X[:, gene_idx], q=98)
                except:
                    max_val = 5

                with col3:
                    st.write(c1)
                    ax1 = sc.pl.umap(
                        adata_dim[adata_dim.obs['sample'].isin(samples1)],
                        color=[gene],
                        frameon=False,
                        sort_order=False,
                        wspace=1,
                        return_fig=True,
                        vmin=0.1,
                        vmax=max_val,
                        cmap='viridis'

                    )

                    st.pyplot(ax1)

                    st.write(gene)

                with col4:
                    st.write(c2)
                    ax2 = sc.pl.umap(
                        adata_dim[adata_dim.obs['sample'].isin(samples2)],
                        color=[gene],
                        frameon=False,
                        sort_order=False,
                        wspace=1,
                        return_fig=True,
                        vmin=0.1,
                        vmax=max_val,
                        cmap='viridis'

                    )

                    st.pyplot(ax2)

                    st.write(gene) 


            string_button = st.button('Generate STRING network diagrams')

            colx, coly, colz = st.columns(3)

            
            if string_button:
                with colx:
                    st.write('STRING network of top up-regulated genes')
                    network_image = plotting.plot_string_network(list(sig_df[sig_df['logFC'] > 0]['gene'])[0:100])
                    st.image(network_image)

                with coly:
                    st.write('STRING network of top down-regulated genes')
                    network_image = plotting.plot_string_network(list(sig_df[sig_df['logFC'] < 0]['gene'])[0:100])
                    st.image(network_image)

                with colz:
                    st.write('STRING network of top regulated genes')
                    network_image = plotting.plot_string_network(list(sig_df['gene'])[0:100])
                    st.image(network_image)



with st.expander(label='Pathway analysis', expanded=True):
    if cell_type != 'All' and contrast != 'None - None':
        cell_type = cell_type.replace(' ', '')
        contrast_map = {f'{x[0]} - {x[1]}': f'{x[0]}-{x[1]}' for x in focus_contrasts[dataset]}
        
        gsva_key = f'{adata_paths[dataset].stem}__gsva__{contrast_map[contrast]}__{cell_type.replace("_", "")}'
        
        st.write(edger_key)
        if gsva_key in gsva_de_dfs.keys():

            plot_gsva_df = gsva_de_dfs[gsva_key]
            plot_gsva_df['point_size'] = 10
                        
            col5, col6 = st.columns(2)

            with col5:
                st.subheader('Pathway Volcano Plot')
                fig = px.scatter(plot_gsva_df, x='logFC', y='-Log10(FDR)', opacity=0.5,
                                 color="Significant", size='point_size', size_max=5, template='plotly_white',
                                 labels={"logFC": "Log2(FoldChange)"},
                                 hover_data=['Pathway'])

                fig.update_layout(legend_font=dict(size=18))
                selected_pathway = plotly_events(fig)

                #st.plotly_chart(fig, theme="streamlit")


            with col6:
                st.subheader('Pathway DE Table')
                pathway_cols = ['Pathway', 'FDR', 'leading_edge', 'DE Genes In Pathway', 
                                'Total Genes In Pathway', '-Log10(FDR)', 'logFC', 'msigdbURL']
                sig_df = plot_gsva_df[plot_gsva_df['Significant']].sort_values(by='FDR')
                st.dataframe(sig_df[pathway_cols].style.format({"FDR": "{:.2E}"}))
                csv = convert_df(sig_df[pathway_cols])
                st.download_button(
                    "Download Table",
                    csv,
                    f"pathways_{gsva_key}.csv",
                    "text/csv",
                    key='download-csv-pathway'
                )

        else:
            st.write('Not enough cells to perform differential analysis.')