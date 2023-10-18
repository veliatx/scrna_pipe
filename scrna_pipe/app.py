import altair as alt
import decoupler as dc
import streamlit as st
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc

from pathlib import Path
from streamlit_echarts import st_echarts
from streamlit_plotly_events import plotly_events

import plotly.graph_objects as go

from scrna_pipe import plotting
from veliadb.base import Session, Gene

st.set_page_config(layout="wide")


@st.cache_data
def load_anndata(adata_path):
    """
    """
    adata = sc.read(adata_path)
    adata_dim = sc.read(adata_path.parent.joinpath(f'{adata_path.stem}_umap.h5ad'))

    return adata, adata_dim#, adata_pb


@st.cache_data
def load_differential_dfs(adata_path):
    """
    """
    gene_de_dfs = {}
    gsva_de_dfs = {}

    for df_path in adata_path.parent.joinpath('edgeR').glob('*.csv'):
        gene_de_dfs[df_path.stem] = pd.read_csv(df_path, index_col=0) 

    for df_path in adata_path.parent.joinpath('gsva').glob('*.csv'):
        gsva_de_dfs[df_path.stem] = pd.read_csv(df_path, index_col=0) 

    return gene_de_dfs, gsva_de_dfs


@st.cache_data
def convert_df(df):
   return df.to_csv(index=False).encode('utf-8')


@st.cache_data
def load_msigdb(pathway_db='reactome'):
    """
    """
    msigdb = dc.get_resource('MSigDB')
    msigdb_subset = msigdb[msigdb['collection']==pathway_db]
    msigdb_subset = msigdb_subset[~msigdb_subset.duplicated(['geneset', 'genesymbol'])]

    return msigdb_subset


def load_gencode_map():
    """
    """
    session = Session()
    gencode_map = {g.hgnc_name: g.ensembl_id \
                    for g in session.query(Gene).\
                                     filter(Gene.ensembl_id != '').all()}
    session.close()
    return gencode_map


#output_dir = Path('/home/ubuntu/scrna_pipe/data')
output_dir = Path('/home/ec2-user/scrna_pipe/data')
output_dir_plate1 = Path('/home/ec2-user/velia-analyses-dev/VAP_20230711_single_cell_moa/outputs/run_6')
output_dir_plate2 = Path('/home/ec2-user/velia-analyses-dev/VAP_20230919_single_cell_pbmc_hits/outputs/run_2')


adata_paths = {
    'PBMC - plate 1': output_dir.joinpath('analysis_plate_1', 'PBMC_coarse.h5ad'),
    #'PBMC 24hr': output_dir.joinpath('analysis', '24hr_coarse.h5ad'),
    'PBMC 5hr - plate 2': output_dir.joinpath('analysis_plate_2', '5hr_coarse_subset.h5ad'),
    #'24hr': output_dir.joinpath('analysis', '24hr_coarse.h5ad'),
}

focus_contrasts = {

    'PBMC - plate 1': [
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
    'PBMC 5hr - plate 2': [
        ('None', 'None'),
        ('Mock__R848_5hr_', 'Mock__Fc_1uM_5hr_'),
        ('VTX0851359__5hr_', 'Mock__Fc_1uM_5hr_'),
        ('VTX0851359__R848_5hr_', 'Mock__R848_5hr_'),
        
        ('VTX0852555__5hr_', 'Mock__Fc_1uM_5hr_'),
        ('VTX0852555__R848_5hr_', 'Mock__R848_5hr_'),
        
        ('VTX0852488__5hr_', 'Mock__Fc_1uM_5hr_'),
        ('VTX0852488__R848_5hr_', 'Mock__R848_5hr_'),
        
        ('IL10__5hr_', 'Mock__Fc_1uM_5hr_'),
        ('IL10__R848_5hr_', 'Mock__R848_5hr_'),
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


datasets = ['PBMC - plate 1', 'PBMC 5hr - plate 2']#'PBMC', 'PBMC 5hr', 'PBMC 24hr']#, 'A549', 'HCT116']


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


st.header('Single Cell Browser')

with st.expander(label='Overview', expanded=True):

    col0, _ = st.columns(2)

    with col0:
        dataset = st.selectbox(
            'Choose a dataset',
            datasets, index=1,
        )
        adata, adata_dim = load_anndata(adata_paths[dataset])

        gene_de_dfs, gsva_de_dfs = load_differential_dfs(adata_paths[dataset])

        cell_types = ['Macrophages', 'T cells', 'B cells', 'Monocytes', 'ILC', 'DC']
        #cell_types = list(set(adata.obs['cell_type']))
        cell_types.sort()
        cell_types.insert(0, 'All')

        cell_type = st.selectbox(
            'Choose a cell type',
            cell_types, index=4
        )


        contrast_map = {f'{x[0]} - {x[1]}': f'{x[0]}-{x[1]}' for x in focus_contrasts[dataset]}
        dataset_abbr = adata_paths[dataset].stem

    if cell_type and cell_type != 'All':
        set_color = color_map[cell_type]
        color_map = {c: "#D3D3D3" for c in color_map.keys()}
        color_map[cell_type] = set_color

    cell_type = cell_type.replace(' ', '')

    plot_df = pd.DataFrame(adata_dim.obsm['X_umap'])
    cell_df = pd.DataFrame(adata_dim.obs['cell_type'])
    cell_df.reset_index(inplace=True)
    plot_df = pd.concat([plot_df, cell_df], axis=1)

    umap = plotting.plot_umap(plot_df, color_map, f'UMAP of {adata.shape[0]} cells')

    st.plotly_chart(umap, theme="streamlit")

    overview_de_df = plotting.overview_count_table(focus_contrasts[dataset],
        gene_de_dfs, dataset_abbr)
    st.markdown('##### Number of DE genes per contrast')
    st.dataframe(overview_de_df)


with st.expander(label='Differential Expression', expanded=True):

    col1, _ = st.columns(2)

    with col1:
        contrast = st.selectbox(
            'Choose a contrast',
            [f'{x[0]} - {x[1]}' for x in focus_contrasts[dataset]]
        )

    if cell_type != 'All' and contrast != 'None - None':

        gene_tab, pathway_gsva_tab, pathway_ora_tab = st.tabs(["Genes", "Pathways - GSVA", "Pathways - ORA"])
        
        gsva_key = f'{adata_paths[dataset].stem}__gsva__{contrast_map[contrast]}__{cell_type.replace("_", "")}'
        edger_key = f'{adata_paths[dataset].stem}__edgeR__{contrast_map[contrast]}__{cell_type.replace("_", "")}'

        if gsva_key in gsva_de_dfs.keys():

            plot_gene_df = gene_de_dfs[edger_key]
            plot_gene_df['point_size'] = 5
            plot_gene_df['Significant'] = plot_gene_df.apply(lambda x: abs(x.logFC) > 1 and x.FDR < .05, axis=1)

            plot_gsva_df = gsva_de_dfs[gsva_key]
            plot_gsva_df['point_size'] = 5

            with gene_tab:

                col2, col3 = st.columns(2)

                with col2:
                    st.subheader('Gene Volcano Plot')
                    select_gene_df = plot_gene_df.copy()

                    fig = plotting.plot_gene_volcano(select_gene_df)
                    
                    selected_points = plotly_events(fig)

                with col3:
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

                try:
                    contrast_tuple = contrast.split(' - ')

                    plot_df = plotting.extract_cluster_map(plot_gene_df.copy(), adata, cell_type, contrast_tuple)

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
                    col4, col5 = st.columns(2)

                    with col4:
                        st.write(gene, c1)
                        ax1 = plotting.plot_static_umap(adata_dim, samples1, gene)
                        st.pyplot(ax1)

                    with col5:
                        st.write(gene, c2)
                        ax2 = plotting.plot_static_umap(adata_dim, samples2, gene)
                        st.pyplot(ax2)

                string_button = st.button('Generate STRING network diagrams')

                col6, col7, col8 = st.columns(3)

                
                if string_button:
                    with col6:
                        st.write('STRING network of top up-regulated genes')
                        network_image = plotting.plot_string_network(list(sig_df[sig_df['logFC'] > 0]['gene'])[0:100])
                        st.image(network_image)

                    with col7:
                        st.write('STRING network of top down-regulated genes')
                        network_image = plotting.plot_string_network(list(sig_df[sig_df['logFC'] < 0]['gene'])[0:100])
                        st.image(network_image)

                    with col8:
                        st.write('STRING network of top regulated genes')
                        network_image = plotting.plot_string_network(list(sig_df['gene'])[0:100])
                        st.image(network_image)


            with pathway_tab:

                col9, col10 = st.columns(2)

                with col9:
                    st.subheader('Pathway Volcano Plot')
                    fig = plotting.plot_pathway_volcano(plot_gsva_df)
                    selected_pathway = plotly_events(fig)

                with col10:
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

                col11, col12 = st.columns(2)

                with col11:
                    st.subheader('Top pathways')
                    fig = plotting.plot_pathway_volcano(plot_gsva_df)
                    selected_pathway = plotly_events(fig)
            
            with pathway_ora_tab:
                st.subheader('Pathway Enrichment')
                pathway_db = st.selectbox(
                    'Choose Pathway Source',
                    ('reactome', 'hallmark', 'go_molecular_function', 'immunesigdb'), index=0)

                col_map = {'gene': 'gene_symbol', 
                           'logCPM': 'baseMean', 
                           'logFC': 'log2FoldChange',
                           '-Log10(FDR)': 'lfcSE',
                           'PValue': 'stat',
                           'FDR': 'padj'}
                
                results_df = plot_gene_df.copy()
                results_df.rename(columns=col_map, inplace=True)
                results_df.set_index('gene_symbol', inplace=True)

                msigdb = load_msigdb(pathway_db)
                
                top_genes = results_df[results_df['padj'] < 0.05]

                enr_pvals = dc.get_ora_df(
                    df=top_genes,
                    net=msigdb_subset,
                    source='geneset',
                    target='genesymbol'
                )
                
                ax = dc.plot_dotplot(top_path_df, x='Combined score', y = 'Term', s='Odds ratio', c = 'FDR p-value', scale = 0.5, figsize=(7,10))

                st.pyplot(ax)
                

        else:
            st.write('Not enough cells to perform differential analysis.')


with st.expander(label='Differential Comparison', expanded=True):

    col13, col14 = st.columns(2)

    with col13:
        contrast1 = st.selectbox(
            'Choose contrast #1',
            [f'{x[0]} - {x[1]}' for x in focus_contrasts[dataset]],
            index=3,
        )

    with col14:
        contrast2 = st.selectbox(
            'Choose contrast #2',
            [f'{x[0]} - {x[1]}' for x in focus_contrasts[dataset]],
            index=9
        )

    with st.container():
        if contrast1 != 'None - None' and contrast2 != 'None - None':
            gencode_map = load_gencode_map()
            
            tab3, tab4 = st.tabs(["Genes", "Pathways"])

            with tab3:
                scatter_gene, scatter_gene_df = plotting.plot_diff_scatter(gene_de_dfs, contrast_map, 
                    cell_type, dataset_abbr, contrast1, contrast2, diff_type='gene',
                    algorithm='edgeR')
                
                selected_gene = plotly_events(scatter_gene, override_height=1000)

                try:
                    entry = scatter_gene_df.iloc[int(selected_gene[0]["pointIndex"])]['gene']
                    st.write(entry)
                except:
                    st.write('')

                col_map = {col: "{:.2E}" for col in scatter_gene_df.columns if col[0:3] == 'FDR'}

                st.dataframe(scatter_gene_df.style.format(col_map))
                csv = convert_df(scatter_gene_df)
                st.download_button(
                        "Download Table",
                        csv,
                        f"genes-diff_{contrast1}.csv",
                        "text/csv",
                        key='download-csv-diff-gene'
                )
                
            with tab4:
                scatter_path, scatter_path_df = plotting.plot_diff_scatter(gsva_de_dfs, contrast_map, 
                    cell_type, dataset_abbr, contrast1, contrast2, diff_type='Pathway',
                    algorithm='gsva')

                selected_pathway = plotly_events(scatter_path, override_height=1000, override_width=1400)
                
                col_map = {col: "{:.2E}" for col in scatter_path_df.columns if col[0:3] == 'FDR'}

                st.dataframe(scatter_path_df.style.format(col_map))
                csv = convert_df(scatter_path_df)
                st.download_button(
                        "Download Table",
                        csv,
                        f"pathways-diff_{contrast1}.csv",
                        "text/csv",
                        key='download-csv-diff-pathway'
                )

                try:
                    entry = scatter_path_df.iloc[int(selected_pathway[0]["pointIndex"])]
                    st.write(entry)
                except:
                    st.write('')




