import altair as alt
import decoupler as dc
import streamlit as st
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc

import pickle

from pathlib import Path
from streamlit_echarts import st_echarts
from streamlit_plotly_events import plotly_events

import plotly.graph_objects as go

from scrna_pipe import plotting, utils

st.set_page_config(layout="wide")


@st.cache_data
def load_anndata(adata_path):
    """
    """
    adata = sc.read(adata_path.parent.joinpath(f'{adata_path.stem}_processed.h5ad'))
    adata_dim = sc.read(adata_path.parent.joinpath(f'{adata_path.stem}_processed_umap.h5ad'))
    
    with open(adata_path.parent.joinpath(f'{adata_path.stem}_cpm.pkl'), 'rb') as f:
        cpm_df = pickle.load(f)
    
    return adata, adata_dim, cpm_df


@st.cache_data
def load_differential_dfs(adata_path):
    """
    """
    gene_de_dfs = {}
    gsva_de_dfs = {}
    ora_de_dfs = {}

    for df_path in adata_path.parent.joinpath('edgeR').glob('*.csv'):
        gene_de_dfs[df_path.stem] = pd.read_csv(df_path, index_col=0) 

    for df_path in adata_path.parent.joinpath('gsva').glob('*.csv'):
        gsva_de_dfs[df_path.stem] = pd.read_csv(df_path, index_col=0) 

    for df_path in adata_path.parent.joinpath('ora').glob('*.csv'):
        ora_de_dfs[df_path.stem] = pd.read_csv(df_path, index_col=0) 

    return gene_de_dfs, gsva_de_dfs, ora_de_dfs


@st.cache_data
def convert_df(df):
   return df.to_csv(index=False).encode('utf-8')


#output_dir = Path('/home/ubuntu/scrna_pipe/data')
output_dir = Path('/home/ec2-user/scrna_pipe/data')
output_dir_plate1 = Path('/home/ec2-user/velia-analyses-dev/VAP_20230711_single_cell_moa/outputs/run_6')
output_dir_plate2 = Path('/home/ec2-user/velia-analyses-dev/VAP_20230919_single_cell_pbmc_hits/outputs/run_2')
output_dir_plate3 = Path('/home/ec2-user/velia-analyses-dev/VAP_20231010_single_cell_plate_3/outputs/run_2')

adata_paths = {
    #'PBMC - plate 1': output_dir.joinpath('analysis_plate_1', 'PBMC_coarse.h5ad'),
    #'PBMC 24hr': output_dir.joinpath('analysis', '24hr_coarse.h5ad'),
    #'PBMC 5hr - plate 2': output_dir.joinpath('analysis_plate_2', '5hr_coarse_subset.h5ad'),
    'PBMC 5hr R848 - aggregate': output_dir.joinpath('analysis_plate_3', '5hr_pbmc_R848_scvi'),
    'PBMC 24hr R848 - aggregate': output_dir.joinpath('analysis_plate_3', '24hr_pbmc_R848_scvi'),
    'PBMC 5hr PAM3 - aggregate': output_dir.joinpath('analysis_plate_3', '5hr_pbmc_PAM3_scvi'),

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
    'PBMC 5hr R848 - aggregate': [
        ('None', 'None'),
        
        ('Mock__R848_5hr_', 'Mock__5hr_'),
        ('IL10__5hr_', 'Mock__5hr_'),
        ('IL10__R848_5hr_', 'Mock__R848_5hr_'),
        
        ('VTX0851359__5hr_', 'Mock__5hr_'),
        ('VTX0851359__R848_5hr_', 'Mock__R848_5hr_'),
        ('VTX0851359__R848_5hr_', 'Mock__5hr_'),
        
        ('VTX0852555__5hr_', 'Mock__5hr_'),
        ('VTX0852555__R848_5hr_', 'Mock__R848_5hr_'),
        ('VTX0852555__R848_5hr_', 'Mock__5hr_'),
    ],
    'PBMC 5hr PAM3 - aggregate': [
        ('None', 'None'),

        ('Mock__PAM3_5hr_', 'Mock__5hr_'),
        ('IL10__5hr_', 'Mock__5hr_'),
        ('IL10__PAM3_5hr_', 'Mock__PAM3_5hr_'),

        ('VTX0850629__5hr_', 'Mock__5hr_'),
        ('VTX0850629__PAM3_5hr_', 'Mock__PAM3_5hr_'),
        ('VTX0850629__PAM3_5hr_', 'Mock__5hr_'),

        ('VTX0850815__5hr_', 'Mock__5hr_'),
        ('VTX0850815__PAM3_5hr_', 'Mock__PAM3_5hr_'),
        ('VTX0850815__PAM3_5hr_', 'Mock__5hr_'),

        ('VTX0850805__5hr_', 'Mock__SLC6A8_buffer_5hr_'),
        ('VTX0850805__LPS_5hr_', 'Mock__LPS_SLC6A8_buffer_5hr_'),
        ('VTX0850805__5hr_', 'Mock__5hr_'),

    ],
    'PBMC 24hr R848 - aggregate': [
        ('None', 'None'),
        
        ('Mock__R848_24hr_', 'Mock__24hr_'),
        ('IL10__24hr_', 'Mock__24hr_'),
        ('IL10__R848_24hr_', 'Mock__R848_24hr_'),
        
        ('VTX0851359__24hr_', 'Mock__24hr_'),
        ('VTX0851359__R848_24hr_', 'Mock__R848_24hr_'),
        ('VTX0851359__R848_24hr_', 'Mock__24hr_'),
        
        ('VTX0852555__24hr_', 'Mock__24hr_'),
        ('VTX0852555__R848_24hr_', 'Mock__R848_24hr_'),
        ('VTX0852555__R848_24hr_', 'Mock__24hr_'),
        
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


datasets = ['PBMC - plate 1', 'PBMC 5hr - plate 2', 
            'PBMC 5hr R848 - aggregate', 'PBMC 24hr R848 - aggregate', 'PBMC 5hr PAM3 - aggregate',
            ]


color_map =  {
    "Macrophages": "#1f77b4",  
    "B cells": "#ff7f0e",
    "Naive B cells": "#ff7f0e",
    "T cells": "#2ca02c",
    "Naive T cells": "#1f77b4",
    "Regulatory T cells": "#1f77b4",
    "Cycling T cells": "#1f77b4",
    "Effector Memory T cells": "#2ca02c",
    "Monocytes": "#d62728",  
    "pDC": "#9467bd",
    "DC1": "#9467bd",
    "Migratory DCs": "#9467bd",
    "ILC": "#8c564b",
    "ILC3": "#8c564b",
    "DC": "#e377c2",
    "Endothelial cells": "#7f7f7f",
    "Plasma cells": "#bcbd22",
    "ETP": "#1f77b4",
    "NK cells": "#1f77b4",
}


st.header('Single Cell Browser')

with st.expander(label='Overview', expanded=True):
    col0, _ = st.columns(2)

    dataset = st.selectbox(
        'Choose a dataset',
        datasets, index=2,
    )
    adata, adata_dim, cpm_dfs = load_anndata(adata_paths[dataset])

    gene_de_dfs, gsva_de_dfs, ora_de_dfs = load_differential_dfs(adata_paths[dataset])

    cell_types = list(set(adata.obs['cell_type_display']))
    cell_types.sort()
    cell_types.insert(0, 'All')

    cell_type = st.selectbox(
        'Choose a cell type',
        cell_types, index=cell_types.index('T cells'),
        key='cell_type_select'
    )

#umap_tab, diff_tab, diff_compare_tab = st.tabs(["Cell types", "Differential", "Differential comparison"])

with st.expander(label='Cell Type Annotation', expanded=True):

    contrast_map = {f'{x[0]} - {x[1]}': f'{x[0]}-{x[1]}' for x in focus_contrasts[dataset]}
    dataset_abbr = adata_paths[dataset].stem

    if cell_type and cell_type != 'All':
        set_color = color_map[cell_type]
        color_map = {c: "#D3D3D3" for c in color_map.keys()}
        color_map[cell_type] = set_color

    cell_type = cell_type.replace(' ', '')

    plot_df = pd.DataFrame(adata_dim.obsm['X_umap'])
    cell_df = pd.DataFrame(adata_dim.obs['cell_type_display'])
    cell_df.reset_index(inplace=True)
    plot_df = pd.concat([plot_df, cell_df], axis=1)

    umap = plotting.plot_umap(plot_df, color_map, f'UMAP of {adata.shape[0]} cells')

    st.plotly_chart(umap, theme="streamlit", use_container_width=True)

    st.markdown('##### Cell type percentages')
    cell_type_df = pd.DataFrame(adata_dim.obs.groupby('cell_type').size()/adata_dim.obs.shape[0], columns=['Cell Type %']).T
    cell_type_df.sort_values(by='Cell Type %', axis=1, ascending=False, inplace=True)
    st.dataframe(cell_type_df.style.format("{:.2%}"))

    st.markdown('##### Cell type counts')
    cell_type_detail_df = utils.overview_cell_type_table(adata_dim)
    st.dataframe(cell_type_detail_df)

    st.markdown('##### Number of DE genes per contrast')
    overview_de_df = utils.overview_de_count_table(focus_contrasts[dataset],
        gene_de_dfs, dataset_abbr)

    overview_de_df.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
    overview_de_df = overview_de_df.astype(int)
    overview_de_df.columns = [f'{x[0]}: {x[1]}' for x in overview_de_df.columns]

    st.dataframe(overview_de_df.style.background_gradient(cmap='Blues'))


with st.expander(label='Differential Expression', expanded=True):

    col1, _ = st.columns(2)

    with col1:
        contrast = st.selectbox(
            'Choose a contrast',
            [f'{x[0]} - {x[1]}' for x in focus_contrasts[dataset]],
            index=5
        )

    if cell_type != 'All' and contrast != 'None - None':

        gene_tab, pathway_gsva_tab, pathway_ora_tab = st.tabs(["Genes", "Pathways - GSVA", "Pathways - ORA"])
        
        gsva_key = f'{adata_paths[dataset].stem}__gsva__{contrast_map[contrast]}__{cell_type.replace("_", "")}'
        edger_key = f'{adata_paths[dataset].stem}__edgeR__{contrast_map[contrast]}__{cell_type.replace("_", "")}'

        if gsva_key in gsva_de_dfs.keys():

            plot_gene_df = gene_de_dfs[edger_key]
            plot_gene_df['point_size'] = 10
            plot_gene_df['Significant'] = plot_gene_df.apply(lambda x: abs(x.logFC) > 1 and x.FDR < .05, axis=1)
            plot_gene_df = plot_gene_df[plot_gene_df['logCPM'] > 5]
            
            plot_gsva_df = gsva_de_dfs[gsva_key]
            plot_gsva_df = plot_gsva_df[plot_gsva_df['DE Genes In Pathway'] > 3]
            plot_gsva_df['point_size'] = 10

            with gene_tab:

                with st.container():
                    st.subheader('Gene Volcano Plot')
                    select_gene_df = plot_gene_df.copy()

                    fig = plotting.plot_gene_volcano(select_gene_df)
                    
                    selected_points = plotly_events(fig)
                
                col2, col3 = st.columns(2)
                
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
                
                with col3:
                    st.subheader('Gene specific data')

                    cpm_df = cpm_dfs[cell_type]
                    genes = list(cpm_df.index)
                    genes.append('')

                    select_gene_df = select_gene_df[select_gene_df['Significant']]

                    if selected_points:
                        gene = select_gene_df.iloc[selected_points[0]["pointIndex"]]['gene']
                    else:
                        gene = ''

                    selected_gene = st.selectbox(
                        'Select a gene',
                        genes, index=genes.index(gene),
                    )

                    c1, c2 = contrast.split(' - ')

                    if selected_gene != '':
                        st.markdown('#### Boxplot')

                        ax1 = plotting.plot_gene_boxplot(cpm_df, selected_gene, cell_type, (c1, c2))
                        st.pyplot(ax1)

                        st.markdown('#### UMAP projection')


                        samples1 = [f'{c1}_{i}' for i in range(1, 4)]
                        samples2 = [f'{c2}_{i}' for i in range(1, 4)]

                        col4, col5 = st.columns(2)

                        with col4:
                            st.write(selected_gene, c1)
                            ax2 = plotting.plot_static_umap(adata_dim, samples1, selected_gene)
                            st.pyplot(ax2)

                        with col5:
                            st.write(selected_gene, c2)
                            ax3 = plotting.plot_static_umap(adata_dim, samples2, selected_gene)
                            st.pyplot(ax3)


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


            with pathway_gsva_tab:

                col9, col10 = st.columns(2)

                st.subheader('Pathway Volcano Plot')
                fig = plotting.plot_pathway_volcano(plot_gsva_df)
                st.plotly_chart(fig, use_container_width=True)
                #selected_pathway = plotly_events(fig, select_event=True)#, override_width="100%")
                #st.write(selected_pathway)

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
            
            with pathway_ora_tab:
                
                colx, coly = st.columns(2)
                col11, col12 = st.columns(2)

                st.subheader('Pathway Enrichment')
                
                with colx:

                    pathway_db = st.selectbox(
                        'Choose Pathway Source',
                        ('reactome_pathways', 'hallmark', 'pid_pathways', 'kegg_pathways',
                         'go_biological_process', 'go_molecular_function', 'immunesigdb',
                         'chemical_and_genetic_perturbations', 
                         'vaccine_response', 'wikipathways'), index=0)
                                        
                    ora_name = edger_key.replace('_edgeR_', f'_ora-{pathway_db}_')
                    ora_flag = True
                    try:
                        top_path_df = ora_de_dfs[ora_name]
                    except:
                        ora_flag = False

                if ora_flag:
                    with col11:
                        
                        scale_param = 3/np.max(top_path_df['Odds ratio'])
                        
                        ax = dc.plot_dotplot(top_path_df[0:20], x='Combined score', y = 'Term', s='Odds ratio',
                                            c='FDR p-value', scale=scale_param, return_fig=True)
                        ax.set_size_inches((5, 8))

                        st.pyplot(ax, use_container_width=True, clear_figure=True)

                    with col12:
                        st.dataframe(top_path_df.style.format({"FDR": "{:.2E}"}))
                        csv = convert_df(top_path_df)
                        st.download_button(
                            "Download Table",
                            csv,
                            f"pathways-ora_{gsva_key}.csv",
                            "text/csv",
                            key='download-csv-pathway-ora'
                        )

        else:
            st.write('Not enough cells to perform differential analysis.')

with st.expander(label='Differential Comparison', expanded=True):

    col13, col14 = st.columns(2)

    with col13:
        contrast1 = st.selectbox(
            'Choose contrast #1',
            [f'{x[0]} - {x[1]}' for x in focus_contrasts[dataset]],
            index=5,
        )

    with col14:
        contrast2 = st.selectbox(
            'Choose contrast #2',
            [f'{x[0]} - {x[1]}' for x in focus_contrasts[dataset]],
            index=3,
        )

    with st.container():
        if contrast1 != 'None - None' and contrast2 != 'None - None':
            
            tab3, tab4 = st.tabs(["Genes", "Pathways"])

            with tab3:
                scatter_gene, scatter_gene_df = plotting.plot_diff_scatter(gene_de_dfs, contrast_map, 
                    cell_type, dataset_abbr, contrast1, contrast2, diff_type='gene',
                    algorithm='edgeR')
                
                scatter_gene.update_layout(legend_font=dict(size=18), height=800)

                st.plotly_chart(scatter_gene, use_container_width=True)
                #selected_gene = plotly_events(scatter_gene, override_height=1000)

                col15, col16 = st.columns(2)

                with col15:
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
                
                with col16:
                    cpm_df = cpm_dfs[cell_type]
                    genes = list(cpm_df.index)

                    selected_contrast_gene = st.selectbox(
                        'Select a gene',
                        genes,  key='contrast_gene'
                    )

                    if selected_contrast_gene:
                        fig_contrast_gene = plotting.plot_gene_boxplot(cpm_df, selected_contrast_gene, cell_type)
                        st.pyplot(fig_contrast_gene, clear_figure=True)
                
            with tab4:

                scatter_path, scatter_path_df = plotting.plot_diff_scatter(gsva_de_dfs, contrast_map, 
                    cell_type, dataset_abbr, contrast1, contrast2, diff_type='Pathway',
                    algorithm='gsva')

                scatter_path.update_layout(legend_font=dict(size=18), height=800)

                st.plotly_chart(scatter_path, use_container_width=True)
                #selected_pathway = plotly_events(scatter_path, override_height=1000)
                
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





