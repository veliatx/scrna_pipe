import streamlit as st
import pandas as pd
import numpy as np
from streamlit_echarts import st_echarts
from streamlit_plotly_events import plotly_events
import plotly.express as px

df = px.data.iris()

cell_type = st.selectbox(
    'Choose a cell type',
    ['setosa', 'virginica', 'versicolor'])

st.write('You selected:', cell_type)

color_map = {'setosa': "#ff0000", 
             'virginica': "#00ff00", 
             'versicolor': "#0000ff"}

if cell_type:
    set_color = color_map[cell_type]
    color_map = {c: "#D3D3D3" for c in color_map.keys()}
    color_map[cell_type] = set_color

fig = px.scatter(df, x="sepal_width", y="sepal_length", 
                color="species", color_discrete_map=color_map)

st.plotly_chart(fig, theme="streamlit")


