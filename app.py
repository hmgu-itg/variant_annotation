#!/usr/bin/env python3

import argparse
import gzip
import json
import pandas as pd
import re

import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
from dash.dependencies import Input, Output
import pandas as pd

def generateTable(dataframe):
    return html.Table([
        html.Thead(
            html.Tr([html.Th(col) for col in dataframe.columns])
        ),
        html.Tbody([
            html.Tr([
                html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
            ]) for i in range(len(dataframe))
        ])
    ])

def table2fig(t):
    return go.Figure(data=[go.Table(header=dict(values=[["Field"],["Value"]],fill_color='paleturquoise',align='left'),
    cells=dict(values=[t.index,t.Value],fill_color='lavender',align='left'))])

# Parsing command line arguments:
parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument("--input", "-i", help="Required: input data", required=True)
args=parser.parse_args()
fname=args.input

data=None
with gzip.GzipFile(fname,"r") as fin:
    data=json.loads(fin.read().decode('utf-8'))

tables=[pd.read_json(data[x]) for x in data if re.search("variant_table",x)]
for t in tables:
    t["Field"]=t.index
print(tables[0])

figs=dict()
for t in tables:
    figs[t.loc["Location","Value"]]=table2fig(t)

app=dash.Dash(__name__)

elements=[html.H1(children='Variant info'),html.Label("Mapping"),dcc.Dropdown(clearable=False,id="select",options=[{"label":x.loc["Location","Value"],"value":x.loc["Location","Value"]} for x in tables],value=tables[0].loc["Location","Value"])]
elements.append(dcc.Graph(id="info_table",figure=figs[tables[0].loc["Location","Value"]]))

app.layout = html.Div(children=elements)

@app.callback(Output("info_table","figure"),[Input("select","value")])
def update_figure(location_str):
    return figs[location_str]

if __name__ == '__main__':
    app.run_server(debug=True)
