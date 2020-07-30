#!/usr/bin/env python3

import argparse
import gzip
import json
import pandas as pd
import re

import plotly.graph_objects as go
import dash
import dash_table
import dash_core_components as dcc
import plotly.figure_factory as ff
import dash_html_components as html
import plotly.express as px
from dash.dependencies import Input, Output
import pandas as pd

def table2fig(t):
    return go.Figure(data=[go.Table(header=dict(values=[["Field"],["Value"]],fill_color='paleturquoise',align='left'),
    cells=dict(values=[t.index,t.Value],fill_color='lavender',align='left'))])

# -----------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument("--input", "-i", help="Required: input data", required=True)
args=parser.parse_args()
fname=args.input

# -----------------------------------------------------------------------------------------------------------------------

data=None
with gzip.GzipFile(fname,"r") as fin:
    data=json.loads(fin.read().decode('utf-8'))

# -----------------------------------------------------------------------------------------------------------------------

tables=list()
mappings=dict()
for x in data:
    m=re.search("variant_table(\d+)",x)
    if m:
        t=pd.read_json(data[x])
        tables.append(t)
        mappings[m.group(1)]=t.loc["Location","Value"]

gnomADfigs=dict()
for x in data:
    m=re.search("gnomad_table(\d+)",x)
    if m:
        gnomADfigs[mappings[m.group(1)]]=ff.create_table(pd.read_json(data[x]))

regfigs=dict()
for x in data:
    m=re.search("regulation_table(\d+)",x)
    if m:
        t=pd.read_json(data[x])
        if len(t)>0:
            regfigs[mappings[m.group(1)]]=dash_table.DataTable(style_cell={'textAlign': 'left','whiteSpace': 'normal','height': 'auto'},style_header={'backgroundColor': 'white','fontWeight': 'bold'},columns=[{"name": i, "id": i} for i in t.columns],data=t.to_dict("records"))
#ff.create_table(t)
        else:
            regfigs[mappings[m.group(1)]]=None

gwasfigs=dict()
for x in data:
    m=re.search("gwas_table(\d+)",x)
    if m:
        t=pd.read_json(data[x])
        if len(t)>0:
            gwasfigs[mappings[m.group(1)]]=dash_table.DataTable(style_cell={'textAlign': 'left','whiteSpace': 'normal','height': 'auto'},style_cell_conditional=[{'if': {'column_id': 'Trait'},'width': '25%'}],style_header={'backgroundColor': 'white','fontWeight': 'bold'},columns=[{"name": i, "id": i} for i in t.columns],data=t.to_dict("records"))
#            gwasfigs[mappings[m.group(1)]]=ff.create_table(t)
        else:
            gwasfigs[mappings[m.group(1)]]=None

# -----------------------------------------------------------------------------------------------------------------------

figs=dict()
for t in tables:
    figs[t.loc["Location","Value"]]=table2fig(t)

# -----------------------------------------------------------------------------------------------------------------------

app=dash.Dash(__name__)

elements=[html.H1(children='Variant info',style={"textAlign":"center"}),dcc.Dropdown(placeholder="Select a mapping",clearable=False,id="select",options=[{"label":x.loc["Location","Value"],"value":x.loc["Location","Value"]} for x in tables],value=tables[0].loc["Location","Value"])]
elements.append(html.Button('Details',id='details-btn',className="collapsible",n_clicks=0))
elements.append(dcc.Graph(id="info_table",figure=figs[tables[0].loc["Location","Value"]],className="content"))
elements.append(html.Button('GnomAD Allele Frequencies',id='gnomad-btn',className="collapsible",n_clicks=0))
elements.append(dcc.Graph(id="gnomad_table",figure=gnomADfigs[tables[0].loc["Location","Value"]],className="content"))

elements.append(html.Button('ENSEMBL Regulation',id='reg-btn',className="collapsible",n_clicks=0,disabled=regfigs[tables[0].loc["Location","Value"]] is None))
#elements.append(dcc.Graph(id="reg_table",figure=regfigs[tables[0].loc["Location","Value"]] if regfigs[tables[0].loc["Location","Value"]] is not None else go.Figure(),className="content"))
elements.append(html.Div(children=regfigs[tables[0].loc["Location","Value"]],className="content",id="reg_table") if regfigs[tables[0].loc["Location","Value"]] is not None else html.Div(children=[],className="content",id="reg_table"))

elements.append(html.Button('GWAS signals nearby',id='gwas-btn',className="collapsible",n_clicks=0,disabled=gwasfigs[tables[0].loc["Location","Value"]] is None))
#elements.append(dcc.Graph(id="gwas_table",figure=gwasfigs[tables[0].loc["Location","Value"]] if gwasfigs[tables[0].loc["Location","Value"]] is not None else go.Figure(),className="content"))
elements.append(html.Div(children=gwasfigs[tables[0].loc["Location","Value"]],className="content",id="gwas_table") if gwasfigs[tables[0].loc["Location","Value"]] is not None else html.Div(children=[],className="content",id="gwas_table"))

app.layout = html.Div(children=elements)

# -----------------------------------------------------------------------------------------------------------------------

@app.callback(Output("info_table","figure"),[Input("select","value")])
def update_figure(location_str):
    return figs[location_str]

@app.callback(Output("gnomad_table","figure"),[Input("select","value")])
def update_gnomad(location_str):
    return gnomADfigs[location_str]

@app.callback(Output("info_table","style"),[Input("details-btn","n_clicks")])
def collapseInfo(val):
    if val%2:
        return {"display":"block"}
    else:
        return {"display":"none"}

@app.callback(Output("gnomad_table","style"),[Input("gnomad-btn","n_clicks")])
def collapseGnomAD(val):
    if val%2:
        return {"display":"block"}
    else:
        return {"display":"none"}

@app.callback([Output("reg_table","children"),Output("reg-btn","disabled")],[Input("select","value")])
def update_reg(location_str):
    if regfigs[location_str] is None:
        return [],True
    else:
        return [regfigs[location_str]],False

@app.callback(Output("reg_table","style"),[Input("reg-btn","n_clicks")])
def collapseReg(val):
    if val%2:
        return {"display":"block"}
    else:
        return {"display":"none"}

@app.callback([Output("gwas_table","children"),Output("gwas-btn","disabled")],[Input("select","value")])
def update_gwas(location_str):
    if gwasfigs[location_str] is None:
        return [],True
    else:
        return [gwasfigs[location_str]],False

@app.callback(Output("gwas_table","style"),[Input("gwas-btn","n_clicks")])
def collapseGwas(val):
    if val%2:
        return {"display":"block"}
    else:
        return {"display":"none"}

# -----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    app.run_server(debug=True)
