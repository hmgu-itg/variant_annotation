#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
import json, requests
from bokeh.io import output_file, show

from bokeh.palettes import Spectral10
from bokeh.plotting import *
from bokeh.models import *

from varannot import query

infile=sys.argv[1]
chrcol=sys.argv[2]
pvalcol=sys.argv[3]
pscol=sys.argv[4]
rscol=sys.argv[5]
mafcol=sys.argv[6]
outfile=sys.argv[7]

pd.options.mode.chained_assignment = None  
df=pd.read_csv(infile, sep="\t",index_col=False)
output_file(outfile)
hover= HoverTool(tooltips = [
        ("==============", "=============="),
        ("name", "   @"+rscol),
        ("ld","@ld"),
        ("MAF", "@"+mafcol),
        ("overlaps gene ", "@gene"),
        ("consequence ", "@consequence"),
        ("known associations", "@traits")
])
p=figure(tools=[hover, 'box_zoom,wheel_zoom,reset,tap'], width=1500)

e=df
e['col']=pd.cut(e['ld'], 9, labels=range(1,10))
collol=pd.DataFrame({'pal':Spectral10})
e['col_traits']=pd.notnull(e['traits']).astype(int)
e['col']=np.asarray(collol.loc[e['col']]).flatten()
e[pvalcol]=e[pvalcol].astype(float)
e["traits"]=e["traits"].astype(str)
e['logp']=-np.log10(e[pvalcol])

print(e)

for _,row in e.iterrows():
    pos=row[pscol]
    lp=row["logp"]
    traits=row["traits"]
    if traits=="nan":
        continue
    p.line(x=[pos,pos],y=[lp,e['logp'].max()],line_color="firebrick",line_width=1,line_dash="dashed")
    label=Label(x=pos,y=e['logp'].max(),text=traits,angle=90,angle_units="deg",text_align="right",text_color="firebrick",text_font_size="11pt",text_font_style="italic")
    p.add_layout(label)

overlapping_genes=query.restQuery(query.makeGeneOverlapQueryURL(str(e['#chr'][0]),e['ps'].min(),e['ps'].max(),build="38"))
print(json.dumps(overlapping_genes,indent=4,sort_keys=True))
genes_df=pd.DataFrame(json.loads(json.dumps(overlapping_genes)))
#genes_df.to_csv(sys.stdout,sep="\t",index=False)

# overlapping_GWASCAT_vars=query.restQuery(query.makeOverlapVarGWASCATQueryURL(str(e['#chr'][0]),e['ps'].min(),e['ps'].max(),build="38"))
# cat=pd.DataFrame(json.loads(json.dumps(overlapping_GWASCAT_vars)))
# print(cat)
#print("")
#cat.to_csv(sys.stdout,sep="\t",index=False)

# TODO: max POST size
# pheno_vars=query.restQuery(query.makeRSPhenotypeQueryURL(build="38"),data=utils.list2string(cat["id"].tolist()),qtype="post")
# for rsid in pheno_vars:
#     p.ray(x=pheno_vars[rsid]['mappings'][0]['end'],color="firebrick",length=0,angle=90,angle_units="deg")
#     label=Label(x=pheno_vars[rsid]['mappings'][0]['end'],y=e['logp'].max(),text=pheno_vars[rsid]['phenotypes'][0]['trait'],angle=90,angle_units="deg",text_align="right",text_color="firebrick",text_font_size="11pt",text_font_style="italic")
#     p.add_layout(label)

e=ColumnDataSource(e)

url="http://ensembl.org/Homo_sapiens/Variation/Explore?v=@"+rscol
taptool=p.select(type=TapTool)
taptool.callback=OpenURL(url=url)
p.circle(pscol,'logp',line_width=2,source=e,size=9,fill_color='col',line_color="black", line_alpha='col_traits')

p2=figure(width=1500,height=300,x_range=p.x_range,y_range=(0,1),tools=['tap'])
if len(genes_df)>0:
    ys=np.linspace(start=1/(len(genes_df)+1),stop=len(genes_df)/(len(genes_df)+1),num=len(genes_df))
    #ys=np.random.rand(len(genes_df))
    print(ys)
    genes_df['y']=ys
    genes_df['color']="cornflowerblue"
    genes_df['color'].loc[genes_df['biotype']=="protein_coding"]="goldenrod"
    genes_df['sens']="<"
    genes_df['sens'].loc[genes_df['strand']>0]=">"
    genes_df['name']=genes_df['sens']+genes_df['external_name']
    p2.segment(x0=genes_df['start'],x1=genes_df['end'],y0=ys,y1=ys,line_width=4,color=genes_df['color'])
    labels=LabelSet(x='start',y='y',text='name',source=ColumnDataSource(genes_df))
    p2.add_layout(labels)

q=gridplot([[p], [p2]])
save(q)

