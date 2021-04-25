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
from varannot import utils

infile=sys.argv[1]

pd.options.mode.chained_assignment = None  
d=pd.read_csv(infile, sep="\t",index_col=False)
output_file(infile+".html")
hover= HoverTool(tooltips = [
        ("==============", "=============="),
        ("name", "   @"+sys.argv[4]),
#        ("RS-id", "@id"),
        ("ld","@ld"),
        ("MAF", "@"+sys.argv[5]),
        ("overlaps gene ", "@gene"),
        ("consequence ", "@consequence"),
        ("known associations", "@traits")
])
p=figure(tools=[hover, 'box_zoom,wheel_zoom,reset,tap'], width=1500)

e=d
e['col']=pd.cut(e['ld'], 9, labels=range(1,10))
collol=pd.DataFrame({'pal':Spectral10})

#e['traits']=e['traits'].astype(str)
e['col_traits']=pd.notnull(e['traits']).astype(int)
#e['col_traits']=e['traits'].astype(str)

#e['col_traits'].loc[pd.isnull(e['traits'])]=0
#e['col_traits'].loc[not pd.isnull(e['traits'])]=1
#e['col_traits'].loc[e['col_traits']=="none"]=0
#e['col_traits'].loc[e['col_traits']!=0]=1

e['col']=np.asarray(collol.loc[e['col']]).flatten()
e[sys.argv[2]]=e[sys.argv[2]].astype(float)
e['logp']=-np.log10(e[sys.argv[2]])

print(e)
#e.to_csv(sys.stdout,sep="\t",index=False)
#print(collol)

overlapping_genes=query.restQuery(query.makeGeneOverlapQueryURL(str(e['#chr'][0]),e['ps'].min(),e['ps'].max(),build="38"))
print(json.dumps(overlapping_genes,indent=4,sort_keys=True))
d=pd.DataFrame(json.loads(json.dumps(overlapping_genes)))
#d.to_csv(sys.stdout,sep="\t",index=False)

overlapping_GWASCAT_vars=query.restQuery(query.makeOverlapVarGWASCATQueryURL(str(e['#chr'][0]),e['ps'].min(),e['ps'].max(),build="38"))
cat=pd.DataFrame(json.loads(json.dumps(overlapping_GWASCAT_vars)))
print(cat)
#print("")
#cat.to_csv(sys.stdout,sep="\t",index=False)

# TODO: max POST size
pheno_vars=query.restQuery(query.makeRSPhenotypeQueryURL(build="38"),data=utils.list2string(cat["id"].tolist()),qtype="post")
for rsid in pheno_vars:
    p.ray(x=pheno_vars[rsid]['mappings'][0]['end'],color="firebrick",length=0,angle=90,angle_units="deg")
    label=Label(x=pheno_vars[rsid]['mappings'][0]['end'],y=e['logp'].max(),text=pheno_vars[rsid]['phenotypes'][0]['trait'],angle=90,angle_units="deg",text_align="right",text_color="firebrick",text_font_size="11pt",text_font_style="italic")
    p.add_layout(label)

e=ColumnDataSource(e)

url="http://ensembl.org/Homo_sapiens/Variation/Explore?v=@"+sys.argv[4]
taptool=p.select(type=TapTool)
taptool.callback=OpenURL(url=url)
p.circle(sys.argv[3],'logp',line_width=2,source=e,size=9,fill_color='col',line_color="black", line_alpha='col_traits')

p2=figure(width=1500, height=300, x_range=p.x_range, tools=['tap'])
ys=np.random.rand(len(d['end']))
d['y']=ys

d['color']="cornflowerblue"
d['color'].loc[d['biotype']=="protein_coding"]="goldenrod"

d['sens']="<"
d['sens'].loc[d['strand']>0]=">"
d['name']=d['sens']+d['external_name']
p2.segment(x0=d['start'],x1=d['end'],y0=ys,y1=ys,line_width=4,color=d['color'])
labels=LabelSet(x='start',y='y',text='name',source=ColumnDataSource(d))
p2.add_layout(labels)
q=gridplot([[p], [p2]])
save(q)

