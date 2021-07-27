import logging
import time
import sys
import pandas as pd
import re
import os
import plotly.express as px
import tempfile as tf

from varannot import config
from varannot import utils

LOGGER=logging.getLogger(__name__)

# ==============================================================================================================================

def df2svg(df):
    if df is None:
        return None
    if len(df)==0:
        return None
    df=df.set_index("Experiment")
    fig=px.imshow(df)
    out=tf.NamedTemporaryFile(dir=config.OUTPUT_DIR_FIG,suffix=".svg")
    fname=out.name
    out.close()
    fig.write_image(fname)
    return fname

# ==============================================================================================================================

def df2heatmap(df):
    if len(df)==0:
        return None
    df=df.set_index("Experiment")
    fig=px.imshow(df)
    out=tf.NamedTemporaryFile(delete=False,mode="w")
    fig.write_html(out.name,full_html=False)
    out.close()
    with open (out.name, "r") as f:
        data=f.readlines()
    f.close()
    if os.path.isfile(out.name):
        os.remove(out.name)
    return "".join(data)

# ==============================================================================================================================

def getGxaDFLocal(ID):
    '''
    Given a gene ID, get data from GXA (baseline expression)

    Input: gene ID
    Output: dataframe
    '''

    LOGGER.debug("Input: %s" % ID)
    df=config.GXA_DF
    df2=df.loc[df["Gene ID"]==ID].drop(["Gene ID","Gene Name"],axis=1)
    isn=df2.drop("Experiment",axis=1).isnull()
    df3=df2[~isn.all(axis=1)]
    LOGGER.debug("Output: %d records" % len(df3))
    return df3

