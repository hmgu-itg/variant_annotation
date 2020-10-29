import logging
import time
import sys
import shutil
import pandas as pd
import re
import os
import plotly.express as px
import tempfile as tf

from selenium.webdriver.firefox.options import Options
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

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

    fig = px.imshow(df)
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

def getGxaHeatmap(ID):
    return df2heatmap(getGxaDF(ID))

# ==============================================================================================================================

def getGxaDF(ID):
    '''
    Given a gene ID, get data from GXA (baseline expression)

    Input: gene ID
    Output: none, data is saved in temporary files (their names saved in "config.GXA_FILES")
    '''

    LOGGER.debug("Input: %s" % ID)

    fn=saveGxaData(ID)
    if fn is None:
        LOGGER.warning("Could not retreive GXA data for %s" % ID)
        return pd.DataFrame(columns=["Empty"])

    df=pd.read_table(fn,comment="#",header=0).rename(columns={"Unnamed: 0":"Experiment"})
    df=df.set_index("Experiment")
    # tissues with highest values
    tissues=set()
    idxs=list()
    for idx in df.index:
        idxs.append(idx)
        z=df.loc[idx,:].nlargest(n=int(config.GXA_HIGHEST))
        tissues.update(set(z.index.format()))
            
    df2=df.loc[idxs,list(tissues)]
    df2["Experiment"]=df2.index
    return df2

# ==============================================================================================================================

def getGxaDFLocal(ID):
    '''
    Given a gene ID, get data from GXA (baseline expression)

    Input: gene ID
    Output: none, data is saved in temporary files (their names saved in "config.GXA_FILES")
    '''

    LOGGER.debug("Input: %s" % ID)

    df=pd.read_table(config.GXA_FILE,header=0)
    df2=df.loc[df["Gene ID"]==ID].drop(["Gene ID","Gene Name"],axis=1)

    # # tissues with highest values
    # tissues=set()
    # idxs=list()
    # for idx in df.index:
    #     idxs.append(idx)
    #     z=df.loc[idx,:].nlargest(n=int(config.GXA_HIGHEST))
    #     tissues.update(set(z.index.format()))
            
    # df2=df.loc[idxs,list(tissues)]
    # df2["Experiment"]=df2.index

    return df2

# ======================================================================================================================

def saveGxaData(ID):
    """
    Retreives and saves GXA baseline expression values for a given gene ID

    Input  : gene ID
    Output : name of the file where the data are stored
    """

    try:
        os.remove("/tmp/expression_atlas-homo_sapiens.tsv")
    except OSError:
        pass

    profile = webdriver.FirefoxProfile()
    profile.set_preference("browser.download.panel.shown", False)
    profile.set_preference("browser.helperApps.neverAsk.openFile","text/plain")
    profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "text/plain")
    profile.set_preference("browser.download.folderList", 2);
    profile.set_preference('browser.download.dir', '/tmp')
    profile.set_preference('browser.download.manager.showWhenStarting', False)

    options = Options()
    options.add_argument("--headless")

    browser = webdriver.Firefox(firefox_profile=profile,firefox_options=options)
    browser.get(utils.getGxaURL(ID))
    element=None
    try:
        # 60 seconds to load
        element = WebDriverWait(browser, 60).until(EC.element_to_be_clickable((By.XPATH, '//button[text()="Download"]')))
    except Exception as e:
        LOGGER.error(type(e).__name__)
        browser.quit()
        return None

    element.click()
    time.sleep(5)
    browser.quit()

    try:
        os.remove(config.OUTPUT_DIR+"/expression_atlas-homo_sapiens.tsv")
    except OSError:
        pass

    shutil.move("/tmp/expression_atlas-homo_sapiens.tsv",config.OUTPUT_DIR)

    return config.OUTPUT_DIR+"/expression_atlas-homo_sapiens.tsv"
