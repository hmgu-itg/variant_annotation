import config
import logging
import time
import sys
import shutil
import pandas as pd
import re

from selenium.webdriver.firefox.options import Options
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

from . import config
from . import utils

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ==============================================================================================================================

def getGxa(ID):
    '''
    Geiven a gene ID, get data from GXA (baseline expression)

    Input: gene ID
    Output: none, data is saved in temporary files (their names saved in "config.GXA_FILES")
    '''

    LOGGER.debug("Input: %s" % ID)

    fn=saveGxaData(ID)
    if fn is None:
        LOGGER.warning("Could not retreive GXA data for %s" % ID)
        return None

    df=pd.read_table(fn,comment="#",header=0)
    df=df.set_index("Unnamed: 0")
    # tissues with highest values
    tissues=set()
    idxs=list()
    for idx in df.index:
        # only selected experiments
        if re.search("GTEx",idx.format(),re.I) or re.search("FANTOM",idx.format(),re.I):
            idxs.append(idx)
            z=df.loc[idx,:].nlargest(n=int(config.GXA_HIGHEST))
            tissues.update(set(z.index.format()))
            
    df.loc[idxs,list(tissues)].to_csv(config.OUTPUT_DIR+"/"+ID+".tsv",sep="\t")

# ======================================================================================================================

def saveGxaData(ID):
    """
    Retreives and saves GXA baseline expression values for a given gene ID

    Input  : gene ID
    Output : name of the file where the data are stored
    """

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
        # give it up to 60 seconds to load
        element = WebDriverWait(browser, 60).until(EC.element_to_be_clickable((By.XPATH, '//button[text()="Download"]')))
    except:
        LOGGER.error("Exception occured")
        return None

    element.click()
    time.sleep(5)
    browser.quit()
    shutil.move("/tmp/expression_atlas-homo_sapiens.tsv",config.OUTPUT_DIR)
    return config.OUTPUT_DIR+"/expression_atlas-homo_sapiens.tsv"
