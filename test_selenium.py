#!/usr/bin/env python3

import time
import sys
from selenium.webdriver.firefox.options import Options
from html5print import HTMLBeautifier
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

profile = webdriver.FirefoxProfile()
profile.set_preference("browser.download.panel.shown", False)
profile.set_preference("browser.helperApps.neverAsk.openFile","text/plain")
profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "text/plain")
profile.set_preference("browser.download.folderList", 2);
profile.set_preference('browser.download.dir', '/home/andrei')
profile.set_preference('browser.download.manager.showWhenStarting', False)

options = Options()
options.add_argument("--headless")

browser = webdriver.Firefox(firefox_profile=profile,firefox_options=options)
#browser.get("https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5214/Results?ref=aebrowse")
browser.get("https://www.ebi.ac.uk/gxa/genes/ENSG00000229183?bs=%7B%22homo%20sapiens%22%3A%5B%22ORGANISM_PART%22%5D%7D&ds=%7B%22kingdom%22%3A%5B%22animals%22%5D%7D#baseline")
element=None
try:
    element = WebDriverWait(browser, 60).until(EC.element_to_be_clickable((By.XPATH, '//button[text()="Download"]')))
except:
    print("Exception occured",file=sys.stderr)
    #print(browser.page_source.encode("utf-8"))

#print(HTMLBeautifier.beautify(browser.page_source.encode("utf-8"),4))
element.click()
#print("Element clicked",file=sys.stderr)
time.sleep(5)
browser.quit()
