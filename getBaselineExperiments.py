#!/usr/bin/env python3

import time
import sys

from selenium.webdriver.firefox.options import Options
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select

sys.stdout=open(sys.stdout.fileno(),mode='w',encoding='utf8',buffering=1)

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
browser.get("https://www.ebi.ac.uk/gxa/experiments?species=homo%20sapiens&experimentType=baseline")
element=None
try:
    element=WebDriverWait(browser,15).until(EC.presence_of_element_located((By.XPATH, "//label[text()='Entries per page:']")))
    if not element is None:
        #print("Found element")
        children=element.find_elements_by_xpath(".//select")
        #print("Found %d children" % len(children))
        sel=Select(children[0])
        #print("Created Select")
        sel.select_by_visible_text('All')
        #print("Selected")
        time.sleep(15)
        print(browser.page_source)
    browser.quit()
except Exception as e:
    print("Exception occured: "+type(e).__name__,file=sys.stderr)
    browser.quit()


