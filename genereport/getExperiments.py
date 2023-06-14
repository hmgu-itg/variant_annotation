#!/usr/bin/env python3

import time
import sys
import argparse

from selenium.webdriver.firefox.options import Options
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select

def main():
    parser=argparse.ArgumentParser(description="Download a list of GXA experiments")
    parser.add_argument('--type','-t',action="store",help="Experiment type (default: %(default)s)",choices=["baseline","differential"],default="baseline",required=False)
    parser.add_argument('--output','-o',action="store",help="Output file",required=True)

    if len(sys.argv[1:])==0:
        parser.print_help()
        sys.exit(0)

    try:
        args=parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    tp=args.type
    outfile=args.output

#---------------------------------------------------------------------------------------------------------------------------

    sys.stdout=open(sys.stdout.fileno(),mode='w',encoding='utf8',buffering=1)

    profile = webdriver.FirefoxProfile()
    profile.set_preference("browser.download.panel.shown", False)
    profile.set_preference("browser.helperApps.neverAsk.openFile","text/plain")
    profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "text/plain")
    profile.set_preference("browser.download.folderList", 2);
    profile.set_preference('browser.download.dir', '/tmp')
    profile.set_preference('browser.download.manager.showWhenStarting', False)

    options=Options()
    options.add_argument("--headless")

    browser=webdriver.Firefox(firefox_profile=profile,firefox_options=options)
    browser.get("https://www.ebi.ac.uk/gxa/experiments?species=homo%20sapiens&experimentType="+tp)
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
            html=browser.page_source
            with open(outfile,"w",encoding="utf-8") as f:
                print(html,file=f)
            browser.quit()
            sys.exit(0)
        else:
            browser.quit()
            sys.exit(1)
    except Exception as e:
        print("Exception occured: "+type(e).__name__+" : "+str(e),file=sys.stderr)
        browser.quit()
        sys.exit(1)

if __name__=="__main__":
    main()
