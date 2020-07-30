#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
import logging
import plotly.express as px
import tempfile as tf

from modules import gxa
from modules import utils
from modules import config

LOGGER=logging.getLogger("test_gxa")
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("gxa").setLevel(logging.DEBUG)
logging.getLogger("utils").setLevel(logging.DEBUG)

#----------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Test GXA module")
#parser.add_argument('--output','-o', action="store",help="Output directory")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--id','-i', action="store",help="gene ID",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

ID=args.id
outdir="./"+ID

#---------------------------------------------------------------------------------------------------------------------------

if not utils.createDir(outdir):
    LOGGER.error("Could not create output dir %s" % outdir)

config.OUTPUT_DIR=outdir

df=gxa.getGxaDF(ID)
print(df)
# out=tf.NamedTemporaryFile(delete=False,mode="w")
# fig = px.imshow(df)
# print(fig.write_html(out.name,full_html=False))
# out.close()
# with open (out.name, "r") as f:
#     data=f.readlines()
# f.close()
# print("".join(data))
# if os.path.isfile(out.name):
#     os.remove(out.name)

