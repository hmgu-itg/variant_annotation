#!/usr/bin/env python3

import json
import logging

from varannot import utils

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)

logging.getLogger("varannot.utils").addHandler(ch)
logging.getLogger("varannot.utils").setLevel(logging.DEBUG)

#------------------------------------------------------------------------------------------------------------------

input_data=[{"chr":"chr1","start":"10034785","end":"10039857","id":"ID"}]
output=utils.runLiftOver(input_data)

# output should be: chr1:10094843-10099915

print("Expected output: chr1:10094843-10099915\n")

for x in output:
    print("%s\t%s\t%s\t%s\n" % (x["chr"],x["start"],x["end"],x["id"]))
