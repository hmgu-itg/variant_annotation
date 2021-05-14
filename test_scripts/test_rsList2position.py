#!/usr/bin/env python3

import json
import logging

from varannot import variant

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)

# logging.getLogger("python.varannot.query").addHandler(ch)
# logging.getLogger("python.varannot.query").setLevel(logging.DEBUG)
logging.getLogger("varannot.variant").addHandler(ch)
logging.getLogger("varannot.variant").setLevel(logging.DEBUG)

#------------------------------------------------------------------------------------------------------------------

L=["rs1261047350","rs777575161","rs1250100640"]

x=variant.rs2position("rs1261047350",build="38")
print(json.dumps(x,indent=4,sort_keys=True))

x=variant.rs2position("rs1261047350",build="38",alleles=True)
print(json.dumps(x,indent=4,sort_keys=True))

x=variant.rsList2position(L,build="38")
print(json.dumps(x,indent=4,sort_keys=True))
