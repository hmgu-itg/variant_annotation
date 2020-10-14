#!/usr/bin/env python3

import json

from python.varannot import variant

L=["rs1261047350","rs777575161","rs1250100640"]

x=variant.rs2position("rs1261047350",build="38")
print(json.dumps(x,indent=4,sort_keys=True))

x=variant.rsList2position(L,build="38")
print(json.dumps(x,indent=4,sort_keys=True))
