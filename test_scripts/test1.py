#!/usr/bin/env python3

import argparse
import sys

from modules import variant
from modules import gene

sys.stdout=open(sys.stdout.fileno(),mode='w',encoding='utf8',buffering=1)

# ------------------------------------------------------------------------------------------------------------------------

build="38"

parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument("--id", "-i", help="Required: input variation, rsID",required=True)

args=parser.parse_args()
VAR_ID=args.id

variant_data=variant.getVariantInfo(VAR_ID,build)
if variant_data is None:
    sys.exit(1)
chrpos=variant.getChrPosList(variant_data["mappings"])

all_genes=list()
for i in range(0,len(chrpos)):
    gene_list=gene.getGeneList(chrpos[i][0],chrpos[i][1],build=build)
    if gene_list:
        all_genes.extend(gene_list)

print("%s\t%d" %(VAR_ID,len(all_genes)))
