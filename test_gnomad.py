#!/usr/bin/env python3

import os
import gnomad

#import requests
#import pprint
#prettyprint = pprint.PrettyPrinter(indent=2).pprint

# def get_variant_list(gene_id, dataset="gnomad_r2_1"):
#     # Note that this is GraphQL, not JSON.
#     fmt_graphql = """
#     {
#         gene(gene_id: "%s") {
#           variants(dataset: %s) {
#             consequence
#             pos
#             rsid
#             variant_id: variantId
#           }
#         }
#       }
#     """
#     # This part will be JSON encoded, but with the GraphQL part left as a
#     # glob of text.
#     req_variantlist = {
#         "query": fmt_graphql % (gene_id, dataset)
#         }
#     print(req_variantlist)
#     response = fetch(req_variantlist)
#     return response["data"]["gene"]["variants"]


#prettyprint(get_variant_list("ENSG00000010610"))

print(os.environ.copy())

x=gnomad.getPopulationAF("rs782819098")
if x is not None:
    print(x)

x=gnomad.getPopulationAF("12-6818479-C-T")
if x is not None:
    print(x)

x=gnomad.getPopulationAF("12-6818480-C-T")
if x is not None:
    print(x)
