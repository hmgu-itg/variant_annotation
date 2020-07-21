#!/usr/bin/env python

import requests
import pprint
prettyprint = pprint.PrettyPrinter(indent=2).pprint
import re

def fetch(jsondata, url="https://gnomad.broadinstitute.org/api"):
    # The server gives a generic error message if the content type isn't
    # explicitly set
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json

def get_variant_list(gene_id, dataset="gnomad_r2_1"):
    # Note that this is GraphQL, not JSON.
    fmt_graphql = """
    {
        gene(gene_id: "%s") {
          variants(dataset: %s) {
            consequence
            pos
            rsid
            variant_id: variantId
          }
        }
      }
    """
    # This part will be JSON encoded, but with the GraphQL part left as a
    # glob of text.
    req_variantlist = {
        "query": fmt_graphql % (gene_id, dataset)
        }
    print(req_variantlist)
    response = fetch(req_variantlist)
    return response["data"]["gene"]["variants"]

def getVariantAF(rs, dataset="gnomad_r3"):
    fmt_graphql = """
    {
        variant (rsid: "%s",dataset: %s) {
            chrom
            pos
            ref
            alt
            genome{populations {
                id
                ac
                an
            }
           }
        }
    }
    """
    req_af = {
        "query": fmt_graphql % (rs, dataset)
        }

    response = fetch(req_af)
    return response["data"]

#prettyprint(get_variant_list("ENSG00000010610"))
afs=getVariantAF("12-6818479-C-T")
prettyprint(afs)
for p in afs["variant"]["genome"]["populations"]:
    if not re.search("_",p["id"]) and not re.search("MALE",p["id"],re.I):
        print("%s: %.2e" %(p["id"],float(p["ac"])/float(p["an"])))

