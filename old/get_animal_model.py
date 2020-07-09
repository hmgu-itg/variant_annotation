'''
importing libraries that are independent from other functions.
This first block don't have to be executed each time:
'''

# library for simple web queries:
import requests, sys

# libray for parsing htlm files:
from lxml import html

# library for generate python data structure from string:
import ast

# Library to execute shell commands:
import subprocess

# Importing date library:
import datetime

# LIbrary for importing options:
import argparse

# Importing the templating library that creates html files:
egg_path='/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/Django-1.8.6-py2.7.egg'
sys.path.append(egg_path)

# Load the jinja library's namespace into the current module.
import jinja2

# Reading saved functions:
from variations import *
from variations_draw import *
from genes import *
from genes_draw import *

# OK, Let's say we have a series of stable gene IDs:
gene_IDs =["ENSG00000168826", "ENSG00000168824", "ENSG00000168818", "ENSG00000145220"]

# OK, let's see what we have:
for id in gene_IDs:




def get_mouse_ID (human_ID):
    '''
    This function looks up the mouse ID of a given gene ID
    Then checks for
    '''
    URL = "http://grch37.rest.ensembl.org/homology/id/%s?content-type=application/json&target_species=mouse&aligned=0&sequence=none&type=orthologues" % human_ID
    data = submit_REST(URL)

    # Now from the returned data we have to pick all possible homolgue IDs:
    mouse_IDs = []
    for homolog in data["data"][0]["homologies"]:
        try:
            mouse_IDs.append(homolog["target"]["id"])
        except:
            continue

    return mouse_IDs