# VarAnnot - a simple variant annotation tool

A variant annotation pipeline that synthetises information from various sources to get a better understanding of the discovered associations.

In its current form, the script is probably not portable, as it was heavily embeded into the Sanger infrastructure and relies on a handful of local files. In the forseeable future, I'm planning to re-factor these dependencies to make sure we extract all required information from online sources.

### Currently the following sources are implemented:

* Ensembl 
* OMIM (requires fresh API key)
* Pubmed
* Uniprot
* GTEx (needs to be fixed - currently uses local file)
* GWAS Catalog (needs to be fixed - currently uses local file)
* ExAC
* MGI mouse database
* Expression Atlas (needs to be fixed, as recently the API has changed)

### The following functional scores are implemented:

* CADD - relies on local file
* GWAVA - relies on executable GWAVA code
* Eigen - relies on local file

### Issues needs to be fixed:

* **Genome builds**: the current version of the script only supports GRCh37, but it needs to be fixed. Probably the only supported build will be GRCh38, however won't be an easy migration as sources are on different builds, we have to harmonize them.
* **local files**: some of the features require the presence of local files. The script has to run without those files with less functionality without failing.
* **up-to-date API**: some of the API calls of the script is broken, needs to be fixed.

### Todos

1) Assemble a "minimal viable product" that provides a frame that can be extended in the future by isolating local dependencies.
2) Re-factor code: split functional units into modules.
3) Add all functionality

### Dependencies:

The following Python modules are required: requests, sys, lxml, ast, subprocess, datetime, jinja2, pandas, json

Good luck running it! 
