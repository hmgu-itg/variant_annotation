# VarAnnot - a simple variant annotation tool

A variant annotation pipeline that combines information from various sources to get a better understanding of the discovered associations.


### Currently the following sources are queried:

* ENSEMBL
* PubMed
* UniProt
* GTEx
* GWAS Catalog
* GnomAD
* MGI mouse database
* Gene Expression Atlas

### The following functional scores are implemented:

* GWAVA
* SIFT
* PolyPhen

### Basic help
```
singularity /storage/hmgu/containers/worker_2.5 <command> -h
```

where command is one of:
1. ```prepareData```
2. ```annotate```

### Usage
1. Download and process data from ENSEMBL Regulation, GWAS Catalog and GTEx
   1. 
