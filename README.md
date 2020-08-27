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

### The following functional scores are used for annotation:

* GWAVA
* SIFT
* PolyPhen

### Basic help

```
singularity /storage/hmgu/containers/worker_2.5 <command> -h
```

where ```command``` is one of:
1. ```prepareData```
2. ```annotate```

### Workflow

1. Download GWAVA scipt and necessary data from [Sanger](https://www.sanger.ac.uk/sanger/StatGen_Gwava). **On AK GWAVA is installed in ```/storage/hmgu/software/gwava/```**

1. Download and process data from ENSEMBL Regulation, GWAS Catalog and GTEx using ```prepareData``` command. **On AK prepared data are in ```/storage/hmgu/referenceData/variant_annotator_data/```**

2. Perform variant annotation using ```annotate``` command.

Example: 
```_singularity exec -B /tmp -B ~ -B /storage/hmgu/ /storage/hmgu/containers/worker_2.5 annotate -v debug -i rs12570947 -o ~/output -d /storage/hmgu/referenceData/variant_annotator_data/ --html -g /storage/hmgu/software/gwava/```

The above command will create ```rs12570947.html``` in ```~/output``` directory
