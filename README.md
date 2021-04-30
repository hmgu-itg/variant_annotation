# Association peak plotting

### Basic help

```
singularity exec <singularity options> <annotator container name> plotpeaks -h
```

```
Usage: plotpeaks -t <p-value threshold>
                 -i <input association file>
                 -c <chromosome column name>
                 -p <p-value column name>
                 -m <MAF column name>
                 --id <variant ID column name>
                 --a1 <allele 1 column name>
                 --a2 <allele 2 column name>
                 --pos <variant position column name>
                 --plink <comma separated PLINK prefixes>
                 --flank <flank size (bp)>: optional, default: 500000
                 --dbsnp <dbSNP VCF>: optional
                 -h : print this help and exit
```

# VarAnnot - a simple variant annotation tool

A variant annotation pipeline that combines information from various sources to get a better understanding of the discovered associations.


### Currently following sources are queried:

* ENSEMBL
* PubMed
* UniProt
* GTEx
* GWAS Catalog
* GnomAD
* MGI mouse database
* Gene Expression Atlas

### Following functional scores are used for annotation:

* GWAVA
* SIFT
* PolyPhen

### Basic help

```
singularity exec variantAnnotator <command> -h
```

where ```command``` is one of:
1. ```prepareData```
2. ```annotate```

### Workflow

1. Download GWAVA script from [Sanger](https://www.sanger.ac.uk/sanger/StatGen_Gwava). **On AK GWAVA is installed in ```/storage/hmgu/software/gwava/```**

1. Download and process data from ENSEMBL Regulation, GWAS Catalog and GTEx using ```prepareData``` command. **On AK prepared data are in ```/storage/hmgu/referenceData/variant_annotator_data/```**

2. Perform variant annotation using ```annotate``` command.

Example: 
```singularity exec -B /tmp -B ~ -B /storage/hmgu/ /storage/hmgu/containers/worker_2.5 annotate -v debug -i rs12570947 -o ~/output -d /storage/hmgu/referenceData/variant_annotator_data/ --html -g /storage/hmgu/software/gwava/```

The above command will create ```rs12570947.html``` and ```rs12570947.log```  in ```~/output``` directory
