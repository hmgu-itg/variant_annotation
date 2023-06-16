#!/bin/bash

if [[ $# -lt 1 ]]; then
    echo "Arguments: chr:start-end {output file; default: STDOUT}"
    exit 0
fi

if [[ $# -gt 1 ]]; then
    curl -s -X 'GET' 'https://bioindex.hugeamp.org/api/bio/query/top-associations?q='$(echo $1|sed "s/:/%3A/"|sed "s/,//g")'&fmt=row' -H 'accept: application/json' |./json2table.py > $2
else
    curl -s -X 'GET' 'https://bioindex.hugeamp.org/api/bio/query/top-associations?q='$(echo $1|sed "s/:/%3A/"|sed "s/,//g")'&fmt=row' -H 'accept: application/json'|./json2table.py
fi
