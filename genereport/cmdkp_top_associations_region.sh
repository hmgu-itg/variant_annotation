#!/bin/bash

curl -s -X 'GET' 'https://bioindex.hugeamp.org/api/bio/query/top-associations?q='$(echo $1|sed "s/:/%3A/"|sed "s/,//g")'&fmt=row'   -H 'accept: application/json' -o $2
