#!/bin/bash

curl -s -X 'GET' 'https://bioindex.hugeamp.org/api/bio/query/clumped-variants?q='$1'%2C'$2'&fmt=row' -H 'accept: application/json' -o $3

