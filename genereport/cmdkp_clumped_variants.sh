#!/bin/bash

if [[ $# -lt 1 ]]; then
    echo "Arguments: phenotype clump"
    exit 0
fi

rm -f .token
curl -s -X 'GET' 'https://bioindex.hugeamp.org/api/bio/query/clumped-variants?q='$1'%2C'$2'&fmt=row' -H 'accept: application/json'|./json2table.py 2>.token

token=$(cat .token)
>&2 echo "Token: $token"
while [[ "$token" != "None" ]];do
    curl -s -X 'GET' 'https://bioindex.hugeamp.org/api/bio/query/clumped-variants?q='$1'%2C'$2'&fmt=row' -H 'accept: application/json'|./json2table.py 2>.token | grep -v "varId"
    token=$(cat .token)
    >&2 echo "Token: $token"
done


