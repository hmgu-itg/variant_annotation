#!/bin/bash

if [[ $# -lt 2 ]]; then
    echo "Arguments: pheno chr:start-end"
    exit 0
fi

rm -f .token

curl -s -X 'GET' 'https://bioindex.hugeamp.org/api/bio/query/associations?q='$1'%2C'$(echo $2|sed "s/:/%3A/"|sed "s/,//g")'&fmt=row' -H 'accept: application/json'|./json2table.py 2>.token

token=$(cat .token)
>&2 echo "Token: $token"
while [[ "$token" != "None" ]];do
    curl -sX 'GET' 'https://bioindex.hugeamp.org/api/bio/cont?token='$token -H 'accept: application/json'|./json2table.py 2>.token | grep -v "varId"
    token=$(cat .token)
    >&2 echo "Token: $token"
done

