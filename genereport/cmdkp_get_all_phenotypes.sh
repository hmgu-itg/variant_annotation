#!/bin/bash

curl -sX 'GET' 'https://bioindex.hugeamp.org/api/portal/phenotypes' -H 'accept: application/json' | ./json2table.py 
