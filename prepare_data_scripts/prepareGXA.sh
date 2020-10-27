#!/bin/bash

getBaselineExperiments.py | sed 's/ class="[^"]*"//g' | sed 's/<\/div><\/div>/\n/g'|sed 's/<div><div><span>/\n<div><div><span>/g' |perl -lne 'print $1 if /^\<div\>\<div\>\<span\>(.*)\<\/div\>\<\/label\>$/;'

