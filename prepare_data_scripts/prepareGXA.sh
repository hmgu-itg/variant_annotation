#!/bin/bash

getBaselineExperiments.py | sed 's/ class="[^"]*"//g' | sed 's/<\/div><\/div>/\n/g'|sed 's/<div><div><span>/\n<div><div><span>/g' |perl -lne 'print $1 if /^\<div\>\<div\>\<span\>(.*)\<\/div\>\<\/label\>$/;' |perl -lne 'print $1."\t".$2."\t".$3 if /.*\/experiments\/(.*)?\/Results.*Design\"\>\<ul\>(.*)?\<\/ul\>\<\/a\>.*\<div\>\<span\>\<ul\>(.*)?\<\/ul\>\<\/span\>\<\/div\>/;' | sed 's/<\/li><li>/;/g' | sed 's/<li>//g' | sed 's/<\/li>//g' | grep "organism part" | grep -v disease | grep -v "cell line" 
