#!/bin/bash
sed -E 's/^.*CRV=//g' $1 | grep -v "#"\
    |cut -d "|" -f3-7,14,19,21,12,24,25,26\
    | tr '|' '\t' | sed -E 's/\"//g' > ${1}_Clean.tsv
