#!/usr/bin/env bash
# go get -v gitlab.com/tyabuta/go-excel-sheets/cmd/...

# get file and sheet names of xlsx files

# input
# $1: path to xlsx file folder
# $2: path to output folder

ls $1 > $2file_name.csv
excel-sheets -d "," $1 > $2sheet_name.csv
