#!/usr/bin/env bash

# get contents after the header rows in csv files
# check total row numbers before execution
# wc -l output/*.csv

# args
# $1: maximum_row_numbers
# $2: word_to_search
# $3: file_name_w/o_extension
# $4: raw_file_folder
# $5: output_file_folder

grep -A $1 $2 "$4"/"$3".csv > "$5"/"$3".csv
