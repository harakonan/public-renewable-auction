#!/usr/bin/env bash

# US spot price data processing step 1
# @author: Konan Hara
# @date: November 18, 2024

# variables:
# $1 = ${DATA}
# $2 = ${CODE}

set -euo pipefail

DATA_DIRECTORY=$1
CODE_DIRECTORY=$2

# make directory for ICE data processing
# include US CPI data
mkdir -p "$DATA_DIRECTORY/processed/step1/ICE"

# # get file and xlsx sheet names for convenience
# $2/tools/getsheetname.sh $1/original/US/ICE $1/processed/step1/ICE/

# xls to csv
# requires gnumeric on your environment
# e.g., brew install gnumeric
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/Indiana Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_Indiana.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/Mass Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_Mass.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/MID-C Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_MidC.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/NP 15 Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_NP15.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/Palo Verde Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_PaloVerde.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/PJM West Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_PJMWest.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/SP 15_1 Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_SP15_1.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-historical/SP 15_2 Hub.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2013before_SP15_2.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-2014final.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2014.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-2015final.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2015.csv"
ssconvert "$DATA_DIRECTORY/original/US/ICE/ice_electric-2016final.xls" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2016.csv"

# convert xlsx to csv
# xlsx2csv.py: https://github.com/dilshod/xlsx2csv
python3 "$CODE_DIRECTORY/tools/xlsx2csv.py" -n '2017' -q 'all' "$DATA_DIRECTORY/original/US/ICE/ice_electric-2017final.xlsx" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2017.csv"
python3 "$CODE_DIRECTORY/tools/xlsx2csv.py" -n '2018' -q 'all' "$DATA_DIRECTORY/original/US/ICE/ice_electric-2018final.xlsx" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2018.csv"
python3 "$CODE_DIRECTORY/tools/xlsx2csv.py" -n '2019' -q 'all' "$DATA_DIRECTORY/original/US/ICE/ice_electric-2019final.xlsx" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2019.csv"
python3 "$CODE_DIRECTORY/tools/xlsx2csv.py" -n '2020' -q 'all' "$DATA_DIRECTORY/original/US/ICE/ice_electric-2020final.xlsx" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2020.csv"
python3 "$CODE_DIRECTORY/tools/xlsx2csv.py" -n '2021' -q 'all' "$DATA_DIRECTORY/original/US/ICE/ice_electric-2021final.xlsx" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2021.csv"
python3 "$CODE_DIRECTORY/tools/xlsx2csv.py" -n '2022' -q 'all' "$DATA_DIRECTORY/original/US/ICE/ice_electric-2022final.xlsx" "$DATA_DIRECTORY/processed/step1/ICE/ice_electric_2022.csv"

# US CPI data processing
# # get file and xlsx sheet names for convenience
# $2/tools/getsheetname.sh $1/original/US $1/processed/step1/cpi_
python3 "$CODE_DIRECTORY/tools/xlsx2csv.py" -n 'BLS Data Series' -q 'all' "$DATA_DIRECTORY/original/US/CPI_2001_2022.xlsx" "$DATA_DIRECTORY/processed/step1/cpi_2001_2022.csv"
# get contents after the header rows in csv files
# check maximum row numbers before execution
# wc -l $1/processed/step1/cpi_2001_2022.csv
"$CODE_DIRECTORY/tools/trim-before-header.sh" 50 -e\ '"Year",' cpi_2001_2022 "$DATA_DIRECTORY/processed/step1" "$DATA_DIRECTORY/processed"
