#!/usr/bin/env bash

# US spot price data processing
# download ICE's electricity wholesale daily spot price from
# https://www.eia.gov/electricity/wholesale/
# @author: Konan Hara
# @date: November 18, 2024

# variables:
# $1 = ${DATA}

# set current directory to original ICE folder
mkdir $1/original/US/ICE/
cd $1/original/US/ICE/

# download ICE's electricity wholesale daily spot price
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-historical.zip
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2014final.xls
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2015final.xls
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2016final.xls
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2017final.xlsx
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2018final.xlsx
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2019final.xlsx
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2020final.xlsx
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2021final.xlsx
wget https://www.eia.gov/electricity/wholesale/xls/archive/ice_electric-2022final.xlsx

# unzip data before 2013
mkdir $1/original/US/ICE/ice_electric-historical/
cd $1/original/US/ICE/ice_electric-historical/
unzip "$1/original/US/ICE/ice_electric-historical.zip"
