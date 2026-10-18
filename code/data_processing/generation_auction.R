# Brazil renewable auction data processing
# 1. select variables
# 2. translate variable names to English
# 3. tidy up

# @author: Konan Hara

# @date: April 20, 2025

# @software: R version 4.4.1

# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)
library(readxl)
library(lubridate)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# load generation auction results until 2022/10
# Resultado Vendedor (seller side)
col_types <- rep("numeric", 48)
col_types[c(3:16,35,46)] <- "text"
col_types[44:45] <- "date"
generation_auction_vendor <- as.data.table(read_excel(
	  file.path(pathtodata, "original", "Brazil", "generation_auction_results_202210.xlsx")
	, sheet = "Resultado Vendedor"
	, skip = 7
	, col_types = col_types))
# Resultado Consolidado (seller and buyer combined)
col_types <- rep("blank", 48)
col_types[c(3,5,32,37)] <- "text"
col_types[37] <- "date"
generation_auction_full <- as.data.table(read_excel(
	  file.path(pathtodata, "original", "Brazil", "generation_auction_results_202210.xlsx")
	, sheet = "Resultado Consolidado"
	, skip = 9
	, col_types = col_types))

# check and translate variable names
var_auction <- c(
	  "Tipo Leilão"
	, "Leilão"
	, "Edital do Leilão ANEEL"
	, "Produto"
	, "Sigla Vendedor"
	, "CNPJ"
	, "CEG"
	, "Nome da Usina"
	, "UF"
	, "Fonte"
	, "Rio/Combustível"
	, "Submercado"
	, "Investimento na Data do Leilão (R$)"
	, "Investimento Atualizado (R$)"
	, "Potência\r\n(MW)"
	, "Garantia Física\r\n(MW médio)"
	, "Fator Lote"
	, "Energia\r\nContratada no ano A\r\n(MW médio)"
	, "Energia Contratada \r\nnos demais anos\r\n(MW médio)"
	, "Total de energia contratada\r\n(MWh)"
	, "Preço de\r\nVenda na Data do Leilão\r\n(R$/MWh)"
	, "Unidade de Preço"
	, "Preço de\r\nVenda Atualizado\r\n(R$/MWh)"
	, "Início Suprimento"
	, "Fim Suprimento"
    )
var_auction_eng <- c(
	  "auction_type"
	, "auction_name"
	, "aneel_id"
	, "auction_category"
	, "bidder_name"
	, "provider_id"
	, "plant_id"
	, "plant_name"
	, "state"
	, "source_type"
	, "fuel_type"
	, "region"
	, "investment_raw"
	, "investment"
	, "capacity"
	, "physical_assurance"
	, "quantity_minimum_unit"
	, "contract_quantity_first"
	, "contract_quantity_year"
	, "contract_quantity"
	, "bid_price_raw"
	, "bid_price_unit"
	, "bid_price"
	, "contract_start_date"
	, "contract_end_date"
	)
setnames(generation_auction_vendor, var_auction, var_auction_eng)

# select variables
generation_auction_vendor <- generation_auction_vendor[, var_auction_eng, with = FALSE]

# convert non-numeric provider IDs
# e.g., "15.087.610/0001-41" -> "15087610000141"
generation_auction_vendor[, provider_id := gsub("\\.|/|-","",provider_id)]
# check whether all of the provider IDs are numeric
generation_auction_vendor[, provider_id := as.numeric(provider_id)]
# -> ok

# replace "-" with NA
generation_auction_vendor[plant_id == "-", plant_id := NA]
generation_auction_vendor[plant_name == "-", plant_name := NA]
generation_auction_vendor[state == "-", state := NA]
generation_auction_vendor[source_type == "-", source_type := NA]
generation_auction_vendor[fuel_type == "-", fuel_type := NA]

# one round of auction is determined by (auction_name,auction_category)
distinct(generation_auction_vendor[,.(auction_name,auction_category)])
auction_id_var <- c("auction_name","auction_category")

# translate source_type into English
sort(unique(generation_auction_vendor$source_type))
source_type <- c(
	  "Biomassa"
	, "Carvão"
	, "Eólica"
	, "Gás de Processo"
	, "Gás Natural"
	, "Hidrelétrica (CGH)"
	, "Hidrelétrica (PCH)"
	, "Hidrelétrica (UHE)"
	, "Óleo Combustível"
	, "Resíduo Sólido Urbano"
	, "Solar Fotovoltaica"
)
source_type_eng <- c(
	  "Biomass"
	, "Coal"
	, "Wind Energy"
	, "Process Gas"
	, "Natural Gas"
	, "Mini Hydro"
	, "Small Hydro"
	, "Large Hydro"
	, "Fuel Oil"
	, "Municipal Solid Waste"
	, "Solar Photovoltaic"
	)
source_type_translate <- data.table(source_type = source_type, source_type_eng = source_type_eng)
generation_auction_vendor <- merge(generation_auction_vendor, source_type_translate, by = "source_type", all.x = TRUE)
generation_auction_vendor[, source_type := source_type_eng]
generation_auction_vendor[, source_type_eng := NULL]

# merge with auction characteristics extracted from generation_auction_full
# check and change variable names
var_auction <- c(
	  "Leilão"
	, "Produto"
	, "Data de\r\nRealização\r\ndo leilão"
	, "Tipo de contrato\r\n(QTD/DIS)"
    )
var_auction_eng <- c(
	  "auction_name"
	, "auction_category"
	, "auction_date"
	, "contract_type"
    )
setnames(generation_auction_full, var_auction, var_auction_eng)
# collapse to auction-level
generation_auction_full <- copy(distinct(generation_auction_full))

# remove " " from auction_category for compatibility
generation_auction <- merge(generation_auction_vendor, generation_auction_full, by = auction_id_var, all.x = TRUE)

# create variables
# lead time (decimal years)
generation_auction[, lead_time := time_length(contract_start_date - auction_date, "years")]
# contract length (decimal years)
generation_auction[, contract_length := time_length(contract_end_date - contract_start_date, "years")]

# save
fwrite(generation_auction, file.path(pathtodata, "processed", "generation_auction.csv"))
