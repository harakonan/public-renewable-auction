# Brazil renewable auction data processing
# collect firm names from "ANEEL/AUCTION REPORT/"

# @author: Konan Hara

# @date: January 15, 2026

# @software: R version 4.5.1

# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

################################################################################
# use below for updating the firm data "generation_auction_firm.csv"
################################################################################

# load generation auction data
generation_auction <- fread(file.path(pathtodata, "processed", "generation_auction.csv"))

# focus on wind energy auctions
generation_auction <- generation_auction[source_type == "Wind Energy"]

# save columns to manually fill
generation_auction_firm <- generation_auction[,.(aneel_id,auction_date,plant_name,bidder_name)]
setorder(generation_auction_firm, "auction_date", "plant_name")
generation_auction_firm[, auction_date := NULL]
# firm name default to bidder's name
generation_auction_firm[, firm := bidder_name]
# share default to 100(%)
generation_auction_firm[, share := 100]
fwrite(generation_auction_firm, file.path(pathtodata, "processed", "generation_auction_firm_check.csv"))

# sanity check of firm data
generation_auction_firm <- fread(file.path(pathtodata, "processed", "generation_auction_firm.csv"))
# can connect all?
generation_auction_firm_check <- fread(file.path(pathtodata, "processed", "generation_auction_firm_check.csv"))
bidder_variable <- c("aneel_id","plant_name","bidder_name")
generation_auction_firm_check <- merge(
	  generation_auction_firm_check[, bidder_variable, with = FALSE]
	, generation_auction_firm
	, by = bidder_variable
	, all.x = TRUE)
generation_auction_firm_check[is.na(firm)]
# -> OK, empty

################################################################################
# use below for updating the parent data "generation_auction_parent.csv"
################################################################################

# assign parent company
generation_auction_firm <- fread(file.path(pathtodata, "processed", "generation_auction_firm.csv"))
generation_auction_firm[, firm := toupper(firm)]
generation_auction_parent_check <- generation_auction_firm[, .(share = max(share, na.rm = TRUE)), by = .(firm)]
setorder(generation_auction_parent_check, firm)
# parent default to the firm name
generation_auction_parent_check[, parent := firm]
# save
fwrite(generation_auction_parent_check, file.path(pathtodata, "processed", "generation_auction_parent_check.csv"))

# notes on assigning parents
# - web search + "corporate_structure.csv" (transmission line project)
# - SPEs are created by a parent company for a specific purpose
#	, such as holding assets or isolating financial risks.
# - many minimal share firms (e.g., 0.01%) seem to be a firm made to manage the plant
#	did not conduct parent search in details for these firms because they are
#	unlikely the bidding decision maker
# - use cnpj 16849271000156 for USINA GERADORA EÓLICA SANTA MÔNICA SPE II S.A.
# - use cnpj 12445328000109 for USINA GERADORA EÓLICA SANTA MÔNICA SPE S.A.
# - OMEGA acquired 50% of VDB F2 GERAÇÃO DE ENERGIA S.A. from EDF in 2020
# - OMEGA acquired CEA in 2019
# - CGN acquired ATLANTIC in 2019
# - EDP acquired AVENTURA in 2018
# - CPFL acquired DOBREVE in 2014

# sanity check of parent data
generation_auction_parent <- fread(file.path(pathtodata, "processed", "generation_auction_parent.csv"))
# can connect all?
generation_auction_parent_check <- fread(file.path(pathtodata, "processed", "generation_auction_parent_check.csv"))
firm_variable <- "firm"
generation_auction_parent_check <- merge(
	  generation_auction_parent_check[, firm_variable, with = FALSE]
	, generation_auction_parent
	, by = firm_variable
	, all.x = TRUE)
generation_auction_parent_check[is.na(parent)]
# -> OK, empty

################################################################################
# define representative parent for each bidder
################################################################################

# load generation auction firm data
generation_auction_firm <- fread(file.path(pathtodata, "processed", "generation_auction_firm.csv"))
generation_auction_parent <- fread(file.path(pathtodata, "processed", "generation_auction_parent.csv"))

# merge firm and parent data
generation_auction_firm[, firm := toupper(firm)]
generation_auction_parent <- merge(
	  generation_auction_firm
	, generation_auction_parent[,.(firm,parent)]
	, by = "firm"
)
generation_auction_parent <- generation_auction_parent[, .(share = sum(share))
	, by = .(aneel_id,plant_name,bidder_name,parent)]
setorder(generation_auction_parent, "aneel_id", "plant_name")

# define representative parent
bidder_variable <- c("aneel_id","plant_name","bidder_name") 
generation_auction_bidder <- distinct(generation_auction_parent[, bidder_variable, with = FALSE])
generation_auction_represent <- list()
# 1. only one parent associated with the bidder
generation_auction_parent[, dup := .N, by = bidder_variable]
generation_auction_represent[[1]] <- generation_auction_parent[dup == 1]
# 2. share > 50%
generation_auction_represent[[2]] <- generation_auction_parent[dup > 1 & share > 50]
# 3. others
# assign representative parent
generation_auction_represent_other_check <- rbindlist(generation_auction_represent)
generation_auction_represent_other_check <- merge(generation_auction_bidder
	, generation_auction_represent_other_check
	, by = bidder_variable
	, all.x = TRUE)
generation_auction_represent_other_check <- merge(
	  generation_auction_represent_other_check[is.na(parent), bidder_variable, with = FALSE]
	, generation_auction_parent
	, by = bidder_variable)
setorder(generation_auction_represent_other_check, bidder_name, aneel_id, plant_name)
setcolorder(generation_auction_represent_other_check, c("bidder_name","aneel_id","plant_name"))
# representative parent default to the parent
generation_auction_represent_other_check[, represent := parent]
# save
fwrite(generation_auction_represent_other_check, file.path(pathtodata, "processed", "generation_auction_represent_other_check.csv"))

# notes on assigning representative parents
# - web search + "corporate_structure.csv" (transmission line project)
# - take the maximum share firm if it's very close to 50% (e.g., 49%) and there
#	are no other comparable firms
# - wind farm management firm (the firm with the same name as the plant name)
#	is not the representative parent
# - ATLANTIC & CASA DOS VENTOS: ATLANTIC
#	https://www.power-technology.com/data-insights/power-plant-profile-campo-formoso-consorcio-morrinhos-brazil/?cf-view
#	https://www.gem.wiki/Morrinhos_wind_farm
# - ELETROBRAS & FIP SALUS: ELETROBRAS
#	https://www.gem.wiki/Itagua%C3%A7u_Da_Bahia_wind_farm
# - COMERCIAL MINEIRA & EMPRESA BRASILEIRA DE DESENVOLVIMENTO & VENTOS TECNOLOGIA ELÉTRICA: EMPRESA BRASILEIRA DE DESENVOLVIMENTO
#	https://cdm.unfccc.int/methodologies/Projects/PriorCDM/notifications/index_html?s=3780
# - CHUÍ & ELETROBRAS & FIP RIO BRAVO ENERGIA I: ELETROBRAS
#	https://www.gem.wiki/Chu%C3%AD_wind_farm
# - EUROPEAN ENERGY & EÓLICA ENERGIA & EÓLICA TECNOLOGIA: EUROPEAN ENERGY
#	"corporate_structure.csv"
# - ENEL & another: ENEL
#	"corporate_structure.csv"
# - ELETROBRAS & FIP RIO BRAVO ENERGIA I & (wind farm management firm): ELETROBRAS
#	"corporate_structure.csv"
# - ALUPAR & ELETROBRAS & (wind farm management firm): ALUPAR
#	https://renewablesnow.com/news/brazil-alupar-to-raise-stake-in-100-mw-aracati-wind-project-442772/
# - OMEGA ENERGIA & SEAWEST DO BRASIL & ZETA ENERGIA: OMEGA ENERGIA
#	https://www.power-technology.com/data-insights/power-plant-profile-delta-do-parnaiba-brazil/
# - BRAZIL ENERGY & SEQUOIA: SEQUOIA
#	not clear but SEQUOIA plays a larger role in this industry
# - GERABRAS PARTICIPAÇÕES & SERVTEC: SERVTEC
#	https://www.gem.wiki/Bons_Ventos_Da_Serra_wind_farm
# - DOBREVE & LUDE PARTICIPAÇÕES & (wind farm management firm): DOBREVE
#	not clear but DOBREVE plays a larger role in this industry
# - DOBREVE & WF WIND HOLDING II: DOBREVE
#	https://www.gem.wiki/Morro_Dos_Ventos_wind_farm
# - ELETROBRAS & VOLTALIA: VOLTALIA
#	https://www.thewindpower.net/windfarm_en_23251_vamcruz.php
# - OMEGA ENERGIA & ZETA ENERGIA: OMEGA ENERGIA
#	not clear but OMEGA ENERGIA plays a larger role in this industry
#	(also OMEGA is likely the parent of ZETA)

# sanity check of representative parent data
generation_auction_represent_other <- fread(file.path(pathtodata, "processed", "generation_auction_represent_other.csv"))
# can connect all?
generation_auction_represent_other_check <- fread(file.path(pathtodata, "processed", "generation_auction_represent_other_check.csv"))
bidder_parent_variable <- c(bidder_variable, "parent")
generation_auction_represent_other_check <- merge(
	  generation_auction_represent_other_check[, bidder_parent_variable, with = FALSE]
	, generation_auction_represent_other
	, by = bidder_parent_variable
	, all.x = TRUE)
generation_auction_represent_other_check[is.na(represent)]
# -> OK, empty

# collect representative parents
generation_auction_represent[[3]] <- distinct(generation_auction_represent_other[, c(bidder_variable,"represent"), with = FALSE])
setnames(generation_auction_represent[[3]], "represent", "parent")
generation_auction_represent <- rbindlist(generation_auction_represent, fill = TRUE)
generation_auction_represent[, share := NULL]
generation_auction_represent[, dup := NULL]
# save
fwrite(generation_auction_represent, file.path(pathtodata, "processed", "generation_auction_represent.csv"))

################################################################################
# list of representative parents
################################################################################

# load generation auction representative parent data
generation_auction_represent <- fread(file.path(pathtodata, "processed", "generation_auction_represent.csv"))

distinct(generation_auction_represent[,.(parent)])[order(parent)]

