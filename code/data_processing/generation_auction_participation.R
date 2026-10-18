# data cleaning for Brazil generation auction participation data

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

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# load generation auction participation data
generation_auction_participation <- as.data.table(read_excel(
	  file.path(pathtodata, "original", "Brazil", "generation_auction_participation.xlsx")
	, sheet = "table"))

# fill in missing registered capacities from state-level information
generation_auction_participation[is.na(total_capacity)]
generation_auction_participation_supp <- as.data.table(read_excel(
	  file.path(pathtodata, "original", "Brazil", "generation_auction_participation.xlsx")
	, sheet = "supplement"))

# unify source_type values
unique(generation_auction_participation[is.na(total_capacity)]$source_type)
unique(generation_auction_participation_supp$source_type)
generation_auction_participation_supp[source_type == "Fotovoltaica", source_type := "UFV"]
generation_auction_participation_supp[source_type == "Térmica a Biomassa", source_type := "UTE"]
generation_auction_participation_supp[source_type == "Eólica", source_type := "EOL"]
generation_auction_participation_supp[source_type == "UHE", source_type := "Hidro"]
generation_auction_participation_supp[source_type == "PCH", source_type := "Hidro"]
generation_auction_participation_supp[source_type == "CGH", source_type := "Hidro"]

# aggregate values
generation_auction_participation_supp <- copy(generation_auction_participation_supp[, .(
	  number_supp = sum(number)
	, total_capacity_supp = sum(total_capacity)), by = .(auction_name,source_type)])

# fill in missing registered capacities
generation_auction_participation <- merge(generation_auction_participation
	, generation_auction_participation_supp
	, by = c("auction_name","source_type")
	, all.x = TRUE)
generation_auction_participation[is.na(total_capacity), total_capacity := total_capacity_supp]
generation_auction_participation[, number_supp := NULL]
generation_auction_participation[, total_capacity_supp := NULL]

# translate source_type values
sort(unique(generation_auction_participation$source_type))
source_type <- c(
	  "Bagaço de Cana-de-açúcar"
	, "Biomassa – bagaço de cana"
	, "Biomassa – capim elefante"
	, "Biomassa – resíduos de madeira"
	, "Cavaco de Madeira"
	, "Centrais Geradoras Hidrelétricas - CGH"
	, "EOL"
	, "Eólica"
	, "Fotovoltaica"
	, "Heliotérmica"
	, "Hidrelétrica – Ampliações"
	, "Hidrelétrica – Novas Concessões"
	, "Hidrelétricas - UHE"
	, "Hidrelétricas - UHE - Entre 30 e 50 MW"
	, "Hidrelétricas (abaixo de 50 MW)"
	, "Hidro"
	, "Pequenas Centrais Hidrelétricas - PCH"
	, "RSU"
	, "Solar Heliotérmica"
	, "Termelétricas a Biogás"
	, "Termelétricas a Biomassa"
	, "Termelétricas a Carvão"
	, "Termelétricas a Gás Natural"
	, "UFV"
	, "UTE"
)
source_type_eng <- c(
	  "Biomass Sugarcane Bagasse"
	, "Biomass Sugarcane Bagasse"
	, "Biomass Elephant Grass"
	, "Biomass Wood Chips"
	, "Biomass Wood Chips"
	, "Mini Hydro"
	, "Wind Energy"
	, "Wind Energy"
	, "Solar Photovoltaic"
	, "Solar Thermal"
	, "Hydro Expansion"
	, "Hydro New"
	, "Large Hydro"
	, "Large Hydro below 50 MW"
	, "Hydro below 50 MW"
	, "Hydro"
	, "Small Hydro"
	, "Municipal Solid Waste"
	, "Solar Thermal"
	, "Biogas"
	, "Biomass"
	, "Coal"
	, "Natural Gas"
	, "Solar Photovoltaic"
	, "Thermo"
)
source_type_translate <- data.table(source_type = source_type, source_type_eng = source_type_eng)
generation_auction_participation <- merge(generation_auction_participation, source_type_translate, by = "source_type", all.x = TRUE)
generation_auction_participation[, source_type := source_type_eng]
generation_auction_participation[, source_type_eng := NULL]

# aggregate to auction-level

# merge auction category
# load generation auction data
generation_auction <- fread(file.path(pathtodata, "processed", "generation_auction.csv"))
# restrict to c("LEN","LFA","LER")
generation_auction <- generation_auction[auction_type %in% c("LEN","LFA","LER")]

generation_auction_participation[grep("LFA", auction_name)]

# source type category
generation_auction[source_type == "Wind Energy", source_category := "wind"]
generation_auction[source_type == "Solar Photovoltaic", source_category := "solar"]
source_type_hydro <- c(
	  "Large Hydro"
	, "Small Hydro"
	, "Mini Hydro"
)
generation_auction[source_type %in% source_type_hydro, source_category := "hydro"]
source_type_thermo <- c(
	  "Biomass"
	, "Natural Gas"
	, "Fuel Oil"
	, "Coal"
	, "Process Gas"
)
generation_auction[source_type %in% source_type_thermo, source_category := "thermo"]
generation_auction[source_type == "Municipal Solid Waste", source_category := "waste"]
generation_auction_participation[source_type == "Wind Energy", source_category := "wind"]
generation_auction_participation[source_type == "Solar Photovoltaic", source_category := "solar"]
generation_auction_participation[source_type == "Solar Thermal", source_category := "thermal"]
source_type_hydro_participation <- c(
	  "Hydro"
	, "Large Hydro below 50 MW"
	, "Hydro below 50 MW"
	, "Hydro New"
	, "Hydro Expansion"
	, "Large Hydro"
	, "Small Hydro"
	, "Mini Hydro"
)
generation_auction_participation[source_type %in% source_type_hydro_participation, source_category := "hydro"]
source_type_thermo_participation <- c(
	  "Thermo"
	, "Biogas"
	, "Biomass"
	, "Biomass Elephant Grass"
	, "Biomass Sugarcane Bagasse"
	, "Biomass Wood Chips"
	, "Natural Gas"
	, "Fuel Oil"
	, "Coal"
	, "Process Gas"
)
generation_auction_participation[source_type %in% source_type_thermo_participation, source_category := "thermo"]
generation_auction_participation[source_type == "Municipal Solid Waste", source_category := "waste"]

# merge auction category
generation_auction_summary <- copy(distinct(generation_auction[, .(auction_name,auction_category,source_category)]))
generation_auction_participation <- merge(generation_auction_participation, generation_auction_summary
	, by = c("auction_name","source_category"), all.x = TRUE)

# remainings are cancelled auctions or no winners (in the source category)
generation_auction_participation[is.na(auction_category)]

# save
fwrite(generation_auction_participation, file.path(pathtodata, "processed"
	, "generation_auction_participation.csv"))
