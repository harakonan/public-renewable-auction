# Brazil corporate structure data cleaning

# @author: Konan Hara

# @date: May 26, 2025

# @software: R version 4.4.1

# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)
library(bit64)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# corporate structure data
corporate_structure <- fread(file.path(pathtodata, "original", "Brazil", "ANEEL"
		, "corporate-structure.csv")
	, encoding = "Latin-1")

# corporate structure data cleaning
variable_cs <- c(
	  "AnoExercicio"
	, "IdcTrimestreFormulario"
	, "NumCPFCNPJPaiCadeiaSocietaria"
	, "NomRazSocPaiCadeiaSocietaria"
	, "NumNivelCadeiaSocietaria"
	, "NumCPFCNPJSocio"
	, "NomRazaoSocialSocio"
	, "PctParticipacaoNivelAcima"
)
variable_cs_eng <- c(
	  "year"
	, "quarter"
	, "firm_id"
	, "firm"
	, "chain_level"
	, "share_holder_id"
	, "share_holder"
	, "share"
)
corporate_structure <- distinct(corporate_structure[, variable_cs, with = FALSE])
setnames(corporate_structure, variable_cs_eng)
# firm id as int64
corporate_structure[, share_holder_id := as.integer64(share_holder_id)]
# firm id zero as NA
corporate_structure[share_holder_id == 0, share_holder_id := NA]
# only use non NA firm id
corporate_structure <- corporate_structure[!is.na(firm_id) & !is.na(share_holder_id)]
# capitalize firm name
corporate_structure[, firm := toupper(firm)]
corporate_structure[, share_holder := toupper(share_holder)]
# save
fwrite(corporate_structure, file.path(pathtodata, "processed", "corporate_structure.csv"))
