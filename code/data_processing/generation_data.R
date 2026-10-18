# Brazil generation data cleaning

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

# generation data
generation_data <- fread(file.path(pathtodata, "original", "Brazil", "ANEEL"
		, "siga-comercial-generation.csv")
	, encoding = "Latin-1")

# generation agent data
generation_agent <- fread(file.path(pathtodata, "original", "Brazil", "ANEEL"
		, "electricity-generation-agents.csv")
	, encoding = "Latin-1")

# generation data cleaning
# # note: DscTipoOutorga can have several values for the same CodCEG
# length(unique(generation_data$CodCEG))
# generation_data[, dup := .N, by = .(CodCEG)]
# generation_data[dup > 1]
# variable selection and translation
variable_gen <- c(
	  "CodCEG"
	, "SigUFPrincipal"
	, "SigTipoGeracao"
	, "DscFaseUsina"
	, "DatEntradaOperacao"
	, "MdaPotenciaOutorgadaKw"
	, "DatInicioVigencia"
	, "DatFimVigencia"
)
variable_gen_eng <- c(
	  "plant_id"
	, "state"
	, "type"
	, "phase"
	, "operation_start"
	, "capacity_kw"
	, "grant_start"
	, "grant_end"
)
generation_data <- distinct(generation_data[, variable_gen, with = FALSE])
setnames(generation_data, variable_gen_eng)
# translate values into English
# generation type
sort(unique(generation_data$type))
type <- c(
	  "CGH"
	, "EOL"
	, "PCH"
	, "UFV"
	, "UHE"
	, "UTE"
	, "UTN"
)
type_eng <- c(
	  "hydro_mini"
	, "wind"
	, "hydro_small"
	, "solar"
	, "hydro_large"
	, "thermal"
	, "nuclear"
)
type_translate <- data.table(type = type, type_eng = type_eng)
generation_data <- merge(generation_data, type_translate, by = "type", all.x = TRUE)
generation_data[, type := type_eng]
generation_data[, type_eng := NULL]
# project phase
sort(unique(generation_data$phase))
phase <- c(
	  "Construção"
	, "Construção não iniciada"
	, "Operação"
)
phase_eng <- c(
	  "construction"
	, "construction_not_started"
	, "operation"
)
phase_translate <- data.table(phase = phase, phase_eng = phase_eng)
generation_data <- merge(generation_data, phase_translate, by = "phase", all.x = TRUE)
generation_data[, phase := phase_eng]
generation_data[, phase_eng := NULL]
setcolorder(generation_data, variable_gen_eng)

# generation agent data cleaning
# variable selection and translation
variable_agent <- c(
	  "NumCPFCNPJ"
	, "CodCEG"
	, "PctParticipacao"
)
variable_agent_eng <- c(
	  "firm_id"
	, "plant_id"
	, "share"
)
generation_agent <- distinct(generation_agent[, variable_agent, with = FALSE])
setnames(generation_agent, variable_agent_eng)
# firm id as int64
generation_agent[, firm_id := as.integer64(firm_id)]

# merge generation and agent data
generation_data <- merge(generation_data, generation_agent, by = "plant_id", all.x = TRUE)
# save
fwrite(generation_data, file.path(pathtodata, "processed", "generation_data.csv"))
