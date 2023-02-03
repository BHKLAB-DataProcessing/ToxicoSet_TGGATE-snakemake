library(ToxicoGx)
library(gdata)
library(affy)
library(Biobase)
library(xml2)
library(abind)
library("hgu133plus2hsensgcdf")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- paste0(args[1], 'download')
output_dir <- paste0(args[1], 'processed')

cdf <- "hgu133plus2hsensgcdf"
untar(file.path(input_dir, "TGGATES_human_CEL.tar.gz"), exdir = input_dir)
xx <- list.files(file.path(input_dir, "CELfiles_Human"), full.names = T, "\\.CEL$")
eset <- just.rma(filenames = xx, verbose = TRUE, cdfname = cdf)

saveRDS(eset, file.path(output_dir, "eset_Human_2382.rds"))

unlink(file.path(input_dir, "CELfiles_Human"), recursive = TRUE)