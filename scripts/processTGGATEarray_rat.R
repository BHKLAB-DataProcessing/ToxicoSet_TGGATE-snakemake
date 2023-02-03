library(ToxicoGx)
library(gdata)
library(affy)
library(Biobase)
library(xml2)
library(abind)
library("rat2302rnensgcdf")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- paste0(args[1], 'download')
output_dir <- paste0(args[1], 'processed')
data_dir <- paste0(args[1], 'data')

cdf <- "rat2302rnensgcdf"
untar(file.path(input_dir, "TGGATES_rat_CEL.tar.gz"), exdir = input_dir)
rat_cels <- read.csv(file.path(data_dir, "rat_cel.csv"), sep = "\t")
celfn <- paste0(file.path(input_dir, "TGGATES_rat_CELfiles"), '/', rat_cels$x)
# xx <- list.files("/pfs/TGRatArray/TGGATES_rat_CELfiles", full.names = T, "\\.CEL$")
eset <- just.rma(filenames = celfn, verbose = TRUE, cdfname = cdf)

saveRDS(eset, file.path(output_dir, "eset_Rat_3276.rds"))

unlink(file.path(input_dir, "TGGATES_rat_CELfiles"), recursive = TRUE)