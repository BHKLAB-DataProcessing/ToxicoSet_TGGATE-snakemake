library(readxl)
library(ToxicoGx)
library(Biobase)
library(dplyr)
library(gdata)
library(affy)
library(PharmacoGx)
library(xml2)
library(biomaRt)
library(abind)
library(SummarizedExperiment)
library(biocompute)

args = commandArgs(trailingOnly=TRUE)

arg_select <- grep(pattern = 'human|rat', x = args)
arg_select <- args[arg_select]


create_phenoData <- function(species=c("Human","Rat"), verbose = TRUE){
  if (verbose) {message("Creating phenoData object...")}
  
  #load master phenoData file from TG-GATEs
  #Master phenoData file from TG-GATEs: #14 from https://dbarchive.biosciencedbc.jp/en/open-tggates/download.html
  all_attribute <- read.csv("/pfs/getTGGATE/Open-tggates_AllAttribute.tsv", sep = "\t")
  all_attribute$SPECIES <- gsub("Rat", "R.norvegicus", all_attribute$SPECIES)
  all_attribute[grep("TNF", all_attribute[,8], useBytes = TRUE),8] <- "TNFa"
  all_attribute[grepl("\x83\xcag/kg", all_attribute[,20], useBytes = TRUE, fixed = T),20] <- "µg/kg"
  all_attribute[grepl("\x83\xcag/mL", all_attribute[,20], useBytes = TRUE, fixed = T),20] <- "µg/mL"
  all_attribute[grepl("\x83\xcaM", all_attribute[,20], useBytes = TRUE, fixed = T),20] <- "µM"
  
  #all_attribute_old$SPECIES <- gsub("Rat", "R.norvegicus", all_attribute_old$SPECIES)
  #saveRDS(all_attribute_old, "../data/TG/Open-tggates_AllAttribute.rds")
  #subset out unwanted samples (rows)
  all_attribute <- subset(all_attribute, all_attribute$BARCODE != "No ChipData" &
                            all_attribute$TEST_TYPE == "in vitro" &
                            (all_attribute$DOSE_UNIT == "µM" | all_attribute$COMPOUND_NAME == "gentamicin"))
  
  #subset out unwanted columns (all NA; for in vivo studies)
  all_attribute <- Filter(function(x) !all(is.na(x)), all_attribute)
  
  
  ## Manual mapping for drug curation
  drug_curation <- unique(all_attribute[, "COMPOUND_NAME", drop=F]) #get all unique drugs from TG-GATEs in vitro experiments
  drug_curation$unique.drugid <- drug_curation$COMPOUND_NAME #add column for lab drugid mapping; start off with all unique.drugid == tggates.drugid's
  
  lab_curation <- read.csv("/pfs/getTGGATE/drugswithids.csv", stringsAsFactors = F) #load lab annotation file, pubchem updated
  
  #exact/case match with lab annotations
  drug_curation$unique.drugid <- sapply(drug_curation$COMPOUND_NAME, function(x) if (tolower(x) %in% tolower(lab_curation$unique.drugid)) (lab_curation[which(tolower(x) == tolower(lab_curation$unique.drugid)),"unique.drugid"]) else x)
  #synonyms - manually curate
  drug_curation[drug_curation$COMPOUND_NAME == "coumarin", "unique.drugid"] <- "Coumarin"
  drug_curation[drug_curation$COMPOUND_NAME == "cyclosporine A", "unique.drugid"] <- "Cyclosporin A"
  drug_curation[drug_curation$COMPOUND_NAME == "ibuprofen", "unique.drugid"] <- "Ibuprofen"
  drug_curation[drug_curation$COMPOUND_NAME == "phenylanthranilic acid", "unique.drugid"] <- "N-phenylanthranilic acid"
  drug_curation[drug_curation$COMPOUND_NAME == "WY-14643", "unique.drugid"] <- "Pirinixic acid"
  drug_curation[drug_curation$COMPOUND_NAME == "valproic acid", "unique.drugid"] <- "Valproic acid"
  drug_curation[drug_curation$COMPOUND_NAME == "glibenclamide", "unique.drugid"] <- "Glyburide"
  drug_curation[drug_curation$COMPOUND_NAME == "adapin", "unique.drugid"] <- "Doxepin Hydrochloride"
  drug_curation[drug_curation$COMPOUND_NAME == "ketoconazole", "unique.drugid"] <- "Xolegel"
  drug_curation[drug_curation$COMPOUND_NAME == "vitamin A", "unique.drugid"] <- "Retinol"
  drug_curation[drug_curation$COMPOUND_NAME == "imipramine", "unique.drugid"] <- "Trimipramine"
  drug_curation[drug_curation$COMPOUND_NAME == "famotidine", "unique.drugid"] <- "Pepcid"
  drug_curation[drug_curation$COMPOUND_NAME == "penicillamine", "unique.drugid"] <- "Penicillamine"
  drug_curation[drug_curation$COMPOUND_NAME == "ajmaline", "unique.drugid"] <- "Ajmalin"
  drug_curation[drug_curation$COMPOUND_NAME == "nitrosodiethylamine", "unique.drugid"] <- "N-Nitrosodiethylamine"
  drug_curation[drug_curation$COMPOUND_NAME == "bromoethylamine", "unique.drugid"] <- "2-Bromoethylamine"
  drug_curation[drug_curation$COMPOUND_NAME == "galactosamine", "unique.drugid"] <- "D-Galactosamine"
  drug_curation[drug_curation$COMPOUND_NAME == "methylene dianiline", "unique.drugid"] <- "N,N'-Diphenylmethylenediamine"
  drug_curation[drug_curation$COMPOUND_NAME == "amphotericin B", "unique.drugid"] <- "Fungizone"
  drug_curation[drug_curation$COMPOUND_NAME == "naphthyl isothiocyanate", "unique.drugid"] <- "1-Naphthyl isothiocyanate"
  drug_curation[drug_curation$COMPOUND_NAME == "ethinylestradiol", "unique.drugid"] <- "Ethinyl estradiol"
  drug_curation[drug_curation$COMPOUND_NAME == "phenylanthranilic acid", "unique.drugid"] <- "N-Phenylanthranilic acid"
  drug_curation[drug_curation$COMPOUND_NAME == "carboplatin", "unique.drugid"] <- "Carboplatinum"
  drug_curation[drug_curation$COMPOUND_NAME == "ethambutol", "unique.drugid"] <- "Ethambutol Hydrochloride"
  drug_curation[drug_curation$COMPOUND_NAME == "puromycin aminonucleoside", "unique.drugid"] <- "Puromycin aminonucleoside"
  
  drug_curation <- drug_curation[,c(2,1)] #reorder columns
  names(drug_curation)[2] <- "dataset_drugid" #rename column
  rownames(drug_curation) <- drug_curation$unique.drugid #rename rows
  
  
  ################# BACK TO PHENODATA ####################
  ## Add necessary columns
  ##Batchid & conversions
  
  #Human
  if (species == "Human"){
    batch <- read_xlsx("/pfs/getTGGATE/nar-02356-data-e-2014-File006.xlsx", sheet = "Sup_table.2")
    batch <- batch[-1,]
    batch <- batch[,c(1,4)]
    names(batch) <- c("BARCODE", "CELL_NAME_TYPE_ID")
    batch <- subset(batch, batch$BARCODE != "No ChipData")
    all_attribute <- merge(all_attribute, batch, by.x = "BARCODE", by.y = "BARCODE")
    names(all_attribute)[names(all_attribute) == "CELL_NAME_TYPE_ID"] <- "batchid"
    
  } else if (species == "Rat"){ #Rat
    batch <- read_xlsx("/pfs/getTGGATE/nar-02356-data-e-2014-File006.xlsx", sheet = "Sup_table.3")
    batch <- batch[-c(1:3),]
    batch <- batch[,c(1,9)]
    colnames(batch) <- c("drugid", "batchid")
    colnames(batch) <- as.factor(colnames(batch))
    batch$drugid[batch$drugid=="TNFalpha"]<- "TNFa"
    
    all_attribute <- subset(all_attribute, all_attribute$ARR_DESIGN == "Rat230_2" & all_attribute$TEST_TYPE == "in vitro")
    all_attribute <- merge(all_attribute, batch, by.x = "COMPOUND_NAME", by.y = "drugid")
    
    conv <- readRDS("/pfs/getTGGATE/conversions_gentamicin.rds")
    #Include converted doses
    need_conversion <- subset(all_attribute,all_attribute$DOSE_UNIT != "µM",select = -c(DOSE,DOSE_UNIT))
    converted <- merge(need_conversion,conv,by.x="BARCODE",by.y="BARCODE")
    all_attribute <- subset(all_attribute,all_attribute$DOSE_UNIT == "µM")
    all_attribute <- rbind(all_attribute,converted)
  }
  all_attribute <- all_attribute[, names(all_attribute) != "DOSE_UNIT"]
  ##Merge drug mapping to lab annotations to phenoData
  all_attribute <- merge(all_attribute, drug_curation, by.x = "COMPOUND_NAME", by.y = "dataset_drugid")
  #CEL file name
  all_attribute$celfilename <- paste0(as.character(all_attribute$BARCODE),".CEL")
  # Label control or perturbation xptype
  all_attribute$xptype <- NA
  all_attribute$xptype[all_attribute$DOSE_LEVEL != "Control"] <- "perturbation"
  all_attribute$xptype[all_attribute$DOSE_LEVEL == "Control"] <- "control"
  #Substring barcode & assign cellids
  all_attribute$BARCODE<-substr(all_attribute$BARCODE,3,12)
  all_attribute$cellid <- "Hepatocyte"
  # all_attribute$cellid <- all_attribute$BARCODE
  #fix
  #all_attribute$unique.drugid[all_attribute$unique.drugid=="TNFƒ¿"]<- "TNFa"
  #take out " hr" from duration column entries
  all_attribute$SACRI_PERIOD <- gsub(" hr","",all_attribute$SACRI_PERIOD)
  #UID
  all_attribute$UID<-paste("drugid_",all_attribute$COMPOUND.Abbr.,"_",all_attribute$unique.drugid,"_",all_attribute$cellid,"_",all_attribute$SACRI_PERIOD,"hr_rep",all_attribute$INDIVIDUAL_ID, sep="")
  
  #reorder columns
  all_attribute <- all_attribute[,c("BARCODE","ARR_DESIGN","EXP_ID","GROUP_ID","INDIVIDUAL_ID",
                                    "batchid","DOSE","DOSE_LEVEL","SACRI_PERIOD","cellid",
                                    "unique.drugid","COMPOUND_NAME","COMPOUND.Abbr.","COMPOUND_NO",
                                    "DNA...","LDH...","UID","SPECIES","TEST_TYPE","SEX_TYPE",
                                    "ORGAN_ID","MATERIAL_ID","celfilename","xptype","STRAIN_TYPE",
                                    "ADM_ROUTE_TYPE")]
  
  #Rename columns
  colnames(all_attribute) <- c("samplename","chiptype","exp_id","group_id","individual_id","batchid",
                               "concentration","dose_level","duration","cellid","drugid",
                               "dataset_drugid","drugid_abbr","drugid_no","DNA","LDH",
                               "UID","species","test_type","sex_type","organ_id","material_id",
                               "celfilename","xptype","STRAIN_TYPE","ADM_ROUTE_TYPE")
  
  #Rearrange order of rows by ascending order of barcodes
  all_attribute <- dplyr::arrange(all_attribute, samplename)
  #Assign rownames as samplename
  rownames(all_attribute)<-all_attribute$samplename
  
  if (verbose) {message("phenoData object created!")}
  
  return(all_attribute)
}
create_exprsData <- function(species=c("Human","Rat"), phenoData, verbose = TRUE){
  if (verbose) {message("Creating eset object...")}
  
  if (species == "Human"){
    library("hgu133plus2hsensgcdf")
    cdf <- "hgu133plus2hsensgcdf"
    eset <- readRDS("/pfs/processTGHumArray/eset_Human_2382.rds")
  } else {
    library("rat2302rnensgcdf")
    cdf <- "rat2302rnensgcdf"
    eset <- readRDS("/pfs/processTGRatArray/eset_Rat_3276.rds")
  }
  
  
  #human eset
  #eset <- readRDS(paste("../data/TG/eset_",species,"_",nrow(phenoData),".rds", sep = ""))
  
  #rat eset 
  #eset <- readRDS("../data/TG/eset_Rat_3276.rds")
  #eset<-readRDS("../data/TG/eset_Human_2382.rds")
  storageMode(eset)<-"environment"
  #subsetting probes
  eset<-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  
  #rename??
  colnames(eset@assayData$exprs)<-substr(colnames(eset@assayData$exprs),3,12)
  colnames(eset@assayData$se.exprs)<-substr(colnames(eset@assayData$se.exprs),3,12)
  rownames(eset@protocolData@data)<-substr(rownames(eset@protocolData@data),3,12)
  rownames(eset@phenoData@data)<-substr(rownames(eset@protocolData@data),3,12)
  # #lock eset@assayData environment again
  storageMode(eset)<-"lockedEnvironment"
  #
  annotation(eset)<-"rna"
  
  #
  if (verbose) {message("eset object created!")}
  #
  return(eset)
}  
########################################################################################################


create_featureData <- function(species=c("Human","Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating featureData object...")}
  
  if (species == "Human"){
    ensembl_data <- "hsapiens_gene_ensembl"
  } else if (species == "Rat"){
    ensembl_data <- "rnorvegicus_gene_ensembl"
  }
  CELgenes <- rownames(eset@assayData$exprs)
  
  ensembl<-useMart("ensembl", dataset = ensembl_data, host="useast.ensembl.org")
  results <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=gsub("_at","",CELgenes),mart=ensembl)
  uniqueBiomaRt<-results[!duplicated(results$ensembl_gene_id),]
  names(uniqueBiomaRt)<-c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")
  
  if(species == "Rat"){
    names(uniqueBiomaRt)[1] <- "Symbol"
    uniqueBiomaRt$BEST <- NA
    uniqueBiomaRt<-arrange(uniqueBiomaRt,uniqueBiomaRt$gene_id)
    uniqueBiomaRt$gene_id <- paste(uniqueBiomaRt$gene_id,"_at", sep = "")
    rownames(uniqueBiomaRt)<-uniqueBiomaRt$gene_id
    
    if (verbose) {message("featureData object created!")}
    
    return(uniqueBiomaRt)
  }
  
  finalFeature <- uniqueBiomaRt
  
  names(finalFeature)[1] <- "Symbol"
  finalFeature$BEST <- NA
  finalFeature<-arrange(finalFeature,finalFeature$gene_id)
  finalFeature$gene_id <- paste(finalFeature$gene_id,"_at", sep = "")
  rownames(finalFeature)<-finalFeature$gene_id
  
  if (verbose) {message("featureData object created!")}
  return(finalFeature)
}

create_sensitivityProfiles <- function(phenoData, verbose = TRUE){
  if (verbose) {message("Creating Sensitivity Profiles...")}
  
  # create empty data frame of unique UID's
  sensProf <- data.frame(row.names = unique(phenoData$UID),
                         slope_recomputed = double(length(unique(phenoData$UID))),
                         auc_recomputed = double(length(unique(phenoData$UID))),
                         stringsAsFactors = F)
  sensProf[sensProf == 0] <- NA
  
  #subset out phenoData
  for (i in unique(phenoData$UID)){
    #samples of an experimental condition that are not control
    samples <- phenoData[phenoData$UID == i & phenoData$dose_level != "Control", , drop = F]
    if (nrow(samples) >= 3){
      sensProf[i,"slope_recomputed"] <- PharmacoGx::computeSlope(samples$concentration,samples$Viability)
      sensProf[i,"auc_recomputed"] <- PharmacoGx::computeAUC(samples$concentration, samples$Viability,
                                                             conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
    }
  }
  
  if (verbose) {message("Sensitivity Profiles object created!")}
  return(sensProf)
}

create_sensitivityRaw <- function(phenoData, verbose = TRUE){
  if (verbose) {message("Creating Sensitivity Raw...")}
  conc_tested<-c("Control","Low","Middle","High")
  #Create temportary data.frames to use for processing
  #Dose Info from sensInfo & reformatting
  allDose <- subset(phenoData, select=c(UID, concentration,dose_level))
  #Viability Info (& dose info) from phenoData & reformatting
  allViability <- subset(phenoData, select=c(UID, concentration, Viability,dose_level))
  
  #Create empty dose array & format
  doseArray <- array(NA, dim = c(NROW(unique(phenoData$UID)),NROW(conc_tested)))
  rownames(doseArray)<-unique(phenoData$UID)
  colnames(doseArray)<-conc_tested #for now, colnames will be the dosages instead of "doses.1", "doses.2", etc.
  #Create empty viability array & format
  viabilityArray <- array(NA, dim=c(NROW(unique(phenoData$UID)),NROW(conc_tested)))
  rownames(viabilityArray)<-unique(phenoData$UID)
  colnames(viabilityArray)<-conc_tested
  
  for (i in 1:NROW(allDose)){ #for loop iterates as long as the number of rows (#UID's)
    #go down the row of doses by iteration
    dose_row_name <- allDose[i,"UID"] #dose_row_name is the UID at that iteration
    dosage <- allDose[i,"concentration"] #dosage is the dose (concentration) at that iteration
    dose_level <- allDose[i,"dose_level"] #dose_level is the dose level (control, low, middle, high) at that iteration
    doseArray[dose_row_name,dose_level]<-dosage #index into the correct spot using the above variables, and assign 'dosage' there
  }
  for (i in 1:NROW(allViability)){ #for loop iterates as long as the number of rows (#UID's)
    dose_row_name <- allViability[i,"UID"] #UID
    dosage <- allViability[i,"concentration"] #concentration
    viability<-allViability[i,"Viability"] #viability is on the same row, one column over
    dose_level <- allViability[i,"dose_level"]
    viabilityArray[dose_row_name,dose_level]<-viability #index into the correct spot using the above variables, and assign 'viability' there
  }
  
  colnames(doseArray)<-c("Control","doses1","doses2","doses3")
  colnames(viabilityArray)<-c("Control","doses1","doses2","doses3")
  
  #Combine doseArray and viabilityArray -> 3D array
  doseViabilityArray<-abind(doseArray, viabilityArray, along=3)
  dimnames(doseViabilityArray)[[3]]<-c("Dose","Viability")
  
  sensRaw <- doseViabilityArray
  
  if (verbose) {message("Sensitivity Raw object created!")}
  
  return(sensRaw)
}

create_sensitivityInfo <- function(phenoData, doseArray, verbose = TRUE){
  if (verbose) {message("Creating Sensitivity Info...")}
  
  sensInfo <- subset(phenoData,select=c(UID, cellid, drugid, duration, individual_id))
  rownames(sensInfo) <- c()
  sensInfo <- unique(sensInfo)
  #set rownames to UIDs
  rownames(sensInfo) <- sensInfo$UID
  sensInfo <- subset(sensInfo,select=-c(UID))
  #rename columns
  names(sensInfo)<-c("cellid","drugid","duration_h","replicate")
  
  temp <- data.frame(row.names = unique(phenoData$UID),
                     Control = character(length(unique(phenoData$UID))),
                     Low = character(length(unique(phenoData$UID))),
                     Middle = character(length(unique(phenoData$UID))),
                     High = character(length(unique(phenoData$UID))),
                     stringsAsFactors = F)
  temp[temp == ""] <- NA
  
  for (i in unique(phenoData$UID)){
    c <- phenoData[phenoData$UID == i, "dose_level"]
    if ("Control" %in% c){ temp[i, "Control"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "Control"] }
    if ("Low" %in% c) { temp[i, "Low"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "Low"] }
    if ("Middle" %in% c) { temp[i, "Middle"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "Middle"] }
    if ("High" %in% c) { temp[i, "High"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "High"] }
  }
  
  sensInfo <- merge(sensInfo, temp, by.x = 0, by.y = 0, sort = F)
  rownames(sensInfo)<-sensInfo$Row.names
  sensInfo<-subset(sensInfo,select=-c(Row.names))
  
  if (verbose) {message("Sensitivity Info object created!")}
  
  return (sensInfo)
}

create_curationDrug <- function(phenoData, verbose = TRUE){
  curationDrug <- unique(subset(phenoData, select=c(drugid, dataset_drugid)))
  rownames(curationDrug) <- curationDrug$drugid
  names(curationDrug) <- c("unique.drugid", "dataset_drugid")
  
  return(curationDrug)
}

create_curationCell <- function(phenoData, verbose = TRUE){
  curationCell <- unique(subset(phenoData, select=c(cellid)))
  curationCell$dataset_cellid <- curationCell$cellid
  names(curationCell) <- c("unique.cellid", "dataset_cellid")
  rownames(curationCell) <- curationCell$unique.cellid
  
  return(curationCell)
}

create_curationTissue <- function(phenoData, verbose = TRUE){
  curationTissue <- unique(subset(phenoData, select=c(organ_id)))
  curationTissue$dataset_tissueid <- "Liver"
  names(curationTissue)[1] <- "unique.tissueid"
  rownames(curationTissue) <- "Hepatocyte"
  
  return(curationTissue)
}

create_drug <- function(phenoData, verbose = TRUE){
  drug <- unique(subset(phenoData, select=c(drugid, drugid_abbr, drugid_no)))
  rownames(drug) <- drug$drugid
  
  return(drug)
}

create_cell <- function(phenoData, verbose = TRUE){
  cell <- unique(subset(phenoData,select=c(cellid, organ_id, material_id, species, test_type)))
  names(cell)<-c("cellid","tissueid","materialid", "species","testType")
  cell$tissueid<-"Liver"
  rownames(cell) <- cell$cellid
  
  return(cell)
}

getTGGATEs <- function(species=c("Human","Rat"),
                       type=c("DNA", "LDH"),
                       verbose = TRUE){
  
  #get all elements of eset
  
  phenoData <- create_phenoData(species, verbose)
  eset <- create_exprsData(species, phenoData, verbose)
  featureData <- create_featureData(species, eset, verbose)
  
  #adding the 6 additional genes from normalized matrix
  additional_genes <- setdiff(rownames(exprs(eset)), featureData$gene_id)
  vector.is.empty <- function(x) {return(length(x) ==0 )}
  
  if(vector.is.empty(additional_genes) == FALSE){
    additional_genes_mat <- matrix(NA, nrow = length(setdiff(rownames(exprs(eset)), featureData$gene_id)), ncol = ncol(featureData))
    additional_genes_mat[,2] <- additional_genes
    rownames(additional_genes_mat) <- additional_genes_mat[,2]
    colnames(additional_genes_mat) <- colnames(featureData)  
    featureData <- rbind(featureData, additional_genes_mat)
  }else{
    print("no additional genes present")
  }
  if (verbose) {message("Putting the eset together...")}
  pData(eset) <- phenoData
  fData(eset) <- featureData
  #sorting rownames to maintain feature data mapping that m=is otherwise shuffled after converting to SE
  fData(eset) <- fData(eset)[sort(rownames(fData(eset))),]
  stopifnot(all(rownames(fData(eset)) == rownames(exprs(eset))))
  stopifnot(all(rownames(pData(eset)) == colnames(exprs(eset))))
  saveRDS(eset, "/pfs/out/eset_processsed.rds")
  if (verbose) {message("Creating summarized experiment object...")}
      
    #new_SE_TG <-  methods::as(eset, Class ="SummarizedExperiment")
    # new_SE_TG_test <-SummarizedExperiment::SummarizedExperiment(assays=list(rna=exprs(eset)),
    #                                                        rowData = fData(eset), colData=pData(eset), metadata = list(annotation = "rna"))
    # 
    new_SE_TG <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(eset)
    stopifnot(all(rownames(colData(new_SE_TG)) == rownames(pData(eset))))
    stopifnot(all(rownames(rowData(new_SE_TG)) == rownames(fData(eset))))
    
      saveRDS(new_SE_TG, "/pfs/out/new_SE_TG.rds")
    
    if (verbose) {message("Done!")}
  
  if (verbose) {message(paste("Requested ",type, sep = ""))}
  if (type == "DNA"){
    phenoData <- subset(phenoData, select=-c(LDH))
  } else {
    phenoData <- subset(phenoData, select=-c(DNA))
  }
  names(phenoData)[which(names(phenoData) == type)] <- "Viability"
  
  sensitivityProfiles <- create_sensitivityProfiles(phenoData)
  sensitivityRaw <- create_sensitivityRaw(phenoData)
  sensitivityInfo <- create_sensitivityInfo(phenoData, doseArray=as.array(sensitivityRaw[,,1]))
  
  if (verbose) {message("Creating curation objects...")}
  curationDrug <- create_curationDrug(phenoData)
  curationCell <- create_curationCell(phenoData)
  curationTissue <- create_curationTissue(phenoData)
  if (verbose) {message("Done!")}
  if (verbose) {message("Creating cell, drug objects...")}
  drug <- create_drug(phenoData)
  cell <- create_cell(phenoData)
  if(verbose) {message("Done!")}
  
  
  #the conversion function might be incorporated to the package later. This step needs to be updated then
  #source("R_CodeOcean/eSetToSE.R")
  #new_SE_TG <- eSetToSE(eset)
  # 
  # #validate the conversion
  #source("R_CodeOcean/validateESetToSEConversions.R")
  
  #validateESetToSEConversions(eSet = eset, SE = new_SE_TG)
  #if (verbose) {message("eset to SE validated")}
  # 
  if (verbose) {message("Putting ToxicoSet together...")}
  
  TGGATES <- ToxicoSet(paste("TGGATES ",species," ",type, sep = ""),
                       molecularProfiles=list("rna"= new_SE_TG),
                       cell=cell,
                       drug=drug,
                       sensitivityInfo=sensitivityInfo,
                       sensitivityRaw=sensitivityRaw,
                       sensitivityProfiles=sensitivityProfiles,
                       curationDrug=curationDrug,
                       curationCell=curationCell,
                       curationTissue=curationTissue,
                       datasetType=c("both"),
                       verify = TRUE)
  if (verbose) {message("Done!")}
  return(TGGATES)
}

if (arg_select == "human"){
TGGATES_humanldh <- getTGGATEs(species = "Human", type = "LDH")
#tggates_human_dna <- getTGGATEs(species = "Human", type = "DNA")

saveRDS(TGGATES_humanldh, "/pfs/out/TGGATES_humanldh.rds")
#saveRDS(tggates_human_dna, "QC/TGGATES_humandna.rds")

} else {

TGGATES_ratldh <- getTGGATEs(species = "Rat", type = "LDH")
#tggates_rat_dna <- getTGGATEs(species = "Rat", type = "DNA")

saveRDS(TGGATES_ratldh, "/pfs/out/TGGATES_ratldh.rds")
#saveRDS(tggates_rat_dna, "QC/TGGATES_ratdna.rds")

}
                                        
                                        
####CREATE BIOCOMPUTE OBJECT####
                                        
if (arg_select == "human"){
  
  name <- "Open TG-GATEs Human"
  pipeline_input_uri <- c(
    "ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Human/in_vitro/",
    "/pfs/getTGGATE/Open-tggates_AllAttribute.tsv",
    "/pfs/getTGGATE/drugswithids.csv",
    "https://github.com/BHKLAB-Pachyderm/getTGGATE/getTG.R"
  )
  
  pipeline_output_uri <- c(
    "/pfs/processTGHumArray/eset_Human_2382.rds",
    "/pfs/getTGGATE/Open-tggates_AllAttribute.tsv",
    "/pfs/getTGGATE/drugswithids.csv",
    "/pfs/out/TGGATES_humanldh.rds"
  )
  
} else {
  
  name <- "Open TG-GATEs Rat"
  pipeline_input_uri <- c(
    "ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vitro/",
    "/pfs/getTGGATE/Open-tggates_AllAttribute.tsv",
    "/pfs/getTGGATE/drugswithids.csv",
    "https://github.com/BHKLAB-Pachyderm/getTGGATE/getTG.R"
  )
  
  pipeline_output_uri <- c(
    "pfs/processTGRatArray/eset_Rat_3276.rds",
    "/pfs/getTGGATE/Open-tggates_AllAttribute.tsv",
    "/pfs/getTGGATE/drugswithids.csv",
    "/pfs/out/TGGATES_ratldh.rds"
  )
  
}

###########################
#####Provenance Domain#####
###########################

#Created and modified dates
#Sys.setenv(TZ = "EST")
created <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
modified <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")

#Contributions
contributors <- data.frame(
  "name" = c("Anthony Mammoliti", "Petr Smirnov", "Benjamin Haibe-Kains"),
  "affiliation" = c(rep("University Health Network", 3)),
  "email" = c("anthony.mammoliti@uhnresearch.ca", "petr.smirnov@utoronto.ca", "Benjamin.Haibe-Kains@uhnresearch.ca"),
  "contribution" = c("createdBy","createdBy","authoredBy"),
  "orcid" = c(NA,NA,"https://orcid.org/0000-0002-7684-0079"),
  stringsAsFactors = FALSE
)

#License
license <- "https://opensource.org/licenses/Apache-2.0"

#Version of biocompute object
bio_version <- "1.0.0"

#Embargo (none)
embargo <- c()

#Derived from and obsolete after (none)
derived_from <- c()
obsolete_after <- c()

#reviewers (none)
review <- c()

#compile domain
provenance <- compose_provenance_v1.3.0(
  name, bio_version, review, derived_from, obsolete_after,
  embargo, created, modified, contributors, license
)
provenance %>% convert_json()


############################
#####Description Domain#####
############################
times_rnaseq <- as.POSIXct("2020-01-20T1:05:13", format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
#Keywords and platform info
keywords <- c("Biomedical", "Toxicogenomics", "Cellline", "Drug")
platform <- c("Pachyderm", "ORCESTRA (orcestra.ca)", "Linux/Ubuntu")

#Metadata for each pipeline step
pipeline_meta <- data.frame(
  "step_number" = c("1","2","3"),
  "name" = c("Microarray processing",
             "Curated sample and treatment identifier compilation",
             "Build data object"),
  "description" = c("Normalization of microarray data", 
                    "Download of appropriate sample and treatment identifiers from GitHub (curations performed by BHK Lab - http://bhklab.ca)",
                    "Building of ORCESTRA data object"),
  "version" = c(1.0,1.0,1.0),
  stringsAsFactors = FALSE
)

#Inputs for each pipeline step
pipeline_input <- data.frame(
  "step_number" = c("1","2","2","3"),
  "filename" = c("Raw microarray data",
                 "Sample annotation data",
                 "Treatment annotations",
                 "Script for data object generation"),
  "uri" = pipeline_input_uri,
  "access_time" = c(times_rnaseq,created,created,created),
  stringsAsFactors = FALSE
)


#Outputs for each pipeline step
pipeline_output <- data.frame(
  "step_number" = c("1","2","2","3"),
  "filename" = c("Processed microarray data", 
                 "Downloaded sample annotations",
                 "Downloaded treatment annotations",
                 "Data object"),
  "uri" = pipeline_output_uri,
  "access_time" = c(times_rnaseq,created,created,created),
  stringsAsFactors = FALSE
)

#xref (none)
xref <- c()

#pipeline prereq (none)
pipeline_prerequisite <- c()

#compile domain
description <- compose_description_v1.3.0(
  keywords, xref, platform,
  pipeline_meta, pipeline_prerequisite, pipeline_input, pipeline_output
)
description %>% convert_json()


############################
######Execution Domain######
############################

script <- c()
script_driver <- c()

#software/tools and its versions used for data object creation
software_prerequisites <- data.frame(
  "name" = c("Pachyderm", "Docker Image"),
  "version" = c("1.9.3", "v2"),
  "uri" = c(
    "https://www.pachyderm.com", "https://hub.docker.com/r/bhklab/toxicogx"
  ),
  stringsAsFactors = FALSE
)

software_prerequisites[,"access_time"] <- rep(NA, length(software_prerequisites$name))
software_prerequisites[,"sha1_chksum"] <- rep(NA, length(software_prerequisites$name))

external_data_endpoints <- c()
environment_variables <- c()

execution <- compose_execution_v1.3.0(
  script, script_driver, software_prerequisites, external_data_endpoints, environment_variables
)
execution %>% convert_json()


############################
######Extension Domain######
############################

#repo of scripts/data used
scm_repository <- data.frame("extension_schema"= c("https://github.com/BHKLAB-Pachyderm"))
scm_type <- "git"
scm_commit <- c()
scm_path <- c()
scm_preview <- c()

scm <- compose_scm(scm_repository, scm_type, scm_commit, scm_path, scm_preview)
scm %>% convert_json()

extension <- compose_extension_v1.3.0(scm)
extension %>% convert_json()

############################
######Parametric Domain#####
############################

df_parametric <- data.frame(
  "param" = c(),
  "value" = c(),
  "step" = c(),
  stringsAsFactors = FALSE
)

parametric <- compose_parametric_v1.3.0(df_parametric)
parametric %>% convert_json()



############################
######Usability Domain######
############################

#usability of our data objects
text <- c(
  paste0("Pipeline for creating ", name, " data object through ORCESTRA (orcestra.ca), a platform for the reproducible and transparent processing, sharing, and analysis of biomedical data.")
)

usability <- compose_usability_v1.3.0(text)
usability %>% convert_json()


######################
######I/O Domain######
######################

input_subdomain <- data.frame(
  "filename" = c("Raw Microarray data",
                 "Sample annotation data",
                 "Treatment annotations",
                 "Script for data object generation"),
  "uri" = pipeline_input_uri,
  "access_time" = c(times_rnaseq,created,created,created),
  stringsAsFactors = FALSE
)

output_subdomain <- data.frame(
  "mediatype" = c("RDS", "tsv", "csv", "RDS"),
  "uri" = pipeline_output_uri,
  "access_time" = c(times_rnaseq,created,created,created),
  stringsAsFactors = FALSE
)

io <- compose_io_v1.3.0(input_subdomain, output_subdomain)
io %>% convert_json()


########################
######Error Domain######
########################

empirical <- c()
algorithmic <- c()

error <- compose_error(empirical, algorithmic)
error %>% convert_json()


####Retrieve Top Level Fields####
tlf <- compose_tlf_v1.3.0(
  provenance, usability, extension, description,
  execution, parametric, io, error
)
tlf %>% convert_json()


####Complete BCO####

bco <- biocompute::compose_v1.3.0(
  tlf, provenance, usability, extension, description,
  execution, parametric, io, error
)
bco %>% convert_json() %>% export_json(paste0("/pfs/out/",name,"_BCO.json")) %>% validate_checksum()                                        
                                        
                                        
