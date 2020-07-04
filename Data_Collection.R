
setwd('~/bigdata/EarlyCancerDetection/CCMA/')

library(Biobase)
library(GEOquery)
library(readxl)
library(stringr)
library(plyr)

getPhenoFun <- function(phenoData) {
  
  info <- c('geo_accession','title','source_name_ch1','organism_ch1', 'platform_id', 'data_processing',
            'contact_institute','contact_name')
  
  idx1 <- match(info, colnames(phenoData))
  idx2 <- grep(':ch1', colnames(phenoData))
  idx <- c(idx1, idx2)
  
  phenoData <- phenoData[,idx]
  colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
  colnames(phenoData) <- gsub(' ', '.', colnames(phenoData))
  colnames(phenoData)[1:8] <- c('Accession','Title','Source','Organism','Platform','Data.Processing','Contact.Institute','Contact.Name')
  
  return(phenoData)
  
}

getMIRAnnoFun <- function(version=22, species='hsa') {
  
  mir <- readLines(paste0('data/miRBase/', version, '/mature.fa.gz'))
  
  id.idx <- grep(species, mir)
  seq.idx <- id.idx+1
  
  mir.name <- gsub('>', '', sapply(mir[id.idx], function(x) strsplit(x, ' ')[[1]][1]))
  mir.id <- sapply(mir[id.idx], function(x) strsplit(x, ' ')[[1]][2])
  mir.seq <- sapply(mir[seq.idx], function(x) strsplit(x, ' ')[[1]][1])
  
  mir.df <- data.frame(ID=mir.id, Name=mir.name, Sequence=mir.seq, row.names = mir.id, stringsAsFactors = F)
  
  return (mir.df)
  
}



for (version in c('11.0','12.0','13.0', 14:22)) {
  
  dir.name <- file.path('data/miRBase', version)
  
  if (! dir.exists(dir.name)) {
    dir.create(dir.name)
  }
  
  system(paste0('wget ftp://mirbase.org/pub/mirbase/', version, '/mature.fa.gz -P ', dir.name))
  
  mir <- getMIRAnnoFun(version = version, species = 'hsa')
  
  saveRDS(mir, paste0('data/miRBase/mir', version, '.RDS'))
  
}


mirIDNameConversionFun <- function(id, to='name2id', version=22) {
  
  mir <- readRDS(paste0('data/miRBase/mir', version, '.RDS'))
  id <- strsplit(id, ',\\s*|/|\\+')[[1]]
  
  if (to=='id2name') {
    idx <- match(id, mir$ID)
    mir.name <- paste(mir$Name[idx], collapse = ',')
    
  } else if (to=='name2id') {
    idx <- match(id, mir$Name)
    mir.name <- paste(mir$ID[idx], collapse = ',')
    
  }
  
  return(mir.name)
}


####

mir22 <- readRDS(paste0('data/miRBase/mir', 22, '.RDS'))
mir22


for (version in c(21,20,19,18,17,16,15,14,'13.0','12.0','11.0','10.0')) {
  
  mir <- readRDS(paste0('data/miRBase/mir', version, '.RDS'))
  
  idx <- intersect(rownames(mir22), rownames(mir))
  
  mir22[idx,version] <- mir[idx,'Name']
  
}

previous.id <- apply(mir22[,-c(1:4)], 1, function(x) unique(x)[!is.na(unique(x))])
previous.id

for (i in 1:nrow(mir22)) {
  
  pre <- previous.id[[i]]
  current <- mir22$Name[i]
  
  previous.id[[i]] <- paste(pre[! pre %in% current], collapse = '; ')
  
}

mir22$Previous_ID <- unlist(previous.id)
mir22$Previous_ID[mir22$Previous_ID==''] <- NA

saveRDS(mir22, file = 'data/miRBase/miRBase_10.0_22.RDS')


###############################################################################################################

######################################################################################
### 3D-Gene Human miRNA V21_1.0.0

# GSE122497
# GSE106817
# GSE137140
# GSE112264
# GSE113486
# GSE139031
# GSE110651
# GSE119159
# GSE119892
# GSE85677/GSE85680

gse <- 'GSE119159'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)

# GSE122497
colnames(phenoData)[9:13] <- c('Age','Clinical.Stage','Disease.Status','Sex','Tissue') 

table(phenoData$Disease.Status)

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='Esophageal cancer', 
                                   'Esophageal squamous cell carcinoma', 'Healthy')

# GSE106817
colnames(phenoData)[9:10] <- c('Age','Clinical.Stage')
table(phenoData$Disease.Status)

phenoData$Disease.Status <- gsub(' \\w+\\d+', '', phenoData$Title)
phenoData$Disease.Status[phenoData$Disease.Status=='non-Cancer'] <- 'Healthy'

# GSE137140
colnames(phenoData)[9:12] <- c('Age','Disease.Status','Sex','Tissue') 

table(phenoData$Disease.Status)
phenoData$Disease.Status[phenoData$Disease.Status=='Non-cancer control'] <- 'Healthy'

# GSE112264
colnames(phenoData)[9:14] <- c('Age','Clinical.M.Stage','Clinical.N.Stage','Clinical.T.Stage', 'Disease.Status','Sex') 
table(phenoData$Disease.Status)
phenoData$Disease.Status[phenoData$Disease.Status=='non-Cancer'] <- 'Healthy'

# GSE113486
colnames(phenoData)[9:14] <- c('Age','Disease.Status','Pathological.Grade','Pathological.T.Stage','Sex','Tissue') 
table(phenoData$Disease.Status)
phenoData$Disease.Status[phenoData$Disease.Status=='Non-cancer control'] <- 'Healthy'

# GSE139031
colnames(phenoData)[9:10] <- c('Age','Sex') 
table(phenoData$Disease.Status)

phenoData$Disease.Status <- gsub(' \\w+\\d+', '', phenoData$Title)
phenoData$Disease.Status[phenoData$Disease.Status=='Non-cancer control'] <- 'Healthy'

# GSE110651
colnames(phenoData)[9:12] <- c('Age','Disease.Status','New.Distant.Metastasis','Tissue') 

table(phenoData$Disease.Status)

phenoData$Disease.Status <- 'Metastatic breast cancer'

# GSE119159
colnames(phenoData)[9] <- c('Disease.Status')

table(phenoData$Disease.Status)

phenoData$Disease.Status <- 'Chronic hepatitis C'

patients <- gsub(' chronic hepatitis C', '', phenoData$Title)

phenoData$HCC.Curative.Treatment <- c(rep('No HCC Curative Treatment',70),rep('HCC Curative Treatment',69))

hcc_occ <- c('AD-1','AD-7','AD-26','AD-59','AD-116','SL-35','SL-44','SL-48','SL-153','SL-173',
             'SL-192','SL-293','SL-341','SR-3','SR-16')
  
hcc_rec <- c('AD-22','AD-24','AD-47','AD-54','AD-57','AD-75','AD-106','AD-138','AD-148','AD-151',
             'SL-15','SL-46','SL-121','SL-137','SL-151','SL-158','SL-167','SL-248','SL-252','SR-6',
             'SR-14','SR-23','SR-80','SR-103','SR-158')

phenoData$Recurrence.Occurrence <- NA

idx <- match(hcc_occ, patients)
phenoData$Recurrence.Occurrence[idx] <- 'HCC Occurrence'

idx <- match(hcc_rec, patients)
phenoData$Recurrence.Occurrence[idx] <- 'HCC Recurrence'

idx <- which(is.na(phenoData$Recurrence.Occurrence) & phenoData$HCC.Curative.Treatment=='No HCC Curative Treatment')
phenoData$Recurrence.Occurrence[idx] <- 'No HCC Occurrence'

idx <- which(is.na(phenoData$Recurrence.Occurrence) & phenoData$HCC.Curative.Treatment=='HCC Curative Treatment')
phenoData$Recurrence.Occurrence[idx] <- 'No HCC Recurrence'

phenoData$Group <- paste0(phenoData$HCC.Curative.Treatment, ', ', phenoData$Recurrence.Occurrence)

# GSE119892
colnames(phenoData)[9:14] <- c('Age','Disease.Status','Individual','Sex','Time','Tissue') 

table(phenoData$Disease.Status)
phenoData$Disease.Status <- 'Biliary Tract Cancer'

phenoData$Group <- paste0('Biliary Tract Cancer, ', phenoData$Time)
table(phenoData$Group)

# GSE85677/GSE85680
colnames(phenoData)[9:12] <- c('Disease.Status','Individual','Tissue','Treatment') 

table(phenoData$Disease.Status)
phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='CHC','Chronic Hepatitis C','Hepatocellular Carcinoma')

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Treatment)
table(phenoData$Group)


###

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData)
max(exprData)



### Annotation
#annoData <- seriesMatrix@featureData@data
#dim(annoData)

#annoData[1:10,]

#rownames(annoData) == rownames(exprData)

rownames(phenoData) == colnames(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################3
### 3D-Gene Human miRNA V20_1.0.0; 3D-Gene Human miRNA V21_1.0.0

# GSE124158
# GSE134108

gse <- 'GSE124158'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[2]] # 2
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE124158
colnames(phenoData)[9:13] <- c('Age','Sex','Disease.Status','Tissue','Clinical.Stage') 
table(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='HCC'] <- 'Hepatocellular Carcinoma'

# GSE134108
phenoData$Disease.Status <- 'Breast Cancer'

colnames(phenoData)[9:15] <- c('Age','ER.PS','HER2.IHC','KI67','Metastasis.Status','PGR.PS','Tissue') 
table(phenoData$Disease.Status)

phenoData$Metastasis.Status[phenoData$Metastasis.Status=='no-Metastasis'] <- 'No Metastasis'

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Metastasis.Status)
table(phenoData$Group)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData, na.rm = T)
max(exprData, na.rm = T)

### Annotation
#annoData <- seriesMatrix@featureData@data
#dim(annoData)

#View(annoData)

#annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(x)[[1]][1])

#rownames(annoData) == rownames(exprData)

rownames(phenoData) == colnames(exprData)


saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################3
### 3D-Gene Human miRNA V20_1.0.0

# GSE73002
# GSE59856
# GSE85679/GSE85680

gse <- 'GSE85679'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = F, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE73002
colnames(phenoData)[9:10] <- c('Disease.Status','Tissue') 
table(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='benign breast disease'] <- 'Benign Breast Disease'
phenoData$Disease.Status[phenoData$Disease.Status=='breast cancer'] <- 'Breast Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='non-cancer'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='prostate disease'] <- 'Prostate Disease'


# GSE59856
colnames(phenoData)[9:14] <- c('Age','Disease.Status','Sex','Operability','Tissue','Tumor.Stage') 
table(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='benign pancreatic or biliary tract diseases'] <- 'Benign Pancreatic or Biliary Tract Diseases'
phenoData$Disease.Status[phenoData$Disease.Status=='biliary tract cancer'] <- 'Biliary Tract Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='colon cancer'] <- 'Colon Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='esophagus cancer'] <- 'Esophagus Cancer'

phenoData$Disease.Status[phenoData$Disease.Status=='healthy control'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='liver cancer'] <- 'Liver Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='pancreatic cancer'] <- 'Pancreatic Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='stomach cancer'] <- 'Stomach Cancer'


# GSE85679/GSE85680
colnames(phenoData)[9:13] <- c('Disease.Status','Individual','Order.In.Which.Blood.Was.Collected.During.Treatment','Sample','Treatment') 
table(phenoData$Disease.Status)

colnames(phenoData)
phenoData$Disease.Status <- 'Hepatocellular Carcinoma'
phenoData$Treatment

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData, na.rm = T)
max(exprData, na.rm = T)


### Annotation
#annoData <- seriesMatrix@featureData@data
#dim(annoData)

#View(annoData)

#annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(x)[[1]][1])

#rownames(annoData) == rownames(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



##############################################################################################################

######################################################################################
### Agilent-031181 Unrestricted_Human_miRNA_V16.0_Microarray (miRBase release 16.0 miRNA ID version)

# GSE68951
# GSE48137
# GSE59993
# GSE68373
# GSE79943
# GSE54156

gse <- 'GSE54156'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE68951
colnames(phenoData)[9:12] <- c('Disease.Status','Patient.ID','Time.Point','Tissue') 
table(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='lung cancer'] <- 'Lung Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='non-cancerous lung disease (control)'] <- 'Non-Cancerous Lung Disease'

# GSE48137
colnames(phenoData)[9:11] <- c('Individual','Tissue','Treatment') 
table(phenoData$Disease.Status)

phenoData$Disease.Status <- phenoData$Group <- phenoData$Source

phenoData$Disease.Status[phenoData$Disease.Status=='serum-control'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='serum-Wilms-preCT'] <- 'Wilms Tumor'
phenoData$Disease.Status[phenoData$Disease.Status=='serum-wilms-postCT'] <- 'Wilms Tumor'

phenoData$Group[phenoData$Group=='serum-control'] <- 'Healthy'
phenoData$Group[phenoData$Group=='serum-Wilms-preCT'] <- 'Wilms Tumor, PreCT'
phenoData$Group[phenoData$Group=='serum-wilms-postCT'] <- 'Wilms Tumor, ProCT'

table(phenoData$Group)

colnames(phenoData)

# GSE59993
colnames(phenoData)[9:15] <- c('Age','Drawing.Year','Hemolysis.Score','Hemolysis.Score.Class','Plasmaid',
                               'Time.From.Surgery','Tissue') 
table(phenoData$Disease.Status)

phenoData$Disease.Status <- 'Breast Cancer'

table(phenoData$Hemolysis.Score.Class)

phenoData$Hemolysis.Score.Class <- ifelse(phenoData$Hemolysis.Score.Class=='hemolyzed','Hemolyzed','No Hemolysis')

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Hemolysis.Score.Class)
table(phenoData$Group)


# GSE68373
colnames(phenoData)[9:17] <- c('Age','Drawing.Year','ER.Class','Hemolysis.Score','Hemolysis.Score.Class',
                               'Metastasis.Class','Plasmaid','Time.From.Surgery','Tissue') 
table(phenoData$Disease.Status)

phenoData$Disease.Status <- 'Breast Cancer'

table(phenoData$Hemolysis.Score.Class)

phenoData$Hemolysis.Score.Class <- ifelse(phenoData$Hemolysis.Score.Class=='hemolyzed','Hemolyzed','No Hemolysis')
phenoData$Metastasis.Class <- ifelse(phenoData$Metastasis.Class=='metastasis','Metastasis','No Metastasis')

phenoData$Group <- NA
idx <- which(!is.na(phenoData$Hemolysis.Score.Class))
phenoData$Group[idx] <- paste0(phenoData$Disease.Status[idx], ', ', phenoData$Hemolysis.Score.Class[idx])

idx <- which(!is.na(phenoData$Metastasis.Class))
phenoData$Group[idx] <- paste0(phenoData$Disease.Status[idx], ', ', phenoData$Metastasis.Class[idx])

# GSE79943
colnames(phenoData)[9:14] <- c('Age','Time.Point','Disease.Status','Histology','Tumor.Stage','Tissue') 
table(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='	Ovarian cancer'] <- '	Ovarian Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='Paraovarian cyst'] <- 'Paraovarian Cyst'

phenoData$Time.Point <- ifelse(phenoData$Time.Point=='preoperative', 'Preoperative', 'Postoperative')

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Time.Point)

# GSE54156
colnames(phenoData)[9:12] <- c('Age','Disease.Status','Sex','Tissue') 
table(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='MM'] <- 'Multiple Myeloma'
phenoData$Disease.Status[phenoData$Disease.Status=='healthy control'] <- 'Healthy'

phenoData$Group <- phenoData$Disease.Status
phenoData$Group


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData)
max(exprData)


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 16)[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

#saveRDS(annoData, file='data/Annotation/Agilent_031181_Unrestricted_Human_miRNA_V16.0_Microarray.RDS')


idx <- grep('hsa', rownames(exprData))
idx

exprData <- exprData[idx,]

rownames(phenoData) == colnames(exprData)

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




##############################################################################################################

######################################################################################
### Agilent-031181 Unrestricted_Human_miRNA_V16.0_Microarray (miRNA_107_Sep09)

# GSE100508


gse <- 'GSE100508'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)

# GSE100508
colnames(phenoData)[9:15] <- c('Age','Sex','Melanoma.At.Sampling.Date','Clinical.Stage','Plasma.Input.ml','Primary.Melanoma.Thickness.mm','Tissue') 
table(phenoData$Disease.Status)

table(phenoData$Source)

phenoData$Disease.Status <- phenoData$Title

phenoData$Disease.Status[grep('healthy',phenoData$Disease.Status)] <- 'Healthy'
phenoData$Disease.Status[grep('tumor|risk|pre|post',phenoData$Disease.Status)] <- 'Melanoma'
phenoData$Disease.Status[grep('Kaposi',phenoData$Disease.Status)] <- 'HIV-associated Kaposi'
phenoData$Disease.Status[grep('MS',phenoData$Disease.Status)] <- 'Multiple Sclerosis'
phenoData$Disease.Status[grep('sEV',phenoData$Disease.Status)] <- 'Cell Line'



phenoData$Group <- phenoData$Source

phenoData$Group[grep('healthy',phenoData$Group)] <- 'Healthy'
phenoData$Group[grep('melanoma bearer',phenoData$Group)] <- 'Melanoma Bearer'
phenoData$Group[grep('melanoma survivor with low relapse risk',phenoData$Group)] <- 'Melanoma Survivor, Low Relapse Risk'
phenoData$Group[grep('melanoma survivor with high relapse risk',phenoData$Group)] <- 'Melanoma Survivor, High Relapse Risk'
phenoData$Group[grep('after',phenoData$Group)] <- 'Melanoma Survivor, After Resection'
phenoData$Group[grep('Kaposi',phenoData$Group)] <- 'HIV-associated Kaposi'
phenoData$Group[grep('Sclerosis',phenoData$Group)] <- 'Multiple Sclerosis'
phenoData$Group[grep('sEV',phenoData$Group)] <- c('Immature Monocyte-dervied Dendritic Cells','Liver Cell Line Huh7',
                                                  'Liver Cell Line Lx-2','SK-Hep1','Mature Monocyte-dervied Dendritic Cells',
                                                  'Monocyte-derived Macrophages')


table(phenoData$Group)

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData)
max(exprData)


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 18)[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


#saveRDS(annoData, file='data/Annotation/Agilent_031181_Unrestricted_Human_miRNA_V16.0_Microarray.RDS')
View(exprData)

idx <- grep('hsa', rownames(exprData))
idx

exprData <- exprData[rownames(annoData),]

rownames(phenoData) == colnames(exprData)

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################3
### Agilent-070156 Human_miRNA_V21.0_Microarray 046064 (miRNA ID version)

# GSE134266
# GSE139164

gse <- 'GSE139164'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)
#View(phenoData)
phenoData$Data.Processing[1]
phenoData$Source

View(phenoData)


# GSE134266
colnames(phenoData)[9:11] <- c('Age','Disease.Status','Tissue')
phenoData$Age <- as.numeric(gsub('y', '', phenoData$Age))

phenoData$Group <- phenoData$Disease.Status

phenoData$Group[phenoData$Group=='benign prostatic hyperplasia'] <- 'Benign Prostatic Hyperplasia'
phenoData$Group[phenoData$Group=='non-bone -metastatic prostate cancer'] <- 'Non-bone Metastatic Prostate Cancer '
phenoData$Group[phenoData$Group=='bone -metastatic prostate cancer'] <- 'Bone Metastatic Prostate Cancer'

phenoData$Disease.Status <- ifelse(phenoData$Group=='Benign Prostatic Hyperplasia','Benign Prostatic Hyperplasia','Prostate Cancer')

# GSE139164
colnames(phenoData)[9:14] <- c('Age','Gender','Subject.Group','Disease.Status','Tissue','Tumor.Stage')
phenoData$Age <- as.numeric(gsub('y', '', phenoData$Age))

phenoData$Disease.Status <- 'Nasopharyngeal Carcinoma'
phenoData$Subject.Group <- ifelse(phenoData$Subject.Group=='radiosensitive', 'Radiosensitive','Radioresistant')

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Subject.Group)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData, na.rm = T)
max(exprData, na.rm = T)


### Annotation
annoData <- seriesMatrix@featureData@data
dim(annoData)

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 21)[[1]][1])
annoData$Name

saveRDS(annoData, file='data/Annotation/Agilent_070156_Human_miRNA_V21.0_Microarray_046064.RDS')


#idx <- grep('hsa', rownames(exprData))
#idx

exprData <- exprData[rownames(annoData),]

rownames(phenoData) == colnames(exprData)

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


#####################################################################################

######################################################################################
### Agilent-021827 Human miRNA Microarray [miRNA_107_Sep09_2_105]

# Release 15, some v11.0, 9.1

# GSE43160
# GSE32099

gse <- 'GSE32099'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE43160
colnames(phenoData)[9:11] <- c('Age','Gender','Survival.Time.Group')

phenoData$Disease.Status <- 'Nasopharyngeal Carcinoma'
phenoData$Survival.Time.Group <- ifelse(phenoData$Survival.Time.Group=='shorter-survival time', 'Shorter Survival Time','Longer Survival Time')


phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Survival.Time.Group)


# GSE32099
colnames(phenoData)[9:18] <- c('Cell.Type','Disease.Status','Legend','MSKCC','OS.Status','OS.Time','PSF.Status','PSF.Time','Subject','Treatment')

filter <- is.na(phenoData$Disease.Status)

phenoData <- phenoData[-filter,]
phenoData$Disease.Status <- 'Renal Cell Carcinoma'
phenoData$Treatment <- ifelse(phenoData$Treatment=='basal','Basal','Sunitinib')

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Treatment)

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))

### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

#annoData$ID <- gsub('_st', '', annoData$ID)
#annoData$ID <- gsub('-star', '*', annoData$ID)

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


saveRDS(annoData, file='data/Annotation/Affymetrix_Multispecies_miRNA-1_Array.RDS')

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### Agilent-021827 Human miRNA Microarray G4470C (Feature Number version)

# GSE39833/GSE40247
# GSE40246/GSE40247
# GSE55139


gse <- 'GSE55139'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE39833/GSE40247
colnames(phenoData)[9:14] <- c('Age','Fraction','Gender','Presence.of.Data.After.Operation.1.Present.0.Not.Tested','Tissue','Tumor.Stage')
phenoData$Age <- as.numeric(gsub('y', '', phenoData$Age))

colnames(phenoData)

phenoData$Disease.Status <- ifelse(grepl('Healthy Control', phenoData$Source), 'Healthy','Colorectal Cancer')
phenoData$Group <- ifelse(grepl('Healthy Control', phenoData$Source), 'Healthy',paste0('Colorectal Cancer, TNM Stage ', phenoData$Tumor.Stage))


unique(phenoData$Group)

# GSE40246/GSE40247
colnames(phenoData)[9:15] <- c('Age','Biological.Source.of.Exosomes','Disease.Status','Sex','Microrna.Source','Time','Tumor.Stage')
phenoData$Age <- as.numeric(gsub('y', '', phenoData$Age))

phenoData$Disease.Status <- 'Colorectal Cancer'
phenoData$Group <- paste0('Colorectal Cancer, TNM Stage ', phenoData$Tumor.Stage)


#GSE55139

phenoData <- phenoData[,-c(9:10,13)]

colnames(phenoData)[9:12] <- c('Age','N','Stage','T')

phenoData$Disease.Status <- 'Colorectal Cancer'
phenoData$Group <- ifelse(grepl('Pre-operative', phenoData$Title), paste0('Colorectal Cancer, Pre-operative'),
                          paste0('Colorectal Cancer, Post-operative'))


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- data.frame(exprData,stringsAsFactors = F)



### Annotation
annoData <- seriesMatrix@featureData@data

#idx <- which(annoData$miRNA_ID=='hsa-miR-720')
#idx


#View(exprData[idx,])

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- as.character(sapply(annoData$ACCESSION_STRING, function(x) strsplit(x, '|', fixed = T)[[1]][4]))
annoData$Name

saveRDS(annoData, file='data/Annotation/Agilent-021827_Human_miRNA_Microarray_G4470C_Feature_Number_version.RDS')

exprData <- exprData[rownames(annoData),]
dim(exprData)

exprData$ID_REF <- annoData$Name
exprData$ID_REF


### select most informative probe, MAX IQR
iqr <- apply(exprData[,-ncol(exprData)], 1, function(v) IQR(v, na.rm = T))
iqr

exprData$IQR <- iqr

probeIndex <- ddply(exprData, .(ID_REF), summarise, probe=which.max(IQR))
probeIndex[1:5,]

probes <- c()
for (i in 1:nrow(probeIndex)) {
  probe <- rownames(exprData)[which(exprData$ID_REF==probeIndex$ID_REF[i])]
  print (probe)
  
  probes <- c(probes, probe[probeIndex$probe[i]])
}

filter <- which(is.na(probes))
filter

#probes <- probes[-filter]
#probes

exprData <- exprData[probes,]
exprData

rownames(exprData) <- exprData$ID_REF

exprData <- exprData[,-c(ncol(exprData)-1,ncol(exprData))]

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


######################################################################################
### Agilent-029297 Human miRNA Microarray (Feature Number version)

# GSE34052
# GSE48485

gse <- 'GSE34052'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE39833/GSE40247
colnames(phenoData)[9:14] <- c('Age','Fraction','Gender','Presence.of.Data.After.Operation.1.Present.0.Not.Tested','Tissue','Tumor.Stage')
phenoData$Age <- as.numeric(gsub('y', '', phenoData$Age))

colnames(phenoData)

phenoData$Disease.Status <- ifelse(grepl('Healthy Control', phenoData$Source), 'Healthy','Colorectal Cancer')
phenoData$Group <- ifelse(grepl('Healthy Control', phenoData$Source), 'Healthy',paste0('Colorectal Cancer, TNM Stage ', phenoData$Tumor.Stage))


unique(phenoData$Group)

# GSE40246/GSE40247
colnames(phenoData)[9:15] <- c('Age','Biological.Source.of.Exosomes','Disease.Status','Sex','Microrna.Source','Time','Tumor.Stage')
phenoData$Age <- as.numeric(gsub('y', '', phenoData$Age))

phenoData$Disease.Status <- 'Colorectal Cancer'
phenoData$Group <- paste0('Colorectal Cancer, TNM Stage ', phenoData$Tumor.Stage)


#GSE55139

phenoData <- phenoData[,-c(9:10,13)]

colnames(phenoData)[9:12] <- c('Age','N','Stage','T')

phenoData$Disease.Status <- 'Colorectal Cancer'
phenoData$Group <- ifelse(grepl('Pre-operative', phenoData$Title), paste0('Colorectal Cancer, Pre-operative'),
                          paste0('Colorectal Cancer, Post-operative'))


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- data.frame(exprData,stringsAsFactors = F)



### Annotation
annoData <- seriesMatrix@featureData@data

#idx <- which(annoData$miRNA_ID=='hsa-miR-720')
#idx


#View(exprData[idx,])

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- as.character(sapply(annoData$ACCESSION_STRING, function(x) strsplit(x, '|', fixed = T)[[1]][4]))
annoData$Name

saveRDS(annoData, file='data/Annotation/Agilent-021827_Human_miRNA_Microarray_G4470C_Feature_Number_version.RDS')

exprData <- exprData[rownames(annoData),]
dim(exprData)

exprData$ID_REF <- annoData$Name
exprData$ID_REF


### select most informative probe, MAX IQR
iqr <- apply(exprData[,-ncol(exprData)], 1, function(v) IQR(v, na.rm = T))
iqr

exprData$IQR <- iqr

probeIndex <- ddply(exprData, .(ID_REF), summarise, probe=which.max(IQR))
probeIndex[1:5,]

probes <- c()
for (i in 1:nrow(probeIndex)) {
  probe <- rownames(exprData)[which(exprData$ID_REF==probeIndex$ID_REF[i])]
  print (probe)
  
  probes <- c(probes, probe[probeIndex$probe[i]])
}

filter <- which(is.na(probes))
filter

#probes <- probes[-filter]
#probes

exprData <- exprData[probes,]
exprData

rownames(exprData) <- exprData$ID_REF

exprData <- exprData[,-c(ncol(exprData)-1,ncol(exprData))]

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))












######################################################################################
### Agilent-046064 Human miRNA Microarray, Release 19.0, 8x60K (miRNA_ID version)

# GSE112840

gse <- 'GSE112840'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix
seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE112840
colnames(phenoData)[9:11] <- c('Age','Sex','Tissue')

unique(phenoData$Source)

phenoData$Disease.Status <- ifelse(grepl('Healthy control', phenoData$Source), 'Healthy','Esophageal Squamous Cell Carcinoma')
phenoData$Group <- ifelse(grepl('Healthy control', phenoData$Source), 'Healthy','Esophageal Squamous Cell Carcinoma, Pre-operative')

unique(phenoData$Group)

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('MIMAT', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- annoData$miRNA_ID

#saveRDS(annoData, file='data/Annotation/Affymetrix_Multispecies_miRNA-1_Array.RDS')

exprData <- exprData[rownames(annoData),]

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


######################################################################################
### Agilent-026867 Human miRNA Microarray_Hsa_Rel14_V5 (Probe Name version)

# GSE53179

gse <- 'GSE53179'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix
seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE53179
colnames(phenoData)[9:10] <- c('Disease.Status','Tissue')

phenoData$Group <- ifelse(phenoData$Disease.Status=='healthy', 'Healthy','Breast Cancer, ER+/HER2-')
phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='healthy', 'Healthy','Breast Cancer')


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$ACCESSION_STRING, function(x) strsplit(x, '|', fixed = T)[[1]][4])
annoData$Name

exprData <- exprData[rownames(annoData),]

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))





######################################################################################
### Agilent-028035 hsa_miRNA14_Virus

# GSE43329

gse <- 'GSE43329'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix
seriesMatrix <- seriesMatrix[[1]]

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source
View(phenoData)

# GSE43329
colnames(phenoData)[9:11] <- c('Age','Sex','Tissue')

phenoData$Disease.Status <- ifelse(grepl('healthy', phenoData$Source), 'Healthy','Nasopharyngeal Carcinoma')
phenoData$Group <- phenoData$Disease.Status


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data
View(annoData)

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




######################################################################################
### [miRNA-1] Affymetrix Multispecies miRNA-1 Array

# GSE46729

gse <- 'GSE46729'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE46729
colnames(phenoData)[9:12] <- c('Age','Batch','Sex','Status')

phenoData$Disease.Status <- ifelse(grepl('control', phenoData$Title), 'Healthy','Non-Small Cell Lung Cancer')
phenoData$Group <- phenoData$Status

phenoData$Group[phenoData$Group=='control'] <- 'Healthy'
phenoData$Group[phenoData$Group=='surgery'] <- 'Non-Small Cell Lung Cancer, Prior to Surgery'
phenoData$Group[phenoData$Group=='1 month'] <- 'Non-Small Cell Lung Cancer, One Month Post Surgery'
phenoData$Group[phenoData$Group=='6 month'] <- 'Non-Small Cell Lung Cancer, Six Months Post Surgery'




platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$ID <- gsub('_st', '', annoData$ID)
annoData$ID <- gsub('-star', '*', annoData$ID)

annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################3
### [miRNA-2] Affymetrix Multispecies miRNA-2 Array

# GSE44281
# GSE39845 - whole blood, tissue

gse <- 'GSE39845'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)
# GSE44281
colnames(phenoData)[9:14] <- c('Age','Batch','Chip.Lot.AB','Pair.Index','Race','Status')

phenoData$Disease.Status <- ifelse(phenoData$Status=='case','Breast Cancer', 'Healthy')
phenoData$Group <- ifelse(phenoData$Status=='case', 'Developed Breast Cancer','Remained Breast Cancer Free')



# GSE39845
colnames(phenoData)[9:11] <- c('Affected.Site','Disease.Status','Tissue')

keep <- grep('Whole blood', phenoData$Source)
keep

phenoData <- phenoData[keep,] 

phenoData$Group <- phenoData$Disease.Status
phenoData$Group[phenoData$Group=='healthy'] <- 'Healthy'
phenoData$Group[phenoData$Group=='Stage 1 colon cancer'] <- 'Colon Cancer, Stage I'
phenoData$Group[phenoData$Group=='Stage 2 colon cancer'] <- 'Colon Cancer, Stage II'
phenoData$Group[phenoData$Group=='Stage 3 colon cancer'] <- 'Colon Cancer, Stage III'
phenoData$Group[phenoData$Group=='Stage 4 colon cancer'] <- 'Colon Cancer, Stage IV'
phenoData$Group[phenoData$Group=='Stage 1 rectal cancer'] <- 'Rectal Cancer, Stage I'
phenoData$Group[phenoData$Group=='Stage 2 rectal cancer'] <- 'Rectal Cancer, Stage II'
phenoData$Group[phenoData$Group=='Stage 3 rectal cancer'] <- 'Rectal Cancer, Stage III'
phenoData$Group[phenoData$Group=='Stage 4 rectal cancer'] <- 'Rectal Cancer, Stage IV'

phenoData$Disease.Status[grep('healthy',phenoData$Disease.Status)] <- 'Healthy'
phenoData$Disease.Status[grep('colon',phenoData$Disease.Status)] <- 'Colon Cancer'
phenoData$Disease.Status[grep('rectal',phenoData$Disease.Status)] <- 'Rectal Cancer'


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$ID)
annoData <- annoData[idx,]

filter <- grep('hp|v11', annoData$ID)
annoData <- annoData[-filter,]

annoData$ID <- gsub('_st', '', annoData$ID)
annoData$ID <- gsub('-star', '*', annoData$ID)

annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name

sum(annoData$Name=='NA')

exprData <- exprData[rownames(annoData),]
dim(exprData)


rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

exprData <- exprData[,rownames(phenoData)]

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


######################################################################################
### [miRNA-4] Affymetrix Multispecies miRNA-4 Array

# GSE85589
# GSE124489


gse <- 'GSE124489'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE85589
colnames(phenoData)[9:11] <- c('Age','CA19-9','Sex')

phenoData$Disease.Status <- gsub('\\(|\\)', '', str_extract(phenoData$Source, '\\(\\w+\\s*\\w+\\)'))
phenoData$Disease.Status

unique(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='pancreatic cancer'] <- 'Pancreatic Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='intrahepatic cholangiocarcinoma'] <- 'Intrahepatic Cholangiocarcinoma'
phenoData$Disease.Status[phenoData$Disease.Status=='stomach cancer'] <- 'Stomach Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='colorectal cancer'] <- 'Colorectal Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='GIST'] <- 'Gastrointestinal Stromal Tumours'
phenoData$Disease.Status[phenoData$Disease.Status=='cholelithiasis'] <- 'Cholelithiasis'
phenoData$Disease.Status[phenoData$Disease.Status=='healthy control'] <- 'Healthy'

phenoData$Group <- phenoData$Disease.Status

# GSE124489
colnames(phenoData)[9:11] <- c('Disease.Status','Ethnicity','Tissue')

phenoData$Disease.Status[phenoData$Disease.Status=='healthy'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='newly diagnosed multiple myeloma'] <- 'Multiple Myeloma'

phenoData$Group <- phenoData$Disease.Status

unique(phenoData$Group)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data
annoData$Name <- annoData$Accession


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




######################################################################################
### [miRNA-4] Affymetrix Multispecies miRNA-4 Array [ProbeSet ID version]

# GSE98181
# GSE141208

gse <- 'GSE98181'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE98181
colnames(phenoData)[9:11] <- c('Sex','Group','Tissue')

phenoData$Disease.Status <- ifelse(phenoData$Group=='case', 'Breast Cancer', 'Healthy')
phenoData$Group <- ifelse(phenoData$Group=='case', 'Developed Breast Cancer', 'Remained Breast Cancer Free')

# GSE141208
colnames(phenoData)[9:10] <- c('Disease.Status','Tissue')
unique(phenoData$Disease.Status)

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='Retinoblastoma', 'Retinoblastoma', 'Healthy')
phenoData$Tissue[phenoData$Tissue=='EV'] <- 'Extracellular Vesicles'

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Tissue)

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$`Transcript ID(Array Design)`)
annoData <- annoData[idx,]

keep <- which(annoData$`Sequence Type`=='miRNA')
annoData <- annoData[keep,]

annoData$Name <- annoData$Accession


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



#################################################################################
### Applied Biosystems Human Taqman MicroRNA Array v3.0

# GSE59520

gse <- 'GSE59520'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE59520
colnames(phenoData)[9] <- c('Sample.Type')

table(phenoData$Sample.Type)

phenoData$Disease.Status <- NA

phenoData$Disease.Status[grepl('CONTROL', phenoData$Title) & grepl('normal', phenoData$Sample.Type)] <- 'Healthy'
phenoData$Disease.Status[grepl('SE', phenoData$Title)] <- 'Seminomas'
phenoData$Disease.Status[grepl('NS', phenoData$Title)] <- 'Non-seminomas'
phenoData$Disease.Status[is.na(phenoData$Disease.Status)] <- 'Non-Germ Cell Tumor'


phenoData$Group <- phenoData$Disease.Status

phenoData$Group[phenoData$Sample.Type=='serum from patient with epidermoid cyst'] <- 'Non-Germ Cell Tumor, Epidermoid Cyst'
phenoData$Group[phenoData$Sample.Type=='serum from patient with epidermis cyst'] <- 'Non-Germ Cell Tumor, Epidermis Cyst'
phenoData$Group[phenoData$Sample.Type=='serum from patient with low grade liposarcoma'] <- 'Non-Germ Cell Tumor, Low Grade Liposarcoma'
phenoData$Group[phenoData$Sample.Type=='serum from patient with neuro endocrine tumor'] <- 'Non-Germ Cell Tumor, Neuro Endocrine Tumor'
phenoData$Group[phenoData$Sample.Type=='serum from patient with paratesticular hemorrhage'] <- 'Non-Germ Cell Tumor, Paratesticular Hemorrhage'

phenoData$Group[phenoData$Sample.Type=='serum from patient with yolk sac tumor'] <- 'Non-seminomas, Yolk Sac Tumor'
phenoData$Group[phenoData$Sample.Type=='serum from patient with teratoma'] <- 'Non-seminomas, Teratoma'
phenoData$Group[phenoData$Sample.Type=='serum from patient with embryonal carcinoma'] <- 'Non-seminomas, Embryonal Carcinoma'
phenoData$Group[phenoData$Sample.Type=='serum from patient with embryonal carcinoma + yolk sac tumor'] <- 'Non-seminomas, Embryonal Carcinoma and Yolk Sac Tumor'
phenoData$Group[phenoData$Sample.Type=='serum from patient with embryonal carcinoma + teratoma'] <- 'Non-seminomas, Embryonal Carcinoma, Teratoma'

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID, '_', annoData$Panel)

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

exprData[filter,]

idx <- intersect(rownames(annoData), rownames(exprData))
idx

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name

rownames(phenoData) == colnames(exprData)

colnames(exprData) <- rownames(phenoData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



#################################################################################
### Applied Biosystems TaqMan Array Human MicroRNA A+B Cards Set v3.0
# GSE47125
# GSE65708

### Applied Biosystems TaqMan Array Human MicroRNA Cards (A+B Card Set v3)
# GSE50013

### Applied Biosystems Taqman Low Density Array Human microRNA Card A
# GSE64591

gse <- 'GSE64591'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE47125
colnames(phenoData)[9] <- 'Group'
phenoData$Disease.Status <- 'Prostate Cancer'

phenoData$Group <- gsub('biochemical recurrence', 'Biochemical Recurrence', phenoData$Group)
phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Group)

# GSE65708
colnames(phenoData)[9:11] <- c('Age', 'Disease.Status', 'Sex')

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='Chronic hepatitis B (CHB)', 'Chronic Hepatitis B', 'Hepatocellular Carcinoma')
phenoData$Group <- phenoData$Disease.Status

# GSE50013
colnames(phenoData)[9:12] <- c('Age', 'Sex', 'Disease.Status', 'Tissue')

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='control group', 'Healthy', 'Hepatocellular Carcinoma')
phenoData$Group <- phenoData$Disease.Status

# GSE64591
colnames(phenoData)[9:13] <- c('Age', 'Disease.Status', 'Sex', 'Morphology', 'Smoking.Status')
unique(phenoData$Group)

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='Non-cancer', 'Healthy', 'Non-Small Cell Lung Cancer')
phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Smoking.Status)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Annottion

annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')

idx <- grep('hsa', annoData$ID)
annoData <- annoData[idx,]


### Expression

### GSE47125
system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47125/suppl/GSE47125.xls.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE47125.xls.gz')

exprData <- read_excel('data/fromGEO/GSE47125.xls')


### GSE65708
system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65708/suppl/GSE65708_Matrix_normalized_by_MCR.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE65708_Matrix_normalized_by_MCR.txt.gz')

exprData <- read.delim('data/fromGEO/GSE65708_Matrix_normalized_by_MCR.txt', header = T, sep = '\t', stringsAsFactors = F)


### GSE50013
exprData <- exprs(seriesMatrix)
dim(exprData)

idx <- intersect(rownames(annoData), rownames(exprData))
idx
dim(exprData)

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name


# GSE64591

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64591/suppl/GSE64591_normalized.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE64591_normalized.txt.gz')

exprData <- read.delim('data/fromGEO/GSE64591_normalized.txt', sep = '\t', header = T, 
                       skip = 7, stringsAsFactors = F)

exprData <- data.frame(exprData, stringsAsFactors = F)
filter <- which(duplicated(exprData$ID_REF))
exprData <- exprData[-filter,]

rownames(exprData) <- exprData$ID_REF
exprData <- exprData[,-1]

rownames(annoData) <- annoData$ID.2

###
filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

rownames(exprData)[filter]

idx <- intersect(rownames(annoData), rownames(exprData))
idx
dim(exprData)

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name

rownames(phenoData) == colnames(exprData)
gsub(' ', '.', phenoData$Title) == colnames(exprData)

colnames(exprData) <- rownames(phenoData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




#################################################################################
### Applied Biosystems Taqman Low Density Array Human microRNA Card A, Applied Biosystems Taqman Low Density Array Human microRNA Card B

# GSE46355

# GSE46355-GPL17039, Applied Biosystems Taqman Low Density Array Human microRNA Card A
# GSE46355-GPL17040, Applied Biosystems Taqman Low Density Array Human microRNA Card B

gse <- 'GSE46355'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

# GSE46355-GPL17039
seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

# GSE46355-GPL17040
seriesMatrix <- seriesMatrix[[2]]
seriesMatrix


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE46355-GPL17039, GSE46355-GPL17040
colnames(phenoData)
colnames(phenoData)[9:19] <- c('Age', 'Associated.DCIS', 'Sex', 'Grade', 'Histolgical.Tumour.Type','Invasive.Tumor.Size',
                               'Lymph.Node.Status','Lymphovascular.Invasion','Perineural.Invasion','Tumor.Stage','Whole.Tumor.Size')
unique(phenoData$Group)

phenoData$Disease.Status <- ifelse(phenoData$Source=='Control Healthy', 'Healthy', 'Breast Cancer')
phenoData$Group <- ifelse(phenoData$Source=='Control Healthy', 'Healthy', 'Luminal A Breast Cancer')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

idx <- grep('hsa', annoData$ID)
annoData <- annoData[idx,]


### Expression
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

exprData[filter,]

idx <- intersect(rownames(annoData), rownames(exprData))
idx
dim(exprData)

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name

rownames(phenoData) == colnames(exprData)
colnames(exprData) <- rownames(phenoData)


saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


#========================================================================================================================================
### Annotation
# 
system('wget https://assets.thermofisher.com/TFS-Assets/LSG/brochures/megaplex-pools-array-card-content.xlsx -P data/fromGEO/')

# GPL16850
gpl <- read.delim('data/fromGEO/GPL16850.txt', header = T, sep = '\t', stringsAsFactors = F)
gpl


system('wget https://assets.thermofisher.com/TFS-Assets/LSG/brochures/megaplex-pools-array-card-content.xlsx -P data/fromGEO/')

annoData.A <- read_excel('data/fromGEO/megaplex-pools-array-card-content.xlsx', sheet = 'Human A')
annoData.A

annoData.A$ID <- gpl$ID[1:nrow(annoData.A)]
annoData.A$Panel <- 'A'

View(annoData.B)

annoData.B <- read_excel('data/fromGEO/megaplex-pools-array-card-content.xlsx', sheet = 'Human B')
annoData.B
View(annoData.B)

annoData.B$ID <- gpl$ID[(nrow(annoData.A)+1):nrow(gpl)]

annoData.B$Panel <- 'B'

annoData <- rbind(annoData.A, annoData.B)
annoData

annoData <- data.frame(annoData, stringsAsFactors = F)

idx <- grep('hsa', annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls.)
idx
length(idx)

annoData <- annoData[idx,]

annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls. <- gsub('-5$', '-5p', annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls.)
annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls. <- gsub('-3$', '-3p', annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls.)

annoData$Name <- sapply(annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls., 
                        function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 22)[[1]][1])
annoData$Name

sum(annoData$Name=='NA')

View(annoData[annoData$Name=='NA',])

idx <- which(annoData$Name=='NA')
idx

annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls.[idx] <- 
  c('hsa-miR-1-3p','hsa-miR-137-3p','hsa-miR-147b-3p','hsa-miR-198','hsa-miR-203a-3p','hsa-miR-217-5p',
    'hsa-miR-301b-3p','hsa-miR-320a-3p','hsa-miR-375-3p','hsa-miR-636','hsa-miR-520b-3p','hsa-miR-520e-3p',
    'hsa-miR-549a-3p','hsa-miR-566','hsa-miR-190b-5p','hsa-miR-1254','hsa-miR-1249-3p')

annoData$Name <- sapply(annoData$miRBase.ID..v22..or.NCBI.Name..for.Controls., 
                        function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 22)[[1]][1])
annoData$Name

sum(annoData$Name=='NA')

View(annoData[annoData$Name=='NA',])

idx <- which(annoData$Name=='NA')
idx

annoData$Name[idx] <- c('MIMAT0003230', 'MIMAT0005905')

rownames(annoData) <- annoData$ID

#rownames(annoData) <- ifelse(annoData$Panel=='A', annoData$Assay.Name, paste(annoData$Assay.Name, annoData$Assay.ID, sep = '-'))

saveRDS(annoData, file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
#========================================================================================================================================



######################################################################################
### Exiqon Human miRCURY LNA Universal RT miRNA PCR Human Serum/Plasma focus panel (V3.0)

# Release 19

# GSE117063

gse <- 'GSE117063'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE117063
colnames(phenoData)[9:10] <- c('Disease.Status','Tissue')

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='Control', 'Healthy', 'Diffuse Large B Cell Lymphoma')
phenoData$Group <- phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]



### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '19')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### Exiqon Human miRCURY LNA Universal RT miRNA PCR Human Serum/Plasma focus panel (V3.R)

# Release 15

# GSE51410

gse <- 'GSE51410'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)

colnames(phenoData)[9:12] <- c('Barretts.Esophagus.Size', 'Disease.Status','Sample.Type','Tumor.Stage')


phenoData$Disease.Status[phenoData$Disease.Status=='Healthy control'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='Barrett\'s esophagus'] <- 'Barretts Esophagus'

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### Exiqon LNA RT-PCR Human panels (1 & 2)

# Release 15

# GSE41922/GSE42128

gse <- 'GSE41922'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE41922/GSE42128
colnames(phenoData)
colnames(phenoData)[9:13] <- c('Age','Disease.Status','Estrogen.Receptor','HER2','Node.Positivity')

phenoData$dDisease.Status[phenoData$Disease.Status=='Healthy volunteer'] <- 'Healthy'
phenoData$Group <- phenoData$Disease.Status


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


saveRDS(annoData, file='data/Annotation/Exiqon_LNA_RT-PCR_Human_panels_1_2.RDS')


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




######################################################################################
### Exiqon miRCURY LNA microRNA array, 7th generation [miRBase v18, condensed Probe_ID version]

# Release 18

# GSE113956

gse <- 'GSE113956'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE113956
colnames(phenoData)
colnames(phenoData)[9:18] <- c('Age','Disease.Status','Drinking','Lymph.Node.Metastasis','Pathological.Grade',
                               'Sex','Smoking','Survival.Status.After.5.Years','Tissue','TNM.Stage')

phenoData$Disease.Status <- ifelse(grepl('healthy', phenoData$Source), 'Healthy', 'Oral Squamous Cell Carcinoma')
phenoData$Group <- phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$Human_miRNA)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '18')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

#saveRDS(annoData, file='data/Annotation/Exiqon_LNA_RT-PCR_Human_panels_1_2.RDS')


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### Exiqon miRCURY LNA™ Universal RT miRNA PCR Human Serum/Plasma focus panel

# Release 15

# GSE57661

gse <- 'GSE57661'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE57661
colnames(phenoData)
colnames(phenoData)[9:18] <- c('Age','ER.Status','HER2.Status','Menopausal.Status','Nodal.Status','GPL','Tissue','Tumor.Grade',
                               'Tumor.Size.mm','WHO.Diagnosis')

phenoData$Disease.Status <- 'Breast Cancer'
phenoData$Group <- ifelse(grepl('pre-operative', phenoData$Title), 'Breast Cancer, Pre-operative', 'Breast Cancer, Pro-operative')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57661/suppl/GSE57661_Normalized_data.txt.gz -P data/fromGEO/')

exprData <-read.table('data/fromGEO/GSE57661_Normalized_data.txt.gz', header = T, stringsAsFactors = F, sep = '\t')
exprData[1:5,1:5]

idx <- grep('hsa', exprData$Sample)
exprData <- exprData[idx,]

mir.name <- exprData$Sample
mir.name

mir.name <- sapply(mir.name, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
mir.name

sum(mir.name=='NA')

filter <- which(mir.name=='NA')
mir.name <- mir.name[-filter]

exprData <- exprData[-filter,]
exprData <- exprData[,-1]

samples <- gsub('-operative serum-', '.', phenoData$Title)
idx <- match(samples, colnames(exprData))
idx


samples == colnames(exprData[,idx])

colnames(exprData)[idx] <- phenoData$Accession

rownames(exprData) <- mir.name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


######################################################################################
### miRCURY LNA™ Universal RT miRNA PCR Human panel I and II

# Release 20

# GSE143489/GSE143491
# GSE143490/GSE143491

gse <- 'GSE143490'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE143489/GSE143491
colnames(phenoData)[9] <- 'Protocols'
phenoData$Disease.Status <- 'Head and Neck Tumors'
phenoData$Group <- '20Gy irradiation'

# GSE143490/GSE143491
colnames(phenoData)[9] <- 'Protocols'
phenoData$Disease.Status <- 'Head and Neck Tumors'
phenoData$Group <- 'No irradiation'

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))



### Expression

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143489/suppl/GSE143489_non_normalized_matrix.txt.gz -P data/fromGEO/')
system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143490/suppl/GSE143490_non_normalized_matrix.txt.gz -P data/fromGEO/')


exprData <-read.table('data/fromGEO/GSE143490_non_normalized_matrix.txt.gz', header = T, stringsAsFactors = F, sep = '\t')
exprData[1:5,1:5]

idx <- grep('hsa', exprData$ID_REF)
exprData <- exprData[idx,]

mir.name <- exprData$ID_REF
mir.name

mir.name <- sapply(mir.name, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '20')[[1]][1])
mir.name

sum(mir.name=='NA')

# GSE143489/GSE143491
exprData <- exprData[,-c(1:2)]
samples <- gsub('circulating miRNA \\[|\\]', '', phenoData$Title)
samples


# GSE143490/GSE143491
exprData <- exprData[,-c(1)]
samples <- paste0('X', gsub('circulating miRNA \\[|\\]', '', phenoData$Title))
samples

idx <- match(samples, colnames(exprData))
idx


samples == colnames(exprData[,idx])

colnames(exprData)[idx] <- phenoData$Accession

rownames(exprData) <- mir.name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


### GSE143491
expr1 <- readRDS('data/rData/GSE143489_GPL27988_Expression.RDS')
expr2 <- readRDS('data/rData/GSE143490')




######################################################################################
### Exiqon human V3 microRNA PCR panel I+II

# GSE65071

gse <- 'GSE65071'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE65071
colnames(phenoData)[9:10] <- c('Disease.Status', 'OS.Type')
phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='healthy', 'Healthy', 'Osteosarcoma')

phenoData$Group <- phenoData$Disease.Status
phenoData$Group[phenoData$OS.Type=='metastatic'] <- 'Osteosarcoma, Metastatic'
phenoData$Group[phenoData$OS.Type=='localized'] <- 'Osteosarcoma, Localized'

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))

### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '19')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### Exiqon microRNA Ready-to-Use PCR, Human panel I+II, V2.M

# GSE37472

gse <- 'GSE37472'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE37472
colnames(phenoData)

colnames(phenoData)[9:14] <- c('Activation.Agent', 'Cell.Type', 'Disease.State', 'Disease.Status', 'Time', 'Tissue')
phenoData$Disease.Status[phenoData$Disease.Status=='non-cancer'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='carcinoma in situ'] <- 'Oral Carcinoma In Situ'
phenoData$Disease.Status[grepl('Squamous cell carcinoma', phenoData$Disease.Status)] <- 'Oral Squamous Cell Carcinoma'

phenoData$Group <- ifelse(phenoData$Disease.Status=='Healthy', 'Healthy', paste0(phenoData$Disease.Status, ', ', phenoData$Time)) 

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))

### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### Exiqon miRCURY LNA microRNA array; 7th generation; batch 208500-2, 208510; lot 35004 - hsa, mmu & rno (miRBase 18.0)

# GSE40738

gse <- 'GSE40738'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE40738
colnames(phenoData)
table(phenoData$Disease.Status)

colnames(phenoData)[9:15] <- c('Age', 'Disease.Status', 'Sex', 'Individual', 'Lung.Disease', 'Time', 'Tissue')
phenoData$Disease.Status[phenoData$Disease.Status=='lung without lung cancer or non-cancerous nodule'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='adenocarcinoma of lung'] <- 'Lung Adenocarcinoma'
phenoData$Disease.Status[phenoData$Disease.Status=='squamous cell carcinoma of lung'] <- 'Lung Squamous Cell Carcinoma'
phenoData$Disease.Status[phenoData$Disease.Status=='poorly differentiated non-small cell cancer of lung'] <- 'Poorly Differentiated Lung Non-small Cell Cancer of Lung'

phenoData$Disease.Status[phenoData$Disease.Status=='neuroendocrine carcinoma of lung'] <- 'Neuroendocrine Carcinoma of Lung'
phenoData$Disease.Status[phenoData$Disease.Status=='large cell carcinoma of lung'] <- 'Large Cell Carcinoma of Lung'
phenoData$Disease.Status[phenoData$Disease.Status=='pleural plaque of lung'] <- 'Pleural Plaque of Lung'

phenoData$Disease.Status[phenoData$Disease.Status=='granuloma of lung'] <- 'Granuloma of Lung'
phenoData$Disease.Status[phenoData$Disease.Status=='hamartoma of lung'] <- 'Hamartoma of Lung'
phenoData$Disease.Status[phenoData$Disease.Status=='amyloid of lung'] <- 'Amyloid of Lung'

phenoData$Disease.Status[phenoData$Disease.Status=='pneumonia of lung'] <- 'Pneumonia of Lung'
phenoData$Disease.Status[phenoData$Disease.Status=='lung after resection of lung cancer'] <- 'Lung After Resection of Lung Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='sarcomatoid carcinoma of lung'] <- 'Sarcomatoid Carcinoma of Lung'
phenoData$Disease.Status[phenoData$Disease.Status=='subpleural fibrosis of lung'] <- 'Subpleural Fibrosis of Lung'

phenoData$Group <- phenoData$Disease.Status

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='Healthy', 'Healthy', 'Lung Cancer')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))

### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

annoData$miRNA_ID <- sapply(annoData$miRNA_ID_LIST, function(x) strsplit(x, ',')[[1]][1])

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '18')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name


### Later
exprData <- exprData[-which(duplicated(rownames(exprData))),]

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


######################################################################################
### Exiqon miRCURY LNA Universal RT microRNA PCR Human Panel I

# GSE102166

gse <- 'GSE102166'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE40738
colnames(phenoData)
table(phenoData$Disease.Status)

colnames(phenoData)[9:13] <- c('Age', 'Metastatic.Stage', 'Progression.Free.Survival', 'Response', 'Sex')

phenoData$Disease.Status <- 'Melanoma'

phenoData$Group <- ifelse(grepl('during', phenoData$Source), paste0('Metastatic Melanoma, During MAPKi Treatment, ', phenoData$Response),
                          paste0('Metastatic Melanoma, Before MAPKi Treatment, ', phenoData$Response))

unique(phenoData$Group)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))

### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '20')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


#################################################################################
### Taqman Array Human MicroRNA A Card v2.0

# GSE75392
# GSE47652

gse <- 'GSE75392'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE75392
colnames(phenoData)[9] <- 'Disease.Status'

phenoData$Disease.Status[phenoData$Disease.Status=='healty'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='helthy'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='chronic myeloid leukemia'] <- 'Chronic Myeloid Leukemia'

phenoData$Source[phenoData$Source=='exosome'] <- 'Exosome'

phenoData$Group <- paste0(phenoData$Disease.Status, ', ', phenoData$Source)


# GSE47652
colnames(phenoData)[9:11] <- c('Disease.Status','Sample.Group','Tissue')

phenoData$Disease.Status <- ifelse(grepl('healthy volunteers', phenoData$Disease.Status), 'Healthy', 'Chronic Myeloid Leukemia')

phenoData$Group <- phenoData$Sample.Group

phenoData$Group[phenoData$Group=='sustained CMR for more than 6 months after discontinuation of imatinib (STOP-IM group)'] <- 'Chronic Myeloid Leukemia, STOP-IM'
phenoData$Group[phenoData$Group=='receiving imatinib with CMR; undetermined minimal disease (UMD)'] <- 'Chronic Myeloid Leukemia, UMD'
phenoData$Group[phenoData$Group=='healthy volunteers (controls)'] <- 'Healthy'


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

rownames(exprData)==rownames(annoData)

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(exprData) <- annoData$Probe_ID


### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

#rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID)

#idx <- which(!rownames(annoData) %in% rownames(exprData))
#idx
#annoData[idx,]

#filter <- which(!rownames(exprData) %in% rownames(annoData))
#filter

#exprData[filter,]

#exprData <- exprData[-filter,]

idx <- intersect(rownames(annoData), rownames(exprData))
idx

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name

rownames(phenoData) == colnames(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



#################################################################################
### TaqMan Low Density Arrays Card A and B v2

# GSE70080

gse <- 'GSE70080'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


# GSE70080
colnames(phenoData)[9] <- c('Disease.Status')

phenoData$Disease.Status[phenoData$Disease.Status=='COPD'] <- 'Chronic Obstructive Pulmonary Disease'
phenoData$Disease.Status[phenoData$Disease.Status=='lung'] <- 'Lung Cancer'
phenoData$Disease.Status[phenoData$Disease.Status=='normal'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='kalibrator'] <- 'Kalibrator'

phenoData$Group <- phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID)

idx <- grep('hsa', annoData$ID)
annoData <- annoData[idx,]


### Expression
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
dim(exprData)

#idx <- which(!rownames(annoData) %in% rownames(exprData))
#idx
#annoData[idx,]

#filter <- which(!rownames(exprData) %in% rownames(annoData))
#filter

#exprData[filter,]

#exprData <- exprData[-filter,]

idx <- intersect(rownames(annoData), rownames(exprData))
idx
dim(exprData)

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name
View(exprData)

rownames(phenoData) == colnames(exprData)
#exprData <- exprData[,rownames(phenoData)]

colnames(exprData) <- rownames(phenoData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




#################################################################################
### TaqMan microRNA Custom Low-Density Array

# GSE76462

gse <- 'GSE76462'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9:11] <- c('Screening.Group','Sample.Type','Tumor.Type')

unique(phenoData$Sample.Type)
unique(phenoData$Tumor.Type)

phenoData$Disease.Status <- phenoData$Tumor.Type
phenoData$Disease.Status[is.na(phenoData$Disease.Status)] <- 'Healthy'

phenoData$Group <- phenoData$Disease.Status


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
dim(exprData)

### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

annoData$ID.2 <- paste0(gsub('\\d+$', '', annoData$ID), annoData$Assay.ID)
annoData$ID.2

#saveRDS(annoData, file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')

rownames(annoData) <- annoData$ID.2

rownames(exprData) <- gsub('mir','miR',rownames(exprData))

#rownames(exprData) <- str_extract(rownames(exprData), '\\d+$')

#idx <- which(!rownames(annoData) %in% rownames(exprData))
#idx
#annoData[idx,]

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

rownames(exprData)[filter]

idx <- intersect(rownames(annoData), rownames(exprData))
idx
dim(exprData)

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name

rownames(phenoData) == colnames(exprData)
#exprData <- exprData[,rownames(phenoData)]

colnames(exprData) <- rownames(phenoData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


#################################################################################
### TaqMan® Array Human MicroRNA Cards Set v2.0 A/B

# GSE67075

gse <- 'GSE67075'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9:10] <- c('Disease.Status','Tissue')

phenoData$Group <- phenoData$Disease.Status

phenoData$Disease.Status[phenoData$Disease.Status=='Control'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='Polyp'] <- 'Colon Polyp'
phenoData$Disease.Status[grepl('CRC', phenoData$Disease.Status)] <- 'Colorectal Cancer'

phenoData$Group[phenoData$Group=='Control'] <- 'Healthy'
phenoData$Group[phenoData$Group=='Polyp'] <- 'Colon Polyp'

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67075/suppl/GSE67075_Normalized_data.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE67075_Normalized_data.txt.gz')

exprData <- read.delim('data/fromGEO/GSE67075_Normalized_data.txt', sep = '\t', header = T, stringsAsFactors = F)

### select most informative probe, MAX IQR
iqr <- apply(exprData[,-1], 1, function(v) IQR(v, na.rm = T))
iqr

exprData$IQR <- iqr

probeIndex <- ddply(exprData, .(ID_REF), summarise, probe=which.max(IQR))
probeIndex[1:5,]

View(probeIndex)

probes <- c()
for (i in 1:nrow(probeIndex)) {
  probe <- rownames(exprData)[which(exprData$ID_REF==probeIndex$ID_REF[i])]
  print (probe)
  
  probes <- c(probes, probe[probeIndex$probe[i]])
}

filter <- which(is.na(probes))
filter

#probes <- probes[-filter]
#probes

exprData <- exprData[probes,]
exprData

rownames(exprData) <- exprData$ID_REF

exprData <- exprData[,-c(1,ncol(exprData))]

idx <- grep('hsa', rownames(exprData))
exprData <- exprData[idx,]

mir <- gsub('-\\d+$', '', rownames(exprData))
mir

mir <- sapply(mir, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
mir

sum(mir=='NA')

filter <- which(mir=='NA')
exprData <- exprData[-filter,]
mir <- mir[-filter]


rownames(exprData) <- mir
phenoData$Title==colnames(exprData)
rownames(phenoData) == colnames(exprData)

idx <- match(phenoData$Title, colnames(exprData))
idx

exprData <- exprData[,idx]

colnames(exprData) <- phenoData$Accession

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))

### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



#################################################################################
### TaqMan® OpenArray® Human MicroRNA Panel (miRBase v14)

# GSE63108

gse <- 'GSE63108'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9:10] <- c('Disease.Status', 'Tumor.Stage')

phenoData$Disease.Status[phenoData$Disease.Status=='none'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='non-dysplastic Barrett\'s Oesophagus'] <- 'Non-dysplastic Barretts Oesophagus'

phenoData$Group <- phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- paste0(annoData$Assay.ID, '_', annoData$Assay.Name)

idx <- grep('hsa', annoData$ID)
annoData <- annoData[idx,]
View(annoData)

#rownames(exprData) <- gsub('mir','miR',rownames(exprData))

#rownames(exprData) <- str_extract(rownames(exprData), '^\\d+')

#idx <- which(!rownames(annoData) %in% rownames(exprData))
#idx
#annoData[idx,]

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

rownames(exprData)[filter]

rownames(exprData) <- gsub('HSA-MIR', 'hsa-miR', rownames(exprData))


idx <- intersect(rownames(annoData), rownames(exprData))
idx
dim(exprData)

exprData <- exprData[idx,]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData[idx,]$Name
View(exprData)

rownames(phenoData) == colnames(exprData)
#exprData <- exprData[,rownames(phenoData)]

colnames(exprData) <- rownames(phenoData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))





######################################################################################
### State Key Laboratory Human microRNA array 1858

# Release 15

# GSE118613
# GSE93850


gse <- 'GSE118613'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE118613
colnames(phenoData)[9:10] <- c('Disease.Status','Tissue')

unique(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='normal control'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='head-neck tumor'] <- 'Head and Neck Tumor'
phenoData$Disease.Status[phenoData$Disease.Status=='nasopharyngeal carcinoma'] <- 'Nasopharyngeal Carcinoma'

phenoData$Group <- phenoData$Disease.Status



# GSE93850
colnames(phenoData)[9] <- c('Disease.Status')

unique(phenoData$Disease.Status)

phenoData$Disease.Status[phenoData$Disease.Status=='healthy volunteer'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='glioblastoma multiforme patient'] <- 'Glioblastoma Multiforme'

phenoData$Group <- phenoData$Disease.Status


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

### Annotation
annoData <- seriesMatrix@featureData@data
annoData$Name <- annoData$Mature_Acc

rownames(exprData) == annoData$ID

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

filter <- which(duplicated(rownames(exprData)))
filter

exprData <- exprData[-filter,]

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### CombiMatrix_Human_1.6K_miRNA array [condensed version]

# Release 11.0 (10.0 reported)

# GSE16512

gse <- 'GSE16512'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9:11] <- c('Sex','Tumor.Site','Tumor.Stage')

phenoData$Disease.Status <- ifelse(phenoData$Tumor.Site=='Normal', 'Healthy', paste0(phenoData$Tumor.Site, ' Cancer'))
phenoData$Group <- phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '11.0')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### BGISEQ-500 (Homo sapiens)

# GSE128004

gse <- 'GSE128004'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9:11] <- c('Patient.Diagnosis','Tissue.Compartment','Tissue')

phenoData$Disease.Status <- phenoData$Patient.Diagnosis
phenoData$Disease.Status[phenoData$Disease.Status=='Normal control'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='GNBi'] <- 'Ganglioneuroblastoma Intermixed'
phenoData$Disease.Status[phenoData$Disease.Status=='NB'] <- 'Neuroblastoma'

phenoData$Group <- phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128004/suppl/GSE128004_RPKM-NB-exomalmiRNA.xls.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE128004_RPKM-NB-exomalmiRNA.xls.gz')

exprData <- read_excel('data/fromGEO/GSE128004_RPKM-NB-exomalmiRNA.xls')

exprData[1:5,1:5]

exprData <- data.frame(exprData, stringsAsFactors = F)

rownames(exprData) <- exprData$gene

exprData <- exprData[,-c(1:2)]


mir <- sapply(rownames(exprData), function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '21')[[1]][1])
mir

sum(mir=='NA')

rownames(exprData) <- mir

rownames(phenoData) == colnames(exprData)
colnames(exprData) <- phenoData$Accession

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




######################################################################################
### Nanostring human miRNA panel (NS_H_MIR_V3A)

# Release 21

# GSE112462

gse <- 'GSE112462'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9:10] <- c('Cell.Type','Tissue')

unique(phenoData$Cell.Type)

phenoData$Disease.Status <- phenoData$Cell.Type
phenoData$Disease.Status[phenoData$Disease.Status=='NORMAL/NONDISEASED'] <- 'Healthy'
phenoData$Disease.Status[phenoData$Disease.Status=='ASTROCYTOMAS'] <- 'Astrocytomas'
phenoData$Disease.Status[phenoData$Disease.Status=='OLIGODENDROGLIOMA'] <- 'Oligodendroglioma'
phenoData$Disease.Status[phenoData$Disease.Status=='GLIOBLASTOMA'] <- 'Glioblastoma'

phenoData$Group <-phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '21')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


######################################################################################
### Serum_microRMA

# GSE101841

gse <- 'GSE101841'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9:11] <- c('Tissue.Type', 'Tissue', 'Trastuzumab.Therapy')
phenoData$Disease.Status <- 'Breast Cancer'

phenoData$Group <- ifelse(phenoData$Trastuzumab.Therapy=='sensitive', 'HER2 Positive Metastatic Breast Cancer, Sensitive to Trastuzumab Therapy', 
                          'HER2 Positive Metastatic Breast Cancer, Resistant to Trastuzumab Therapy')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]


annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '18')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

filter <- which(duplicated(rownames(exprData)))
filter

exprData <- exprData[-filter,]

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################
### febit Homo Sapiens miRBase 13.0

# Release 15

# GSE20994

gse <- 'GSE20994'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

colnames(phenoData)[9] <- 'Disease.Status'

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='normal control', 'Healthy', 'Melanoma')
phenoData$Group <- phenoData$Disease.Status


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




######################################################################################
### febit Homo sapiens miRBase 15.0

# Release 15

# GSE31309

gse <- 'GSE31309'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


colnames(phenoData)

colnames(phenoData)[9:14] <- c('Sex','TNM.HER2.Status','TNM.N.Stage','TNM.T.Stage','TNM.Tumor.Grading','Year.of.Birth')
phenoData$Disease.Status <- ifelse(phenoData$TNM.HER2.Status=='healthy', 'Healthy', 'Breast Cancer')

phenoData$Group <- phenoData$Disease.Status
phenoData$Group[phenoData$TNM.HER2.Status=='negative'] <- 'Breast Cancer, HER2 Negative'
phenoData$Group[phenoData$TNM.HER2.Status=='positive'] <- 'Breast Cancer, HER2 Positive'


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')


exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


######################################################################################
### Illumina Human v2 MicroRNA expression beadchip

# Release 15, some v11.0, 9.1

# GSE25609
# GSE41526
# GSE22981

gse <- 'GSE22981'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

# GSE25609
colnames(phenoData)[9:10] <- c('Group', 'Tissue')

unique(phenoData$Group)

phenoData$Disease.Status <- phenoData$Group

phenoData$Disease.Status[grepl('Control', phenoData$Disease.Status)] <- 'Healthy'
phenoData$Disease.Status[grepl('Adenoma', phenoData$Disease.Status)] <- 'Advanced Adenoma'
phenoData$Disease.Status[grepl('colorectal', phenoData$Disease.Status)] <- 'Colorectal Cancer'

phenoData$Group[phenoData$Group=='Control'] <- 'Healthy'
phenoData$Group[phenoData$Group=='Advanced Adenoma after removal'] <- 'Advanced Adenoma After Removal'
phenoData$Group[phenoData$Group=='colorectal cancer'] <- 'Colorectal Cancer'
phenoData$Group[phenoData$Group=='colorectal cancer after surgical intervention'] <- 'Colorectal Cancer After Surgical Intervention'

# GSE41526
colnames(phenoData)[9:11] <- c('Sex', 'Sample.Type', 'Tissue')

unique(phenoData$Sample.Type)
phenoData$Disease.Status <- phenoData$Sample.Type

phenoData$Disease.Status[grepl('control', phenoData$Disease.Status)] <- 'Healthy'
phenoData$Disease.Status[grepl('breast cancer', phenoData$Disease.Status)] <- 'Breast Cancer'
phenoData$Disease.Status[grepl('lung cancer', phenoData$Disease.Status)] <- 'Lung Cancer'
phenoData$Disease.Status[grepl('colon cancer', phenoData$Disease.Status)] <- 'Colon Cancer'

phenoData$Group <- phenoData$Disease.Status

phenoData$Group[phenoData$Sample.Type=='pre-resection breast cancer'] <- 'Breast Cancer, Pre-resection'
phenoData$Group[phenoData$Sample.Type=='post-resection breast cancer'] <- 'Breast Cancer, Post-resection'


# GSE22981
colnames(phenoData)[9:11] <- c('Disease.Status', 'Race', 'Tissue')

phenoData$Disease.Status <- ifelse(phenoData$Disease.Status=='Control', 'Healthy', 'Breast Cancer')
phenoData$Group <- phenoData$Disease.Status

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

filter <- grep(':', annoData$miRNA_ID)
filter
annoData <- annoData[-filter,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


saveRDS(annoData, file='data/Annotation/Illumina_Human_v2_MicroRNA_expression_beadchip.RDS')

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))


####################################################################################################################################################


####################################################################################################################################################




######################################################################################
### Agilent Human 8x60k MicroRNA Array (miRBase 16.0) with non human spike-ins (Agilent design ID 047209)

# Release 20

# E-MTAB-4667


gse <- 'E-MTAB-4667'

system(paste0('mkdir data/fromArrayExpress/', gse))

### Phenotype

system(paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', gse, '/', gse, '.sdrf.txt -P data/fromArrayExpress/', gse, '/'))

phenoData <- read.table(paste0('data/fromArrayExpress/', gse, '/', gse, '.sdrf.txt'), header = T, stringsAsFactors = F, sep = '\t')
View(phenoData)

# rownames(phenoData) <- gsub(' ', '.', phenoData$Source.Name)

#rownames(phenoData) <- gsub('_|-| ', '.', phenoData$Assay.Name)


idx <- c(1, grep('Characteristics', colnames(phenoData)))
idx

phenoData <- phenoData[,idx]
colnames(phenoData) <- c('Title','Organism','Disease','Source')

rownames(phenoData) <- phenoData$Title


platform <- 'A-MTAB-585'
platform


### annotation
cmd <- paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', platform, '/', platform, '.adf.txt -P data/fromArrayExpress/', gse, '/')
system(cmd)

files <- list.files(paste0('data/fromArrayExpress/', gse))
files

fl <- files[grep('A-MTAB|A-GEOD', files)]
fl

#annoData <- readLines(paste0('data/fromArrayExpress/', gse, '/', fl))
#annoData
#View(annoData)

annoData <- read.table(paste0('data/fromArrayExpress/', gse, '/', fl), header = T, sep = '\t', skip = 19, stringsAsFactors = F)
annoData
View(annoData)

idx <- grep('hsa', annoData$X.2)
idx

annoData <- annoData[idx,]

idx <- grep('MIMAT', annoData)
idx
annoData <- annoData[idx]

annoData <- lapply(annoData, function(x) strsplit(x, '\t')[[1]])
annoData <- do.call(rbind, annoData)

annoData <- data.frame(annoData, stringsAsFactors = F)
rownames(annoData) <- annoData$X4

annoData$Name <- annoData$X3


saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


#### 

### Expression
cmd <- paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', gse, '/', gse, '.processed.1.zip -P data/fromArrayExpress/', gse, '/')
system(cmd)
cmd <- paste0('unzip data/fromArrayExpress/', gse, '/', gse, '.processed.1.zip -d data/fromArrayExpress/', gse, '/')
system(cmd)

files <- list.files(paste0('data/fromArrayExpress/', gse))
files

fl <- files[which(!grepl('MTAB|GEOD', files))]
fl

exprData <- read.table(paste0('data/fromArrayExpress/', gse, '/', fl), header = T, stringsAsFactors = F, sep = '\t')
View(exprData)
dim(exprData)
dim(phenoData)

idx <- grep('hsa', exprData$Hybdridization.REF)
exprData <- exprData[idx,]

rownames(exprData) <- exprData$Hybdridization.REF
exprData <- exprData[,-1]

samples <- intersect(rownames(phenoData), colnames(exprData))
samples

exprData <- exprData[,samples]
phenoData <- phenoData[samples,]

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))

#exprData <- data.frame(t(exprData), stringsAsFactors = F)

#rownames(exprData) <- gsub('.', '-', rownames(exprData), fixed = T)

rownames(phenoData) == colnames(exprData)

#exprData <- exprData[,rownames(phenoData)]

mir.name <- sapply(rownames(exprData), function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 20)[[1]][1])
mir.name

sum(mir.name=='NA')

rownames(exprData)[which(mir.name=='NA')]

idx <- which(mir.name=='NA')

mir.name[idx] <- sapply(rownames(exprData)[idx], function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 19)[[1]][1])

idx <- which(mir.name=='NA')
idx

rownames(exprData)[idx]
mir.name[idx] <- c(mirIDNameConversionFun(id = 'hsa-miR-1274a', to = 'name2id', version = 16),
                   mirIDNameConversionFun(id = 'hsa-miR-1274b', to = 'name2id', version = 16),
                   mirIDNameConversionFun(id = 'hsa-miR-1280', to = 'name2id', version = 18),
                   mirIDNameConversionFun(id = 'hsa-miR-3647-3p', to = 'name2id', version = 17),
                   mirIDNameConversionFun(id = 'hsa-miR-3647-5p', to = 'name2id', version = 17),
                   mirIDNameConversionFun(id = 'hsa-miR-517b', to = 'name2id', version = 17),
                   mirIDNameConversionFun(id = 'hsa-miR-720', to = 'name2id', version = 18))

View(data.frame(mir.name, rownames(exprData)))

rownames(exprData) <- mir.name


saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))





######################################################################################
### Agilent-026094 GT Human miRNA Rel14 v1.1

# Release 15, some v11.0, 9.1

# GSE42072

gse <- 'GSE42072'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix
seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

#annoData$ID <- gsub('_st', '', annoData$ID)
#annoData$ID <- gsub('-star', '*', annoData$ID)

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


saveRDS(annoData, file='data/Annotation/Affymetrix_Multispecies_miRNA-1_Array.RDS')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))





######################################################################################
### Agilent-046064 Human miRNA Microarray, Release 19.0, 8x60K (miRNA_ID version)

# Release 15, some v11.0, 9.1

# GSE112840

gse <- 'GSE112840'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix
seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('MIMAT', annoData$miRNA_ID)
annoData <- annoData[idx,]


#annoData$ID <- gsub('_st', '', annoData$ID)
#annoData$ID <- gsub('-star', '*', annoData$ID)

annoData$Name <- annoData$miRNA_ID

saveRDS(annoData, file='data/Annotation/Affymetrix_Multispecies_miRNA-1_Array.RDS')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))




######################################################################################
### Illumina HiSeq 2000 (Homo sapiens)

# Release 15, some v11.0, 9.1

# GSE130654
# GSE122488
# GSE58410 # release 19
# GSE94564
# GSE104251


gse <- 'GSE104251'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix
seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130654/suppl/GSE130654_Readcount_TPM.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE130654_Readcount_TPM.txt.gz')

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122488/suppl/GSE122488_normalized_microRNA_counts.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE122488_normalized_microRNA_counts.txt.gz')

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58410/suppl/GSE58410_predicted_and_know_miRNA_species_identified_in_this_study.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE58410_predicted_and_know_miRNA_species_identified_in_this_study.txt.gz')

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104251/suppl/GSE104251_miRNA_Expression_Normalized.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE104251_miRNA_Expression_Normalized.txt.gz')



exprData <- read.table('data/fromGEO/GSE130654_Readcount_TPM.txt', header = T, sep = '\t', stringsAsFactors = F)
exprData

rownames(exprData) <- exprData$sRNA.readcount

idx <- grep('tpm', colnames(exprData))
exprData <- exprData[,idx]


exprData <- read.table('data/fromGEO/GSE122488_normalized_microRNA_counts.txt', header = T, sep = '\t', stringsAsFactors = F, skip = 3)
exprData

rownames(exprData) <- exprData$miRNA
exprData <- exprData[,-1]

exprData <- read.table('data/fromGEO/GSE58410_predicted_and_know_miRNA_species_identified_in_this_study.txt', 
                       header = T, sep = '\t', stringsAsFactors = F, skip = 1)
exprData

rownames(exprData) <- exprData$X
exprData <- exprData[,-c(1:3)]

exprData <- read.table('data/fromGEO/GSE104251_miRNA_Expression_Normalized.txt', 
                       header = T, sep = '\t', stringsAsFactors = F, skip = 1)
exprData




idx <- grep('hsa', rownames(exprData))
exprData <- exprData[idx,]

mir.name <- sapply(rownames(exprData), function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 19)[[1]][1])
mir.name
sum(mir.name=='NA')
idx <- which(mir.name=='NA')
rownames(exprData)[idx]

mir.name[idx] <- sapply(rownames(exprData)[idx], function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 22)[[1]][1])

rownames(exprData) <- mir.name

colnames(exprData)
rownames(phenoData)

phenoData$Title <- gsub('-', '.', phenoData$Title)

idx <- match(phenoData$Title, colnames(exprData))
idx

exprData <- exprData[,idx]

colnames(exprData) <- rownames(phenoData)


# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))





######################################################################################
### Illumina NextSeq 500 (Homo sapiens)

# Release 15, some v11.0, 9.1

# GSE94533

gse <- 'GSE94533'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix
seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94533/suppl/GSE94533_Processed_file_Cohort1.xlsx.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE94533_Processed_file_Cohort1.xlsx.gz')


exprData <- read_excel('data/fromGEO/GSE94533_Processed_file_Cohort1.xlsx')
exprData

rownames(exprData) <- exprData$sRNA.readcount

idx <- grep('tpm', colnames(exprData))
exprData <- exprData[,idx]


exprData <- read.table('data/fromGEO/GSE122488_normalized_microRNA_counts.txt', header = T, sep = '\t', stringsAsFactors = F, skip = 3)
exprData

rownames(exprData) <- exprData$miRNA
exprData <- exprData[,-1]

exprData <- read.table('data/fromGEO/GSE58410_predicted_and_know_miRNA_species_identified_in_this_study.txt', 
                       header = T, sep = '\t', stringsAsFactors = F, skip = 1)
exprData

rownames(exprData) <- exprData$X
exprData <- exprData[,-c(1:3)]

exprData <- read.table('data/fromGEO/GSE104251_miRNA_Expression_Normalized.txt', 
                       header = T, sep = '\t', stringsAsFactors = F, skip = 1)
exprData




idx <- grep('hsa', rownames(exprData))
exprData <- exprData[idx,]

mir.name <- sapply(rownames(exprData), function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 19)[[1]][1])
mir.name
sum(mir.name=='NA')
idx <- which(mir.name=='NA')
rownames(exprData)[idx]

mir.name[idx] <- sapply(rownames(exprData)[idx], function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 22)[[1]][1])

rownames(exprData) <- mir.name

colnames(exprData)
rownames(phenoData)

phenoData$Title <- gsub('-', '.', phenoData$Title)

idx <- match(phenoData$Title, colnames(exprData))
idx

exprData <- exprData[,idx]

colnames(exprData) <- rownames(phenoData)


# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))






######################################################################################3
### Agilent-070156 Human miRNA
# Release 21

# E-MTAB-8026

gse <- 'E-MTAB-8026'

system('mkdir data/fromArrayExpress/E-MTAB-8026')

### Phenotype

system('wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8026/E-MTAB-8026.sdrf.txt -P data/fromArrayExpress/E-MTAB-8026/')

phenoData <- read.table('data/fromArrayExpress/E-MTAB-8026/E-MTAB-8026.sdrf.txt', header = T, stringsAsFactors = F, sep = '\t')
View(phenoData)

rownames(phenoData) <- gsub(' ', '.', phenoData$Source.Name)
phenoData <- phenoData[,1:5]
colnames(phenoData) <- c('Title','Organism','Source','Individual','Disease')


platform <- '	A-GEOD-20712'
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
system('wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8026/E-MTAB-8026.processed.1.zip -P data/fromArrayExpress/E-MTAB-8026/')
system('unzip data/fromArrayExpress/E-MTAB-8026/E-MTAB-8026.processed.1.zip -d data/fromArrayExpress/E-MTAB-8026/')

exprData <- read.table('data/fromArrayExpress/E-MTAB-8026/expression_matrix.csv', header = T, stringsAsFactors = F, sep = '\t')
View(exprData)

rownames(exprData) <- exprData$miRNA
exprData <- exprData[,-1]
exprData[1:5,1:5]

rownames(phenoData) == colnames(exprData)

mir.name <- sapply(rownames(exprData), function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 21)[[1]][1])
View(mir.name)
rownames(exprData) <- mir.name

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))















######################################################################################
### GenoExploxer Homo sapiens 1.4K miRBASE version 13.0

# Release 13.0

# GSE25019

gse <- 'GSE25019'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

filter <- grep('Pre', annoData$miRNA_ID, ignore.case = T)
annoData <- annoData[-filter,]

annoData$miRNA_ID <- gsub('-0| ', '-', annoData$miRNA_ID)
annoData$miRNA_ID <- gsub('mir', 'miR', annoData$miRNA_ID)

annoData$Name22 <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '22')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')


#saveRDS(annoData, file='data/Annotation/Exiqon_LNA_RT-PCR_Human_panels_1_2.RDS')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

rownames(exprData) <- annoData$Name
View(exprData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))







######################################################################################
### Affymetrix GeneChip® miRNA 4.0 Array [miRNA-4_0]

# E-MTAB-3888
# E-MTAB-6652


gse <- 'E-MTAB-6652'

system(paste0('mkdir data/fromArrayExpress/', gse))

### Phenotype

system(paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', gse, '/', gse, '.sdrf.txt -P data/fromArrayExpress/', gse, '/'))

phenoData <- read.table(paste0('data/fromArrayExpress/', gse, '/', gse, '.sdrf.txt'), header = T, stringsAsFactors = F, sep = '\t')
View(phenoData)

# rownames(phenoData) <- gsub(' ', '.', phenoData$Source.Name)

rownames(phenoData) <- gsub('_|-| ', '.', phenoData$Assay.Name)


idx <- c(1, grep('Characteristics', colnames(phenoData)))
idx

phenoData <- phenoData[,idx]
colnames(phenoData) <- c('Title','Organism','Source','Disease')


platform <- 'A-MTAB-537'
platform

cmd <- paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', platform, '/', platform, '.adf.txt -P data/fromArrayExpress/', gse, '/')
system(cmd)

files <- list.files(paste0('data/fromArrayExpress/', gse))
files

fl <- files[grep('A-MTAB', files)]
fl

annoData <- readLines(paste0('data/fromArrayExpress/', gse, '/', fl))

View(annoData)

idx <- grep('hsa', annoData)
idx

annoData <- annoData[idx]

idx <- grep('MIMAT', annoData)
idx
annoData <- annoData[idx]

annoData <- lapply(annoData, function(x) strsplit(x, '\t')[[1]])
annoData <- do.call(rbind, annoData)

annoData <- data.frame(annoData, stringsAsFactors = F)
rownames(annoData) <- annoData$X4

annoData$Name <- annoData$X3


saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


#### 

### Expression
cmd <- paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', gse, '/', gse, '.processed.1.zip -P data/fromArrayExpress/', gse, '/')
system(cmd)
cmd <- paste0('unzip data/fromArrayExpress/', gse, '/', gse, '.processed.1.zip -d data/fromArrayExpress/', gse, '/')
system(cmd)

files <- list.files(paste0('data/fromArrayExpress/', gse))
files

fl <- files[which(!grepl('MTAB', files))]
fl


exprData <- read.table(paste0('data/fromArrayExpress/', gse, '/', fl), header = T, stringsAsFactors = F, sep = '\t')
View(exprData)

idx <- which(exprData$Hybridization.REF %in% rownames(annoData))
idx

exprData <- exprData[idx,]

rownames(exprData) <- exprData$Hybridization.REF
exprData <- exprData[rownames(annoData),]

rownames(exprData) <- annoData$Name

exprData <- exprData[,-c(1, ncol(exprData))]
colnames(exprData) <- gsub('.CEL.pimg', '', colnames(exprData))
colnames(exprData)[4] <- 'Synovial.Sarcoma.cell.free.serum.4'

rownames(phenoData) %in% colnames(exprData)
rownames(phenoData)
colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))





######################################################################################
### Exiqon miRCURY LNA™ Universal RT miRNA PCR Human Serum/Plasma focus panel (V3.R)
# Release 19


# E-MTAB-6304

gse <- 'E-MTAB-6304'

system(paste0('mkdir data/fromArrayExpress/', gse))

### Phenotype

system(paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', gse, '/', gse, '.sdrf.txt -P data/fromArrayExpress/', gse, '/'))

phenoData <- read.table(paste0('data/fromArrayExpress/', gse, '/', gse, '.sdrf.txt'), header = T, stringsAsFactors = F, sep = '\t')
View(phenoData)

# rownames(phenoData) <- gsub(' ', '.', phenoData$Source.Name)

rownames(phenoData) <- gsub('_|-| ', '.', phenoData$Assay.Name)


idx <- c(1, grep('Characteristics', colnames(phenoData)))
idx

phenoData <- phenoData[,idx]
colnames(phenoData) <- c('Title','Organism','Disease','Grade','Age','Sex','Source')

rownames(phenoData) <- phenoData$Title

platform <- 'A-GEOD-18693'
platform

cmd <- paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', platform, '/', platform, '.adf.txt -P data/fromArrayExpress/', gse, '/')
system(cmd)

files <- list.files(paste0('data/fromArrayExpress/', gse))
files

fl <- files[grep('A-MTAB|A-GEOD', files)]
fl

annoData <- readLines(paste0('data/fromArrayExpress/', gse, '/', fl))

View(annoData)

idx <- grep('hsa', annoData)
idx

annoData <- annoData[idx]

idx <- grep('MIMAT', annoData)
idx
annoData <- annoData[idx]

annoData <- lapply(annoData, function(x) strsplit(x, '\t')[[1]])
annoData <- do.call(rbind, annoData)

annoData <- data.frame(annoData, stringsAsFactors = F)
rownames(annoData) <- annoData$X4

annoData$Name <- annoData$X3


saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


#### 

### Expression
cmd <- paste0('wget https://www.ebi.ac.uk/arrayexpress/files/', gse, '/', gse, '.processed.1.zip -P data/fromArrayExpress/', gse, '/')
system(cmd)
cmd <- paste0('unzip data/fromArrayExpress/', gse, '/', gse, '.processed.1.zip -d data/fromArrayExpress/', gse, '/')
system(cmd)

files <- list.files(paste0('data/fromArrayExpress/', gse))
files

fl <- files[which(!grepl('MTAB|GEOD', files))]
fl

exprData <- read.table(paste0('data/fromArrayExpress/', gse, '/', fl), header = T, stringsAsFactors = F, sep = '\t')
View(exprData)

rownames(exprData) <- exprData$X
exprData <- exprData[,-1]

exprData <- data.frame(t(exprData), stringsAsFactors = F)

rownames(exprData) <- gsub('.', '-', rownames(exprData), fixed = T)

rownames(phenoData) == colnames(exprData)

exprData <- exprData[,rownames(phenoData)]

mir.name <- sapply(rownames(exprData), function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 19)[[1]][1])
mir.name

sum(mir.name=='NA')

rownames(exprData) <- mir.name


saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))






