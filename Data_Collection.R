
setwd('~/bigdata/EarlyCancerDetection/CCMA/')

library(Biobase)
library(GEOquery)
library(readxl)
library(stringr)

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



######################################################################################3
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
### 3D-Gene Human miRNA V20_1.0.0

# GSE73002
# GSE59856
# GSE85679/GSE85680


gse <- 'GSE73002'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = F, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)
#View(phenoData)
phenoData$Data.Processing[1]
phenoData$Source

### Annotation
annoData <- seriesMatrix@featureData@data
dim(annoData)

View(annoData)

#annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(x)[[1]][1])

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData, na.rm = T)
max(exprData, na.rm = T)


rownames(phenoData) == colnames(exprData)

rownames(annoData) == rownames(exprData)

#View(annoData)

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
### 3D-Gene Human miRNA V20_1.0.0; 3D-Gene Human miRNA V21_1.0.0

# GSE124158
# GSE134108

gse <- 'GSE134108'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)
#View(phenoData)
phenoData$Data.Processing[1]
phenoData$Source

### Annotation
annoData <- seriesMatrix@featureData@data
dim(annoData)

View(annoData)

#annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(x)[[1]][1])

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData, na.rm = T)
max(exprData, na.rm = T)


rownames(phenoData) == colnames(exprData)

rownames(annoData) == rownames(exprData)

#View(annoData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))



######################################################################################3
### Agilent-031181 Unrestricted_Human_miRNA_V16.0_Microarray (miRBase release 16.0 miRNA ID version)

# GSE68951
# GSE48137
# GSE59993
# GSE68373
# GSE79943
# GSE54156

gse <- 'GSE42654'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

### Annotation
#annoData <- seriesMatrix@featureData@data

#View(annoData)

#annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 16)[[1]][1])
#annoData$Name

#saveRDS(annoData, file='data/Annotation/Agilent_031181_Unrestricted_Human_miRNA_V16.0_Microarray.RDS')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData)
max(exprData)

View(exprData)

idx <- grep('hsa', rownames(exprData))
idx

exprData <- exprData[idx,]

rownames(phenoData) == colnames(exprData)

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
### Agilent-021827 Human miRNA Microarray [miRNA_107_Sep09_2_105]

# Release 15, some v11.0, 9.1

# GSE43160

gse <- 'GSE43160'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

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





######################################################################################
### [miRNA-1] Affymetrix Multispecies miRNA-1 Array

# Release 15, some v11.0, 9.1

# GSE46729

gse <- 'GSE46729'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

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

annoData$ID <- gsub('_st', '', annoData$ID)
annoData$ID <- gsub('-star', '*', annoData$ID)

annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
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




######################################################################################3
### [miRNA-2] Affymetrix Multispecies miRNA-2 Array
# Release 15, some v11

# GSE44281
# GSE39845 - whole blood, tissue

gse <- 'GSE39845'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

###
keep <- grep('Whole blood', phenoData$Source)
keep

phenoData <- phenoData[keep,] 


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$ID)
annoData <- annoData[idx,]

filter <- grep('hp|v11', annoData$ID)
annoData <- annoData[-filter,]

annoData$ID <- gsub('_st', '', annoData$ID)
annoData$ID <- gsub('-star', '*', annoData$ID)

annoData$Name <- sapply(annoData$ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name

#idx <- grep('v11', annoData$ID)
#idx

#annoData$ID[idx] <- gsub('v11_', '', annoData$ID[idx])
#annoData$Name[idx] <- sapply(annoData$ID[idx], function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '11.0')[[1]][1])

saveRDS(annoData, file='data/Annotation/Affymetrix_Multispecies_miRNA_2_Array.RDS')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

idx <- grep('hsa', rownames(exprData))
idx

exprData <- exprData[idx,]

filter <- grep('hp|v11', rownames(exprData))
exprData <- exprData[-filter,]

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
### [miRNA-4] Affymetrix Multispecies miRNA-4 Array

# Release 20

# GSE85589

gse <- 'GSE85589'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data
annoData$Name <- annoData$Accession

saveRDS(annoData, file='data/Annotation/Affymetrix_Multispecies_miRNA-4_Array.RDS')

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
### [miRNA-4] Affymetrix Multispecies miRNA-4 Array [ProbeSet ID version]

# Release 20

# GSE98181

gse <- 'GSE98181'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

idx <- grep('hsa', annoData$`Transcript ID(Array Design)`)
annoData <- annoData[idx,]

keep <- which(annoData$`Sequence Type`=='miRNA')
annoData <- annoData[keep,]

annoData$Name <- annoData$Accession

# saveRDS(annoData, file='data/Annotation/Affymetrix_Multispecies_miRNA-4_Array.RDS')

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

### Annotation
annoData <- seriesMatrix@featureData@data
dim(annoData)

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = 21)[[1]][1])
annoData$Name

saveRDS(annoData, file='data/Annotation/Agilent_070156_Human_miRNA_V21.0_Microarray_046064.RDS')



platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]
min(exprData, na.rm = T)
max(exprData, na.rm = T)


rownames(phenoData) == colnames(exprData)

rownames(annoData) == rownames(exprData)

rownames(exprData) <- annoData$Name

#View(annoData)

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


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


#################################################################################
### Applied Biosystems TaqMan Array Human MicroRNA A+B Cards Set v3.0

# GSE47125
# GSE65708
# GSE50013


gse <- 'GSE50013'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

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

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix

annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47125/suppl/GSE47125.xls.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE47125.xls.gz')

system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65708/suppl/GSE65708_Matrix_normalized_by_MCR.txt.gz -P data/fromGEO/')
system('gunzip data/fromGEO/GSE65708_Matrix_normalized_by_MCR.txt.gz')

exprData <- read_excel('data/fromGEO/GSE47125.xls')

exprData <- read.delim('data/fromGEO/GSE65708_Matrix_normalized_by_MCR.txt', header = T, sep = '\t', stringsAsFactors = F)

exprData <- exprs(seriesMatrix)

exprData[1:5,1:5]
min(exprData)
max(exprData)

View(exprData)
dim(exprData)

exprData <- data.frame(exprData, stringsAsFactors = F)
rownames(exprData) <- exprData$ID_REF
exprData <- exprData[,-1]

idx <- which(!rownames(annoData) %in% rownames(exprData))
idx

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

exprData[filter,]

exprData <- exprData[-filter,]

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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID, '_', annoData$Panel)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

exprData <- exprs(seriesMatrix)

exprData[1:5,1:5]
min(exprData)
max(exprData)

View(exprData)
dim(exprData)

idx <- which(!rownames(annoData) %in% rownames(exprData))
idx
annoData[idx,]

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

exprData[filter,]

exprData <- exprData[-filter,]

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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

exprData <- exprs(seriesMatrix)

exprData[1:5,1:5]
min(exprData)
max(exprData)

View(exprData)
dim(exprData)

idx <- which(!rownames(annoData) %in% rownames(exprData))
idx
annoData[idx,]

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

exprData[filter,]

exprData <- exprData[-filter,]

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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- annoData$Assay.ID

rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

exprData <- exprs(seriesMatrix)

exprData[1:5,1:5]
min(exprData)
max(exprData)

View(exprData)
dim(exprData)

rownames(exprData) <- gsub('mir','miR',rownames(exprData))

rownames(exprData) <- str_extract(rownames(exprData), '\\d+$')

idx <- which(!rownames(annoData) %in% rownames(exprData))
idx
annoData[idx,]

filter <- which(!rownames(exprData) %in% rownames(annoData))
filter

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
### Applied Biosystems Taqman Low Density Array Human microRNA Card A, Applied Biosystems Taqman Low Density Array Human microRNA Card B

# GSE46355

gse <- 'GSE46355'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)
seriesMatrix

seriesMatrix <- seriesMatrix[[2]]
seriesMatrix

### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- str_extract(rownames(annoData), '\\d+$')

#rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

exprData <- exprs(seriesMatrix)

exprData[1:5,1:5]
min(exprData)
max(exprData)

View(exprData)
dim(exprData)

rownames(exprData) <- gsub('mir','miR',rownames(exprData))

rownames(exprData) <- str_extract(rownames(exprData), '\\d+$')

idx <- which(!rownames(annoData) %in% rownames(exprData))
idx
annoData[idx,]

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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

### Annotation
annoData <- readRDS(file='data/Annotation/Applied_Biosystems_TaqMan_Array_Human_MicroRNA_A_B_Cards_Set_v3.0.RDS')
annoData

rownames(annoData) <- annoData$Assay.ID

#rownames(annoData) <- paste0(annoData$Assay.Name, '-', annoData$Assay.ID)


platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression

exprData <- exprs(seriesMatrix)

exprData[1:5,1:5]
min(exprData)
max(exprData)

View(exprData)
dim(exprData)

rownames(exprData) <- gsub('mir','miR',rownames(exprData))

rownames(exprData) <- str_extract(rownames(exprData), '^\\d+')

idx <- which(!rownames(annoData) %in% rownames(exprData))
idx
annoData[idx,]

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
### Illumina Human v2 MicroRNA expression beadchip

# Release 15, some v11.0, 9.1

# GSE25609
# GSE41526


gse <- 'GSE41526'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

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

filter <- grep(':', annoData$miRNA_ID)
filter
annoData <- annoData[-filter,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


saveRDS(annoData, file='data/Annotation/Illumina_Human_v2_MicroRNA_expression_beadchip.RDS')

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
### Agilent-021827 Human miRNA Microarray G4470C (Feature Number version)

# Release 15, some v11.0, 9.1

# GSE39833/GSE40247
# GSE40246/GSE40247
# GSE55139


gse <- 'GSE39833'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

#idx <- grep('serum_exo', phenoData$Title)
#idx
#phenoData <- phenoData[idx,]


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$ACCESSION_STRING, function(x) strsplit(x, '|', fixed = T)[[1]][4])
annoData$Name

saveRDS(annoData, file='data/Annotation/Agilent-021827_Human_miRNA_Microarray_G4470C_Feature_Number_version.RDS')

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

exprData <- data.frame(exprData,stringsAsFactors = F)

exprData <- exprData[rownames(annoData),]
dim(exprData)

rownames(phenoData) == colnames(exprData)
exprData <- exprData[,rownames(phenoData)]

rownames(annoData) == rownames(exprData)
sum(rownames(annoData) == rownames(exprData))
dim(exprData)

exprData$ID <- annoData$Name





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
### Exiqon LNA RT-PCR Human panels (1 & 2)

# Release 15

# GSE42128 [3]
# GSE41922

gse <- 'GSE41922'

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

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


saveRDS(annoData, file='data/Annotation/Exiqon_LNA_RT-PCR_Human_panels_1_2.RDS')

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
### Exiqon Human miRCURY LNA Universal RT miRNA PCR Human Serum/Plasma focus panel (V3.0)

# Release 19

# GSE117063

gse <- 'GSE51410'

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

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '19')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


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

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$Human_miRNA)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '18')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


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

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source

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

exprData <- exprData[,-c(1:2)]
exprData <- exprData[,-c(1)]

samples <- gsub('circulating miRNA \\[|\\]', '', phenoData$Title)
samples
samples <- paste0('X', gsub('circulating miRNA \\[|\\]', '', phenoData$Title))
samples

idx <- match(samples, colnames(exprData))
idx


samples == colnames(exprData[,idx])

colnames(exprData)[idx] <- phenoData$Accession

rownames(exprData) <- mir.name
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
### State Key Laboratory Human microRNA array 1858

# Release 15

# GSE118613
# GSE93850


gse <- 'GSE93850'

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
annoData$Name <- annoData$Mature_Acc

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '11.0')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '21')[[1]][1])
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

View(phenoData)

phenoData$Data.Processing[1]
phenoData$Source


### Annotation
annoData <- seriesMatrix@featureData@data

View(annoData)
idx <- grep('hsa', annoData$miRNA_ID)
annoData <- annoData[idx,]

annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '15')[[1]][1])
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
### Serum_microRMA

# Release 15 (13.0 reported)

# GSE101841

gse <- 'GSE101841'

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


annoData$Name <- sapply(annoData$miRNA_ID, function(x) mirIDNameConversionFun(id = x, to = 'name2id', version = '18')[[1]][1])
annoData$Name
sum(annoData$Name=='NA')

filter <- which(annoData$Name=='NA')
annoData <- annoData[-filter,]


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


filter <- which(duplicated(rownames(exprData)))
filter

exprData <- exprData[-filter,]

#exprData <- data.frame(ID=annoData$miRNA,exprData, stringsAsFactors = F)

#View(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_', platform, '_Expression.RDS'))


### eSet
eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData),
                      #featureData = annoData,
                      annotation = platform)

saveRDS(eSet, file=paste0('data/rData/', gse, '_', platform, '_eSet.RDS'))
