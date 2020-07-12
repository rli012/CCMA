
setwd('C:\\Users/rli3/Documents/CCMA/')

ccma.primary <- read.delim('shinyApp/data/CCMA_Datasets_Primary.txt', header = T, sep = '\t', stringsAsFactors = F)
View(ccma.primary)

saveRDS(ccma.primary, file='shinyApp/data/CCMA_Datasets_Primary.RDS')

ccma <- read.delim('shinyApp/data/CCMA_Datasets.txt', header = T, sep = '\t', stringsAsFactors = F)
View(ccma)

saveRDS(ccma, file='shinyApp/data/CCMA_Datasets.RDS')


###
ccma <- readRDS(file='shinyApp/data/CCMA_Datasets.RDS')
colnames(ccma)[7] <- 'N'

saveRDS(ccma, file='shinyApp/data/CCMA_Datasets.RDS')


ccma.primary <- readRDS(file='shinyApp/data/CCMA_Datasets_Primary.RDS')
colnames(ccma.primary)

colnames(ccma.primary)[5] <- 'N'

saveRDS(ccma.primary, file='shinyApp/data/CCMA_Datasets_Primary.RDS')



ccma <- readRDS(file='shinyApp/data/CCMA_Datasets.RDS')
ccma$N[ccma$Dataset=='GSE100508'] <- 54

saveRDS(ccma, file='shinyApp/data/CCMA_Datasets.RDS')


ccma.primary <- readRDS(file='shinyApp/data/CCMA_Datasets_Primary.RDS')
colnames(ccma.primary)

ccma$N[ccma$Dataset=='GSE100508'] <- 54

saveRDS(ccma.primary, file='shinyApp/data/CCMA_Datasets_Primary.RDS')




meta.ccma <- list()
expr.ccma <- list()

for (i in 1:nrow(ccma)) {
  
  accession <- ccma$Accession[i]
  platform <- ccma$Annotation[i]
  dataset <- ccma$Dataset[i]
  
  print (i)
  print (dataset)
  
  exprData <- readRDS(file=paste0('shinyApp/data/Datasets/', accession, '_', platform, '_Expression.RDS'))
  phenoData <- readRDS(file=paste0('shinyApp/data/Datasets/', accession, '_', platform, '_Metadata.RDS'))
  
  expr.ccma[[dataset]] <- exprData
  meta.ccma[[dataset]] <- phenoData
  
  
}

saveRDS(expr.ccma, file='shinyApp/data/CCMA_Expression.RDS')
saveRDS(meta.ccma, file='shinyApp/data/CCMA_Metadata.RDS')


ccma$Dataset[1:2]

meta.ccma[[ccma$Dataset[1]]]


#################

ccma <- readRDS(file='shinyApp/data/CCMA_Datasets.RDS')
idx

i <- 66
ccma$Accession[i]
ccma$Annotation[i]
expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
               ccma$Annotation[i], '_Expression.RDS'))
genes <- rownames(expr)
genes

expr[1:5,1:5]

min(expr, na.rm = T)
max(expr, na.rm = T)


expr[expr<=0.1] <- 0.1

expr <- apply(expr, 2, as.numeric)

expr <- apply(expr, 2, function(v) round(v,3))
expr[1:5,1:5]

expr <- data.frame(expr, stringsAsFactors = F, row.names = genes)

expr <- apply(expr, 2, function(v) round(log2(v),3))
expr[1:5,1:5]

filter <- which(rowSums(is.na(expr))>0.2*ncol(expr))
filter
length(filter)
nrow(expr)

expr <- expr[-filter,]

saveRDS(expr, paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
               ccma$Annotation[i], '_Expression.RDS'))


idx <- c()

for (i in 1:nrow(ccma)) {
  ccma$Accession[i]
  ccma$Annotation[i]
  expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                         ccma$Annotation[i], '_Expression.RDS'))
  
  filter <- which(rowSums(is.na(expr))>0.2*ncol(expr))
  
  if (length(filter) > 0) {
    idx <- c(idx, i)
    print (i)
  }
  
}



i <- 38
ccma$Accession[i]
ccma$Annotation[i]
expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                       ccma$Annotation[i], '_Expression.RDS'))

filter <- which(colSums(is.na(expr))>0.5*nrow(expr))
filter

expr <- expr[,-filter]

pheno <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                        ccma$Annotation[i], '_Metadata.RDS'))
pheno

pheno <- pheno[-filter,]


saveRDS(pheno, paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                     ccma$Annotation[i], '_Metadata.RDS'))

nrow(pheno)




min.val <- c()
max.val <- c()
na.val <- c()
char.val <- c()

for (i in 1:nrow(ccma)) {
  expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                         ccma$Annotation[i], '_Expression.RDS'))
  
  char.val <- c(char.val, sum(!is.numeric(expr)))
  expr <- apply(expr, 2, as.numeric)
  
  min.val <- c(min.val, min(expr, na.rm = T))
  max.val <- c(max.val, max(expr, na.rm = T))
  na.val <- c(na.val, sum(is.na(expr)))
  
}

which(max.val>100)
which(na.val>0)

which(char.val>0)





for (i in 1:nrow(ccma)) {
  message(i)
  ccma$Accession[i]
  ccma$Annotation[i]
  expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                         ccma$Annotation[i], '_Expression.RDS'))
  
  print (nrow(expr))
  idx <- which(rowSums(is.na(expr))>0)
  print (length(idx))
  
  
}
