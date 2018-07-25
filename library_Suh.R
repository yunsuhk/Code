##file name contents as a vector in working directory
list.files()

## read CSV with string data (default <-  as.is = FALSE for factor)
read.csv("data.csv", as.is = TRUE)

## write CSV without row number (default <- row.names = TRUE for row number)
write.csv(data, "data.csv", row.names = FALSE)

## rename of column
install.packages("dplyr")
drug <-  dplyr::rename(df, new_colName=old_colName)

#character in col to factor, levels= possible entire levels of data itself, labels = another optional vector for levels.
df$col <- factor (df$col, levels=c("low", "med", "high"))

#remove row with a certin condition of row
df <- subset(df, subset=(colA > 0)) #select rows only with colA>0

#remove list elements using condition
list[lapply(list, function(x){
  expression... #expression with condition for removal.
})] <- NULL

#reshape
drug_summary <- dplyr::select (df, x,y,z) #select x,y,z columns in df

##chi square, Fisher exact test

chisq.test(df$x, df$y)
fisher.test(df$x, df$y, alternative = "two.sided")


## Merging read count file (tumor, normal) by HTseq.

rm(list=ls())

#set working directory.
setwd("~/Dropbox/AEG classification_WES_JAX_2018 0118/data/readCount/readCountGEJ")

#retreive sample information including sampleID".
GEJdata=read.csv("sampleData_20180427.csv", header=T,as.is=T)
dim(GEJdata)

#Loop for merging readcount of each tumor and normal sample.
readTumor=c()
readNorm=c()

for(i in 1:nrow(GEJdata)){
  sampleID=GEJdata[i,1]
  filenameT=list.files(path=getwd(), pattern=paste(sampleID, "T", sep=""))
  filenameN=list.files(path=getwd(), pattern=paste(sampleID, "N", sep=""))
  print(c(filenameT, filenameN))
  readT <- read.delim(filenameT, sep="\t", header=T)
  readN <- read.delim(filenameN, sep="\t", header=T)
  readTumor <- cbind(readTumor, readT[,2])
  readNorm <- cbind(readNorm, readN[,2])
}
dim(readTumor)
dim(readNorm)

#annotation of geneID and sample name.
rownames(readTumor) <- readT[,1]
colnames(readTumor) <- paste(GEJdata[,1],"_T", sep="")
rownames(readNorm) <- readN[,1]
colnames(readNorm) <- paste(GEJdata[,1],"_N", sep="")

#merging total readcount file.
readTotal <- cbind(readTumor, readNorm)
dim(readTotal)
head(readTotal)
write.csv(readTotal, file="readTotal.csv")

#convert geneID (Ensembl) to gene symbol

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library("biomaRt")

listEnsemblArchives() ## GRCh37.74 http://Dec2013.archive.ensembl.org

# listMarts(host='Dec2013.archive.ensembl.org')  ## identify GRCh37 v.74

ensemb37 <- useMart(host='Dec2013.archive.ensembl.org', 
                    biomart='ENSEMBL_MART_ENSEMBL', 
                    dataset='hsapiens_gene_ensembl')

# ensemb37(df) = useDataset("hsapiens_gene_ensembl",mart=ensembl) ## for GRCh38.

#filters37 = listFilters(ensemb37)  ## identify filter name
#attribute37 = listAttributes(ensemb37)  ## identify attribute name

gene <- rownames(readTotal)
genelist37 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=gene,mart= ensemb37)
length(gene) # [1] 63681
dim(genelist37) #[1] 63761     2

##match HGNC ID to row name.
readTemp <- merge(readTotal, genelist37, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T) #possible duplication of gene name due to isoform.
readDup <- subset(readTemp, hgnc_symbol != "")  ## remove rows with blank gene symbol.

##remove duplication of isoform at readcount.
exp <- readDup[,- c(1,ncol(readDup))]  ## remove column with geneID and symbol for calculation.
row.mean <- apply(exp,1,mean)  ## calcuate row-wise mean value.
o = order(row.mean,decreasing=T)  ## order with decreasing tendency of mean value
readDup.ord <- readDup[o,]  ##  order with decreasing tendency of mean value
readDup.uniq <- readDup.ord[!duplicated(readDup.ord$hgnc_symbol),] ## !duplicated returns logical vector; False from 2nd (1st is true). remove all duplications except for 1st gene with the largest mean expression.


#ldply <- Split list, apply function, and return results in a data frame.For each element of a list, apply function then combine results into a data frame.

data <- ldply (fs, function(fn){
  v <- read.delim(fn, sep="\t", header=F) #”\t” is regular expression for tab
  v <- dplyr::rename(v, probe=V1, count=V2)
  id <- gsub(".sorted.dp_htseqcount.txt", "", fn) #a kind of grep,  gsub replaces all occurences, sub replaces only the first occuarence of a pattern.
  v <- dplyr::mutate(v, ID=id)
})


# mrnaIDs 를 세로로 transplose 시키되, "-"를 기준으로substring으로 나누고 4번째값을 분해한 후 1번째에서 2번째 값만 취하여 list를 만든다
samp <- lapply (as.list(t(mrnaIDs)), function(t){
  substr(unlist(strsplit(t, "-"))[4], 1, 2)}) 
  
# samp 가 10미만이면 cancer 아니면 normal.
lapply(samp,function(t) {if(t<10) return("1") else return ("0")})



#BSS/WSS

bssWssFast <- function (X, givenClassArr, numClass=2)
  # between squares / within square feature selection
{
  classVec <- matrix(0, numClass, length(givenClassArr))
  for (k in 1:numClass) {
    temp <- rep(0, length(givenClassArr))
    temp[givenClassArr == (k - 1)] <- 1
    classVec[k, ] <- temp
  }
  classMeanArr <- rep(0, numClass)
  ratio <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    overallMean <- sum(X[, j]) / length(X[, j])
    for (k in 1:numClass) {
      classMeanArr[k] <- 
        sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
    }
    classMeanVec <- classMeanArr[givenClassArr + 1]
    bss <- sum((classMeanVec - overallMean)^2)
    wss <- sum((X[, j] - classMeanVec)^2)
    ratio[j] <- bss/wss
  }
  sort(ratio, decreasing = TRUE, index = TRUE)
}

#mrnaData: gene name 이 column으로 된 df. mrnaClassNum: mrnaID가 row인df.
bss <- bssWssFast(mrnaData, t(mrnaClassNum), 2)

#1-100위까지의 DEG만 선택하여 샘플-DEG 간의 df생성.
mrnaDataReduced <- mrnaData[,bss$ix[1:100]]

