if (!file.exists(file.path("../raw_data/counts_Calero_20160325.tsv"))){
  library(BiocFileCache)
  bfc <- BiocFileCache("../raw_data", ask = FALSE)
  lun.zip <- bfcrpath(bfc, 
                      file.path("https://www.ebi.ac.uk/arrayexpress/files",
                                "E-MTAB-5522/E-MTAB-5522.processed.1.zip"))
  lun.sdrf <- bfcrpath(bfc, 
                       file.path("https://www.ebi.ac.uk/arrayexpress/files",
                                 "E-MTAB-5522/E-MTAB-5522.sdrf.txt"))
  unzip(lun.zip, exdir = "../raw_data")
  
  metadata <- read.delim(lun.sdrf, check.names=FALSE, header=TRUE)
  saveRDS(metadata, file.path("../raw_data/416B_metadata.rds"))
}

plate1 <- read.delim(file.path("../raw_data/counts_Calero_20160113.tsv"), 
                     header=TRUE, row.names=1, check.names=FALSE)
plate2 <- read.delim(file.path("../raw_data/counts_Calero_20160325.tsv"), 
                     header=TRUE, row.names=1, check.names=FALSE)
gene.lengths <- plate1$Length # First column is the gene length.
plate1 <- as.matrix(plate1[,-1]) # Discarding gene length (as it is not a cell).
plate2 <- as.matrix(plate2[,-1])
rbind(Plate1=dim(plate1), Plate2=dim(plate2))
all.counts <- cbind(plate1, plate2)
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts))
rowData(sce)$GeneLength <- gene.lengths

isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
is.sirv <- grepl("^SIRV", rownames(sce))
sce <- sce[!is.sirv,] 

metadata <- readRDS(file.path("../raw_data/416B_metadata.rds"))
m <- match(colnames(sce), metadata[["Source Name"]]) # Enforcing identical order.
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]

colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "control")
colData(sce)$Oncogene <- pheno

library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb

library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)

library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
                   column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")

mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))

num.cells <- nexprs(sce, byrow=TRUE)
to.keep <- num.cells > 0
sce <- sce[to.keep,]
sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
