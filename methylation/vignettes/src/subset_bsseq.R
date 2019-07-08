library(bsseq)

bs <- readRDS("bsseq.rds")
bs20 <- chrSelectBSseq(bs, "chr20")
sub <- !(grepl("Withdraw", pData(bs20)$dox) &
                 grepl("Methylated", pData(bs20)$condition))
bs20 <- bs20[,sub]

# reconstruct
M <- (realize(getCoverage(bs20, type="M")))
Cov <- (realize(getCoverage(bs20, type="Cov")))
BStmp <- BSseq(chr = as.character(seqnames(bs20)), pos = start(bs20), 
               M = M, Cov = Cov, sampleNames = pData(bs)$basefile[sub])
BStmp <- realize(BStmp)
pData(BStmp) <- pData(bs)[sub,c(1:14,51)]

pData(BStmp)$condition <- as.factor(pData(BStmp)$condition)
pData(BStmp)$Sample <- sampleNames(BStmp) <- pData(bs)$title[sub]

# save
saveRDS(BStmp, "bsseq_chr20.rds")
