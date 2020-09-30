# SVA analysis and correction using the `sva` R package on DESeq2-variance stabilized data
# This implementation was benchmarked in Germain et al. 2020 ( https://doi.org/10.1186/s13059-020-02136-7 )
dosvacor <- function(SE, form=NULL, form0=~1, ...){
  library(sva)
  library(DESeq2)
  CD <- as.data.frame(colData(SE))
  mm <- model.matrix(form, data=CD)
  mm0 <- model.matrix(form0, data=CD)
  dds <- DESeqDataSetFromMatrix(round(assay(SE)), as.data.frame(colData(SE)), form)
  dds <- estimateSizeFactors(dds)
  en <- as.matrix(assay(vst(dds, blind=FALSE)))
  sv <- sva(en, mm, mm0, n.sv=NULL, ...)
  n.sv <- sv$n.sv
  sv <- sv$sv
  
  colnames(sv) <- paste0("SV",1:ncol(sv))
  X <- cbind(mm, sv)
  mm2 <- cbind(mm[,1,drop=F],sv,mm[,-1,drop=F])
  H <- solve(t(X)%*%X)%*%t(X)
  b <- (H%*%t(en))
  cn <- setdiff(colnames(X),setdiff(colnames(mm), colnames(mm0)))  
  cn <- setdiff(cn, "(Intercept)")
  encor <- en - t(as.matrix(X[,cn]) %*% b[cn,])
  SE <- SE[row.names(encor),]
  colData(SE) <- cbind(colData(SE), sv)
  assays(SE)$corrected <- encor
  return(SE)
}


printDEA <- function(x){
  if(is(x,"TopTags")) x <- x$table
  for(f in c("logCPM", grep("logFC",colnames(x), value=TRUE)))
    x[[f]] <- round(x[[f]],2)
  x$F <- NULL
  print(x, digits=2)
}
