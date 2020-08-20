predict.rnaseq.rsi.rank<-function(x) {
  
  coeffs<-c("ENSG00000115415"=0.0254552046046931,
            "ENSG00000116478" = -0.020446868981449,
            "ENSG00000177606" = 0.0128282697929755,
            "ENSG00000173039" = -0.00381713588908642,
            "ENSG00000097007" = 0.107021274358643,
            "ENSG00000125347" = -0.0441682719957224,
            "ENSG00000180370" = -0.00924314364533049,
            "ENSG00000166501" = -0.00175888592095915,
            "ENSG00000116030" = -0.00025093287356907,
            "ENSG00000169083" = -0.00980092741843779)
  
  # Get specific gene identifiers (may have versions)
  gene_ids<-sapply(names(coeffs),function(y){featureNames(ccle_rsem)[which(startsWith(featureNames(ccle_rsem),y))]})
  stopifnot(length(gene_ids)==10)
  subset.data <- x[gene_ids,]
  predicted<-rep(NA,ncol(subset.data))
  
  for(i in 1:ncol(subset.data)){
    temp.ranks <- rank(subset.data[,i])
    predicted[i] <- sum(temp.ranks*coeffs)
  }
  results <- data.frame("Patient.ID"=colnames(x),"RSI.Rank"=predicted)
  return(results)
}