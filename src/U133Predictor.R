# Probeset						GeneSymbol	RankCoefficient		ExpressionCoefficient
# 211110_s_at					AR		-0.00980092741843779	-0.0183234816987221
# 201466_s_at					CJUN	 0.0128282697929755		 0.0139806837676373
# AFFX-HUMISGF3A/M97935_MA_at	STAT1	 0.0254552046046931		 0.0425780903190549
# 207957_s_at					PRKCB	-0.00175888592095915	 0.0114458950449979
# 201783_s_at					RELA	-0.00381713588908642	-0.00167323516795081
# 202123_s_at					CABL	 0.107021274358643		 0.129226680793857
# 208762_at						SUMO1	-0.00025093287356907	-0.00502939666304479
# 205962_at						PAK2	-0.00924314364533049	-0.00852700671867646
# 201209_at						HDAC1	-0.020446868981449		-0.0181165060212908
# 202531_at						IRF1	-0.0441682719957224		-0.0927956729169349



runPredictor<-function() {
   x<-justRMA()
   write.exprs(x,file="exprs-rma.txt")
   calculateRSI()
}

calculateRSI<-function(file="exprs-rma.txt",outputfile="RSI-predictions.txt") {
	gene.data <- as.matrix(read.table(file,header=TRUE,row.names=1,sep="\t",comment.char="",quote=""))
	results.rank<-predict.u133.rsi.rank(gene.data)
	results.expr<-predict.u133.rsi.expression(gene.data)
	all.results<-data.frame("RSI.Rank"=results.rank$RSI.Rank,"RSI.Expression"=results.expr$RSI.Expression)
	rownames(all.results)<-sub(".[Cc][Ee][Ll]","",results.rank[,"Patient.ID"])
	write.csv(all.results,file=outputfile,quote=F)
	return(all.results)
}
predict.u133.rsi.rank<-function(x) {

	coeffs<-c("AFFX-HUMISGF3A/M97935_MA_at"=0.0254552046046931,
			"201209_at" = -0.020446868981449,
			"201466_s_at" = 0.0128282697929755,
			"201783_s_at" = -0.00381713588908642,
			"202123_s_at" = 0.107021274358643,
			"202531_at" = -0.0441682719957224,
			"205962_at" = -0.00924314364533049,
			"207957_s_at" = -0.00175888592095915,
			"208762_at" = -0.00025093287356907,
			"211110_s_at" = -0.00980092741843779)
	
	subset.data <- x[match(names(coeffs),rownames(x)),]
	predicted<-rep(NA,ncol(subset.data))

	for(i in 1:ncol(subset.data)){
		temp.ranks <- rank(subset.data[,i])
		predicted[i] <- sum(temp.ranks*coeffs)
	}
	results <- data.frame("Patient.ID"=colnames(x),"RSI.Rank"=predicted)
	return(results)
}

predict.u133.rsi.expression<-function(x) {
	coeffs<-c("211110_s_at" = -0.0183234816987221,
				"201466_s_at" = 0.0139806837676373,
				"AFFX-HUMISGF3A/M97935_MA_at" = 0.0425780903190549,
				"207957_s_at" = 0.0114458950449979,
				"201783_s_at" = -0.00167323516795081,
				"202123_s_at" =	0.129226680793857,
				"208762_at" = -0.00502939666304479,
				"205962_at" = -0.00852700671867646,
				"201209_at" = -0.0181165060212908,
				"202531_at" = -0.0927956729169349)

	subset.data <- x[match(names(coeffs),rownames(x)),]
	predicted<-rep(NA,ncol(subset.data))

	for(i in 1:ncol(subset.data)){
		predicted[i] <- sum(subset.data[,i]*coeffs)
	}
	results <- data.frame("Patient.ID"=colnames(x),"RSI.Expression"=predicted)
	return(results)
}
