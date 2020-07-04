#' Quality control of coverage cormalization result
#'
#' Takes coverage normalized result and check PCC between ligation frequency and coverage values
#' @export

checkFreqCovPCC <- function(df, outpdfname='NULL')
{
	pcc_before_covnorm_1 <- cor(df$cov_frag1, df$freq)	
	pcc_before_covnorm_2 <- cor(df$cov_frag2, df$freq)
	print(paste("PCC before coverage normalization: ", round(pcc_before_covnorm_1,4), '/', round(pcc_before_covnorm_2,4), sep=' '))
	pcc_after_distnorm_1 <- cor(df$cov_frag1, df$capture_res)
	pcc_after_distnorm_2 <- cor(df$cov_frag2, df$capture_res)
	print(paste("PCC after coverage normalization: ", round(pcc_after_distnorm_1,4), '/',  round(pcc_after_distnorm_2,4), sep=' '))
	
	if(outpdfname != 'NULL')
	{
		if(!requireNamespace("corrplot", quietly = TRUE)){stop("The package 'corrplot' not be found. Please install it.")}
		pdf(outpdfname)
		M <- cor(data.frame(df$freq, df$cov_frag1, df$cov_frag2, df$capture_res))
		corrplot::corrplot(M, method="ellipse")
		dev.off()
	}
}

#' Plot coverage-interaction frequency before/expected/after normalization
#'
#' Takes coverage normalized result and plot 3 heatmaps
#' @export

plotCovNormRes <- function(df, binsize=100, pct_quantile=0.99, outpdfname="plot_coverage_norm_result.pdf")
{
	if(!requireNamespace("reshape2", quietly = TRUE)){stop("The package 'reshape2' not be found. Please install it.")}
	if(!requireNamespace("gplots", quietly = TRUE)){stop("The package 'gplots' not be found. Please install it.")}

	df2 <- df[c('cov_frag1', 'cov_frag2','freq')]

	cov_frag1_int <- as.integer(df2$cov_frag1/binsize)
	cov_frag2_int <- as.integer(df2$cov_frag2/binsize)

	freq <- df$freq
	exp <- df$exp_value_capture
	norm <- df$capture_res

	df3   <- data.frame(cov_frag1_int, cov_frag2_int, freq)
	df3_1 <- data.frame(cov_frag1_int, cov_frag2_int, exp )
	df3_2 <- data.frame(cov_frag1_int, cov_frag2_int, norm)

	df4   <- reshape2::acast(df3,   df3$cov_frag1_int~df3$cov_frag2_int,     fun.aggregate=sum, value.var='freq')
	df4_1 <- reshape2::acast(df3_1, df3_1$cov_frag1_int~df3_1$cov_frag2_int, fun.aggregate=sum, value.var='exp')
	df4_2 <- reshape2::acast(df3_2, df3_2$cov_frag1_int~df3_2$cov_frag2_int, fun.aggregate=sum, value.var='norm')

	pdf(outpdfname)

	col_bby=colorRampPalette(c('black','blue','yellow'))
	gplots::heatmap.2(df4, col = col_bby, trace='none', Colv=FALSE, Rowv=FALSE, labRow=TRUE, labCol=TRUE, dendrogram='none',density.info='none',breaks = seq(0, quantile(df4, pct_quantile), length.out =100), symbreaks = TRUE, symm=FALSE,symkey=FALSE, margins=c(10,10), key.title=NA, lhei=c(1,6))
	gplots::heatmap.2(df4_1, col = col_bby, trace='none', Colv=FALSE, Rowv=FALSE, labRow=TRUE, labCol=TRUE, dendrogram='none',density.info='none',breaks = seq(0, quantile(df4, pct_quantile), length.out =100), symbreaks = TRUE, symm=FALSE,symkey=FALSE, margins=c(10,10), key.title=NA, lhei=c(1,6))
	gplots::heatmap.2(df4_2, col = col_bby, trace='none', Colv=FALSE, Rowv=FALSE, labRow=TRUE, labCol=TRUE, dendrogram='none',density.info='none',breaks = seq(0, quantile(df4, pct_quantile), length.out =100), symbreaks = TRUE, symm=FALSE,symkey=FALSE, margins=c(10,10), key.title=NA, lhei=c(1,6))

	dev.off()
}
