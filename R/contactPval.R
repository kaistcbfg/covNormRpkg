#' Compute p-value/FDR of normalized interaction by 3P-Weibull fitting
#'
#' Takes formatted dataframe and retunrs dataframe with 3P Weibull fitted value
#' @importFrom propagate fitDistr
#' @importFrom FAdist pweibull3
#' @param df A coverage-distance normalized output df.
#' @param outpdfname fitDistr-generated fitting result pdf name.
#' @return df Added p-value/FDR of interactions as columns.
#' @export

contactPval <- function(df, outpdfname)
{
	dist_res <- as.numeric(df$dist_res)
	#dist_res_filter = dist_res
	#if(fit_thresh > 0){ dist_res_filter = dist_res[which(dist_res <= fit_thresh)] }
	dist_res_filter = dist_res[which(dist_res <= 2)]

	pdf(outpdfname)
	result=fitDistr(dist_res_filter)
	dev.off()
	location=result$fit$"3P Weibull"$par$location
	shape=result$fit$"3P Weibull"$par$shape
	scale=result$fit$"3P Weibull"$par$scale
	p_result_dist <- pweibull3(dist_res,shape,scale,location,lower.tail=FALSE,log.p=FALSE)
	FDR_dist_res <- p.adjust(p_result_dist,method="fdr",n=length(p_result_dist))
	final_data=cbind(df,p_result_dist,FDR_dist_res)

	return(final_data)
}

