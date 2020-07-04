#' Compute p-value/FDR of normalized interaction by 3P-Weibull fitting
#'
#' Takes formatted dataframe and retunrs dataframe with 3P Weibull fitted value
#' @importFrom propagate fitDistr
#' @importFrom FAdist pweibull3
#' @param df A coverage-distance normalized output df.
#' @return df Added p-value/FDR of interactions as columns.
#' @export

contactPval <- function(df)
{
	dist_res <- as.vector(df$dist_res)

	result=fitDistr(dist_res,plot=FALSE)
	location=result$fit$"3P Weibull"$par$location
	shape=result$fit$"3P Weibull"$par$shape
	scale=result$fit$"3P Weibull"$par$scale
	p_result_dist <- pweibull3(dist_res,shape,scale,location,lower.tail=FALSE,log.p=FALSE)
	FDR_dist_res <- p.adjust(p_result_dist,method="fdr",n=length(p_result_dist))
	final_data=cbind(df,p_result_dist,FDR_dist_res)

	return(final_data)
}

