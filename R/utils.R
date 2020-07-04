#' Scale vector for fitting
#' @param df_vec is df$'key' values
#' @noRd 

scale_input <- function(df_vec)
{
	df_vec = log(df_vec)
	df_vec[is.infinite(df_vec)] = 0
	df_vec = (df_vec - mean(df_vec))/sd(df_vec)

	return(df_vec)
}

#' B-spline fitting of vector
#' @param data_vec is input vector for scale
#' @importFrom splines bs
#' @noRd 

b_spline_fit <- function(data_vec, n_degree=3)
{
	nknots = round((length(data_vec)^0.2)*2)	
	probseq = seq(0, 1, length=nknots+2)
	knots_vec = quantile(data_vec, probs=probseq[-c(1,nknots+2)])
	basis_vec = bs(data_vec, degree=n_degree, knots=knots_vec, intercept=FALSE)

	return(basis_vec)
}


#' Save normalized df to chromsome-wide split files
#'
#' Takes normalized df as input
#' @importFrom stringr str_split_fixed
#' @export

saveEachChr <- function(df, savepath, savename)
{
	for (chr in unique(str_split_fixed(df$frag1,'\\.',3)[,1]))
	{
		df2 <- df[which(str_split_fixed(df$frag1,'\\.',3)[,1] == chr),]	
		write.table(df2, file=gzfile(paste(savepath, '/', savename, '_', chr, '_normalized.gz', sep='')), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
	}
}
