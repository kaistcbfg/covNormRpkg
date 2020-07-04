#' Quality control of distance cormalization result
#'
#' Takes distance normalized result and check PCC between ligation frequency and distance values
#' @export

checkFreqDistPCC <- function(df, outpdfname='NULL')
{
	pcc_before_distnorm <- cor(df$dist, df$capture_res)	
	print(paste("PCC before distance normalization: ", round(pcc_before_distnorm,4), sep=''))
	pcc_after_distnorm <- cor(df$dist, df$dist_res)
	print(paste("PCC after distance normalization: ", round(pcc_after_distnorm,4), sep=''))

	if(outpdfname != 'NULL')
	{
		if(!requireNamespace("corrplot", quietly = TRUE)){stop("The package 'corrplot' not be found. Please install it.")}
		pdf(outpdfname)
		M <- cor(data.frame(df$capture_res, df$dist, df$dist_res))
		corrplot::corrplot(M, method="ellipse")
		dev.off()
	}	

}

#' Plot distance-interaction frequency before/after normalization
#'
#' Takes distance normalized result and plot two hexagonal heatmap
#' @export

plotDistNormRes <- function(df, maxX=2000000, maxY1=70, maxY2=15, outpdfname="plot_distance_norm_result.pdf")
{
	if(!requireNamespace("ggplot2", quietly = TRUE)){stop("The package 'ggplot2' not be found. Please install it.")}	

	pdf(outpdfname)

	g1=ggplot2::ggplot(df)
	g2=ggplot2::geom_hex(ggplot2::aes(df$dist, df$capture_res))
	g3=ggplot2::scale_fill_gradientn('Density', colours=c('black', 'brown','red','darkorange','yellow'))
	g4=ggplot2::coord_cartesian(xlim=c(0,maxX),ylim=c(0,maxY1))
	g5=ggplot2::geom_smooth(ggplot2::aes(df$dist, df$capture_res), span=0.001)
	g6=ggplot2::theme_classic()
	plot(g1+g2+g3+g4+g5+g6)

	g1=ggplot2::ggplot(df)
	g2=ggplot2::geom_hex(ggplot2::aes(df$dist, df$dist_res))
	g3=ggplot2::scale_fill_gradientn('Density', colours=c('black', 'brown','red','darkorange','yellow'))
	g4=ggplot2::coord_cartesian(xlim=c(0,maxX),ylim=c(0,maxY2))
	g5=ggplot2::geom_smooth(ggplot2::aes(df$dist, df$dist_res), span=0.001)
	g6=ggplot2::theme_classic()
	plot(g1+g2+g3+g4+g5+g6)

	dev.off()
}
