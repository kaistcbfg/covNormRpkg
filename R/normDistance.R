#' Distance normalization function
#'
#' Takes formatted dataframe and retunrs glm.nb fitted coeffs & normalized df
#' @importFrom MASS glm.nb
#' @param df An input dataframe, coverage normalized.
#' @param max_dist If > 0, Rows with larger $dist values will be dropped.
#' @param max_distnorm_value Outliers by division with too small a number are replaced by this value.
#' @param sample_ratio If NOT -1 && between 0.0-1.0, downsample the df rows to be used for fitting by that ratio.
#' @return df Distance normalized.
#' @export

normDistance <- function(df, max_dist=-1, max_distnorm_value=50, sample_ratio=-1)
{
	if(max_dist > 0){df<-df[which(df$dist < max_dist),]}

	u_vec = df$capture_res

	dist_vec = scale_input(df$dist)

	u_vec_fit <- u_vec
	dist_vec_fit <- dist_vec

	if(sample_ratio > 0)
	{
		if(sample_ratio > 1.0){sample_ratio = 1.0}
		df_ds <- df[sample(nrow(df), size = as.integer(nrow(df) * sample_ratio)),]
		
		u_vec_ds = df_ds$freq

		dist_vec_ds = scale_input(df_ds$dist)

		u_vec_fit <- u_vec_ds
		dist_vec_fit <- dist_vec_ds
	}

	#dist_vec_fit <- b_spline_fit(dist_vec_fit)

	fit<-glm.nb(u_vec_fit ~ dist_vec_fit)
	coeff<-round(fit$coeff,4)

	intercept = coeff[1]
	coeff_dist= coeff[2]

	nknots = round((length(u_vec))^(0.2)*2)
	exp_value_dist=round(exp(intercept + coeff_dist*dist_vec),4)
	#exp_value_dist=round(exp(as.vector(intercept+dist_vec_fit%*%coeff[2:(nknots+4)])), 4)
	#dist_res=round((u_vec)/(exp_value_dist), 4)
	dist_res<- round((u_vec+mean(u_vec))/(exp_value_dist+mean(u_vec)), 4)

	dist_res[is.infinite(dist_res)]=max_distnorm_value

	df_result<-cbind(df,exp_value_dist,dist_res)

	result_list <- list("intercept"=intercept, "coeff_dist"=coeff_dist, "result_df"=df_result)	

	return(result_list)
}
