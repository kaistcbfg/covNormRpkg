#' Coverage Normalization function
#'
#' Takes formatted dataframe and retunrs glm.nb fitted coeffs & normalized df
#' @importFrom MASS glm.nb
#' @param df An input dataframe.
#' @param cov1_thresh&cov2_thresh Threshold for filteraing bins with too low coverage. 
#' @param max_covnorm_value Outliers by division with too small a number are replaced by this value.
#' @param sample_ratio If NOT -1 && between 0.0-1.0, downsample the df rows to be used for fitting by that ratio.
#' @return df Coverage normalized.
#' @export

normCoverage <- function(df, do_shuffle=TRUE, cov1_thresh=200, cov2_thresh=200, max_covnorm_value=50, sample_ratio=-1)
{
	df <- df[which(df$cov_frag1 > cov1_thresh & df$cov_frag2 > cov2_thresh),]

	df$rand = sample(100,size=nrow(df), replace=TRUE)
	if(do_shuffle){df[df$rand > 50,c("cov_frag1","cov_frag2")] = df[df$rand > 50 ,c("cov_frag2","cov_frag1")]}	

	u_vec = df$freq
	capture_vec1=scale_input(df$cov_frag1)
	capture_vec2=scale_input(df$cov_frag2)

	u_vec_fit <- u_vec
	capture_vec1_fit <- capture_vec1
	capture_vec2_fit <- capture_vec2

	if(sample_ratio > 0)
	{
		if(sample_ratio > 1.0){sample_ratio = 1.0}
		df_ds <- df[sample(nrow(df), size = as.integer(nrow(df) * sample_ratio)),]
		
		u_vec_ds = df_ds$freq
		capture_vec1_ds=scale_input(df_ds$cov_frag1)
		capture_vec2_ds=scale_input(df_ds$cov_frag2)

		u_vec_fit <- u_vec_ds
		capture_vec1_fit <- capture_vec1_ds
		capture_vec2_fit <- capture_vec2_ds
	}

	#capture_vec1_fit <- b_spline_fit(capture_vec1_fit)
	#capture_vec2_fit <- b_spline_fit(capture_vec2_fit)

	basis <- cbind(capture_vec1_fit,capture_vec2_fit)
	fit<-glm.nb(u_vec_fit ~ basis)
	coeff<-round(fit$coeff,4)

	intercept = coeff[1]
	coeff_cov1= coeff[2]
	coeff_cov2= coeff[3]

	nknots = round((length(u_vec))^(0.2)*2)
	exp_value_capture=round(exp(intercept + coeff_cov1*capture_vec1 + coeff_cov2*capture_vec2),4)
	#exp_value_capture=round(exp(as.vector(intercept+capture_vec1_fit%*%coeff[2:(nknots+4)]+capture_vec2_fit%*%coeff[(nknots+5):(2*nknots+7)])), 4)
	capture_res=round((u_vec)/(exp_value_capture), 4)
	capture_res[is.infinite(capture_res)] = max_covnorm_value
	df_result<-cbind(df,exp_value_capture,capture_res)

	result_list <- list("intercept"=intercept, "coeff_cov1"=coeff_cov1, "coeff_cov2"=coeff_cov2, "result_df"=df_result)	

	return(result_list)
}

