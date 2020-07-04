#' Filter input data frame rows
#'
#' Takes formatted dataframe and filtered df
#' @importFrom stringr str_split_fixed
#' @param df Input data frame.
#' @param check_trans_values If TRUE, rows containing trans interactions dropped.
#' @param remove_zero_contact If TRUE, rows wifh $freq value 0 dropped.
#' @param self_ligation_distance Distance between two fragments regarded as self-ligation. Default 15kb.
#' @return df with rows filtered
#' @export

filterInputDF <- function(df, check_trans_values=TRUE, remove_zero_contact=TRUE, self_ligation_distance=15000)
{
	df <- df[which(str_split_fixed(df$frag1,"\\.",3)[,1] == str_split_fixed(df$frag2,"\\.",3)[,1]),]
	if (remove_zero_contact){df <- df[which(df$freq > 0),] }
	if (self_ligation_distance > 0){ df <- df[which(df$dist > self_ligation_distance),]}

	return(df)
}

