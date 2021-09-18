library("htmlTable")
library("stringr")
save_html_table <- function(df, col_sep = 3, decimals = 3, file = NULL, ...) {
  
  # Round to decimals
  for (i in 1:ncol(df)) {
    if(is.numeric(df[, i])) 
      df[, i] <- round(df[, i], decimals)
  }
  
  # Hack to add more space between the columns:
  # (just adds non-separable white-space units before and after each name)
  white_space <- paste0(rep("&nbsp;", col_sep), collapse = "")
  names(df) <- paste(white_space, names(df), white_space)
  
  df %>% htmlTable(...) %>% 
    as.character %>% 
    str_remove_all('\n') %>% 
    noquote() %>% 
    cat(file = file)
}

