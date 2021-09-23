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



## Package used for reading command line argument
library("R.utils")

## Use extremely low-rank versions of the models to quickly establish that the 
## code works? 
check_quick <- function() {
  if (!exists("quick")) {
    ## Read in low-rank from the command line (i.e., from the makefile)
    args <- R.utils::commandArgs(trailingOnly = TRUE, asValue = TRUE)
    if (length(args) == 0) {
      cat("You have not specified whether or not you want to use quick, low-rank 
       versions of the models: Setting quick = TRUE\n")
      quick <- TRUE
    } else if (length(args)==1) {
      quick <- as.logical(args[1])
    } else {
      stop("Too many arguments to deal with!")
    }  
  }
  
  if (quick) {
    cat("Running low-rank versions of the models: This should be relatively quick\n")
  } else {
    cat("Not running low-rank versions of the models: This may take a while\n")
  }
  
  return(quick)
}
