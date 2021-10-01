suppressMessages({
library("htmlTable")
library("stringr")
library("R.utils")
})

## Helper functions for running the scripts and providing informative output
print_start_msg <- function(section_name, section_number) {
  cat(paste("\n\n######## STARTING", section_name, "OF SECTION", section_number, "#############\n\n"))
}

print_run_time <- function(section_name, section_number, total_time) {
  cat(paste("\nFinished", tolower(section_name), "of Section", section_number, "in", round(total_time["elapsed"] / 60, 4), "minutes.\n"))
}

source_wrapper <- function(script, section_name, section_number) {
  print_start_msg(section_name, section_number)
  total_time <- system.time(source(paste0("scripts/", script)))
  print_run_time(section_name, section_number, total_time)
  rm(list = setdiff(ls(), c("quick", "print_start_msg", "print_run_time")))
  invisible(gc())
}

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
  
  df %>% htmlTable::htmlTable(...) %>% 
    as.character %>% 
    stringr::str_remove_all('\n') %>% 
    noquote() %>% 
    cat(file = file)
}


## Checks if the argument 'quick' is set, either as a variable in the current 
## workspace or provided from the command line. 
check_quick <- function() {
  
  if (!exists("quick")) {

    ## Read in quick from the command line (i.e., from the makefile)
    args <- R.utils::commandArgs(trailingOnly = TRUE, asValue = TRUE)
    if (length(args) == 0) {
      cat("You have not specified whether or not you want to use quick, very-low-dimensional representations of the models: Setting quick = TRUE.\n")
      quick <- TRUE
    } else if (length(args)==1) {
      # quick <- as.integer(args[1])
      quick <- as.logical(as.numeric(args[1]))
    } else {
      stop("check_quick() assumes that at most one command line argument is provided.")
    }  
  }
  
  if (quick) {
    cat("Running very-low-dimensional representations of the models: This should be relatively quick.\n")
  } else {
    cat("Running full models: This may take a while.\n")
  }
  
  return(quick)
}


