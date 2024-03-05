prepare_data_perpage <- function(metab_interest_sorted, metab_interest_contr, nr_plots_perpage, nr_pat=20, nr_contr=30) {
  total_nr_pages <- 0
  metab_perpage <- list()
  metab_category <- c()
  for (metab_class_index in 1:length(metab_interest_sorted)) { # "acyl_carnitines" "amino_acids" "crea_gua"
    # split list into pages, each page containing max nr_plots_perpage (20) compounds
    metab_interest_perclass <- metab_interest_sorted[[metab_class_index]]
    metab_class <- names(metab_interest_sorted)[metab_class_index]
    # add controls
    metab_interest_contr_perclass <- metab_interest_contr[[metab_class_index]]
    # number of pages for this class
    nr_pages <- ceiling(length(unique(metab_interest_perclass$HMDB_name)) / nr_plots_perpage)
    for (page_nr in 1:nr_pages) {
      total_nr_pages <- total_nr_pages + 1
      select_rows_start <- (nr_pat * nr_plots_perpage * (page_nr-1)) + 1
      select_rows_end   <-  nr_pat * nr_plots_perpage * page_nr
      metab_onepage_pat <- metab_interest_perclass[select_rows_start:select_rows_end, ]
      # same for controls
      select_rows_start_contr <- (nr_contr * nr_plots_perpage * (page_nr-1)) + 1
      select_rows_end_contr   <-  nr_contr * nr_plots_perpage * page_nr
      metab_onepage_pcontr <- metab_interest_contr_perclass[select_rows_start_contr:select_rows_end_contr, ]
      # add controls
      metab_onepage <- rbind(metab_onepage_pat, metab_onepage_pcontr)
      # if a page has fewer than nr_plots_perpage plots, fill page with empty plots
      NA_rows <- which(is.na(metab_onepage$HMDB_name))
      if (length(NA_rows) > 0) {
        # repeat the patient and control variables
        metab_onepage$variable[NA_rows] <- metab_onepage$variable[1:(nr_pat + nr_contr)]
        # for HMDB name, substitute a number of spaces
        for (row_nr in NA_rows) {
          metab_onepage$HMDB_name[row_nr] <- paste0(rep("_", ceiling(row_nr/(nr_pat + nr_contr))), collapse = "")
        }
        metab_onepage$HMDB_name <- gsub("_", " ", metab_onepage$HMDB_name)
        # leave the values at NA
      }
      # put data for one page into object with data for all pages
      metab_perpage[[total_nr_pages]] <- metab_onepage
      # create list of page headers
      metab_category <- c(metab_category, paste(metab_class, page_nr, sep="_"))
    }
  }
  # add page headers to list
  names(metab_perpage) <- metab_category
  
  return(metab_perpage)
  
}
