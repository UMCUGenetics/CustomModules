combine_sum_adduct_parts <- function(scanmode) {
  #' Combine all AdductSum parts in 1 dataframe
  #'
  #' @param scanmode: string with the scanmodus, either positive or negative
  #'
  #' @returns: outlist_tot: dataframe with all adducts for one scanmodus
  adductsum_part_files <- list.files("./", pattern = scanmode)
  outlist_tot <- NULL
  for (i in 1:length(adductsum_part_files)) {
    load(adductsum_part_files[i])
    outlist_tot <- rbind(outlist_tot, adductsum)
  }

  return(outlist_tot)
}

combine_scanmodes_intensities <- function(outlist_pos_adducts_hmdb, outlist_neg_adducts_hmdb) {
  #' Combine the scanmodes and when present in both scanmodes add intensities
  #'
  #' @param outlist_pos_adducts_hmdb: dataframe with adducts for the positive scanmodus
  #' @param outlist_neg_adducts_hmdb: dataframe with adducts for the positive scanmodus
  #'
  #' @returns: outlist: dataframe with intensities for all metabolites present in either or both scanmodes

  # Only continue with patients (columns) that are in both pos and neg, so patients that are in both
  samples_both_modes <- intersect(colnames(outlist_neg_adducts_hmdb), colnames(outlist_pos_adducts_hmdb))
  outlist_neg_adducts_hmdb <- outlist_neg_adducts_hmdb[, samples_both_modes]
  outlist_pos_adducts_hmdb <- outlist_pos_adducts_hmdb[, samples_both_modes]

  # Find indexes of neg hmdb code that are also found in pos and vice versa
  index_neg <- which(rownames(outlist_neg_adducts_hmdb) %in% rownames(outlist_pos_adducts_hmdb))
  index_pos <- which(rownames(outlist_pos_adducts_hmdb) %in% rownames(outlist_neg_adducts_hmdb))

  # Get intensities of metabs present in both modes from pos modus
  outlist_combi_pos_ints <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ] %>% select(-c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID))

  # Get intensities of metabs present in both modes from neg modus
  outlist_combi_neg_ints <- outlist_neg_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ] %>% select(-c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID))

  # HMDB info for metabs present in both modes
  outlist_combi_info <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ]  %>% select(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID)

  # Combine positive and negative numbers and paste back HMDB column
  outlist_combi_ints <- apply(outlist_combi_pos_ints, 2, as.numeric) + apply(outlist_combi_neg_ints, 2, as.numeric)
  rownames(outlist_combi_ints) <- rownames(outlist_combi_pos_ints)
  outlist_combi <- cbind(outlist_combi_ints, outlist_combi_info)

  # Get outlist with metabs in either pos or neg
  outlist_pos <- outlist_pos_adducts_hmdb[-index_pos, ]
  outlist_neg <- outlist_neg_adducts_hmdb[-index_neg, ]

  # Combine all outlists and arrange rownames (HMDB codes)
  outlist <- rbind(outlist_combi, outlist_pos, outlist_neg)
  outlist <- outlist %>% arrange(rownames(outlist))
  return(outlist)
}