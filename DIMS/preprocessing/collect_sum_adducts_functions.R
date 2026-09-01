# functions for combining sum adducts

#' Combine all AdductSum parts in 1 dataframe
#'
#' @param scanmode: string with the scan mode, either "positive" or "negative"
#'
#' @returns: adductsums_all: dataframe with all adducts for one scan mode
combine_sum_adduct_parts <- function(scanmode) {
  adductsum_part_files <- list.files("./", pattern = scanmode)
  adductsums_all <- NULL
  for (i in seq_along(adductsum_part_files)) {
    load(adductsum_part_files[i])
    adductsums_all <- rbind(adductsums_all, adductsum)
  }

  return(adductsums_all)
}

#' Combine the scan modes; add intensities if both scan modes are present
#'
#' @param outlist_pos_adducts_hmdb: intensities for adducts in the positive scan mode (data frame)
#' @param outlist_neg_adducts_hmdb: intensities for adducts in the negative scan mode (data frame)
#'
#' @returns: outlist: dataframe with intensities for all metabolites present in either or both scanmodes
combine_scanmodes_intensities <- function(outlist_pos_adducts_hmdb, outlist_neg_adducts_hmdb) {
  # Only continue with samples (columns) that are in both pos and neg
  samples_both_modes <- intersect(colnames(outlist_neg_adducts_hmdb), colnames(outlist_pos_adducts_hmdb))
  outlist_neg_adducts_hmdb <- outlist_neg_adducts_hmdb[, samples_both_modes]
  outlist_pos_adducts_hmdb <- outlist_pos_adducts_hmdb[, samples_both_modes]

  # Find indices of metabolites in neg that are also found in pos and vice versa
  index_neg <- which(rownames(outlist_neg_adducts_hmdb) %in% rownames(outlist_pos_adducts_hmdb))
  index_pos <- which(rownames(outlist_pos_adducts_hmdb) %in% rownames(outlist_neg_adducts_hmdb))

  # Get intensities of metabs present in both modes from pos mode
  outlist_combi_pos_ints <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ] %>% select(-c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID))

  # Get intensities of metabs present in both modes from neg modus
  outlist_combi_neg_ints <- outlist_neg_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ] %>% select(-c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID))

  # HMDB info for metabs present in both modes
  outlist_combi_info <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ]  %>% select(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID)

  # Sum positive and negative intensities and put back HMDB columns
  outlist_combi_ints <- apply(outlist_combi_pos_ints, 2, as.numeric) + apply(outlist_combi_neg_ints, 2, as.numeric)
  rownames(outlist_combi_ints) <- rownames(outlist_combi_pos_ints)
  outlist_combi <- cbind(outlist_combi_ints, outlist_combi_info)

  # Get remaining metabolites which are not present in both scan modes
  outlist_pos <- outlist_pos_adducts_hmdb[-index_pos, ]
  outlist_neg <- outlist_neg_adducts_hmdb[-index_neg, ]

  # Combine all outlists and arrange rownames (HMDB codes)
  outlist <- rbind(outlist_combi, outlist_pos, outlist_neg)
  outlist <- outlist %>% arrange(rownames(outlist))
  return(outlist)
}
