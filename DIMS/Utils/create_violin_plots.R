# remove parameter default value
# add explanation variable to parameters
create_violin_plots <- function(pdf_dir, pt_name, metab_perpage, top_metab_pt = NULL) {
  #' Create violin plots for each patient
  #'
  #' @param pdf_dir: location where to save the pdf file (string)
  #' @param pt_name: patient code (string)
  #' @param metab_perpage: list of dataframe with a dataframe for each page in the violinplot pdf (list)
  #' @param top_metab_pt: dataframe with increased and decrease metabolites for this patient (dataframe)

  # set parameters for plots
  plot_height <- 9.6
  plot_width <- 6
  fontsize <- 1
  circlesize <- 0.8
  colors_4plot <- c("#22E4AC", "#00B0F0", "#504FFF", "#A704FD", "#F36265", "#DA0641")
  #                   green     blue      blue/purple purple    orange     red

  # patient plots, create the PDF device
  pt_name_sub <- pt_name
  suffix <- ""
  if (grepl("Diagnostics", pdf_dir) & is_diagnostic_patient(pt_name)) {
    prefix <- "MB"
    suffix <- "_DIMS_PL_DIAG"
    # substitute P and M in P2020M00001 into right format for Helix
    pt_name_sub <- gsub("[PM]", "", pt_name)
    pt_name_sub <- gsub("\\..*", "", pt_name_sub)
  } else if (grepl("Diagnostics", pdf_dir)) {
    prefix <- "Dx_"
  } else if (grepl("IEM", pdf_dir)) {
    prefix <- "IEM_"
  } else {
    prefix <- "R_"
  }

  pdf(paste0(pdf_dir, "/", prefix, pt_name_sub, suffix, ".pdf"),
      onefile = TRUE,
      width = plot_width,
      height = plot_height)

  # page headers:
  page_headers <- names(metab_perpage)

  # put table into PDF file, if not empty
  if (!is.null(dim(top_metab_pt))) {
    plot.new()
    # get the names and numbers in the table aligned
    table_theme <- ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 6)),
                                  colhead = list(fg_params = list(fontsize = 8, fontface = "bold")))
    grid.table(top_metab_pt, theme = table_theme, rows = NULL)
    # g <- tableGrob(top_metab_pt)
    # grid.draw(g)
    text(x = 0.45, y = 1.02, paste0("Top deviating metabolites for patient: ", pt_name), font = 1, cex = 1)
  }

  # violin plots
  for (page_index in 1:length(metab_perpage)) {
    # extract list of metabolites to plot on a page
    metab_list_2plot <- metab_perpage[[page_index]]
    # extract original data for patient of interest (pt_name) before cut-offs
    pt_list_2plot_orig <- metab_list_2plot[which(metab_list_2plot$variable == pt_name), ]
    # cut off Z-scores higher than 20 or lower than -5 (for nicer plots)
    metab_list_2plot$value[metab_list_2plot$value >  20] <-  20
    metab_list_2plot$value[metab_list_2plot$value <  -5] <-  -5
    # extract data for patient of interest (pt_name)
    pt_list_2plot <- metab_list_2plot[which(metab_list_2plot$variable == pt_name), ]
    # restore original Z-score before cut-off, for showing Z-scores in PDF
    pt_list_2plot$value_orig <- pt_list_2plot_orig$value
    # remove patient of interest (pt_name) from list; violins will be made up of controls and other patients
    metab_list_2plot <- metab_list_2plot[-which(metab_list_2plot$variable == pt_name), ]
    # subtitle per page
    sub_perpage <- gsub("_", " ", page_headers[page_index])
    # for IEM plots, put subtitle on two lines
    sub_perpage <- gsub("probability", "\nprobability", sub_perpage)
    # add size parameter for showing Z-score of patient per metabolite
    z_size <- rep(3, nrow(pt_list_2plot))
    # set size to 0 if row is empty
    z_size[is.na(pt_list_2plot$value)] <- 0

    # draw violin plot. shape=22 gives square for patient of interest
    ggplot_object <- ggplot(metab_list_2plot, aes(x = value, y = HMDB_name)) +
      theme(axis.text.y = element_text(size = rel(fontsize)), plot.caption = element_text(size = rel(fontsize))) +
      xlim(-5, 20) +
      geom_violin(scale = "width") +
      geom_point(data = pt_list_2plot, aes(color = value), size = 3.5 * circlesize, shape = 22, fill = "white") +
      scale_fill_gradientn(
        colors = colors_4plot, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar",
        aesthetics = "colour"
      ) +
      # add Z-score value for patient of interest at x=16
      geom_text(
        data = pt_list_2plot, aes(16, label = paste0("Z=", round(value_orig, 2))), hjust = "left", vjust = +0.2,
        size = z_size
      ) +
      # add labels. Use font Courier to get all the plots in the same location.
      labs(x = "Z-scores", y = "Metabolites", subtitle = sub_perpage, color = "z-score") +
      theme(axis.text.y = element_text(family = "Courier", size = 6)) +
      # do not show legend
      theme(legend.position = "none") +
      # add title
      ggtitle(label = paste0("Results for patient ", pt_name)) +
      # add vertical lines
      geom_vline(xintercept = 2, col = "grey", lwd = 0.5, lty = 2) +
      geom_vline(xintercept = -2, col = "grey", lwd = 0.5, lty = 2)

    suppressWarnings(print(ggplot_object))

  }

  # add explanation of violin plots, version number etc.
  plot(NA, xlim = c(0, 5), ylim = c(0, 5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  if (length(explanation) > 0) {
    text(0.2, 5, explanation[1], pos = 4, cex = 0.8)
    for (line_index in 2:length(explanation)) {
      text_y_position <- 5 - (line_index * 0.2)
      text(-0.2, text_y_position, explanation[line_index], pos = 4, cex = 0.5)
    }
  }

  # close the PDF file
  dev.off()
}
