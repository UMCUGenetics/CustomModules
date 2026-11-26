## adapted from elementInfo.R, which is adapted from Rdisop function .getElement
# refactor: check where library is initialised: library(Rdisop)
get_element_info <- function(name, elements = NULL) {
  #' Get info on m/z and isotopes for all chemical elements
  #'
  #' @param name: Name of adduct, e.g. Na (string)
  #' @param elements: List of all adducts to take into account (list of strings)
  #'
  #' @return element_info: peak group list with filled-in intensities (matrix)

  # get information on all elements
  if (!is.list(elements) || length(elements) == 0 ) {
    elements <- initializePSE()
  }
  # extract information for a particular adduct
  if (name == "CH3OH+H") {
    # regular_expr should be exact match for name, except for methanol
    regular_expr <- "^CH3OH\\+H$"
  } else {
    regular_expr <- paste0("^", name, "$")
  }
  element_info <- elements[[grep(regular_expr, sapply(elements, function(x) { x$name }))]]
  return(element_info)
}