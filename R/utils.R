#' Generate attribute pattern table header
#' 
#' 
#' @param k     Number of Attributes.
#' @param m     Number of Categories. Default 2 or dichotomous response.
#' @param order Order of the table. Default `k` or the full order.
#' 
#' @return
#' Return a matrix containing the class table
#' 
#' @export
#' @examples
#' attribute_pattern_table_header(3)
#' 
#' attribute_pattern_table_header(4)
attribute_pattern_table_header = function(k, m = 2, order = k) {
  n_class = m^k
  d_to_q = t(GenerateAtable(n_class, k, m, order)$DtoQtable)
  # Construct strings
  apply(d_to_q, 1, paste0, collapse = "")
} 
