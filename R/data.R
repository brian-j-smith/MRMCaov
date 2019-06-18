#' Multi-reader multi-case dataset
#' 
#' @format A data frame with 1140 rows and 6 variables:
#' \describe{
#'   \item{reader}{reader identifier}
#'   \item{treatment}{treatment identifier}
#'   \item{case}{case identifier (factorial design)}
#'   \item{case2}{case identifier (cases nested within readers)}
#'   \item{case3}{case identifier (cases nested within treatments)}
#'   \item{truth}{true case status (1 = positive, 0 = negative)}
#'   \item{rating}{ordinal reader ratings of case status}
#' }
#' 
"VanDyke"