#' Multi-reader multi-case dataset
#'
#' Data from a multi-reader multi-case study conducted by Kundel et al. (1997)
#' to compare the clinical diagnostic accuracy of hard-copy readings of
#' screen-film bedside chest radiographs and both hard- and soft-copy readings
#' of computed radiographs from separate groups of patients in a medical
#' intensive care unit.
#'
#' @format A data frame with 760 rows and 5 variables:
#' \describe{
#'   \item{reader}{reader identifier}
#'   \item{treatment}{numeric identifier for treatment modality}
#'   \item{treatment_name}{name of treatment modality}
#'   \item{case}{patient identifier}
#'   \item{truth}{true case status (1 = diseased, 0 = non-diseased)}
#'   \item{rating}{ordinal reader ratings of case status}
#' }
#'
#' @references
#' Kundel HL, Gefter W, Aronchick J, Miller W Jr, Hatabu H, Whitfill CH, Miller
#' W Sr (1997). Accuracy of bedside chest hard-copy screen-film versus hard- and
#' soft-copy computed radiographs in a medical intensive care unit: receiver
#' operating characteristic analysis. Radiology, 205(3): 859-63.
#'
"Kundel"


#' Multi-reader multi-case dataset
#'
#' @format A data frame with 800 rows and 5 variables:
#' \describe{
#'   \item{Reader}{reader identifier}
#'   \item{Treatment}{treatment identifier}
#'   \item{Case}{case identifier}
#'   \item{Truth}{true case status (1 = abnormal, 0 = normal)}
#'   \item{Rating}{ordinal reader ratings of case normal status (1 = definitely
#'   abnormal, 5 = definitely normal)}
#' }
#'
#' @references
#' Franken EA Jr, Berbaum KS, Marley SM, Smith WL, Sato Y, Kao SC, Milam SG
#' (1992). Evaluation of a digital workstation for interpreting neonatal
#' examinations: a receiver operating characteristic study. Investigational
#' Radiology, 27(9): 732-737.
#'
"Franken"


#' Multi-reader multi-case dataset
#'
#' @format A data frame with 1140 rows and 7 variables:
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
#' @references
#' Van Dyke CW, White RD, Obuchowski NA, Geisinger MA, Lorig RJ, Meziane MA
#' (1993). Cine MRI in the diagnosis of thoracic aortic dissection. 79th
#' Radiological Society of North America Meetings, Chicago, IL.
#'
"VanDyke"
