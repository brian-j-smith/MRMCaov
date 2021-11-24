#' Convert Obuchowski-Rockette to Roe & Metz Parameters
#'
#' Determines Roe & Metz (RM) simulation model parameters for simulating
#' multireader multicase likelihood-of-disease rating data based on real-data or
#' conjectured Obuchowski-Rockette (OR) parameter estimates that describe the
#' distribution of the empirical AUC reader performance measure. The algorithm
#' assumes the constrained unequal-variance RM model (Hillis, 2012), which
#' generalizes the original RM model (Roe and Metz, 1997) by allowing the
#' diseased and nondiseased decision-variable distributions to have unequal
#' variances for each reader with the variance components
#' involving diseased cases constrained to differ by a factor of 1/b^2 from
#' corresponding variance components involving nondiseased cases. This algorithm
#' is described in Hillis (2020).
#'
#' A related function is the  RMH_to_OR function, which determines OR parameters
#' that describe the distribution of empirical AUC estimates computed from
#' inputted RM model parameter values, based on the analytical mapping provided
#' by Hillis (2018).
#'
#' @rdname OR_to_RMH
#'
#' @param AUC1,AUC2 test 1 and 2 expected empirical AUCs.
#' @param var_R,var_TR OR reader and test-by-reader variance components.
#' @param corr1,corr2,corr3 OR error correlations.
#' @param var_error OR error variance.
#' @param n0,n1 number of nondiseased and diseased cases.
#' @param b_method method of estimating RM b parameter.
#' @param mean_sig_input mean-to-sigma ratio, required only if
#'   \code{b_method = "mean_to_sigma"}.
#' @param b_input binormal \emph{b} value, required only if
#'   \code{b_method = "specified"}.
#' @param b_le_1 logical indicating whether the algorithm searches first for
#'   b <= 1 and then, if no solution, for b >= 1; if FALSE, the algorithm
#'   searches only for for b >= 1. Required only if
#'   \code{b_method = "unspecified"}.
#' @param params data frame of above OR parameter values in the columns.
#' @param ... arguments passed to the default method.
#'
#' @details
#' Hillis (2012) modified the original RM model (Roe and Metz, 1997) by allowing
#' variance components involving case to depend on truth (diseased/nondiseased),
#' with variance components involving diseased cases set equal to those
#' involving nondiseased cases multiplied by the factor 1/b^2, b>0.  Hillis
#' (2018) derived analytical formulas, which apply to the Hillis (2012) RM
#' model, that express OR parameters describing the distribution of empirical
#' AUC outcomes computed from RM simulated data as functions of the RM
#' parameters. This mapping from the RM parameters to the OR parameters is
#' implemented in \R by the  RMH_to_OR function.
#'
#' For the constrained unequal-variance RM model, the OR-to-RM algorithm
#' provides the reverse transformation that determines the corresponding RM
#' parameters.  This algorithm is described in Hillis (2020).  The OR_to_RMH
#' function implements this algorithm using an iterative search procedure.
#'
#' \code{b_method} indicates the method for estimating the RM binormal \emph{b}
#' parameter.
#' \itemize{
#'   \item \code{b_method = "unspecified"} should be used when the goal is to
#'     determine RM parameters that result in simulated data for which the
#'     empirical AUC distribution is described by the inputted values for the
#'     OR parameter vector \deqn{\beta_OR = (AUC1, AUC2, var_R, var_TR,
#'     var_error, corr1, corr2, corr3).}
#'   \item \code{b_method = "mean_to_sigma"} should be used when the goal is to
#'     determine RM parameters that result in simulated data for which the
#'     empirical AUC distribution is described by the inputted values for the OR
#'     parameter vector \deqn{\beta1_OR = (AUC1, AUC2, var_R, var_TR, corr1,
#'     corr2, corr3),} and such that the median mean-to-sigma ratio across
#'     readers is equal to \code{mean_sig_input} for the test having the lowest
#'     AUC.  Note that \eqn{\beta1_OR} differs from \eqn{\beta_OR} in that it
#'     does not contain the OR error variance.
#'   \item \code{b_method = "specified"} should be used when the goal is to
#'     determine RM parameters that result in simulated data for which the
#'     empirical AUC distribution is described by the inputted values for the OR
#'     parameter vector \eqn{\beta1_OR} (see above) with \emph{b} equal to
#'     \code{mean_sig_input}.  (E.g., set \code{b_input = 1} for symmetric
#'     ROC curve.)
#' }
#'
#' For \code{b_method = "mean_to_sigma"} or \code{"specified"}, the simulated
#' empirical AUC estimate distribution is specified by the parameter values
#' in \code{params}, except for \code{var_error}. Thus for these two options,
#' \code{var_error} can be equal to \code{NA} or excluded from \code{params}.
#'
#' Parameter \code{mean_sig_input} is the inputted mean-to-sigma ratio needed
#' for \code{b_method = "mean_to_sig"}. See Hillis & Berbaum (2011) for more
#' information.
#'
#' Parameter \code{b_input} is the inputted binormal \emph{b} value needed for
#' \code{b_method = "specified"}.
#'
#' There may not be a solution for a set of OR parameters values.  When this
#' occurs, the function will either produce an approximate solution or indicate
#' what OR input needs to be changed.
#'
#' @return
#' The RM model parameters are returned in a data frame with the following
#' elements.
#'
#' \describe{
#'   \item{delta1}{mean separation of nondiseased and diseased decision-variable
#'     distributions for test 1 across reader population.}
#'   \item{delta2}{mean separation of nondiseased and diseased decision-variable
#'     distributions for test 2 across reader population.}
#'   \item{var_R}{RM reader variance  compnent.}
#'   \item{var_TR}{RM text-by-reader variance component.}
#'   \item{var_C}{RM case variance component.}
#'   \item{var_TC}{RM test-by-case variance.}
#'   \item{var_RC}{RM reader-by-case variance.}
#'   \item{var_error}{RM error variance.}
#'   \item{b}{variance components involving diseased cases are constrained to
#'     differ by a factor of 1/b^2 from corresponding variance components
#'     involving nondiseased cases.}
#' }
#'
#' Related quantities that are also returned in the data frame:
#' \describe{
#'   \item{b_method}{method used to estimate b.}
#'   \item{n0}{number of nondiseased cases per simulated sample.}
#'   \item{n1}{number of diseased cases per simulated sample.}
#'   \item{mean_to_sig1}{expected mean-to-sigma ratio across readers for test
#'     1.}
#'   \item{mean_to_sig2}{expected mean-to-sigma ratio across readers for test
#'     2.}
#'   \item{Pr1_improper}{probability that the test 1 ROC curve for a random
#'     reader will be visually improper (i.e, |mean-to-sigma ratio| < 3).}
#'   \item{Pr2_improper}{probability that the test 2 ROC curve for a random
#'     reader will be visually improper (i.e, |mean-to-sigma ratio| < 3).}
#' }
#'
#' @references
#' Hillis, Stephen L. 2012. "Simulation of Unequal-Variance Binormal Multireader
#' ROC Decision Data: An Extension of the Roe and Metz Simulation Model."
#' \emph{Academic Radiology} no. 19 (12):1518-1528.
#' doi: 10.1016/j.acra.2012.09.011.
#'
#' Hillis, Stephen L. 2018. "Relationship between Roe and Metz simulation model
#' for multireader diagnostic data and Obuchowski-Rockette model parameters."
#' \emph{Statistics in Medicine} no. 37 (13):2067-2093. doi: 10.1002/sim.7616.
#'
#' Hillis, Stephen L. 2020. "Determining Roe and Metz model parameters for
#' simulating multireader multicase confidence-of-disease rating data based on
#' read-data or conjectured Obuchowski-Rockette parameter estimates." Vol.
#' 11316, SPIE Medical Imaging: SPIE. doi.org/10.1117/12.2550541
#'
#' Hillis, Stephen L., and Kevin S. Berbaum. 2011. "Using the mean-to-sigma
#' ratio as a measure of the improperness of binormal ROC curves."
#' \emph{Academic Radiology} no. 18 (2):143-154. doi:
#' 10.1016/j.acra.2010.09.002.
#'
#' Roe, Cheryl A., and Charles E. Metz. 1997. "Dorfman-Berbaum-Metz method for
#' statistical analysis of multireader, multimodality receiver operating
#' characteristic data: Validation with computer simulation." \emph{Academic
#' Radiology} no. 4 (4):298-303. doi: 10.1016/S1076-6332(97)80032-3.
#'
#' @author
#' Stephen L. Hillis, Departments of Radiology and Biostatistics,
#' University of Iowa, \email{steve-hillis@uiowa.edu}
#'
#' Brian J. Smith, Department of Biostatistics, University of Iowa,
#' \email{brian-j-smith@uiowa.edu}
#'
#' @seealso \code{\link{RMH_to_OR}}
#'
#' @examples
#' ## Example 1: Computing RM parameters from OR parameters directly
#' ##--------------------------------------------------------------
#' ## Example 1a: Using b_method ="unspecified" (the default)
#' RM <- OR_to_RMH(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                 corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                 var_R = 0.00154, var_TR = 0.000208, var_error = 0.000788)
#' RM
#' ##  We recommend also computing the OR parameter values ("true values")
#' # that describe the distribution of simulated data based on above RM parameters,
#' # using the RMH_to_OR function. Ideally the true values will be the same as the
#' # inputted OR values used for deriving the RM parameter values. We recommend
#' # always performing this check.  This check is carried out below, as shown below.
#' true_values = RMH_to_OR(RM)
#' true_values
#' #   From the output we see, for this example, that the true OR values are identical to the
#' # inputted OR values
#'
#'
#' # Example 1b: Using b_method = "specified" with b_input = 1
#' #   Note that the error variance does not need to be specified since this b_method
#' # does not utilize it.
#' RM <- OR_to_RMH(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                 corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                 var_R = 0.00154, var_TR = 0.000208,
#'                 b_method = "specified", b_input = 1)
#' RM
#' true_values <- RMH_to_OR(RM)
#' true_values
#'
#' #  From the output we see, for this example, that the true values are identical
#' # (within rounding error) to the inputted OR values (but note that var_error was
#' # not inputted)
#'
#' ## Example 1c: Using b_method = "mean_to_sigma" with mean_to_sig_input = 4.5
#' #   Note the error variance does not need to be  specified since this b_method
#' # does not utilize it.
#' RM <- OR_to_RMH(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                 corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                 var_R = 0.00154, var_TR = 0.000208,
#'                 b_method = "mean_to_sigma", mean_sig_input = 4.5)
#' RM
#' true_values <- RMH_to_OR(RM)
#' true_values
#' #   From the output we see for this example that the true OR values are identical
#' # (within rounding error) to the inputted OR values (but note that var_error was
#' # not inputted)
#'
#' ##---------------------------------------------------------------------
#'
#' ## Example 2: Computing RM parameters from a data frame of OR parameters
#' ## ---------------------------------------------------------------------
#' ## Example 2a: One study
#' vandyke_OR <- data.frame(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                          corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                          var_R = 0.00154, var_TR = 0.000208, var_error = 0.000788)
#' vandyke_RM <- OR_to_RMH(vandyke_OR)
#' vandyke_RM
#' true_values <- RMH_to_OR(vandyke_RM)
#' true_values
#'
#' ## Example 2b: Three studies
#' three_studies_OR <- data.frame(
#'   rbind(
#'     vandyke = c(69, 45, 0.89793704, 0.94083736, 0.432, 0.429, 0.298, 0.00154,
#'                 0.0002, 0.00080229),
#'     franken = c(33, 67, .8477498869, 0.8368950701, 0.521430051, 0.319691199,
#'                 0.3386375697, 0.0000433385, 0.0, 0.0014967254),
#'     kundel = c(66, 29, 0.8038793103, 0.8413662487, 0.507695244, 0.3843523762,
#'                0.4035662578, 0.0007340122, 0, 0.002148844)
#'   )
#' )
#' colnames(three_studies_OR) <- c("n0", "n1", "AUC1", "AUC2", "corr1", "corr2",
#'                                 "corr3", "var_R", "var_TR", "var_error")
#' three_studies_OR
#' three_studies_RM <- OR_to_RMH(three_studies_OR)
#' three_studies_RM
#' true_values <- RMH_to_OR(three_studies_RM)
#' true_values
#' ##   Note above that the true values for corr2 and corr3 for the Franken study
#' # differ slightly from the inputted values; this is because corr2 < corr3 for the
#' # inputted OR values, which is not possible for simulated RM model data.
#'
#' ##Example 2c: Examples 1a, 1b and 1c run using one data frame
#' vandyke_OR <- data.frame(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                          corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                          var_R = 0.00154, var_TR = 0.000208, var_error = 0.000788)
#' vandyke_OR_x3 <- vandyke_OR[c(1,1,1),] #has 3 rows, each equal to vandyke_OR
#' b_method = c("unspecified","mean_to_sigma","specified")
#' mean_sig_input = c(NA,4.5,NA)
#' b_input = c(NA,NA,1)
#' vandyke_OR_3ex <- cbind(vandyke_OR_x3,b_method,mean_sig_input,b_input)
#' vandyke_OR_3ex
#' vandyke_OR_3ex_RM <- OR_to_RMH(vandyke_OR_3ex)
#' vandyke_OR_3ex_RM
#' true_values <- RMH_to_OR(vandyke_OR_3ex_RM)
#' true_values
#'
#'
#' ## Example 3: Printing the alternative x1 -- x7 parameter values
#' ## ---------------------------------------------------
#' ## The OR_to_RMH function first finds the solutions using the alternative RM
#' ## parameterization consisting of b and the alternative parameters
#' ## x1, x2, x3, x4, x5, x6, and x7, and then solves for the conventional RM
#' ## parameters in terms of these alternative parameters.  (See Hillis (2020) for details.)
#' ## Although the user generally has no need to know these parameter values, they
#' ## can be printed out using the all = TRUE print option, as shown below
#' ## using Example 1a:
#'
#' RM <- OR_to_RMH(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                var_R = 0.00154, var_TR = 0.000208, var_error = 0.000788)
#' print(RM,all = TRUE)
#'
OR_to_RMH <- function(...) {
  UseMethod("OR_to_RMH")
}


#' @rdname OR_to_RMH
#'
OR_to_RMH.default <- function(...,
  AUC1, AUC2, var_R, var_TR, corr1, corr2, corr3, var_error = NULL, n0, n1,
  b_method = c("unspecified", "mean_to_sigma", "specified"),
  mean_sig_input = NULL, b_input = NULL, b_le_1 = TRUE
) {

  AUC1_OR <- AUC1
  AUC2_OR <- AUC2
  var_R_OR <- var_R
  var_TR_OR <- var_TR
  corr1_OR <- corr1
  corr2_OR <- corr2
  corr3_OR <- corr3
  var_error_OR <- var_error # can be NA if use b_method = "unspecified" or "mean_to_sigma"
  increment1 <- .001 #initial search step increment for b_method = 1;
  increment2 <- .000001 #search step increment for b_method = 1 after using increment1 steps;
  lower_bd <- 0.01 #lower bound for unspecified method b search
  step_no <- 30 #number of times interval is split in 2 when searching for solution
  flag1 <- 0 #becomes 1 when an x value (x1 -- x7) or b is NA
  b <- NA
  x1 <- NA
  x2 <- NA
  x3 <- NA
  x4 <- NA
  x5 <- NA
  x6 <- NA
  x7 <- NA

  b_method <- match.arg(b_method)

  # check for out of bounds values
  inbounds <- c(
    n0 = (n0 > 0),
    n1 = (n1 > 0),
    AUC1 = (.5 <= AUC1_OR  & 1 >= AUC1_OR),
    AUC2 = (.5 <= AUC2_OR  & 1 >= AUC2_OR),
    var_R = (0 <= var_R_OR),
    var_TR = (0 <= var_TR_OR),
    corr1 = (0 <= corr1_OR & corr1_OR <= 1),
    corr2 = (0 <= corr2_OR & corr2_OR <= 1),
    corr3 = (0 <= corr3_OR & corr3_OR <= 1),
    var_error = (0 <= var_error_OR | is.na(var_error_OR))
  )
  invalid_params <- names(inbounds)[!inbounds]
  if (length(invalid_params)) {
    stop("OR_to_RMH parameter values out of bound for ",
         toString(invalid_params), call. = FALSE)
  }

  x1 <- qnorm(AUC1_OR)
  x2 <- qnorm(AUC2_OR)
  c1 <- 1/(n0*n1)
  c2 <- (n1-1)/(n0*n1)
  c3 <- (n0-1)/(n0*n1)
  c4 <- (1-n0-n1)/(n0*n1)

  # STEP 2: Compute x3 (= 2*sigma_R_OR**2/V)
  x <- .5
  for (step in 2:step_no) {
    temp <- multnorm (x1,x2,x) -  pnorm(x1)*pnorm(x2)
    x <- x - sign(temp - var_R_OR)*.5^step
  }
  x3 <- x
  temp1 <- multnorm (x1,x2,1) -  pnorm(x1)*pnorm(x2)
  if (x3 > 1 - .5^(step_no - 1)) {
    x3 <- NA
    flag1 <- 1
  }
  #STEP 3: Compute x4 (= 2(sigma_R_OR**2 + sigma_TR_OR**2)/V)

  if (flag1 == 1) {
    x4 <- NA
  }
  if (flag1 == 0){
    x <- .5
    for (step in 2:step_no) {
      temp <- .5*(multnorm(x1,x1,x) - (pnorm(x1))**2 + multnorm(x2,x2,x) -
                    (pnorm(x2))**2) - var_R_OR
      x <- x - sign(temp - var_TR_OR)*.5^step
    }
    x4 <- max(x, x3)
    if (x4 > 1 - .5**(step_no - 1)) {
      x4 <- NA
      flag1 <- 1
    }
  }


  # STEP 4 Compute b
  if (flag1 == 1) {
    b <- NA
    cov1 <- NA
    cov2 <- NA
    cov3 <- NA
  }
  if (flag1 == 0) {
    if (b_method == "mean_to_sigma") { # Determine b for given mean-to-sigma ratio
      min_AUC_OR <- min(AUC1_OR,AUC2_OR) #desired median AUC across readers
      sigma <- 1 # sigma for DV for normal cases -- this is for a fixed reader
      r <- mean_sig_input # desired mean-to-sigma ratio
      a1_ <- qnorm(min_AUC_OR)^2 * (1/(1-x4))- r^2
      b1_ <- 2*r^2
      c1_ <- qnorm(min_AUC_OR)^2 *(1/(1-x4))- r^2
      discrim <- (b1_^2 - 4*a1_*c1_)^.5 # square root of the discriminant of the
      # quadratic equation for solving for 1/b
      b <- ((-b1_ - discrim)/(2*a1_))^-1
    }

    if (b_method == "unspecified") { # Determine b to result in inputted OR_var
      flag <- 0
      x <- 1
      b <- NA #   this value will indicate b_method 2 does not work for
      # these data
      if (b_le_1)
      { # if (b_le_1) then the algorithm does searches
        # for a value of b < 1. Otherwise, it searches only for b >=1
        while (flag == 0 & x > lower_bd) {
          # first assume increment1 = .001 <= b <= 1.  There may also be a solution for b > 1, but
          # we do not compute it if a solution exists for b <= 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 onsider that
          r1 <- 1
          r2 <- x^2/(1 + x^2)
          r3 <- 1/(1 + x^2)
          r4 <- 0
          a1 <- multnorm(x1,x1,r1*(1 - x4) + x4)
          a2 <- multnorm(x1,x1,r2*(1 - x4) + x4)
          a3 <- multnorm(x1,x1,r3*(1 - x4) + x4)
          a4 <- multnorm(x1,x1,r4*(1 - x4) + x4)
          var_error_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
          a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
          a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
          a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
          a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
          var_error_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
          var_error_formula <- (var_error_mod1 + var_error_mod2)/2
          if (x == 1) {
            sign <- sign(var_error_formula - var_error_OR)
          }
          sign1 <- sign(var_error_formula - var_error_OR)
          if (sign1 != sign) {
            flag <- 1
            b <- x
            temp_b <- b + increment1
          } else {
            x <- x - increment1
          }
        }

        if (!is.na(b)) { #now compute higher precision solution
          flag <- 0
          x <- temp_b
          while (flag == 0 & x >= temp_b - increment1) {
            r1 <- 1
            r2 <- x^2/(1 + x^2)
            r3 <- 1/(1 + x^2)
            r4 <- 0
            a1 <- multnorm(x1,x1,r1*(1 - x4) + x4)
            a2 <- multnorm(x1,x1,r2*(1 - x4) + x4)
            a3 <- multnorm(x1,x1,r3*(1 - x4) + x4)
            a4 <- multnorm(x1,x1,r4*(1 - x4) + x4)
            var_error_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
            a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
            a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
            a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
            a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
            var_error_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
            var_error_formula <- (var_error_mod1 + var_error_mod2)/2
            if (x >= (temp_b - increment2/2)) {
              sign <- sign(var_error_formula - var_error_OR)
            }
            sign1 <- sign(var_error_formula - var_error_OR)
            if (sign1 != sign) {
              flag <- 1
              b <- x
            } else {
              x <- x - increment2
            }
          }
        }
      } # end of if (b_ge_1 != "yes")

      if (is.na(b)) {
        x <- 1
        while (flag == 0 & x <= 4) {
          r1 <- 1
          r2 <- x**2/(1 + x**2)
          r3 <- 1/(1 + x**2)
          r4 <- 0
          a1 <- multnorm(x1,x1,r1*(1 - x4) + x4)
          a2 <- multnorm(x1,x1,r2*(1 - x4) + x4)
          a3 <- multnorm(x1,x1,r3*(1 - x4) + x4)
          a4 <- multnorm(x1,x1,r4*(1 - x4) + x4)
          var_error_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
          a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
          a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
          a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
          a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
          var_error_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
          var_error_formula <- (var_error_mod1 + var_error_mod2)/2

          if (x == 1) {
            sign <- sign(var_error_formula - var_error_OR)
          }
          sign1 <- sign(var_error_formula - var_error_OR)
          if (sign1 != sign) {
            flag <- 1
            b <- x
            temp_b <- b - increment1
          } else {
            x <- x + increment1
          }
        }

        if (!is.na(b)) { #now compute higher precision solution
          flag <- 0
          x <- temp_b
          while (flag == 0 & x <= temp_b + increment1) {
            r1 <- 1
            r2 <- x^2/(1 + x^2)
            r3 <- 1/(1 + x^2)
            r4 <- 0
            a1 <- multnorm(x1,x1,r1*(1 - x4) + x4)
            a2 <- multnorm(x1,x1,r2*(1 - x4) + x4)
            a3 <- multnorm(x1,x1,r3*(1 - x4) + x4)
            a4 <- multnorm(x1,x1,r4*(1 - x4) + x4)
            var_error_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
            a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
            a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
            a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
            a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
            var_error_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
            var_error_formula <- (var_error_mod1 + var_error_mod2)/2
            if (x <= (temp_b + increment2/2)) {
              sign <- sign(var_error_formula - var_error_OR)
            }
            sign1 <- sign(var_error_formula - var_error_OR)
            if (sign1 != sign) {
              flag <- 1
              b <- x
            } else {
              x <- x + increment2
            }
          }
        }




      }
    }


    if (b_method == "specified") { # use b_input as b value
      b <- b_input
    }
    if (is.na(b)) {
      b <- NA
      flag1 <- 1
    }

    if (flag1 == 0) {
      # Now compute covariances from inputted correlations and the previously
      # computed error variance.
      r1 <- 1
      r2 <- b^2/(1 + b^2)
      r3 <- 1/(1 + b^2)
      r4 <- 0
      a1 <- multnorm(x1,x1,r1*(1 - x4) + x4)
      a2 <- multnorm(x1,x1,r2*(1 - x4) + x4)
      a3 <- multnorm(x1,x1,r3*(1 - x4) + x4)
      a4 <- multnorm(x1,x1,r4*(1 - x4) + x4)
      var_error_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
      a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
      a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
      a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
      var_error_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      var_error_formula <- (var_error_mod1 + var_error_mod2)/2
      cov1 <- corr1_OR*var_error_formula
      cov2 <- corr2_OR*var_error_formula
      cov3 <- corr3_OR*var_error_formula

      # now compute mean-to-sigma ratios
      V <- (1 + 1/b^2)/(1-x4)
      delta1 <- x1*sqrt(V)
      delta2 <- x2*sqrt(V)
      meansig1 <- delta1/((1/b) - 1)
      meansig2 <- delta2/((1/b) - 1)
    }
  }
  # STEP 5: Solve for x5 = RM case variance
  if (flag1 == 1) {
    x5 <- NA
  }
  if (flag1 == 0) {
    c1 <- 1/(n0*n1)
    c2 <- (n1-1)/(n0*n1)
    c3 <- (n0-1)/(n0*n1)
    c4 <- (1-n0-n1)/(n0*n1)
    flag <- 0
    x <- .5
    for (step in 2:step_no) {
      r1 <- x
      r2 <- x/(1 + 1/b^2)
      r3 <- x/(1+b^2)
      r4 <- 0
      a1 <- multnorm(x1,x2,r1*(1 - x4))
      a2 <-  multnorm(x1,x2,r2*(1 - x4))
      a3 <-  multnorm(x1,x2,r3*(1 - x4))
      a4 <-  multnorm(x1,x2,r4*(1 - x4))
      cov3_formula <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      x <- x - sign(cov3_formula - cov3)*.5^step
    }
    x5 <- x
    if ((x5 > 1 - .5^(step_no - 1))) {
      x5 <- NA
      flag1 <- 1
    }
  }

  #  STEP 6: Solve for x6 = RM case + treatment-by-case variances

  if (flag1 == 1) {
    x6 <- NA
  }
  if (flag1 == 0) {
    x <- .5
    for (step in 2:step_no) {
      r1 <- x
      r2 <- x/(1 + 1/b^2)
      r3 <- x/(1+b^2)
      r4 <- 0
      a1 <- multnorm(x1,x1,r1*(1 - x4))
      a2 <- multnorm(x1,x1,r2*(1 - x4))
      a3 <- multnorm(x1,x1,r3*(1 - x4))
      a4 <- multnorm(x1,x1,r4*(1 - x4))
      cov2_formula_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      a1 <- multnorm(x2,x2,r1*(1 - x4))
      a2 <- multnorm(x2,x2,r2*(1 - x4))
      a3 <- multnorm(x2,x2,r3*(1 - x4))
      a4 <- multnorm(x2,x2,r4*(1 - x4))
      cov2_formula_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      cov2_formula <- .5*(cov2_formula_mod1 + cov2_formula_mod2)
      x <- x - sign(cov2_formula - cov2)*.5^step
    }
    x6 <- max(x, x5)
    if (x6 > 1 - .5^(step_no - 1))  {
      x6 <- NA
      flag1 <- 1
    }
  }


  # STEP 7: Compute x7 = RM case + reader-by-case variances*

  if (flag1 == 1) {
    x7 <- NA
  }
  if (flag1 == 0) {
    x <- .5
    for (step in 2:step_no) {
      r1 <- x
      r2 <- x/(1 + 1/b^2)
      r3 <- x/(1+b^2)
      r4 <- 0
      a1 <- multnorm(x1,x2,r1*(1 - x4) + x3)
      a2 <- multnorm(x1,x2,r2*(1 - x4) + x3)
      a3 <- multnorm(x1,x2,r3*(1 - x4) + x3)
      a4 <- multnorm(x1,x2,r4*(1 - x4) + x3)
      cov1_formula <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      x <- x - sign(cov1_formula - cov1)*.5**step
    }
    x7 <- min(max(x, x5), 1-x6+x5)
    if (x7 > 1 - .5^(step_no - 1)) {
      x7 <- NA
      flag1 <- 1
    }}

  #### Compute RM parameters from x1-x7 and b****
  V_RM <- (1 + 1/b^2)/(1 - x4)
  delta1_RM <- x1*sqrt(V_RM)
  delta2_RM <- x2*sqrt(V_RM)
  var_R_RM <- .5*x3*V_RM
  var_TR_RM <- .5*x4*V_RM - .5*x3*V_RM
  var_C_RM <- x5
  var_TC_RM <- x6-x5
  var_RC_RM <- x7 - x5
  var_error_RM <- 1 - (x6 + x7 -x5)
  r_temp <- delta1_RM/(1/b - 1)
  AUC1_pred_RM <- pnorm(delta1_RM/sqrt(V_RM))
  AUC2_pred_RM <- pnorm(delta2_RM/sqrt(V_RM))
  mean_to_sig1 <- delta1_RM/(1/b - 1)
  mean_to_sig2 <- delta2_RM/(1/b - 1)
  sigma_tilda = sqrt(2*var_R_RM + 2*var_TR_RM)/(1/b - 1);
  Pr1_improper <- pnorm((2-mean_to_sig1)/sigma_tilda) -
    pnorm((-3-mean_to_sig1)/sigma_tilda)
  Pr2_improper <- pnorm((2-mean_to_sig2)/sigma_tilda) -
    pnorm((-3-mean_to_sig2)/sigma_tilda)
  if (!is.na(b)) {
    if (b == 1) {
      Pr1_improper <- 0
      Pr2_improper <- 0
    }
  }

  res <- data.frame(
    n0 = n0,
    n1 = n1,
    delta1 = delta1_RM,
    delta2 = delta2_RM,
    var_R = var_R_RM,
    var_TR = var_TR_RM,
    var_C = var_C_RM,
    var_TC = var_TC_RM,
    var_RC = var_RC_RM,
    var_error = var_error_RM,
    b = b,
    b_method = b_method,
    mean_to_sig1 = mean_to_sig1,
    mean_to_sig2 = mean_to_sig2,
    Pr1_improper = Pr1_improper,
    Pr2_improper = Pr2_improper,
    x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, x7 = x7
  )

  NA_fixes <- c(
    "x3" = "Try reducing the value of var_R.",
    "x4" = "Try reducing the value of var_TR.",
    "b" = paste0(
      "If using b_method = \"unspecified\", there are two possible solutions:\n",
      "  a) Try changing (reduce or increase) the value of var_error.\n",
      "  b) Try using one of the other two b_method options, which should always work."
    ),
    "x5" = "Try reducing the value of corr3.",
    "x6" = "Try reducing the value of corr2.",
    "x7" = "Try reducing the value of corr1."
  )
  ind <- match(NA, as.numeric(res[names(NA_fixes)]))
  if (!is.na(ind)) {
    warning("Conversion failed.\n", NA_fixes[ind])
  }

  structure(res, class = c("RMparams", class(res)))

}


print.RMparams <- function(x, all = FALSE, ...) {
  if (!all) x <- x[setdiff(names(x), paste0("x", 1:7))]
  print(as.data.frame(x))
}


#' @rdname OR_to_RMH
#'
OR_to_RMH.data.frame <- function(params, ...) {
  f <- function(row) do.call(OR_to_RMH, c(as.list(row), ...))
  res <- by(params, seq_len(nrow(params)), f)
  names(res) <- rownames(params)
  do.call(rbind, res)
}


#' Convert Roe & Metz Parameters to Corresponding Obuchowski-Rockette Parameters
#'
#' Determines Obuchowski-Rockette (OR) model parameter values that describe the
#' distribution of  empirical AUC reader performance outcomes computed from
#' multireader multicase likelihood-of-disease rating data simulated from the
#' Roe & Metz (RM) simulation model, based on the analytical mapping provided by
#' Hillis (2018).  The function assumes the RM model proposed by (Hillis, 2012),
#' which generalizes the orginal RM model (Roe and Metz, 1997) by allowing the
#' latent confidence-of-disease rating distributions to have unequal
#' diseased-case and nondiseased-case variances, with the variance components
#' involving diseased cases constrained to differ by a factor of 1/b^2, b>0,
#' from corresponding variance components involing nondiseased cases.
#'
#' A related function is the OR_to_RMH function, which determines RM parameter
#' values for simulating multireader multicase confidence-of-disease rating data
#' based on real-data or conjectured Obuchowski-Rockette (OR) parameter
#' estimates that describe the distribution of the empirical AUC reader
#' performance measures.
#'
#' @rdname RMH_to_OR
#'
#' @param ... arguments passed to the default method.
#' @param n0,n1 numbers of nondiseased and diseased cases.
#' @param b b>0, with 1/b^2 = ratio of each diseased-case variance component to
#'   the corresponding diseased-case variance component. It follows that b is
#'   also the conventional binormal-curve slope, i.e., the slope of each
#'   reader's true ROC curve plotted in probit space.
#' @param delta1,delta2 test 1 and test 2 separations of the RM-model
#'   nondiseased and diseased latent likelihood-of-disease rating distribution
#'   means.
#' @param var_R,var_TR RM-model reader and test-by-reader variance components.
#' @param var_C,var_TC,var_RC,var_error RM-model case, test-by-case,
#'   reader-by-case and error variance components for nondiseased cases.
#' @param params data frame of above RM parameter values in the columns.
#'
#' @details
#'
#' Hillis (2012) modified the original RM model (Roe and Metz, 1997) by allowing
#' variance components involving case to depend on truth (diseased/nondiseased),
#' with variance components involving diseased cases set equal to those
#' involving nondiseased cases multiplied by the factor 1/b^2, b>0. Assuming
#' this model, Hillis (2018) derived analytical formulas expressing OR
#' parameters that describe the distribution of empirical AUC outcomes computed
#' from RM=model simulated data as functions of the RM parameters. This mapping
#' from the RM parameters to the OR parameters is implemented in R by the
#' RMH_to_OR function.
#'
#' @return
#' The OR model parameters are returned in a data frame with the following
#' elements.
#'
#' \describe{
#'   \item{...}{arguments passed to the default method.}
#'   \item{AUC1, AUC2}{test 1 and 2 expected empirical AUCs.}
#'   \item{var_R, var_TR}{OR reader and test-by-reader variance components.}
#'   \item{corr1, corr2, corr3}{OR error correlations.}
#'   \item{var_error}{OR error variance.}
#'   \item{n0, n1}{number of nondiseased and diseased cases.}
#' }
#'
#' Related quantities describing the true reader ROC curves that are also
#' returned in the data frame:
#'
#' \describe{
#'   \item{b}{b > 0, with 1/b^2 = (RM diseased variance component) /
#'     (corresponding RM nondiseased variance component).}
#'   \item{mean_to_sig1}{expected mean-to-sigma ratio across readers for test
#'     1.}
#'   \item{mean_to_sig2}{expected mean-to-sigma ratio across readers for test
#'     2.}
#'   \item{Pr1_improper}{probability that the test 1 ROC curve for a random
#'     reader will be visually improper (i.e, |mean-to-sigma ratio| < 3).}
#'   \item{Pr2_improper}{probability that the test 2 ROC curve for a random
#'     reader will be visually improper (i.e, |mean-to-sigma ratio| < 3).}
#' }
#'
#' @references
#' Hillis, Stephen L. 2012. "Simulation of Unequal-Variance Binormal Multireader
#' ROC Decision Data: An Extension of the Roe and Metz Simulation Model."
#' \emph{Academic Radiology} no. 19 (12):1518-1528.
#' doi: 10.1016/j.acra.2012.09.011.
#'
#' Hillis, Stephen L. 2018. "Relationship between Roe and Metz simulation model
#' for multireader diagnostic data and Obuchowski-Rockette model parameters."
#' \emph{Statistics in Medicine} no. 37 (13):2067-2093. doi: 10.1002/sim.7616.
#'
#' Roe, Cheryl A., and Charles E. Metz. 1997. "Dorfman-Berbaum-Metz method for
#' statistical analysis of multireader, multimodality receiver operating
#' characteristic data: Validation with computer simulation." \emph{Academic
#' Radiology} no. 4 (4):298-303. doi: 10.1016/S1076-6332(97)80032-3.
#'
#' @author
#' Stephen L. Hillis, Departments of Radiology and Biostatistics,
#' University of Iowa, \email{steve-hillis@uiowa.edu}
#'
#' Brian J. Smith, Department of Biostatistics, University of Iowa,
#' \email{brian-j-smith@uiowa.edu}
#'
#' @seealso \code{\link{OR_to_RMH}}
#'
#' @examples
#' ##  Example 1: Computing OR parameters from RM parameters directly
#' # RM parameters from first line (A_z = 0.702) of Table 1 in Roe & Metz (1997)
#' # with 50 diseased and 50 nondiseased cases.
#' OR <- RMH_to_OR(n0 = 50, n1 = 50, delta1 = 0.75, delta2 = 0.75,
#'                 var_R = 0.0055, var_TR = 0.0055, var_C = 0.3, var_TC = 0.3,
#'                 var_RC = 0.2, var_error = 0.2, b = 1)
#' OR
#'
#' ##  Example 2: Computing RM parameters from a data frame of RM parameters
#' ##  ---------------------------------------------------------------------
#' ## Example 2a:  RM parameters from first line (A_z = 0.702) of Table 1 in
#' # Roe & Metz (1997) with 50 diseased and 50 nondiseased cases
#' RM_parms_line1 <- data.frame(n0 = 50, n1 = 50, delta1 = 0.75, delta2 = 0.75,
#'                              var_R = 0.0055, var_TR = 0.0055, var_C = 0.3, var_TC = 0.3,
#'                              var_RC = 0.2, var_error = 0.2, b = 1)
#' OR <- RMH_to_OR(RM_parms_line1)
#' OR
#' ## Note below that applying the OR_to_RMH function to the above OR parameters
#' # results in the original RM parameters within rounding error:
#' check <- OR_to_RMH(OR)
#' check
#'
#' ## Example 2b: RM parameters from last 3 lines of Table 1 in Roe & Metz (1997)
#' # using 10 diseased and 25 nondiseased cases
#' RM_3_models <- data.frame(
#'   rbind(
#'     line6 = c(25, 10, 0.75, 0.75, 0.011, 0.011, 0.1, 0.1, 0.2, 0.6, 1),
#'     line7 = c(25, 10, 1.50, 1.50, 0.03, 0.03, 0.1, 0.1, 0.2, 0.6, 1),
#'     line8 = c(25, 10, 2.5, 2.5, 0.056, 0.056, 0.1, 0.1, 0.2, 0.6, 1)
#'   )
#' )
#' colnames(RM_3_models) <- c("n0", "n1", "delta1", "delta2", "var_R", "var_TR",
#'                            "var_C", "var_TC", "var_RC", "var_error", "b")
#' RM_3_models
#' OR_3_models <- RMH_to_OR(RM_3_models)
#' OR_3_models
#'
#' ## Example 2c: RM parameters from last 3 lines of Table 1 in Hillis (2012)
#' # using 10 diseased and 25 nondiseased cases
#' RM_3_models_Hillis <- data.frame(
#'   rbind(
#'     line6 = c(25, 25, 0.821, 0.821, 0.0132, 0.0132, 0.1, 0.1, 0.2, 0.6, 0.84566),
#'     line7 = c(25, 25, 1.831, 1.831, 0.0447, 0.0447, 0.1, 0.1, 0.2, 0.6, 0.71082),
#'     line8 = c(25, 25, 3.661, 3.611, 0.1201, 0.1201, 0.1, 0.1, 0.2, 0.6, 0.55140)
#'   )
#' )
#' colnames(RM_3_models_Hillis) <- c("n0", "n1", "delta1", "delta2", "var_R", "var_TR",
#'                                   "var_C", "var_TC", "var_RC", "var_error", "b")
#' RM_3_models_Hillis
#' OR_3_models_Hillis <- RMH_to_OR(RM_3_models_Hillis)
#' OR_3_models_Hillis
#'
RMH_to_OR <- function(...) {
  UseMethod("RMH_to_OR")
}


#' @rdname RMH_to_OR
#'
RMH_to_OR.default <- function(
  ..., n0, n1, b, delta1, delta2, var_R, var_TR, var_C, var_TC, var_RC,
  var_error
) {

  delta1_RM <- delta1
  delta2_RM <- delta2
  var_R_RM <- var_R
  var_TR_RM <- var_TR
  var_C_RM <- var_C
  var_TC_RM <- var_TC
  var_RC_RM <- var_RC
  var_error_RM <- var_error

  c1 <- 1/(n0*n1)
  c2 <- (n1-1)/(n0*n1)
  c3 <- (n0-1)/(n0*n1)
  c4 <- (1-n0-n1)/(n0*n1)

  # compute x1-x7
  V <- 1 + b^-2 + 2*(var_R_RM + var_TR_RM)
  x1 <- delta1_RM/sqrt(V)
  x2 <- delta2_RM/sqrt(V)
  x3 <- 2*var_R_RM/V
  x4 <- 2*(var_R_RM + var_TR_RM)/V
  x5 <- var_C_RM
  x6 <- var_TC_RM + var_C_RM
  x7 <- var_RC_RM + var_TC_RM + var_C_RM
  AUC1_OR <- pnorm(x1)
  AUC2_OR <- pnorm(x2)

  # compute reader and reader-by-test variance components#
  var_R_OR <- multnorm(x1,x2,x3) - AUC1_OR*AUC2_OR
  var_TR_OR <- .5*(multnorm(x1,x1,x4)-AUC1_OR^2 + multnorm(x2,x2,x4) -
                     AUC2_OR^2) - var_R_OR

  V_RM <- (1 + 1/b^2)/(1 - x4)
  delta1_RM <- x1*sqrt(V_RM)
  delta2_RM <- x2*sqrt(V_RM)
  a1 <- delta1_RM/(1/b)
  a2 <- delta2_RM/(1/b)
  mean_to_sig1 <- a1/(1- b)
  mean_to_sig2 <- a2/(1 - b)
  sigma_tilda <- sqrt(2*var_R_RM + 2*var_TR_RM)/(1/b - 1)
  Pr1_improper <- pnorm((2-mean_to_sig1)/sigma_tilda) -
    pnorm((-3-mean_to_sig1)/sigma_tilda)
  Pr2_improper <- pnorm((2-mean_to_sig2)/sigma_tilda) -
    pnorm((-3-mean_to_sig2)/sigma_tilda)
  if (b == 1) {
    Pr1_improper = 0
    Pr2_improper = 0
  }


  # compute Cov3 #
  r1 <- x5
  r2 <- x5/(1 + b^-2)
  r3 <- x5/(1 + b^2)
  r4 <- 0
  a1 <- multnorm(x1,x2,r1*(1 - x4))
  a2 <- multnorm(x1,x2,r2*(1 - x4))
  a3 <- multnorm(x1,x2,r3*(1 - x4))
  a4 <- multnorm(x1,x2,r4*(1 - x4))
  cov3_OR <- c1*a1 + c2*a2 + c3*a3 + c4*a4

  # compute Cov2 #
  r1 <- x6
  r2 <- x6/(1 + b^-2)
  r3 <- x6/(1 + b^2)
  r4 <- 0
  a1 <- multnorm(x1,x1,r1*(1 - x4))
  a2 <- multnorm(x1,x1,r2*(1 - x4))
  a3 <- multnorm(x1,x1,r3*(1 - x4))
  a4 <- multnorm(x1,x1,r4*(1 - x4))
  cov2_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
  a1 <- multnorm(x2,x2,r1*(1 - x4))
  a2 <- multnorm(x2,x2,r2*(1 - x4))
  a3 <- multnorm(x2,x2,r3*(1 - x4))
  a4 <- multnorm(x2,x2,r4*(1 - x4))
  cov2_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
  cov2_OR <- .5*(cov2_mod1 + cov2_mod2)

  # compute Cov1 #
  r1 <- x7 - x6 + x5
  r2 <- (x7 - x6 + x5)/(1 + b^-2)
  r3 <- (x7 - x6 + x5)/(1 + b^2)
  r4 <- 0
  a1 <- multnorm(x1,x2,(r1*(1 - x4) + x3))
  a2 <- multnorm(x1,x2,(r2*(1 - x4) + x3))
  a3 <- multnorm(x1,x2,(r3*(1 - x4) + x3))
  a4 <- multnorm(x1,x2,(r4*(1 - x4) + x3))
  cov1_OR <- c1*a1 + c2*a2 + c3*a3 + c4*a4

  # compute Error variance correlations between scores#
  r1 <- 1
  r2 <- 1/(1 + b^-2)
  r3 <- 1/(1 + b^2)
  r4 <- 0
  a1 <- multnorm(x1,x1,(r1*(1 - x4) + x4))
  a2 <- multnorm(x1,x1,(r2*(1 - x4) + x4))
  a3 <- multnorm(x1,x1,(r3*(1 - x4) + x4))
  a4 <- multnorm(x1,x1,(r4*(1 - x4) + x4))
  var_error_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
  a1 <- multnorm(x2,x2,(r1*(1 - x4) + x4))
  a2 <- multnorm(x2,x2,(r2*(1 - x4) + x4))
  a3 <- multnorm(x2,x2,(r3*(1 - x4) + x4))
  a4 <- multnorm(x2,x2,(r4*(1 - x4) + x4))
  var_error_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
  var_error_OR <- (var_error_mod1 + var_error_mod2)/2

  # compute correlations#
  corr1_OR <- cov1_OR/var_error_OR
  corr2_OR<- cov2_OR/var_error_OR
  corr3_OR <- cov3_OR/var_error_OR

  res <- data.frame(
    n0 = n0,
    n1 = n1,
    AUC1 = AUC1_OR,
    AUC2 = AUC2_OR,
    var_R = var_R_OR,
    var_TR = var_TR_OR,
    var_error = var_error_OR,
    cov1 = cov1_OR,
    cov2 = cov2_OR,
    cov3 = cov3_OR,
    corr1 = corr1_OR,
    corr2 = corr2_OR,
    corr3 = corr3_OR,
    b = b,
    mean_to_sig1 = mean_to_sig1,
    mean_to_sig2 = mean_to_sig2,
    Pr1_improper = Pr1_improper,
    Pr2_improper = Pr2_improper

  )
  structure(res, class = c("ORparams", class(res)))

}


#' @rdname RMH_to_OR
#'
RMH_to_OR.data.frame <- function(params, ...) {
  f <- function(row) do.call(RMH_to_OR, c(as.list(row), ...))
  res <- by(params, seq_len(nrow(params)), f)
  names(res) <- rownames(params)
  do.call(rbind, res)
}


multnorm <- function(x1, x2, x) {
  mvtnorm::pmvnorm(upper = c(x1, x2), corr = matrix(c(1, x, x, 1), 2))[1]
}
