#' Convert Obuchowski-Rockette to Roe & Metz Parameters
#'
#' Determines Roe & Metz (RM) simulation model parameters for simulating
#' multireader multicase likelihood-of-disease rating data based on real-data or
#' conjectured Obuchowski-Rockette (OR) parameter estimates that describe the
#' distribution of the empirical AUC reader performance measure. The algorithm
#' is for the constrained unequal-variance RM model (Hillis, 2012), which
#' generalizes the orginal RM model (Roe and Metz, 1997) by allowing the
#' diseased and nondiseased decision-variable distributions to have unequal
#' variances for each reader.  Parameter values for the the original RM model
#' can be computed by setting \code{b_method = 3} and \code{b_input = 1}.
#'
#' Hillis (2018) derived formulas expressing OR parameters that describe the
#' distribution of epirical AUC outcomes computed from RM simulated data as
#' functions of the RM parameters.  This algorithm provides the reverse
#' transformation that determines the corresponding RM parameters.  The
#' algorithm is described in Hillis (2020).
#'
#' @rdname OR_to_RM
#'
#' @param AUC1,AUC2 test 1 and 2 AUCs.
#' @param var_R,var_TR OR reader and test-by-reader variances.
#' @param corr1,corr2,corr3 OR fixed reader correlations.
#' @param error_var OR fixed reader error variance.
#' @param n0,n1 number of nondiseased and diseased cases.
#' @param precision precision for estimation of RM using numerical methods.
#' @param b_method method of estimating RM binormal b parameter (= 1, 2, or 3).
#' @param mean_to_sig mean-to-sigma ratio (Hillis & Berbaum, 2011), required
#'   only if \code{b_method = 2}.
#' @param b_input binormal b value, required only if \code{b_method = 3}.
#' @param params dataframe of OR parameter values in the columns or a list of
#' OR parameter values.
#' @param ... arguments passed to the default method.
#'
#' @details
#' \code{precision} indicates the numerical accuracy of the algorithm.
#' Typically, values less than 0.0001 will result in little difference.  For
#' final results we recommend using 0.00001, with smaller values (e.g., 0.001,
#' 0.0001) used for preliminary results to speed up the algorithm.
#'
#' \code{b_method} indicates the method for estimating the RM binormal \emph{b}
#' parameter. \code{b_method = 1} results in RM simulated data for which
#' expected OR estimates are equal to those specified in \code{params}.
#' \code{b_method = 2} results in RM parameters such that the median
#' mean-to-sigma ratio across readers for the RM binormal distributions
#' corresponding to the lowest OR AUC is equal to \code{mean-to_sig}.
#' \code{b_method = 3} allows the user to input the desired value of \emph{b}
#' (e.g., set \code{b_input = 1} to obtain original RM parameter estimates.
#' For \code{b_method = 2,3}, the simulated data empirical AUC estimate
#' distribution is specified by the parameter values specified in
#' \code{OR_params}, \emph{except} for error_var_OR. Thus for b_methods 2 & 3,
#' \code{error_var_OR} can be equal to \code{NA} or excluded from \code{params}.
#'
#' \code{mean_to_sig} is the inputted mean-to-sigma ratio needed for
#' \code{b_method = 2}.
#'
#' \code{b_input} is the inputted \emph{b} value needed for \code{b_method = 3}.
#'
#' The OR model includes the following constraints:  \code{corr1_OR >= corr3_OR}
#' and \code{corr2_OR >= corr3_OR}. These constraints are  included in the
#' algorithm by setting \code{corr2_OR = corr3_OR} if \code{corr2_OR < corr3_OR}
#' (and similarly for \code{corr1_OR} if \code{corr1_OR < corr3_OR}).
#'
#' @return
#' The RM model parameters are returned as a data frame with the following
#' elements.
#' \describe{
#'   \item{n0}{}
#'   \item{n1}{}
#'   \item{b}{}
#'   \item{mu1}{}
#'   \item{mu2}{}
#'   \item{var_R}{}
#'   \item{var_TR}{}
#'   \item{var_C}{}
#'   \item{var_TC}{test-by-case variance.}
#'   \item{var_RC}{reader-by-case variance.}
#'   \item{var_error}{}
#'   \item{b_method}{}
#' }
#'
#' @references
#' Hillis, Stephen. 2020. "Determining Roe and Metz model parameters for
#' simulating multireader multicase confidence-of-disease rating data based on
#' read-data or conjectured Obuchowski-Rockette parameter estimates." Vol.
#' 11316, SPIE Medical Imaging: SPIE. doi.org/10.1117/12.2550541
#'
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
#' Radiology} no. 4 (4):298-303. doi: doi: 10.1016/S1076-6332(97)80032-3.
#'
#' @author
#' Stephen L. Hillis, Departments of Radiology and Biostatistics,
#' University of Iowa, steve-hillis@uiowa.edu
#'
#' @examples
#' ## Example 1a: Computing RM parameters from OR parameters directly
#' RM <- OR_to_RM(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                var_R = 0.00154, var_TR = 0.000208,
#'                error_var = 0.00788)
#'                RM
#'
#'#' ## Example 1aa: Computing RM parms --identical test model
#' RM <- OR_to_RM(n0 = 69, n1 = 45, AUC1 = 0.919, AUC2 = 0.919,
#'                corr1 = 0.433, corr2 = (0.430 + .299)/2, corr3 = (0.430 + .299)/2,
#'                var_R = 0.00154, var_TR = 0.0,
#'                error_var = 0.000788)
#'                RM
#'
#'
#'  RM_to_OR(RM)
#'
#'
#' ## Example 1b: Computing RM parameters from a list of OR parameters
#' vandyke_OR <- list(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                    corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                    var_R = 0.00154, var_TR = 0.000208, error_var = 0.000788)
#' vandyke_RM <- OR_to_RM(vandyke_OR)
#' vandyke_RM
#' RM_to_OR(vandyke_RM) ##Brian -- does not work ##
#'
#'
#' vandyke_OR <- data.frame(n0 = 69, n1 = 45, AUC1 = 0.897, AUC2 = 0.941,
#'                    corr1 = 0.433, corr2 = 0.430, corr3 = 0.299,
#'                    var_R = 0.00154, var_TR = 0.000208, error_var = 0.000788)
#' vandyke_RM <- OR_to_RM(vandyke_OR)
#' vandyke_RM
#' RM_to_OR(vandyke_RM)
#'

#'
#' ## Example 1c: Computing RM parameters from a data frame of OR parameters
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
#'                                 "corr3", "var_R", "var_TR", "error_var")
#' three_studies_RM <- OR_to_RM(three_studies_OR)
#' three_studies_RM
#'
#' RM_to_OR(three_studies_RM)  ##BRIAN -- DOES NOT WORK##
OR_to_RM <- function(...) {
  UseMethod("OR_to_RM")
}


#' @rdname OR_to_RM
#'
OR_to_RM.default <- function(
  AUC1, AUC2, var_R, var_TR, corr1, corr2, corr3, error_var, n0, n1,
  b_method = "unspecified", mean_sig_input, b_input,
  ...
) {

  AUC1_OR <- AUC1
  AUC2_OR <- AUC2
  var_R_OR <- var_R
  var_TR_OR <- var_TR
  corr1_OR <- corr1
  corr2_OR <- corr2
  corr3_OR <- corr3
  error_var_OR <- error_var # can be NA if use b_method 2
  # increment <- precision
  increment1 = .001 #initial search step increment for b_method = 1;
  increment2 = .000001 #search step increment for b_method = 1 after using increment1 steps;
  lower_bd = 0.01 #lower bound for unspecified method b search
  step_no = 30 #number of times interval is split in 2 when searching for solution
  flag1 = 0 #becomes 1 when an x value (x1 -- x7) or b is NA
  b = NA; x1 = NA; x2 = NA; x3 = NA; x4 = NA; x5 = NA; x6 =  NA; x7 =NA

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
    error_var = (0 <= error_var_OR | is.na(error_var_OR)) #,
    # precision = (0 < precision & precision < .1)
  )
  invalid_params <- names(inbounds)[!inbounds]
  if (length(invalid_params)) {
    stop("OR_to_RM parameter values out of bound for ",
         toString(invalid_params), call. = FALSE)
  }

  x1 <- qnorm(AUC1_OR)
  x2 <- qnorm(AUC2_OR)
  c1 <- 1/(n0*n1)
  c2 <- (n1-1)/(n0*n1)
  c3 <- (n0-1)/(n0*n1)
  c4 <- (1-n0-n1)/(n0*n1)

  # STEP 2: Compute x3 (= 2*sigma_R_OR**2/V)
  x = .5;
  for (step in 2:step_no) {
    temp <- multnorm (x1,x2,x) -  pnorm(x1)*pnorm(x2)
    x <- x - sign(temp - var_R_OR)*.5^step;
  }
  x3 <-  x
  temp1 <- multnorm (x1,x2,1) -  pnorm(x1)*pnorm(x2)
  if (x3 > 1 - .5^(step_no - 1) )
  {x3 = NA
  flag1 = 1
  }
  #STEP 3: Compute x4 (= 2(sigma_R_OR**2 + sigma_TR_OR**2)/V)

  if (flag1 == 1) {x4 = NA}
  if (flag1 == 0){
    x = .5
    for (step in 2:step_no) {
      temp <- .5*(multnorm(x1,x1,x) - (pnorm(x1))**2 + multnorm(x2,x2,x) - (pnorm(x2))**2) - var_R_OR
      x <- x - sign(temp - var_TR_OR)*.5^step
    }
    x4 <- max(x, x3);
    if (x4 > 1 - .5**(step_no - 1) )
    {x4 = NA
    flag1 = 1
    }}


  # STEP 4 Compute b
  if (flag1 == 1) {b = NA;cov1 =NA; cov2 = NA; cov3 =NA }
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
      while (flag == 0 & x >= lower_bd) {
        # first assume increment1 = .001 <= b <= 1.  There may also be a solution for b > 1, but
        # we do compute it if a solution exists for b <= 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 onsider that
        r1 <- 1
        r2 <- x^2/(1 + x^2)
        r3 <- 1/(1 + x^2)
        r4 <- 0
        a1 <- multnorm(x1,x1,r1*(1 - x4) + x4)
        a2 <- multnorm(x1,x1,r2*(1 - x4) + x4)
        a3 <- multnorm(x1,x1,r3*(1 - x4) + x4)
        a4 <- multnorm(x1,x1,r4*(1 - x4) + x4)
        error_var_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
        a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
        a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
        a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
        a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
        error_var_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
        error_var_formula <- (error_var_mod1 + error_var_mod2)/2
        if (x == 1) {
          sign <- sign(error_var_formula - error_var_OR)
        }
        sign1 <- sign(error_var_formula - error_var_OR)
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
          error_var_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
          a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
          a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
          a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
          a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
          error_var_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
          error_var_formula <- (error_var_mod1 + error_var_mod2)/2
          if (x >= (temp_b - increment2/2)) {
            sign <- sign(error_var_formula - error_var_OR)
          }
          sign1 <- sign(error_var_formula - error_var_OR)
          if (sign1 != sign) {
            flag <- 1
            b <- x
          } else {
            x <- x - increment2
          }
        }
      }



      if (is.na(b)) {
        x <- 1
        while (flag == 0 & x <= 4) {
          r1 = 1
          r2 = x**2/(1 + x**2)
          r3 = 1/(1 + x**2)
          r4 = 0
          a1 = multnorm(x1,x1,r1*(1 - x4) + x4)
          a2 = multnorm(x1,x1,r2*(1 - x4) + x4)
          a3 = multnorm(x1,x1,r3*(1 - x4) + x4)
          a4 = multnorm(x1,x1,r4*(1 - x4) + x4)
          error_var_mod1 = c1*a1 + c2*a2 + c3*a3 + c4*a4
          a1 = multnorm(x2,x2,r1*(1 - x4) + x4)
          a2 = multnorm(x2,x2,r2*(1 - x4) + x4)
          a3 = multnorm(x2,x2,r3*(1 - x4) + x4)
          a4 = multnorm(x2,x2,r4*(1 - x4) + x4)
          error_var_mod2 = c1*a1 + c2*a2 + c3*a3 + c4*a4
          error_var_formula = (error_var_mod1 + error_var_mod2)/2

          if (x == 1) {
            sign = sign(error_var_formula - error_var_OR)
          }
          sign1 = sign(error_var_formula - error_var_OR)
          if (sign1 != sign) {
            flag =1
            b = x
          } else {
            x <- x + increment1
          }
        }
      }
    }


    if (b_method == "specified") { # use b_input as b value
      b <- b_input}
    #   print(b)
    if (is.na(b)) {
      b <- NA
      flag1 = 1}
    print(flag1)
    #   stop()

    if (flag1 ==0){
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
      error_var_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      a1 <- multnorm(x2,x2,r1*(1 - x4) + x4)
      a2 <- multnorm(x2,x2,r2*(1 - x4) + x4)
      a3 <- multnorm(x2,x2,r3*(1 - x4) + x4)
      a4 <- multnorm(x2,x2,r4*(1 - x4) + x4)
      error_var_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
      error_var_formula <- (error_var_mod1 + error_var_mod2)/2
      cov1 <- corr1_OR*error_var_formula
      cov2 <- corr2_OR*error_var_formula
      cov3 <- corr3_OR*error_var_formula

      # now compute mean-to-sigma ratios
      V <- (1 + 1/b^2)/(1-x4)
      delta1 <- x1*sqrt(V)
      delta2 <- x2*sqrt(V)
      meansig1 <- delta1/((1/b) - 1)
      meansig2 <- delta2/((1/b) - 1)
    }
  }
  # STEP 5: Solve for x5 = RM case variance
  if (flag1 == 1) {x5 =NA}
  if (flag1 == 0) {
    c1 <- 1/(n0*n1)
    c2 <- (n1-1)/(n0*n1)
    c3 <- (n0-1)/(n0*n1)
    c4 <- (1-n0-n1)/(n0*n1)
    flag <- 0
    x <- .5;
    for (step in 2:step_no) {
      r1 <- x;   r2 <- x/(1 + 1/b^2);  r3 <- x/(1+b^2);  r4 <- 0
      a1 <- multnorm(x1,x2,r1*(1 - x4))
      a2 <-  multnorm(x1,x2,r2*(1 - x4))
      a3 <-  multnorm(x1,x2,r3*(1 - x4))
      a4 <-  multnorm(x1,x2,r4*(1 - x4))
      cov3_formula <- c1*a1 + c2*a2 + c3*a3 + c4*a4  ;
      x <- x - sign(cov3_formula - cov3)*.5^step
    }
    x5 <- x
    if ((x5 > 1 - .5^(step_no - 1))) {
      x5 = NA
      flag1 == 1
    }}

  #  STEP 6: Solve for x6 = RM case + treatment-by-case variances

  if (flag1 == 1) {x6 = NA}
  if (flag1 == 0) {
    x <= .5
    for (step in 2:step_no) {
      r1 <- x;   r2 <- x/(1 + 1/b^2);  r3 <- x/(1+b^2);  r4 <- 0
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
      x6 = NA
      flag1 = 1
    }}


  # STEP 7: Compute x7 = RM case + reader-by-case variances*

  if (flag1 == 1) {x7 = NA}
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
      x7 = NA
      flag1 = 1
    }}

  #### Compute RM parameters from x1-x7 and b****
  V_RM <- (1 + 1/b^2)/(1 - x4)
  mu1_RM <- x1*sqrt(V_RM)
  mu2_RM <- x2*sqrt(V_RM)
  var_R_RM <- .5*x3*V_RM
  var_TR_RM <- .5*x4*V_RM - .5*x3*V_RM
  var_C_RM <- x5
  var_TC_RM <- x6-x5
  var_RC_RM <- x7 - x5
  var_error_RM <- 1 - (x6 + x7 -x5)
  r_temp <- mu1_RM/(1/b - 1)
  AUC1_pred_RM <- pnorm(mu1_RM/sqrt(V_RM))
  AUC2_pred_RM <- pnorm(mu2_RM/sqrt(V_RM))
  mean_sig1 <- mu1_RM/(1/b - 1)
  mean_sig2 <- mu2_RM/(1/b - 1)
  mean_sig1_025 <- max(0,mu1_RM - 1.96*sqrt(2*(var_R_RM + var_TR_RM)))/(b^-1 - 1)
  mean_sig2_025 <- max(0,mu2_RM - 1.96*sqrt(2*(var_R_RM + var_TR_RM)))/(b^-1 - 1)

  data.frame(
    n0 = n0,
    n1 = n1,
    mu1 = mu1_RM,
    mu2 = mu2_RM,
    var_R = var_R_RM,
    var_TR = var_TR_RM,
    var_C = var_C_RM,
    var_TC = var_TC_RM,
    var_RC = var_RC_RM,
    var_error = var_error_RM,
    b_method = b_method,
    mean_sig1 = mean_sig1,
    mean_sig2 = mean_sig2,
    mean_sig1_025 = mean_sig1_025,
    mean_sig2_025 = mean_sig2_025,
    x1 = x1, x2=x2,x3=x3,x4=x4,b=b, x5=x5,x6=x6,x7=x7
  )

}


#' @rdname OR_to_RM
#'
OR_to_RM.data.frame <- function(params, ...) {
  f <- function(row) do.call(OR_to_RM, c(as.list(row), ...))
  res <- by(params, seq_len(nrow(params)), f)
  names(res) <- rownames(params)
  do.call(rbind, res)
}


#' @rdname OR_to_RM
#'
OR_to_RM.list <- function(params, ...)
{
  do.call(OR_to_RM, c(params, ...))
}


RM_to_OR <- function(...) {
  UseMethod("RM_to_OR")
}


RM_to_OR.default <- function(
  n0, n1, b, mu1, mu2, var_R, var_TR, var_C, var_TC, var_RC, var_error, ...
) {

  mu1_RM <- mu1
  mu2_RM <- mu2
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
  x1 <- mu1_RM/sqrt(V)
  x2 <- mu2_RM/sqrt(V)
  x3 <- 2*var_R_RM/V
  x4 <- 2*(var_R_RM + var_TR_RM)/V
  x5 <- var_C_RM
  x6 <- var_TC_RM + var_C_RM
  x7 <- var_RC_RM + var_TC_RM + var_C_RM
  AUC1_OR <- pnorm(x1)
  AUC2_OR <- pnorm(x2)

  # compute reader and reader-by-test variance components#
  var_R_OR <- multnorm(x1,x2,x3) - AUC1_OR*AUC2_OR
  var_TR_OR <- .5*(multnorm(x1,x1,x4)-AUC1_OR^2 + multnorm(x2,x2,x4)-AUC2_OR^2) -var_R_OR

  V_RM <- (1 + 1/b^2)/(1 - x4)
  mu1_RM <- x1*sqrt(V_RM)
  mu2_RM <- x2*sqrt(V_RM)
  a1 <- mu1_RM/(1/b)
  a2 <- mu2_RM/(1/b)
  mean_to_sig1 <- a1/(1- b)
  mean_to_sig2 <- a2/(1 - b)
  mean_sig1_025 = max(0,mu1_RM - 1.96*sqrt(2*(var_R_RM + var_TR_RM)))/(b^-1 - 1)
  mean_sig2_025 = max(0,mu2_RM - 1.96*sqrt(2*(var_R_RM + var_TR_RM)))/(b^-1 - 1)

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
  error_var_mod1 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
  a1 <- multnorm(x2,x2,(r1*(1 - x4) + x4))
  a2 <- multnorm(x2,x2,(r2*(1 - x4) + x4))
  a3 <- multnorm(x2,x2,(r3*(1 - x4) + x4))
  a4 <- multnorm(x2,x2,(r4*(1 - x4) + x4))
  error_var_mod2 <- c1*a1 + c2*a2 + c3*a3 + c4*a4
  error_var_OR <- (error_var_mod1 + error_var_mod2)/2

  # compute correlations#
  corr1_OR <- cov1_OR/error_var_OR
  corr2_OR<- cov2_OR/error_var_OR
  corr3_OR <- cov3_OR/error_var_OR

  data.frame(
    n0 = n0,
    n1 = n1,
    AUC1 = AUC1_OR,
    AUC2 = AUC2_OR,
    var_R = var_R_OR,
    var_TR = var_TR_OR,
    error_var = error_var_OR,
    cov1 = cov1_OR,
    cov2 = cov2_OR,
    cov3 = cov3_OR,
    corr1 = corr1_OR,
    corr2 = corr2_OR,
    corr3 = corr3_OR,
    mean_to_sig1 = mean_to_sig1,
    mean_to_sig2 = mean_to_sig2,
    mean_sig1_025 = mean_sig1_025,
    mean_sig2_025 = mean_sig2_025
  )

}


RM_to_OR.data.frame <- function(params, ...) {
  f <- function(row) do.call(RM_to_OR, c(as.list(row), ...))
  res <- by(params, seq_len(nrow(params)), f)
  names(res) <- rownames(params)
  do.call(rbind, res)
}


RM_to_OR.list <- function(params, ...)
{
  do.call(RM_to_OR, c(params, ...))
}


multnorm <- function(x1, x2, x) {
  mvtnorm::pmvnorm(upper = c(x1, x2), corr = matrix(c(1, x, x, 1), 2))[1]
}
