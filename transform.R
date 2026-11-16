##
# Transform Functions
# Keith Chamberlain
##


#' @title Data Transformations
#'
#' @description
#' Functions to perform data transformations.
#'
#' @usage
#' `score(y, fun=mean, spread_agg=sum, power=2, var_power=1/2, ladder=I, ...)`
#'
#' @param y numeric vector to transform.
#' @param fun central tendency function.
#' @param spread_agg function used with spread determination.
#' @param power power to apply to deviation definition.
#' @param var_power power to apply to the spread definition.
#' @param ladder initial data transformation to apply.
#' @param ... Other options to control the function fun.
#'
#' @value A vector of the transformed values
#'
#' @references
#' Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions,
#'    Vol. 1, Ch. 14 (Wiley, 1994).
#'
#' @seealso [base::scale()] for scaling and centering matrices.
#'
#' @examples:
#'
#' ## Generate Normal and Zero Bounded Data
#' data <- rnorm(10000, mean=3, sd=5)
#' log_scale <- method_of_moments_log(mean=8, sd=4)
#' data2 <- exp(rnorm(10000, mean=log_scale['mu'], sd=log_scale['sigma']))
#' x <- 1:length(data)
#'
#' ## z-transform of Data
#' z <- score(y=data, fun=mean, spread_agg=sum, power=2, var_power=1/2,
#'            ladder=I, ...)
#'
#' ## Scale using Median Deviation and Variance.
#' m <- score(y=data, fun=median, spread_agg=sum, power=1, var_power=1)
#'
#' ##
#' ## Visualize
#' plot(z~x, type="p")
#' plot(m~x, type="p")
#' hist(data2)
#'
#' @export


score <- function(y, fun=mean, spread_agg=sum, power=2,
                  var_power=1/2, ladder=I, ...){
    '
    calculate a transform of y.
    - z-score: fun=mean, spread_agg=sum, power=2, and var_power is 1/2
    - a median-score: fun=median, spread_agg=sum, power=1, var_power=1
    '
    y<-ladder(y)
    (y-fun(y, ...))/(
        spread_agg(abs(y-fun(y, ...))^power)/(length(y)-1))^var_power
}


method_of_moments_log <- function(mean, sd){
    '
    method of moments transformation between mean and sd to log scale.
    '
    sigma=sqrt(log(1+sd^2/mean^2))
    mu=log(mean) - sigma^2/2
    c("sigma"=sigma, "mu"=mu)
}


non_saturated <- function(fun=tanh, inv=atanh, unique=2^8, ...){
    '
    an estimate of the unsaturated range for the transform
    '
    y <- fun(seq(-12, 12, 12/unique), ...)
    index <- y > min(y)+2/unique & y < max(y)-2/unique
    inv_y <- inv(y[index], ...)
    c("min"=min(inv_y, ...), "max"=max(inv_y, ...))
}


transform.z <- function(y, ...){
    '
    tranform to z-score for y
    '
    (y-mean(y, ...))/sd(y, ...)
}


transform.z.restore <- function(y, mu, sigma, ...){
    '
    inverse transform the z-transform
    '
    y*sigma+mu
}


transform.minmax <- function(y, ...){
    '
    transform to range of [0, 1]
    '
    (y-min(y, ...))/(max(y, ...)-min(y, ...))
}


transform.minmax.restore <- function(y, new_min, new_max, ...){
    '
    restore to range of y
    '
    y * (new_max-new_min)+new_min
}


transform.activation <- function(y, fun=tanh, inv=atanh, unique=2^8, ...){
    '
    transform to range of [new_min, new_max] based on activation limits and
    unique bit depth
    '
    lims <- non_saturated(fun=fun, inv=inv, unique=unique, ...)
    fun(transform.minmax.restore(transform.minmax(y, ...),
                      new_min=lims["min"], new_max=lims["max"], ...))
}


transform.activation.restore <- function(y, fun=atanh, new_min, new_max, ...){
    '
    transform to range of [new_min, new_max] based on fun and new limits.
    used as inverse transform to transform.activation
    '
    transform.minmax.restore(transform.minmax(fun(y), ...),
                             new_min=new_min, new_max=new_max, ...)
}


transform.logit <- function(y, base=exp(1)){
    '
    transform from probability to infinite range
    '
    log(y/(1-y), base=base)
}


transform.sigmoid <- function(y, base=exp(1)){
    '
    transform from infinite range to probabilities
    '
    1/(1+base^(-y))
}
