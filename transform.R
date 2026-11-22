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


#' @title Modified Z-Score
#'
#' @description
#' Function to calculate the modified z-score.
#'
#' @usage
#' `transform.modifiedz(y, na.rm=FALSE, constant=0.6744671045)
#'
#' @param y numeric vector to transform.
#' @param na.rm to remove NAs in calcuation but not in result.
#' @param constant to scale to sigma units.
#'
#' @returns A vector of the transformed values
#'
#' @details
#' Calculates and returns the modified z-score of the median deviated data
#' normalized by the median absolute deviation, and optionally scales by a
#' constant to resolve sigma units (e.g. +/-3sd).
#'
#' @references
#' Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute
#'     deviation. Journal of the American Statistical Association, 88(424),
#'     1273–1283.
#' Hubert, M., & Vandervieren, E. (2008). An adjusted boxplot for skewed
#'     distributions. Computational Statistics & Data Analysis, 52(12),
#'     5186–5201.
#'
#' @seealso [transform.z()] for the z-score. [transform.z.inverse()] for an
#' inverse to the z-score.
#'
#' @examples:
#' ## Generate income data bounded at 22000, with a mean of 70390 and sd of
#' ##    64000.
#' data <- qnorm(
#'     runif(n=10000, min=pnorm(q=22000, mean=70390, sd=64000), max=1),
#'     mean=70390, sd=64000)
#' hist(data, breaks=100, main="Histogram of Surrogate Income", xlab="Data")
#' ## modified z-score of Data
#' z <- transform.modifiedz(y=data)
#' ## Visualize
#' hist(z, breaks=100, main="Histogram of Modified-z Scores", xlab="Modified-z")
#'
#' @export
transform.modifiedz <- function(y, na.rm=FALSE, constant=0.6744671045){
    '
    transform y to the modified z-score.
    '
    ((y-median(y, na.rm=na.rm))/
         median(abs(y-median(y, na.rm=na.rm))) / constant)
}


#' @title Double Modified Z-Score Transform
#'
#' @description
#' Function to calculate the double modified z-score.
#'
#' @usage
#' `transform.2modifiedz(y, na.rm=FALSE, constant=0.6744671045)
#'
#' @param y numeric vector to transform.
#' @param na.rm to remove NAs in calculation but not in result.
#' @param constant to scale to sigma units or another value. Unscaled = 1.
#'
#' @returns A vector of the transformed values
#'
#' @details
#' Calculates and returns the double modified z-score of the median deviated
#' data normalized by the median absolute deviation of the lower and upper half,
#' and optionally scales by a constant to resolve sigma units (e.g. +/-3sd).
#'
#' @references
#' Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute
#'     deviation. Journal of the American Statistical Association, 88(424),
#'     1273–1283.
#' Hubert, M., & Vandervieren, E. (2008). An adjusted boxplot for skewed
#'     distributions. Computational Statistics & Data Analysis, 52(12),
#'     5186–5201.
#'
#' @seealso [transform.z()] for the z-score. [transform.z.inverse()] for an
#' inverse to the z-score.
#'
#' @examples:
#' ## Generate skewed 22000 bounded income data with mean of 70390 and sd of
#' ##     64000.
#' #' data <- qnorm(
#'     runif(n=10000, min=pnorm(q=22000, mean=70390, sd=64000), max=1),
#'     mean=70390, sd=64000)
#' hist(data, breaks=100, main="Histogram of Surrogate Income", xlab="Data")
#' ## modified z-score of Data
#' z <- transform.2modifiedz(y=data)
#' ## Visualize
#' hist(z, breaks=100, main="Histogram of Modified-z Scores", xlab="Modified-z from Each Side")
#'
#' @export
transform.2modifiedz <- function(y, na.rm=FALSE, constant=0.6744671045){
    '
    transform y to an adapting double modified-z score
    '
    deviation <- y-median(y, na.rm=na.rm)
    qs <- quantile(deviation, prob=c(0.25, 0.5, 0.75), na.rm=na.rm)
    names(qs) <- c("1/4", "1/2", "3/4")
    deviation / c(
        median(abs(deviation[deviation < qs["1/2"]])) / constant,
        median(abs(deviation[deviation >= qs["1/2"]])) / constant
    )
}


#' @title Finds the modes of data
#'
#' @description
#' Function to calculate the modes of data.
#'
#' @param values the values to find the peaks or valleys in.
#' @param sign whether to invert the search and find the valleys.
#' @param logic logical; whether to return the logical or numeric indices identifying the modes.
#'
#' @returns an integer or logical vector of indices for the identified values.
#'
#' @details
#' Calculates and returns a logical `logic=TRUE` or numeric `logic=FALSE` vector
#' of indices of the peaks or valleys in data depending on `sign`. May be slow
#' if there are many, so it is recommended to smooth first.
#'
#' @examples
#' ## Generate synthetic income data with a lower bound.
#' data <- qnorm(
#'    runif(n=10000, min=pnorm(q=22000, mean=70390, sd=64000), max=1),
#'    mean=70390, sd=64000)
#' ## Grab the smoothed distribution
#' smoothed_distribution <- density(data, adjust = 2)
#' ## Find the distributions modes
#' utility.modes(smoothed_distribution$y, logic=FALSE)
#'
utility.modes <- function(values, sign=1, logic=TRUE){
    '
    Indices of local maxima (sign = 1) or minima (sign = -1), edges included.
    For use with `density(data)$y`, for example.
    '
    modes <- logical(length=length(values))
    values <- sign*values
    for(i in 2:(length(modes)-1)){
        # Usual case
        if(values[i] > values[i-1] & values[i] > values[i+1]) modes[i] <- TRUE
        # First value larger than next
        if(values[i-1] > values[i] & (i-1) == 1) modes[i] <- TRUE
        # Last value larger than previous
        if(values[i+1] > values[i] & (i+1)==length(modes)) modes[i] <- TRUE
    }
    if(logic) modes else which(modes)
}


#' @title N-Modified Z-Score Transform
#'
#' @description
#' Function to calculate the modified z-score for n modes.
#'
#' @usage
#' `transform.nmodifiedz(y, na.rm=FALSE, constant=0.6744671045)
#'
#' @param y numeric vector to transform.
#' @param na.rm to remove NAs in calculation but not in result.
#' @param constant to scale to sigma units or another value. Unscaled = 1.
#'
#' @returns A vector of the transformed values
#'
#' @details
#' Calculates and returns the n-modified z-score of the median deviated
#' data normalized by the median absolute deviation of the lower and upper half,
#' of each mode, and optionally scales by a constant to resolve sigma units.
#'
#' @seealso [transform.z()] for the z-score. [transform.z.inverse()] for an
#' inverse to the z-score. [transform.modifiedz()] for the modified z score,
#' [transform.2modifiedz()] for the modified z score from each side.
#'
#' @examples:
#' ## Generate skewed 22000 bounded income data with mean of 70390 and sd of
#' ##     64000.
#' data1 <- qnorm(
#'    runif(n=10000, min=pnorm(q=22000, mean=70390, sd=64000), max=1),
#'    mean=70390, sd=64000)
#' data2 <- qnorm(
#'     runif(n=10000, min=0, max=pnorm(q=400000, mean=300000, sd=80000)),
#'     mean=400000, sd=80000)
#' data <- c(data1, data2)
#' hist(data, breaks=100, main="Histogram of Surrogate Income", xlab="Data")
#' ## modified z-score of Data
#' z <- transform.nmodifiedz(y=data)
#' ## Visualize
#' hist(z, breaks=100, main="Histogram of Modified-z Scores", xlab="Modified-z from Each Side of Each Mode")
#'
#' @export
transform.nmodifiedz <- function(y, ..., constant=0.6744671045, na.rm=FALSE){
    '
    Calculate the modified z-score for each mode, if there are fewer maxima
    than minima.
    '
    ## Density
    density_values <- density(y, ..., na.rm=na.rm)
    ## Indexes of Maxima and Minima
    maxima_indexes <- utility.modes(density_values$y, sign=1, logic=FALSE)
    minima_indexes <- utility.modes(density_values$y, sign=-1, logic=FALSE)
    nmode_z <- rep(NA, length=length(density_values$y))
    sections <- findInterval(data, density_values$x[minima_indexes])
    for(i in 1:length(maxima_indexes)){
        center <- density_values$x[maxima_indexes[i]]
        deviation <- data[sections==i] - center
        nmode_z[sections==i] <- deviation * constant / c(
            median(abs(deviation[deviation < 0])),
            median(abs(deviation[deviation >=0]))
        )
    }
    nmode_z
}

