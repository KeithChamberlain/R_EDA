##
# Imputation Methods
# Keith Chamberlain
##


impute.fun <- function(y, fun, ...){
    '
    uses fun to impute missing values
    '
    y[is.na(y)] <- fun(y[!is.na(y)], ...)
    y
}


impute.distribution <- function(y, distribution=rnorm, ...){
    '
    uses distribution to impute missing values.
    '
    distribution_of_data <- distribution(n=length(y), ...)
    ifelse(is.na(y), distribution_of_data, y)
}


impute.walk <- function(y, distribution, ...){
    '
    uses rolling sum of a distribution for a random walk
    '
    ifelse(is.na(y),
           cumsum(c(y[!is.na(y)][1], distribution(n=length(y), ...)[-1])), y)
}


impute.lm.smooth <- function(y, x, degree=6){
    '
    uses lm to impute missing values
    '
    errorCondition(message="use impute.fun.formula() instead")
    df <- data.frame(y, x)
    nas <- is.na(df$y)
    fit <- lm(y~poly(x, degree=degree), data=df, na.action=na.omit)
    df[nas, "y"] <- predict(fit, newdata=df[nas, ])
    df$y
}


impute.categorical <- function(y){
    '
    imputes categorical missingness
    '
    yy <- sort(table(y), decreasing=TRUE)
    yy <- yy[names(yy) != ""]
    probs <- yy / length(yy)
    y[is.na(y)] <- sample(x=names(yy), size=sum(is.na(y)), prob=probs, replace=TRUE)
    y
}


impute.fun.formula <- function(fun, formula, data, dv, ...){
    '
    Use fun and multivariate formula to impute dv missing values. Any missing
    values in independent variables will be silently passed to dv.
    '
    nas <- as.logical(rowSums(
        data.frame(lapply(colnames(data), function(x) is.na(data[, x])))))
    model <- fun(formula=as.formula(formula), data=data[!nas, ], ...)
    data[nas, dv] <- predict(model, newdata=data[nas, ])
    data
}


impute.fun.index <- function(fun, y, index, ...){
    '
    Use index to identify points to impute. NAs handled separately with
    parameter through `...`
    '
    y[index] <- fun(y[!index], ...)
}


impute.fun.formula.index <- function(fun, formula, data, dv, index_col, ...){
    '
    Use fun that has a predict method, and multivariate formula to impute dv
    values found in index_col. Missing independent variables values are
    silently passed to the dependent variable. Send na policy through `...`.
    '
    model <- fun(formula=as.formula(formula), data=data[!data[, index_col], ],
                 ...)
    data[data[, index_col], dv] <- predict(
        model, newdata=data[data[, index_col], ])
    data
}


impute.lin <- function(y, index){
    '
    Imputes numeric missingness by the previous and next non-missing values.
    '
    cls <- class(y) # Used to reset the class if mean changes it.
    # If the first value is set, grab the closest next one
    if(index[1]) y[1] <- y[!is.na(y)][1]
    # If the last value is set, grab the closest one before it
    if(index[length(index)]) y[length(index)] <- y[!is.na(y)][length(!is.na(y))]
    # For the rest find the next valid first and last value and take their mean
    for(i in 2:(length(y)-1)){
        if(index[i]){
            j <- k <- 1
            while(index[i-j] & (i-j) > 0) j <- j - 1
            while(index[i+k] & (i+k) < length(index)) k <- k + 1
            y[i] <- mean(c(y[i-j], y[i+k]))
        }
    }
    class(y) <- cls
    y
}
