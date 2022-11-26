###########################################################################
### Function borrowed from the R package 'ExclusiveLasso'               ###
###                                                                     ###
### Reference: https://github.com/DataSlingers/ExclusiveLasso           ###
###########################################################################

GLM_FAMILIES <- c(gaussian=0,
                  binomial=1,
                  poisson=2)

#' Fit a GLM with Exclusive Lasso Regularization
#'
#' Fit a generalized linear model via maximum penalized likelihood
#' using the exclusive lasso penalty. The regularization path is computed
#' along a grid of values for the regularization parameter (lambda).
#' The interface is intentionally similar to that of \code{\link[glmnet]{glmnet}} in
#' the package of the same name.
#' @references
#' Campbell, Frederick and Genevera I. Allen. "Within Group Variable Selection
#'     with the Exclusive Lasso". Electronic Journal of Statistics 11(2),
#'     pp.4220-4257. 2017. \doi{10.1214/17-EJS1317}

exclusive_lasso <- function(X, y, groups, family=c("gaussian", "binomial", "poisson"),
                            weights, offset, nlambda=100,
                            lambda.min.ratio=ifelse(nobs < nvars, 0.01, 1e-04),
                            lambda, standardize=TRUE, intercept=TRUE,
                            lower.limits = rep(-Inf, nvars),
                            upper.limits = rep(Inf, nvars),
                            thresh=1e-07, thresh_prox=thresh,
                            skip_df=FALSE,
                            algorithm=c("cd", "pg")){

    tic <- Sys.time()

    ####################
    ##
    ## Input Validation
    ##
    ####################

    nobs <- NROW(X);
    nvars <- NCOL(X);

    if(length(y) != nobs){
        stop(sQuote("NROW(X)"), " and ", sQuote("length(y)"), " must match.")
    }

    if(length(groups) != nvars){
        stop(sQuote("NCOL(X)"), " and ", sQuote("length(groups)"), " must match.")
    }

    if(anyNA(X) || anyNA(y)){
        stop(sQuote("exclusive_lasso"), " does not support missing data.")
    }

    if(!all(is.finite(X))){
        stop("All elements of ", sQuote("X"), " must be finite.")
    }

    if(!all(is.finite(y))){
        stop("All elements of ", sQuote("y"), " must be finite.")
    }

    ## Index groups from 0 to `num_unique(groups) - 1` to represent
    ## in a arma::ivec
    groups <- match(groups, unique(groups)) - 1

    family <- match.arg(family)

    if(family == "poisson"){
        if(any(y < 0)){
            stop(sQuote("y"), " must be non-negative for Poisson regression.")
        }
    }

    if(family == "binomial"){
        if(any(y < 0) || any(y > 1)){
            stop(sQuote("y"), " must be in [0, 1] for logistic regression.")
        }
    }

    if(missing(weights)){
        weights <- rep(1, nobs)
    }

    if(length(weights) != nobs){
        stop(sQuote("NROW(X)"), " and ", sQuote("length(weights)"), " must match.")
    }

    if(any(weights <= 0)){
        stop("Observation weights must be strictly positive.")
    }

    if(sum(weights) != nobs){
        weights <- weights * nobs / sum(weights)
        warning(sQuote("sum(weights)"), " is not equal to ", sQuote("NROW(X)."), " Renormalizing...")
    }

    if(missing(offset)){
        offset <- rep(0, nobs)
    }

    if(length(offset) != nobs){
        stop(sQuote("NROW(X)"), " and ", sQuote("length(offset)"), " must match.")
    }

    nlambda <- as.integer(nlambda)

    if((lambda.min.ratio <= 0) || (lambda.min.ratio >= 1)){
        stop(sQuote("lambda.min.ratio"), " must be in the interval (0, 1).")
    }

    if(standardize){
        ## FIXME -- This form of standardizing X isn't quite right with observation weights
        Xsc <- scale(X, center=TRUE, scale=TRUE)
        X_scale <- attr(Xsc, "scaled:scale", exact=TRUE)
        X_center <- attr(Xsc, "scaled:center", exact=TRUE)

        if(!all(is.finite(Xsc))){
            stop("Non-finite ", sQuote("X"), " found after standardization.")
        }
    } else {
        Xsc <- X
        X_scale <- rep(1, nvars)
        X_center <- rep(0, nvars)
    }

    if(missing(lambda)){
        lambda_max <- max(abs(crossprod(Xsc, y - offset - weighted.mean(y, weights) * intercept)/nobs))
        lambda <- logspace(lambda.min.ratio * lambda_max, lambda_max, length.out=nlambda)
    }

    if(length(lambda) < 1){
        stop("Must solve for at least one value of lambda.")
    }

    if(any(lambda <= 0)){
        stop("All values of ", sQuote("lambda"), " must be positive.")
    }

    if(is.unsorted(lambda)){
        warning("User-supplied ", sQuote("lambda"), " is not increasing. Sorting for maximum performance.")
        lambda <- sort(lambda)
    }

    if(thresh_prox <= 0){
        stop(sQuote("thresh_prox"), " must be positive.")
    }

    if(thresh <= 0){
        stop(sQuote("thresh"), " must be positive.")
    }

    if(any(is.na(upper.limits)) || any(is.nan(upper.limits))){
        stop(sQuote("upper.limits"), " should be either finite or +/-Inf.")
    }

    if(any(is.na(upper.limits)) || any(is.nan(lower.limits))){
        stop(sQuote("lower.limits"), " should be either finite or +/-Inf.")
    }

    if(length(lower.limits) == 1L){
        lower.limits <- rep(lower.limits, nvars)
    }

    if(length(upper.limits) == 1L){
        upper.limits <- rep(upper.limits, nvars)
    }

    if(length(upper.limits) != nvars){
        stop(sQuote("upper.limits"), " must be of length ", sQuote("NCOL(X)."))
    }

    if(length(lower.limits) != nvars){
        stop(sQuote("lower.limits"), " must be of length ", sQuote("NCOL(X)."))
    }

    if(any(upper.limits <= lower.limits)){
        stop(sQuote("upper.limits"), " must be strictly greater than ", sQuote("lower.limits."))
    }

    algorithm <- match.arg(algorithm)

    if((family == "gaussian") && getOption("ExclusiveLasso.gaussian_fast_path", TRUE)){
        if(algorithm == "cd"){
            res <- exclusive_lasso_gaussian_cd(X=Xsc, y=y, groups=groups,
                                               lambda=lambda, w=weights, o=offset,
                                               lower_bound=lower.limits, upper_bound=upper.limits,
                                               thresh=thresh, intercept=intercept)
        } else {
            res <- exclusive_lasso_gaussian_pg(X=Xsc, y=y, groups=groups,
                                               lambda=lambda, w=weights, o=offset,
                                               lower_bound=lower.limits, upper_bound=upper.limits,
                                               thresh_prox=thresh_prox, thresh=thresh,
                                               intercept=intercept)
        }
    } else {
        if(algorithm == "cd"){
            res <- exclusive_lasso_glm_cd(X=Xsc, y=y, groups=groups,
                                          lambda=lambda, w=weights, o=offset,
                                          family=GLM_FAMILIES[family],
                                          lower_bound=lower.limits, upper_bound=upper.limits,
                                          thresh=thresh, thresh_prox=thresh_prox,
                                          intercept=intercept)
        } else {
            res <- exclusive_lasso_glm_pg(X=Xsc, y=y, groups=groups,
                                          lambda=lambda, w=weights, o=offset,
                                          family=GLM_FAMILIES[family],
                                          lower_bound=lower.limits, upper_bound=upper.limits,
                                          thresh=thresh, thresh_prox=thresh_prox,
                                          intercept=intercept)
        }
    }

    ## Convert intercept to R vector (arma::vec => R column vector)
    res$intercept <- as.vector(res$intercept)

    ## Convert coefficients and intercept back to original scale
    if(standardize){
        ## To get back to the original X, we multiply by X_scale,
        ## so we divide beta to keep things on the same unit
        res$coef <- res$coef / X_scale
        if(intercept){
            ## Map back to original X (undo scale + center)
            ##
            ## We handled the scaling above, now we adjust for the
            ## centering of X: beta(X - colMeans(X)) = beta * X - beta * colMeans(X)
            ## To uncenter we add back in beta * colMeans(X), summed over all observations
            res$intercept <- res$intercept - colSums(res$coef * X_center)
        }
    }

    ## Degrees of freedom -- calculated using original scale matrix
    ## (though it shouldn't really matter)
    if(!skip_df){
        df <- calculate_exclusive_lasso_df(X, lambda, groups, res$coef)
    } else {
        df <- NULL
    }

    if(!is.null(colnames(X))){
        rownames(res$coef) <- colnames(X)
    }

    result <- list(coef=res$coef,
                   intercept=res$intercept,
                   y=y,
                   X=X,
                   standardize=standardize,
                   groups=groups,
                   lambda=lambda,
                   weights=weights,
                   offset=offset,
                   family=family,
                   df=df,
                   algorithm=algorithm,
                   nnz=apply(res$coef, 2, function(x) sum(x != 0)),
                   time=Sys.time() - tic)

    class(result) <- c("ExclusiveLassoFit", class(result))

    result
}

has_intercept <- function(x){
    !is.null(x$intercept)
}

has_offset <- function(x){
    any(x$offset != 0)
}

#' @export
print.ExclusiveLassoFit <- function(x, ..., indent=0){
    icat("Exclusive Lasso Fit", "\n", indent=indent)
    icat("-------------------", "\n", indent=indent)
    icat("\n", indent=indent)
    icat("N: ", NROW(x$X), ". P: ", NCOL(x$X), ".\n", sep="", indent=indent)
    icat(length(unique(x$groups)), "groups. Median size", median(table(x$groups)), "\n", indent=indent)
    icat("\n", indent=indent)
    icat("Grid:", length(x$lambda), "values of lambda. \n", indent=indent)
    icat("  Miniumum:", min(x$lambda), "\n", indent=indent)
    icat("  Maximum: ", max(x$lambda), "\n", indent=indent)
    if(!is.null(x$df)){
      icat("  Degrees of freedom: ", min(x$df), " --> ", max(x$df), "\n", indent=indent)
    }
    icat("  Number of selected variables:", min(x$nnz), " --> ", max(x$nnz), "\n", indent=indent)
    icat("\n", indent=indent)
    icat("Fit Options:\n", indent=indent)
    icat("  - Family:        ", capitalize_string(x$family), "\n", indent=indent)
    icat("  - Intercept:     ", has_intercept(x), "\n", indent=indent)
    icat("  - Standardize X: ", x$standardize, "\n", indent=indent)
    icat("  - Algorithm:     ", switch(x$algorithm, pg="Proximal Gradient", cd="Coordinate Descent"),
         "\n", indent=indent)
    icat("\n", indent=indent)
    icat("Time: ", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n", indent=indent)
    icat("\n", indent=indent)

    invisible(x)
}

# Refit exclussive lasso on new lambda grid
# Used internally by predict(exact=TRUE)
#' @importFrom utils modifyList


update_fit <- function(object, lambda, ...){
    ARGS <- list(X=object$X,
                 y=object$y,
                 groups=object$groups,
                 weights=object$weights,
                 offset=object$offset,
                 family=object$gamily,
                 standardize=object$standardize,
                 intercept=has_intercept(object),
                 lambda=lambda)

    ARGS <- modifyList(ARGS, list(...))

    do.call(exclusive_lasso, ARGS)
}


coef.ExclusiveLassoFit <- function(object, lambda=s, s=NULL,
                                   exact=FALSE, group_threshold=FALSE, ...){

    predict(object, lambda=lambda, type="coefficients",
            exact=exact, group_threshold=group_threshold, ...)
}

#' Make predictions using the exclusive lasso.
e.
#' @examples
#' n <- 200
#' p <- 500
#' groups <- rep(1:10, times=50)
#' beta <- numeric(p);
#' beta[1:10] <- 3
#'
#' X <- matrix(rnorm(n * p), ncol=p)
#' y <- X %*% beta + rnorm(n)
#'
#' exfit <- exclusive_lasso(X, y, groups)
#' coef(exfit, lambda=1)
#' predict(exfit, lambda=1, newx = -X)


predict.ExclusiveLassoFit <- function(object, newx, lambda=s, s=NULL,
                                      type=c("link", "response", "coefficients"),
                                      group_threshold=FALSE,
                                      exact=FALSE, offset, ...){
    type <- match.arg(type)

    ## Get coefficients first
    if(!is.null(lambda)){
        if(exact){
            object <- update_fit(object, lambda=lambda, ...)

            if(has_intercept(object)){
                int <- Matrix(object$intercept, nrow=1,
                              sparse=TRUE, dimnames=list("(Intercept)", NULL))
            } else {
                int <- Matrix(0, nrow=1, ncol=length(object$lambda),
                              sparse=TRUE, dimnames=list("(Intercept)", NULL))
            }

            coef <- rbind(int, object$coef)
        } else {
            if(has_intercept(object)){
                int <- Matrix(object$intercept, nrow=1,
                              sparse=TRUE, dimnames=list("(Intercept)", NULL))
            } else {
                int <- Matrix(0, nrow=1, ncol=length(object$lambda),
                              sparse=TRUE, dimnames=list("(Intercept)", NULL))
            }

            coef <- rbind(int, object$coef)
            lambda <- clamp(lambda, range=range(object$lambda))

            coef <- lambda_interp(coef,
                                  old_lambda=object$lambda,
                                  new_lambda=lambda)
        }
    } else {
        if(has_intercept(object)){
             int <- Matrix(object$intercept, nrow=1,
                           sparse=TRUE, dimnames=list("(Intercept)", NULL))
        } else {
             int <- Matrix(0, nrow=1, ncol=length(object$lambda),
                           sparse=TRUE, dimnames=list("(Intercept)", NULL))
        }

        coef <- rbind(int, object$coef)
    }

    if(group_threshold){
            coef[-1,,drop=FALSE] <- Matrix(apply(coef[-1,,drop=FALSE], 2, do_group_threshold, object$groups), sparse=TRUE)
    }

    if(type == "coefficients"){
        return(coef) ## Done
    }

    if(missing(newx)){
        link <- object$offset + cbind(1, object$X) %*% coef
    } else {
        if(missing(offset)){
            if(has_offset(object)){
                stop("Original fit had an offset term but", sQuote("offset"), "not supplied.")
            } else {
                offset <- rep(0, NROW(newx))
            }
        }
        link <- offset + cbind(1, newx) %*% coef
    }

    link <- as.matrix(link)

    if(type == "link"){
        link
    } else {
        ## Returning response
        switch(object$family,
               gaussian=link,
               binomial=inv_logit(link),
               poisson=exp(link))
    }
}

do_group_threshold <- function(x, groups){
    for(g in unique(groups)){
        g_ix <- (g == groups)
        x[g_ix] <- x[g_ix] * (abs(x[g_ix]) == max(abs(x[g_ix])))
    }
    x
}
