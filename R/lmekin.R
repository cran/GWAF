lmekin <- function (fixed, data = parent.frame(), random, varlist = NULL,
    variance, sparse = c(20, 0.05), rescale = T, pdcheck = T,
    subset, weight, na.action)
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "data", "weights", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    if (missing(variance))
        theta <- NULL
    else theta <- variance
    reSt <- reStruct(random, REML = F, data = NULL)
    gform <- getGroupsFormula(reSt)
    if (is.null(gform)) {
        temp.fixed <- fixed
        gvars <- NULL
    } else {
        gvars <- all.vars(random)
        fvars <- all.vars(formula)
        gvars <- gvars[is.na(match(gvars, fvars))]
        temp.fixed <- paste(deparse(as.vector(fixed)), collapse = "")
        temp.fixed <- paste(temp.fixed, paste(gvars, collapse = "+"),sep = "+")
        temp.fixed <- as.formula(temp.fixed)
    }
    m$formula <- temp.fixed
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Terms <- terms(fixed)
    X <- model.matrix(Terms, m)
    Y <- model.extract(m, "response")
    n <- length(Y)
    weights <- model.extract(m, "weights")
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0) rep(0, n) else 
              if (tt == 1) m[[offset]] else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    ncluster <- length(gvars)
    if (ncluster == 0) stop("No grouping variables found")
    groups <- getGroups(m, gform)
    temp <- coxme.varcheck(ncluster, varlist, n, gvars, groups, sparse, rescale, pdcheck)
    varlist <- temp$varlist
    kindex <- temp$kindex
    ntheta <- temp$ntheta

    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i]==1) theta.names <- c(theta.names,gvars[i]) else theta.names <- c(theta.names,paste(gvars[i],1:ntheta[i],sep=""))
    }
    if (length(theta)==0) theta <- rep(0,sum(ntheta)) else if (length(theta)!=sum(ntheta)) stop("Wrong length for theta")
    names(theta) <- theta.names

    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i]==1) theta.names <- c(theta.names,gvars[i]) else theta.names <- c(theta.names,paste(gvars[i],1:ntheta[i],sep=""))
    }
    if (length(theta)==0) theta <- rep(0,sum(ntheta)) else if (length(theta)!=sum(ntheta)) stop("Wrong length for theta")
    names(theta) <- theta.names

    tindex <- which(theta == 0)
    if (ncluster > 1)
        stop("function can have only 1 random effect")
    varlist <- varlist[[1]]
    kindex <- kindex[, 1]
    if (max(kindex) != n)
        stop("The random effect must be 1 per subject")
    ntheta <- ntheta[1]
    kindex2 <- integer(n)
    kindex2[kindex] <- 1:n
    logfun <- function(itheta, X, Y, varlist, theta, tindex,
        center) {
        theta[tindex] <- exp(itheta)
        tkmat <- varlist[[1]]
        tkmat@blocks <- tkmat@blocks * theta[1]
        diag(tkmat) <- diag(tkmat + 1)
        if (length(varlist) > 1) {
            for (i in 2:length(varlist)) tkmat@blocks <- varlist[[i]]@blocks *
                theta[i] + tkmat@blocks
        }
        tkmat@blocks <- tkmat@blocks/tkmat@blocks[1]
        gk <- gchol(tkmat)
        newx <- solve(gk, X, full = FALSE)
        newy <- solve(gk, Y, full = FALSE)
        resid <- qr.resid(qr(newx), newy)
        n <- length(Y)
        loglik <- (n/2) * (log(mean(resid^2)) - center) + sum(log(diag(gk)))/2
        loglik
    }
    newX <- X[kindex2, ]
    newY <- as.vector(Y[kindex2])
    dimnames(newX) <- NULL
    if (length(tindex) > 0) {
        center <- log(mean((Y - mean(Y))^2))
        nfit <- optim(par = rep(-1, length(tindex)), logfun,
            method = "L-BFGS-B", lower = log(1e-05), X = newX,
            Y = newY, varlist = varlist, theta = theta, tindex = tindex,
            center = center)
        iter <- nfit$counts
        theta[tindex] <- exp(nfit$par)
    } else iter <- 0
    tkmat <- varlist[[1]]
    tkmat@blocks <- tkmat@blocks * theta[1]
    diag(tkmat) <- diag(tkmat + 1)
    if (length(varlist) > 1) {
        for (i in 2:length(varlist)) tkmat@blocks <- varlist[[i]]@blocks*theta[i]+tkmat@blocks
    }
    gk <- gchol(tkmat)
    xok <- as.matrix(solve(gk, newX, full = F))
    yok <- solve(gk, newY, full = FALSE)
    lfit <- lm(yok ~ 0 + xok)
    names(lfit$coefficients) <- dimnames(X)[[2]]
    ls <- summary(lfit)
    resid.var <- mean(lfit$residuals^2)
    theta <- c(theta * resid.var, resid.var)
    names(theta) <- c(theta.names, "resid")
    fitted <- c(X %*% lfit$coef)
    residuals <- Y - fitted
    frail <- residuals[kindex2]
    names(frail) <- groups
    fcoef <- lfit$coef
    call$fixed <- fixed
    call$random <- random
    fit <- list(coefficients = list(fixed = fcoef, random = frail),
        theta = theta, variance = ls$cov.unscaled * resid.var,
        ctable = ls$coefficients, residuals = residuals, fitted.values = fitted,
        effects = lfit$effects, rank = lfit$rank, assign = lfit$assign,
        df.residual = lfit$df.residual - length(theta), loglik = (-n/2) *
            (log(mean(lfit$residuals^2)) + 1 + log(2 * pi)) -
            sum(log(diag(gk)))/2, iter = iter, n = n, call = call,
        method = "ML")
    na.action <- attr(m, "na.action")
    if (length(na.action))
        fit$na.action <- na.action
    oldClass(fit) <- c("lmekin")
    fit
}

