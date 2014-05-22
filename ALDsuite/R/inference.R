# inference.R
# functions to infer statistical significance to admixture estimates
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created December 18, 2012
# Last Modified May 22, 2014

mald <- function(adm, data, formula = NULL, pop = 1, family = 'binomial', keep.models = FALSE)
{
    # default formula
    if(is.null(formula))
        formula <- formula('case ~ g')

    # check that names match
    if(is.null(rownames(data)))
    {
        warning("Rownames of data undefined. Assuming data is properly sorted")
    }else{
        if(any(rownames(data) != dimnames(adm$gammas)[[1]]))
            stop("Row names of data seem to be incorrectly sorted or labeled")
    }

    # this is how many markers we are working with
    J <- dim(adm$gammas)[2]

    # if we are keeping all model objects, run with lapply
    models <- lapply(1:J, mald.j, formula = formula, pop = pop, adm = adm,
                     data = data, family = family, keep.models = keep.models)

    # return results
    return(models)
}

mald.j <- function(j, adm, data, formula, pop, family, keep.models)
{
    data$g <- apply(adm$gammas[,j,,pop], 1, sum)

    model <- glm(formula, data = data, family = family)

    if(keep.models)
        return(summary(model))

    return(list(coefficients = summary(model)$coefficients))
}

lgs.case.only <- function(adm, phi1, phi2 = phi1^2, phi0 = 1, pop = 1, cases = NULL, dev = FALSE,
                          phased = FALSE)
{
    warning('lgs.case.only() is not currently working properly. All output is suspect.')
    # checks
    if(is.null(cases))
        cases <- rownames(adm$A0)

    # get number of populations we are dealing with
    npops <- dim(adm$P)[2]

    # convert pop to number if label given
    if(is.character(pop)) ######### check that we have a valid population!
        pop <- which(dimnames(adm$P)[[2]] == pop)

    if(!phased)
    {
        stop('Unphased analysis not yet implemented')
        # get names of Aj that we want
        probs2 <- which(dimnames(adm$gammas)[[3]] == paste('g', pop, pop, sep = ''))

        probs1 <- grep(pop, dimnames(adm$gammas)[[3]])
        probs1 <- probs1[probs1 != probs2] # probs2 is of length 1

        probs0 <- which(!1:length(dimnames(adm$gammas)[[3]]) %in% c(probs1, probs2))

        ### numerator ###
        if(length(probs0) == 1)
        {
            p0 <- adm$gammas[cases,,probs0] # probability of 0 risk alleles
        }else{
            p0 <- apply(adm$gammas[cases,,probs0], 1:2, sum)
        }

        if(length(probs1) == 1)
        {
            p1 <- adm$gammas[cases,,probs1] # probability of 1 risk allele
        }else{
            p1 <- apply(adm$gammas[cases,,probs1], 1:2, sum)
        }

        # always of length 1
        p2 <- adm$gammas[cases,,probs2] # probability of 2 risk alleles
    }else{
        p2 <- adm$gammas[cases,,1,pop] * adm$gammas[cases,,2,pop]
        p1 <- adm$gammas[cases,,1,pop] * (1 - adm$gammas[cases,,2,pop]) +
              (1 - adm$gammas[cases,,1,pop]) * adm$gammas[cases,,2,pop]
        p0 <- (1 - adm$gammas[cases,,1,pop]) * (1 - adm$gammas[cases,,2,pop])
    }

    # product of the wighted, summed probabilities
    num <- apply(log(p0*phi0 + p1*phi1 + p2*phi2), 2, sum)

    ### denominator ###
    den <- sum(log((1 - adm$A0[cases,pop])^2 * phi0 +
                   2 * adm$A0[cases,pop] * (1 - adm$A0[cases,pop]) * phi1 +
                   adm$A0[cases,pop]^2 * phi2))

    ### likelihood ###
    lod <- num - den

    return(lod)
}

lgs <- function(formula, adm, phi1, phi0 = 1, phi2 = phi1^2, pop = 1, covars = NULL, family = binomial)
{
    warning('lgs() has not been properly tested. All output is suspect.')
    ##### checks & setup #####

    # to do:
    # check that adm and covars are ordered the same (by individual)

    ## formula checks & parsing ##

    # if no left hand side, run lgs.case.only --- must no have covars! i.e. formula must be "~ 1"
    ## if(length(formula) == 2)
    ## {
    ##     if(formula == formula(~ 1))
    ##         return(lgs.case.only(adm, phi0 = phi0, phi1 = phi1, phi2 = phi2, pop = pop))

    ##     print(formula)
    ##     stop("The intercept-only formula (i.e. '~ 1') is the only acceptible one sided option.")
    ## }

    # family check - most of this section is taken from glm()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())

    if (is.function(family))
        family <- family()

    if (is.null(family$family))
    {
        print(family)
        stop("'family' not recognized")
    }

    # formula parsing
    if(is.null(covars))
        covars <- environment(formula)

    covars$.A <- adm$final$A0[, pop] # this assumes indivdiuals are ordered the same as in covars!!!!!!!
    covars$.Ak <- apply(adm$final$Ak[,,pop], 1, prod)

    tmp <- attributes(terms(formula))
    ## formula <- formula(paste(attributes(tmp$factors)$dimnames[[1]][1], ' ~ ',
    ##                          paste(c('.A', tmp$term.labels), collapse = '+')))

    if(tmp$response)
    {
        formula <- reformulate(termlabels = c('.A', '.Ak', tmp$term.labels),
                               response = as.character(tmp$variables)[2],
                               intercept = tmp$intercept)
    }else{
        covars$case.only <- rep(1, length(covars$.A))
        formula <- reformulate(termlabels = c('.A', '.Ak', tmp$term.labels),
                               response = "case.only",
                               intercept = 0)
    }


    ## population parsing ##

    # total number of populations
    npops <- dim(adm$final$P)[2]

    # convert pop to number if label given
    if(is.character(pop))
        pop <- which(dimnames(adm$final$P)[[2]] == pop)

    # sort names of Aj into appropriate groups
    probs2 <- which(dimnames(adm$final$Aj)[[3]] == paste('g', pop, pop, sep = ''))

    probs1 <- grep(pop, dimnames(adm$final$Aj)[[3]])
    probs1 <- probs1[probs1 != probs2] # probs2 is of length 1

    probs0 <- which(!1:length(dimnames(adm$final$Aj)[[3]]) %in% c(probs1, probs2))


    # get probability totals
    if(length(probs0) == 1)
    {
        p0 <- adm$final$Aj[,,probs0] # probability of 0 risk alleles
    }else{
        p0 <- apply(adm$final$Aj[,,probs0], 1:2, sum)
    }

    if(length(probs1) == 1)
    {
        p1 <- adm$final$Aj[,,probs1] # probability of 1 risk allele
    }else{
        p1 <- apply(adm$final$Aj[,,probs1], 1:2, sum)
    }

    # always of length 1
    p2 <- adm$final$Aj[,,probs2] # probability of 2 risk alleles


    ##### null model #####

    Ho <- glm(formula, family = family, data = covars)
    if(family$family != 'binomial')
    {
        e <- Ho$residuals
        e.sd <- sd(e)

        # P(y | Ho) --- denominator
        logPy.Ho <- sum(log(dnorm(e, sd = e.sd)))
    }else{ # binomial
        logPy.Ho <- sum(log(Ho$fitted.values))
    }

    ##### Alternate model #####

    phi <- log(p0*phi0 + p1*phi1 + p2*phi2)
    phibar <- apply(phi, 2, mean)

    if(family$family != 'binomial')
    {
        # put e and phibar into matrices with similar dimensions as phi
        e.mat <- matrix(rep(e, length(phibar)), nrow = length(e)) # each row is identical
        phibar.mat <- matrix(rep(phibar, each = length(e)), nrow = length(e)) # each column is identical

        e.alternate <- e.mat - phibar.mat + phi

        logPy.Ha <- apply(log(dnorm(e.alternate, sd = e.sd)), 2, sum)
    }else{
        # fitted probabilities / odds
        p <- Ho$fitted.values
        lod <- log(p) - log(1 - p)

        # set up matricies for addition
        lod.mat <- matrix(rep(lod, length(phibar)), nrow = length(lod)) # each row is identical
        phibar.mat <- matrix(rep(phibar, each = length(lod)), nrow = length(lod)) # each column is identical

        lod.alternate <- lod - phibar.mat + phi
        p.alternate <- exp(lod.alternate) / (1 + exp(lod.alternate))

        logPy.Ha <- apply(log(p.alternate), 2, sum)
    }

    ### Bayes factor ###
    Bf <- logPy.Ha - logPy.Ho

    return(Bf)
}
