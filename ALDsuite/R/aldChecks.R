# aldChecks.R
# Quality control/error checks for admixture()
# Randall Johnson and Emily Mulhern
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# Leidos Biomedical Research, Inc
# Created June 23, 2014
# Last Modified October 7, 2014


###################################
# Original R functions for checks #
###################################

##### be sure all variables are formatted properly #####
######## check for quality of genetic data, too ########

ald.qc <- function(Pm.prior, haps, geno, gender, chr, pos, burn, iter, every, indiv.id,
                   marker.id, pop.id, lambda, tau, omega, rand.seed, cores, cl, dev,
                   verbose, debug, sex.chr, indiv.geno.check, marker.geno.check,
                   hwe.thresh, bad.indiv, bad.marker, male.het.X, fast, A0, Ak)
{
######### make sure we have things defined #########
    if(!is.null(haps) & is.null(dimnames(haps)[[2]]))
        stop("Missing marker names for haps.")

    if(is.null(dimnames(geno)[[2]]))
        stop("Missing marker names for geno.")

    if(is.null(indiv.id))
        indiv.id <- dimnames(geno)[[1]]

    if(is.null(marker.id))
        marker.id <- names(Pm.prior)

    if(is.null(pop.id))
        stop("pop.id is missing")


######### Type checks #########

    # character checks
    if(is.factor(indiv.id))
        indiv.id <- as.character(indiv.id)

    if(is.factor(marker.id))
        marker.id <- as.character(marker.id)

    if(is.factor(pop.id))
        pop.id <- as.character(pop.id)

    # boolean checks
    if(!is.logical(dev))
    {
        warning("Incorrect data type for dev, ignoring input.")
        dev <- FALSE
    }

    if(!is.logical(verbose))
    {
        warning("Incorrect data type for verbose, ignoring input.")
        verbose <- FALSE
    }

    if(!is.logical(debug))
    {
        warning("Incorrect data type for debug, ignoring input.")
        debug <- FALSE
    }

    # numeric checks
    if(!is.numeric(burn) | length(burn) != 1 | burn < 0)
    {
        warning("Invalid value for burn, ignoring input.")
        burn <- 100
    }

    if(!is.numeric(iter) | length(iter) != 1 | iter < 0)
    {
        warning("Invalid value for iter, ignoring input.")
        iter <- 200
    }

    if(!is.numeric(every) | length(every) != 1 | every < 1)
    {
        warning("Invalid value for every, ignoring input.")
        every <- 5
    }

    if(!is.numeric(lambda) | length(lambda) != 1 | lambda <= 0)
    {
        warning("Invalid value for lambda, ignoring input.")
        lambda <- 6
    }

    if(!is.numeric(tau) | !length(tau) %in% c(1, length(pop.id))  | tau <= 0)
    {
        warning("Invalid value for tau, ignoring input.")
        tau <- 300
    }

    if(!is.null(omega) & (!is.numeric(omega) | any(omega < 0) |
                          (length(omega) != length(pop.id) & !is.null(pop.id))))
    {
        warning("Invalid value for omega, ignoring input.")
        omega <- NULL
    }

    if(!is.null(rand.seed) & !is.numeric(rand.seed))
    {
        warning("Invalid value for rand.seed, ignoring input.")
        rand.seed <- NULL
    }else{
        if(!is.null(rand.seed))
            set.seed(rand.seed)
    }

    if(!is.numeric(cores) | cores < 1)
    {
        warning("Invalid value for cores, ignoring input.")
        cores <- detectCores()
    }

    if(!is.numeric(indiv.geno.check) | indiv.geno.check < 0 | indiv.geno.check > 1)
    {
        warning("Bad value for indiv.geno.check, ignoring input.")
        indiv.geno.check <- 0.98
    }

    if(!is.numeric(marker.geno.check) | marker.geno.check < 0 | marker.geno.check > 1)
    {
        warning("Bad value for marker.geno.check, ignoring input.")
        marker.geno.check <- 0.98
    }

    if(!is.null(pop.id) & length(pop.id) != length(Pm.prior[[1]]$freq))
        stop("Incorrect length of pop.id.")

    # other data types
    if(!is.null(cl) & !"cluster" %in% class(cl))
        stop("Invalid value for cl.")


########## check geno/haplo-types ##########

    if(is.null(geno) & is.null(haps))
        stop("Either geno or haps must be specified")

    # check that this matches with haps if given
    if(!is.null(haps))
        if(!all(apply(haps, 1:2, sum) == geno))
            stop("Both haps and geno given, but genotypes don't match")

    # geno format check
    if(any(dimnames(geno)[[1]] != indiv.id))
        warning(sum(dimnames(geno)[[1]] != indiv.id), " individual names in geno do not match indiv.id.")

    # gender format check
    if(!is.null(gender) & length(gender) != dim(geno)[1])
        stop("Length of gender should be ", dim(geno)[1], '.')

    # names of Pm.prior should be in column names of haps (and geno?)-needs to be updated
    if(!all(names(Pm.prior) %in% dimnames(geno)[[2]]))
        stop("Each element of Pm.prior must be associated with a column of geno.")

    if(!is.null(haps))
    {
        # haps format check
        if(any(dimnames(haps)[[1]] != indiv.id))
            warning(sum(dimnames(haps)[[1]] != indiv.id), " individual names in haps do not match indiv.id.")

        # names of Pm.prior should be in column names of haps (and geno?)-needs to be updated
        if(!all(names(Pm.prior) %in% dimnames(haps)[[2]]))
            stop("Each element of Pm.prior must be associated with a column of geno.")
    }

    if(any(names(gender) != indiv.id))
        warning(sum(names(gender) != indiv.id), " names of the gender variable do not match indiv.id.")

    # chr format check
    if(length(chr) != length(Pm.prior))
        stop("Length of chr should be ", length(Pm.prior), '.')

    if(any(names(chr) != marker.id))
        warning(sum(names(chr) != marker.id), " names of the chr variable do not match marker.id.")

    # Pm.prior format check
    if(!all(names(Pm.prior[[1]][[3]]) == pop.id) & !is.null(pop.id))
        stop("Pm.prior must have the same names as pop.id when pop.id is not null.")

    # pos/d format check
    if(length(names(Pm.prior)) != length(pos))
        stop("Length of pos (", length(pos),  ") should be match length of Pm.prior (", length(names(Pm.prior)), ').')

    if(any(names(Pm.prior) != marker.id))
        stop("Each marker in Pm.prior needs to be present in marker.id when not null.")

    if(any(unlist(lapply(Pm.prior, length)) != 3) | # each marker should have 3 elements
       length(table(unlist(lapply(Pm.prior, lapply, length)))) != 1) # each element should have length K
        stop("Formatting error in Pm.prior detected.")


########## be sure gender is specified if anything in chr == sex.chr ##########
    if(is.null(gender) & any(chr == sex.chr))
        stop("Sex chromosomes detected but gender is null")


########## gender check ##########
    if(any(chr == sex.chr))
    {
         # check males
        if(any(na.omit(unlist(geno[gender == 'M', chr == sex.chr])) == 1))
        {
             # get number of genotypes with more than one allele
            sex.check <- apply(geno[gender == 'M', chr == sex.chr] == 1, 1, sum, na.rm = TRUE)
            sex.check <- sex.check[sex.check > male.het.X]

             # make heterozygous loci missing
            geno[gender == 'M', chr == sex.chr] <- apply(geno[gender == 'M', chr == sex.chr], 1:2,
                                                         function(x) ifelse(x == 1, NA, x))

             # throw a warning about how many we have above the threshold
            if(length(sex.check) > 0)
                warning(length(sex.check), " males have between ", min(sex.check), " and ",
                        max(sex.check), " heterozygote X-chromosome genotypes.")
        }else{
            sex.check <- NULL
        }

        # check females -- given the non-random selection of SNPs (based on frequencies in races) and
        #                  the fact that some individuals might have very high admixture percentages
        #                  for some races, I'm not sure how exactly to do this for females...
        # PLINK uses the inbreeding coefficient: F ... might consider wether this would work here.
        autosomes <- apply(geno[gender == 'F', chr != sex.chr] > 1, 1, sum) / sum(chr != sex.chr)
        sex.check <- apply(geno[gender == 'F', chr == sex.chr] > 1, 1, sum) / sum(chr == sex.chr)
    }else{
        sex.check <- NULL
    }

########## check for complete genotyping (by individual / by marker) ##########
    indiv.check <- apply(!is.na(haps), 1, sum) / dim(haps)[2]
    marker.check <- apply(!is.na(haps), 2, sum) / dim(haps)[1]

    if(any(indiv.check < indiv.geno.check))
        warning(sum(indiv.check < indiv.geno.check), " individuals with less than ",
                round(indiv.geno.check * 100),
                "% complete genotyping detected.")

    if(any(marker.check < marker.geno.check))
        warning(sum(marker.check < marker.geno.check), " markers with less than ",
                round(marker.geno.check * 100),
                "% complete genotyping detected.")


########## check HWE ##########
    hwe <- NULL

    # have had some recent changes to the arguments...recheck this section!

    if(is.null(geno))
    {
        tmp <- array(0, dim = c(dim(haps)[1] / 2, dim(haps)[2], 2),
                     dimnames = list(dimnames(haps)[[1]][1:(dim(haps)[1]/2) * 2 - 1],
                                     dimnames(haps)[[2]], c('Mother', 'Father')))

        tmp[,,1] <- haps[1:dim(tmp)[1] * 2 - 1,]
        tmp[,,2] <- haps[1:dim(tmp)[1] * 2,]
        haps <- tmp

        geno <- apply(haps, 1:2, sum)
    }


    hwe <- apply(geno[,chr != sex.chr], 2, table)
    if(is.list(hwe))
    {
        tmp <- matrix(0, nrow = 3, ncol = length(hwe),
                      dimnames = list(c('0', '1', '2'), colnames(geno)[chr != sex.chr]))

        tmp['0',] <- sapply(hwe, `[`, '0')
        tmp['1',] <- sapply(hwe, `[`, '1')
        tmp['2',] <- sapply(hwe, `[`, '2')

        tmp[is.na(tmp)] <- 0

        hwe <- tmp
    }

    hwe <- hwexact(hwe['0',], hwe['1',], hwe['2',])
    names(hwe) <- marker.id[!chr == sex.chr]

    # only do females on this set...can't really say anything with males
    if(any(chr == sex.chr) & any(gender == 'F'))
    {
        hweX <- apply(geno[gender == 'F', chr == sex.chr], 2, table)
        if(is.list(hweX))
        {
            tmp <- matrix(0, nrow = 3, ncol = length(hweX),
                          dimnames = list(c('0', '1', '2'), colnames(geno)[chr == sex.chr]))

            tmp['0',] <- sapply(hweX, `[`, '0')
            tmp['1',] <- sapply(hweX, `[`, '1')
            tmp['2',] <- sapply(hweX, `[`, '2')

            tmp[is.na(tmp)] <- 0

            hweX <- tmp
        }

        hweX <- hwexact(hweX['0',], hweX['1',], hweX['2',])
        names(hweX) <- marker.id[chr == sex.chr]

        hwe <- c(hwe, hweX)
    }

    if(any(hwe < hwe.thresh))
        warning(sum(hwe < hwe.thresh), " markers failed HWE test.")

##### check for allele flips later on in the process #####
}
