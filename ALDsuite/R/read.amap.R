# read.amap.R (ancestrymap)
# Read ANCESTRYMAP output into R
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created August 2005
# Last Modified May 25, 2011


read.amap <- function(file.name, remove.files=TRUE)
{
  # read the lines in and extract data one table at a time
  x <- readLines(file.name)
  x <- ifelse(nchar(x) > 100 & substr(x, 1, 1) == '#', ' ', x)


##############
# Parameters #
##############

  # start by getting risk(s) ... when multiple risks are used the get.amap.section gets confused ...
  to.get <- c(1:length(x))[substr(x, 1, 5) == 'risk:'][1]
  risk <- x[to.get]
  risk <- as.numeric(unlist(strsplit(substr(risk, 7, nchar(risk)), ' ')))

  # take out header line and clean up to eliminate problems with multiple risks
  x[to.get-1] <- paste('#', x[to.get-1])
  rm(to.get)

  param <- get.amap.section(x, '### THE INPUT PARAMETERS', file.out=file.name, header=FALSE)
  if(is.null(param))
  {
    stop('Is this output from ANCESTRYMAP? It is missing the parameter section.')
  }

  # get parameter names and values
  param.names <- as.character(param[[1]])
  param <- as.list(param[[2]])
  names(param) <- param.names

  # switch yes's and no's to T's and F's
  param <- ifelse(param == 'YES', TRUE, ifelse(param == 'NO', FALSE, param))

  # each parameter name should end with a ':'
  n <- nchar(names(param))
  param <- param[substr(names(param), n, n) == ':'] # sometimes something extra sneeks in

  # drop ':'
  names(param) <- gsub(':', '', names(param))

  # make a named list of all parameters
  param <- lapply(param, function(x){
    if(is.na(suppressWarnings(as.numeric(x))) | is.logical(x))
      return(x)
    else
    return(as.numeric(x))})

  # clean up and add any lost risk values
  param$risk <- risk

  rm(param.names, risk)


####################
# Genetic Distance #
####################

  gen.dist <- get.amap.section(x, '###GENETIC DISTANCE FOR ALL CHROMOSOMES', file.out=file.name, messed.up.header=TRUE)
  if(!is.null(gen.dist))
  {
    gen.dist <- subset(gen.dist, gen.dist$v1 == 'chrom:', select=c(2,4,6,8))
    names(gen.dist) <- c('chr.num', 'first.snp', 'last.snp', 'gen.dist')
    gen.dist$first.snp <- as.integer(as.character(gen.dist$first.snp))
    gen.dist$chr.num <- as.integer(as.character(gen.dist$chr.num))
  }

###########################
# Males Heterozygous at X #
###########################

  het.x <- get.amap.section(x, '###HETXCHECK RESULTS BEGIN', file.out=file.name, header=FALSE)
  if(!is.null(het.x))
  {
    het.x$v1 <- NULL
    names(het.x) <- c('snp.id', 'num.het', 'num.homozy')
    het.x$snp.id <- as.character(het.x$snp.id)
  }

##########
# Counts #
##########

  tmp <- get.amap.section(x, '###COUNTS', file.out=file.name, header=FALSE)
  counts <- list(fake.markers = as.numeric(as.character(tmp$v5[1])),
                 real.markers = as.numeric(as.character(tmp$v10[1])),
                 fake.space = as.numeric(as.character(tmp$v15[1])),
#                ignore.x = as.numeric(as.character(v11[2])),
                 num.markers = as.numeric(as.character(tmp$v4[2])),
                 num.samples = as.numeric(as.character(tmp$v8[2])),
                 num.case = as.numeric(as.character(tmp$v4[3])),
                 num.cont = as.numeric(as.character(tmp$v8[3])),
                 ignore.sample = as.numeric(as.character(tmp$v13[3])))

  temp <- get.amap.section(x, '### CHECKGENO RESULTS FOLLOW:', file.out=file.name, header=FALSE)
  n <- length(temp$v4)
  counts$good.geno <- temp$v4[n] # these will be NULL if temp is NULL
  counts$bad.geno <- temp$v8[n]

######################
# Physical Map Check #
######################

  phys.check <- get.amap.section(x, '###PHYSCHECK RESULTS FOLLOW:', file.out=file.name, header=FALSE)
  if(!is.null(phys.check))
  {
    phys.check$v1 <- NULL
    names(phys.check) <- c('snp1.id', 'snp2.id', 'snp1.gen.pos', 'snp2.gen.pos', 'snp1.phys.pos', 'snp2.phys.pos')
    phys.check$snp1.id <- as.character(phys.check$snp1.id)
    phys.check$snp2.id <- as.character(phys.check$snp2.id)
  }

#########################
# Hardey-Weinburg Check #
#########################

  hw <- get.amap.section(x, '###HWCHECK RESULTS FOLLOW:', file.out=file.name, header=FALSE)
  if(!is.null(hw))
  {
    hw <- subset(hw, hw$v1 == 'hwcheck', select=c(2:5))
    names(hw) <- c('snp.id', 'snp.index', 'chr.num', 'hw.score')
    hw$snp.id <- as.character(hw$snp.id)
    hw$snp.index <- as.integer(as.character(hw$snp.index))
    hw$chr.num <- as.numeric(as.character(hw$chr.num))
  }

###############################
# Duplicate Individuals Check #
###############################

  dups <- get.amap.section(x, '###CHECKDUP RESULTS FOLLOW:', file.out=file.name, header=FALSE)
  if(!is.null(dups))
  {
    dups.sub <- subset(dups, dups$v1 == 'match:', select=c(2,4))

    dups <- subset(dups, dups$v1 == 'dup?', select=c(2,3))
    names(dups) <- c('indiv1', 'indiv2')

    dups$indiv1 <- as.character(dups$indiv1)
    dups$indiv2 <- as.character(dups$indiv2)
    dups$match <- as.numeric(as.character(dups.sub[[1]]))
    dups$mismatch <- as.numeric(as.character(dups.sub[[2]]))
  }

#####################
# Individual Checks #
#####################

  if(any(substr(x, 1, 10) == 'checkindiv'))
  {
    to.skip <- min( (1:length(x))[substr(x, 1, 10) == 'checkindiv']) - 1 # check how many rows to skip
    nrows <- length(x[substr(x, 1, 10) == 'checkindiv']) #check howmany we want

    check.indiv <- read.table(file.name, skip=to.skip, nrows=nrows, comment.char='*') # get data
    check.indiv$V1 <- NULL

    names(check.indiv) <- c('indiv.id', 'theta', 'v2', 'score', 'v4', 'v5', 'v6')
                                        # see warnings for individuals with no autosomal data
    check.indiv <- subset(check.indiv, select=c('indiv.id', 'theta', 'score'))
  }else{
    check.indiv <- NULL
  }

############
# Warnings #
############

  if(any(substr(x, 1, 11) == '*** warning'))
  {
    warnings <- x[substr(x, 1, 11) == '*** warning']

    tmp <- lapply(warnings, warning)
  }else{
    warnings <- NULL
  }

######################
# Individual Details #
######################

  individuals <- get.amap.section(x, '###DETAILS ABOUTTHE INDIVIDUALS', file.out=file.name)
  if(!is.null(individuals))
  {
    individuals$status <- factor(individuals$status)
  }

##################
# Marker Details #
##################

  markers <- get.amap.section(x, '###DETAILS ABOUT THE MARKERS', file.out=file.name, messed.up.header=TRUE)
  if(!is.null(markers))
  {
    names(markers)[1:8] <- c('snp.id', 'chr.num', 'gen.pos', 'phys.pos', 'afr.vart', 'afr.ref', 'eur.vart', 'eur.ref')
  }else{
    markers <- try(
                   get.amap.section(readLines(param$snpoutfilename), '###DETAILS ABOUT THE MARKERS',
                                    file.out=param$snpoutfilename, messed.up.header=TRUE),
                  silent=TRUE)

    if(class(markers) != 'try-error' & !is.null(markers)){
      names(markers)[1:8] <- c('snp.id', 'chr.num', 'gen.pos', 'phys.pos', 'afr.vart', 'afr.ref',
                               'eur.vart', 'eur.ref')
    }else{
      markers <- NULL
    }
  }



#############################
# Scores from EM Iterations #
#############################

  em.scores <- get.amap.section(x, '##SCORES FROM EXPECTATION_MAXIMIZATION ALGORITHM ITERATIONS', file.out=file.name, header=FALSE)
  if(!is.null(em.scores))
  {
    em.scores <- subset(em.scores, select=c(3,4))
    names(em.scores) <- c('iteration', 'score')
  }

################
# MCMC Results #
################

  mcmc.results <- get.amap.section(x, '###RESULTS FOR EACH MARKOV CHAIN MONTE CARLO ITERATION', file.out=file.name, header=FALSE)
  if(!is.null(mcmc.results))
  {
    mcmc.results <- subset(mcmc.results, mcmc.results$v1 == 'estglob', select=c(2:8))
    names(mcmc.results) <- c('parameter', 'iteration', 'p1', 'p2', 'xp1', 'xp2', 'average')
    mcmc.results$parameter <- factor(mcmc.results$parameter)
  }

#######################
# Posterior Estimates #
#######################

  estimates <- get.amap.section(x, '###POSTERIOR ESTIMATES', file.out=file.name, header=FALSE)
  if(!is.null(estimates))
  {
    names <- paste(estimates$v1, estimates$v2, sep='.')
    estimates <- as.list(estimates$v3)
    names(estimates) <- names
  }

######################
# Genome-Wide Scores #
######################

  genome.score <- get.amap.section(x, '###GENOME_WIDE SCORE FOR ALL THE MODELS', file.out=file.name, messed.up.header=TRUE)
  if(!is.null(genome.score))
  {
    genome.score[[1]] <- NULL
    names(genome.score) <- c('risk1', 'risk2', 'crisk', 'score')
    rownames(genome.score) <- as.numeric(rownames(genome.score)) - 1

    genome.score$risk1 <- as.numeric(as.character(genome.score$risk1))
    genome.score$risk2 <- as.numeric(as.character(genome.score$risk2))
    genome.score$crisk <- as.numeric(as.character(genome.score$crisk))
  }

##################################
# Individual Theta/Lambda Values #
##################################

  temp <- get.amap.section(x, '###THETA or M, LAMBDA VALUES FOR ALL INDIVIDUALS', file.out=file.name)
  if(is.null(individuals))
  {
    individuals <- temp
  }else{
    if(!is.null(temp))
      individuals <- merge(individuals, temp, all=TRUE)
  }

##############################
# Allele Frequency Estimates #
##############################

  temp <- get.amap.section(x, '###ALLELE FREQUENCY ESTIMATES WITH STANDARD ERROR', file.out=file.name)
  if(is.null(markers))
  {
    markers <- temp
  }else{
    if(!is.null(temp))
      markers <- merge(markers, temp, all=TRUE)
  }

########################
# Lag and Correlations #
########################

# do I really need this??? It is all messed up since version 0.0.3
  lag.corr <- NULL
#  lag.corr <- get.amap.section(x, '###LAG AND CORRELATIONS', file.out=file.name, header=FALSE)
#
#  if(!is.null(lag.corr))
#  {
#    lag.corr <- subset(lag.corr, !is.na(v4), select=c(v2, v4))
#    names(lag.corr) <- c('lag', 'corr')
#
#    n <- length(lag.corr$lag)
#
#    lag.corr$label <- c(rep('llike', n/4), rep('log10fac', n/4), rep('factor', n/4), rep('tauscal', n/4))
#  }

##########################
# Scores for Each Marker #
##########################

  marker.stats <- get.amap.section(x, '###SCORES FOR EACH MARKER (fakes used for global score)', file.out=file.name, header=FALSE)
  if(!is.null(marker.stats))
    names(marker.stats) <- c('snp.index', 'chr', 'snp.id', 'phys.pos', 'gen.pos', 'lgs', 'ccs', 'g.case', 'g.control', 'rpower')

##############################
# Scores for Each Chromosome #
##############################

  chr.stats <- get.amap.section(x, '###SCORES FOR EACH CHROMOSOME', file.out=file.name)
  if(!is.null(chr.stats))
  {
    chr.stats <- subset(chr.stats, chr.stats$chr %in% c(1:23))
    chr.stats$chr.num <- as.numeric(as.character(chr.stats$chr))
    chr.stats$lgs.max <- as.numeric(as.character(chr.stats$lgs.max))
  }


#################
# Genome Scores #
#################

# this is not working ...
  best.scores <- NULL
#  best.scores <- get.amap.section(x, '###BESTSCORES: Maximum genome-wide score for the locus-genome statistic (LGS_MAX), and the maximum and minimum genome-wide scores for the case-control statistic (CCS_MAX and CCS_MIN)', file.out=file.name, header=FALSE)
#  if(!is.null(best.scores))
#    best.scores <- list(lgs = best.scores$v2,
#                        ccs.max = best.scores$v3,
#                        ccs.min = best.scores$v4)
#
#  best.scores$genome.log.factor <- get.amap.section(x, '###GENOME LOG FACTOR: log-likelihood of the locus genome statistic averaged over all the markers in the genome', file.out=file.name, header=FALSE)$v3

#############
# Map Check #
#############

  temp <- get.amap.section(x, '###MAPCHECK RESULTS FOLLOW:', file.out=file.name, header=FALSE)
  if(is.null(markers))
  {
    markers <- temp
  }else{
    if(!is.null(temp))
    {
      names(temp) <- c('drop', 'snp.id', 'snp.index', 'ancestry.diff')
      temp$drop <- NULL

      markers <- merge(markers, temp, all=TRUE)
    }
  }

###################
# Frequency Check #
###################

  temp <- get.amap.section(x, '###FREQCHECK RESULTS FOLLOW:', file.out=file.name, header=FALSE)
  if(!is.null(temp))
  {
    names(temp) <- c('drop', 'snp.id', 'drop', 'score.all', 'score.controls',
                     'drop', 'drop', 'cal.af.fq', 'cal.eur.fq')
    temp <- subset(temp, drop == 'freqcheck', select=names(temp)[names(temp) != 'drop'])

    temp$score.all <- as.numeric(as.character(temp$score.all))
    temp$score.controls <- as.numeric(as.character(temp$score.controls))

    if(is.null(markers))
      markers <- temp
    else
      markers <- merge(markers, temp, all=TRUE)
  }

##########################
# Get Ancestry Estimates #
##########################

  temp <- grep('### Gammas for marker:', x)
  if(length(temp) == 1)
  {
    temp <- get.amap.section(x, x[temp+1], file.out=file.name)
    names(temp) <- c('g0', 'g1', 'g2', 'logscore', 'drop')
    temp$indiv.id <- row.names(temp)

    gammas <- subset(temp, select = c('indiv.id', 'g0', 'g1', 'g2', 'logscore'))
  }else{
    gammas <- NULL
  }

########################
# Read in Output Files #
########################

  # Individuals:
  temp <- try(read.table(param$indoutfilename, sep='', header=FALSE), silent=TRUE)
  if(class(temp) != 'try-error'){
    names(temp) <- c('indiv.id', 'gender', 'status', 'obs')
  }else{
    temp <- NULL
  }

  thetas <- try(read.table(param$thetafilename, sep='', header=TRUE), silent=TRUE)

  if(class(thetas) != 'try-error')
  {
    names(thetas) <- c('indiv.index', 'indiv.id', 'tmean', 'tsdev', 'txmean',
                       'txsdev', 'status')

    if(!is.null(temp)){
      temp <- merge(temp, thetas, all=TRUE)
    }else{
      temp <- thetas}
  }else{
    thetas <- NULL
  }

  lambdas <- try(read.table(param$lambdafilename, sep='', header=TRUE), silent=TRUE)

  if(class(lambdas) != 'try-error')
  {
    names(lambdas) <- gsub('_', '', gsub('am', '', tolower(names(lambdas))))
    names(lambdas) <- gsub('vi', 'v.i', names(lambdas)) # add '.' after indiv
    temp <- merge(temp, lambdas, all.x=TRUE)
  }else{
    lambdas <- NULL
  }

  if(is.null(individuals)){
    individuals <- temp
  }else{
    if(!is.null(temp))
      individuals <- merge(individuals, temp, all.x=TRUE)
  }

  # Ethnicity
  ethnic <- try(read.table(param$ethnicfilename, sep='', header=TRUE), silent=TRUE)

  if(class(ethnic) != 'try-error')
  {
    names(ethnic) <- gsub('_', '.', tolower(names(ethnic)))
    names(ethnic)[names(ethnic) == 'avg.ethnicity'] <- 'case.ethnicity'
  }else{
    ethnic <- NULL
  }

  temp <- try(read.table(param$ethnic.control, sep='', header=TRUE), silent=TRUE)
  if(class(temp) != 'try-error')
  {
    names(temp) <- gsub('_', '.', tolower(names(temp)))
    names(temp)[names(temp) == 'avg.ethnicity'] <- 'cont.ethnicitiy'

    ethnic <- merge(ethnic, temp)
  }else{
    temp <- NULL
  }

  # More Marker info
  freq <- try(read.table(param$freqfilename, sep='', header=TRUE), silent=TRUE)
  if(class(freq) != 'try-error')
  {
    names(freq) <- gsub('_', '.', tolower(names(freq)))

    if(!is.null(temp)){
      temp <- merge(temp, freq, all=TRUE)
    }else{
      temp <- freq}
  }else{
    freq <- NULL
  }

  if(is.null(markers)){
    markers <- temp
  }else{
    if(!is.null(temp))
      markers <- merge(markers, temp, all.x=TRUE)}

  # Other output
  output <- try(read.table(param$output, sep='', header=TRUE), silent=TRUE)
  if(class(output) != 'try-error'){
    names(output) <- gsub('_', '.', tolower(names(output)))
  }else{
    output <- NULL
  }

  pubx <- try(read.table(param$pubxname, sep='', header=TRUE), silent=TRUE)
  if(class(pubx) != 'try-error'){
    names(pubx) <- gsub('_', '.', tolower(names(pubx)))
  }else{
    pubx <- NULL
  }


###############################
# Create S3 Object and Return #
###############################

  x <- list(file.name = file.name,
            param = param,
            genome.score = genome.score,
            markers = markers,
            individuals = individuals,
            gammas = gammas,
            checks = list(counts = counts,
              dups = dups,
              hw = hw,
              het.x = het.x,
              phys.check = phys.check,
              gen.dist = gen.dist,
              check.indiv = check.indiv,
              pubx = pubx),
            stats = list(best.scores = best.scores,
              chr.stats = chr.stats,
              marker.stats = marker.stats,
              estimates = estimates,
              em.scores = em.scores,
              mcmc.results = mcmc.results,
              lag.corr = lag.corr),
            ethnic = ethnic,
            output = output,
            warnings = warnings)

  class(x) <- 'amap'
  return(x)
}
