# amap.R (ancestrymap)
# Run ancestrymap (or atleast set everything up)
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick

# this file is deprecated...should remove at some point

amap <- function(param.file.name = 'param', preset = NULL, risk = 2.0, indivname = 'indiv.dat',
                 snpname = 'snps.dat', genotypename = 'geno.dat', badsnpname = NA,
                 tlreest = TRUE, seed = NA, splittau = TRUE, fancyxtheta = TRUE, output = NA,
                 trashdir = NA, checkit = FALSE, details = FALSE, numburn = 1, numiters = 5,
                 emiter = 10, dotoysim = FALSE, cleanit = TRUE, reestiter = 5,
                 thetafilename = NA, lambdafilename = NA, freqfilename = NA,
                 indoutfilename = NA, snpoutfilename = NA, ethnicfilename = 'ethnic.out',
                 noxdata = FALSE, fakespacing = NA, thxpars = NA, thpars = NA, lampars = NA,
                 lamxpars = NA, markersim = NA, simnumindivs = NA, risksim = NA,
                 tauscal = NA, wrisk = NA, lrisk = NA, controlrisk = NA, risk2 = NA,
                 taulsdev = NA, taulmean = NA, allmale = NA, allcases = NA,
                 usecontrols = NA, pbfmodern = NA, pbxname = NA, genotoyoutfilename = NA,
                 indtoyoutfilename = NA, unknowngender = NA,
                 ethnic.control = NA, parameter.name = NULL, value = NULL,
                 run = TRUE, in.prefix = NULL, out.prefix = NULL, remove.files = TRUE)
{
  # check input
  if(is.null(param.file.name))
    stop("Must specify a parameter file name")
  if(length(parameter.name) != length(value))
    stop("parameter.name and value must have same length")

  if(is.na(seed))
    seed <- floor(runif(1)*10000)

  # change variables from a vectors to ' ' delimeted strings
  risk <- paste(risk, collapse=' ')
  if(!is.na(thxpars))
    thxpars <- paste(thxpars, collapse=' ')
  if(!is.na(thpars))
    thpars <- paste(thpars, collapse=' ')
  if(!is.na(lampars))
    lampars <- paste(lampars, collapse=' ')
  if(!is.na(lamxpars))
    lamxpars <- paste(lamxpars, collapse=' ')
  if(!is.na(tauscal))
    tauscal <- paste(tauscal, collapse=' ')

  # if we aren't checking, set defaults for unspecefied output file names
  if(!checkit)
  {
    if(is.na(thetafilename))
      thetafilename <- 'theta.out'

    if(is.na(output))
      output <- 'out.out'

    if(is.na(snpoutfilename))
      snpoutfilename <- 'snps.out'

    if(is.na(indoutfilename))
      indoutfilename <- 'ind.out'

    if(is.na(freqfilename))
      freqfilename <- 'freq.out'

    if(is.na(lambdafilename))
      lambdafilename <- 'lambda.out'
  }

  # add output prefix to output file names
  if(!is.null(out.prefix))
  {
    # make sure we get a '/' between directory and file
    if(substr(out.prefix, nchar(out.prefix), nchar(out.prefix)) != '/')
      out.prefix <- paste(out.prefix, '/', sep='')

    if(!is.na(thetafilename))
      thetafilename <- paste(out.prefix, thetafilename, sep='')

    if(!is.na(output))
      output <- paste(out.prefix, output, sep='')

    if(!is.na(snpoutfilename))
      snpoutfilename <- paste(out.prefix, snpoutfilename, sep='')

    if(!is.na(indoutfilename))
      indoutfilename <- paste(out.prefix, indoutfilename, sep='')

    if(!is.na(freqfilename))
      freqfilename <- paste(out.prefix, freqfilename, sep='')

    if(!is.na(lambdafilename))
      lambdafilename <- paste(out.prefix, lambdafilename, sep='')

    if(!is.na(ethnic.control))
      ethnic.control <- paste(out.prefix, ethnic.control, sep='')
  }

  # add prefixes to input file names
  if(!is.null(in.prefix))
  {
    # make sure we get a '/' between directory and file
    if(substr(in.prefix, nchar(in.prefix), nchar(in.prefix)) != '/')
      in.prefix <- paste(in.prefix, '/', sep='')

    indivname <- paste(in.prefix, indivname, sep='')

    genotypename <- paste(in.prefix, genotypename, sep='')

    snpname <- paste(in.prefix, snpname, sep='')

    if(!is.na(badsnpname))
      badsnpname <- paste(in.prefix, badsnpname, sep='')
  }

  # set up param *** NOTE: if you change this list don't forget to change .preset
  parameter.name <- c('risk', 'indivname', 'snpname', 'genotypename',
                      'badsnpname', 'tlreest', 'seed', 'splittau', 'fancyxtheta',
                      'output', 'trashdir', 'checkit', 'details', 'numburn', 'numiters',
                      'emiter', 'dotoysim', 'cleanit', 'reestiter', 'thetafilename',
                      'lambdafilename', 'freqfilename', 'indoutfilename',
                      'snpoutfilename', 'ethnicfilename', 'noxdata', 'fakespacing',
                      'thxpars', 'thpars', 'lampars', 'lamxpars', 'markersim',
                      'simnumindivs', 'risksim', 'tauscal', 'wrisk', 'lrisk',
                      'controlrisk', 'risk2', 'taulsdev', 'taulmean', 'allmale',
                      'allcases', 'usecontrols', 'pbfmodern', 'pbxname',
                      'genotoyoutfilename', 'indtoyoutfilename', 'unknowngender',
                      'ethnic.control', parameter.name)
  value <- c(risk, indivname, snpname, genotypename, badsnpname,
             tlreest, seed, splittau, fancyxtheta, output, trashdir, checkit,
             details, numburn, numiters, emiter, dotoysim, cleanit, reestiter,
             thetafilename, lambdafilename, freqfilename, indoutfilename,
             snpoutfilename, ethnicfilename, noxdata, fakespacing, thxpars,
             thpars, lampars, lamxpars, markersim, simnumindivs, risksim,
             tauscal, wrisk, lrisk, controlrisk, risk2, taulsdev, taulmean,
             allmale, allcases, usecontrols, pbfmodern, pbxname,
             genotoyoutfilename, indtoyoutfilename, unknowngender,
             ethnic.control, value)
  names(value) <- parameter.name

  param <- .preset(parameter.name, value, preset)

  param <- subset(param, !is.na(value))
  param$parameter.name <- paste(param$parameter.name, ':', sep='')

  param$value <- ifelse(param$value==TRUE, 'YES',
                        ifelse(param$value==FALSE, 'NO', as.character(param$value)))

  # write file
  write.table(param, file=param.file.name, col.names=FALSE, quote=FALSE, row.names=FALSE)

  # Do ethicity for controls? ... Certainly there is a better way to do this!
  if(!is.null(ethnic.control))
  {
    x <- read.table(indivname)
    x[[3]] <- ifelse(x[[3]] == 'Case', 'Control', ifelse(x[[3]] == 'Control', 'Case', x[[3]]))
    write.table(x, file='case.cont.reversed', col.names=FALSE,  quote=FALSE, row.names=FALSE)

    param$value[param$parameter.name == 'indivname:'] <- 'case.cont.reversed'
    param$value[param$parameter.name == 'ethnicfilename:'] <- ethnic.control

    write.table(param, file='param.reverse', col.names=FALSE, quote=FALSE, row.names=FALSE)
  }

  # Run ANCESTRYMAP
  if(run)
  {
    if(version$arch %in% c('i386', 'powerpc'))
    {
      warning("Skipping call to ANCESTRYMAP -- Not supported on this distribution")
      return()
    }

    if(!is.null(ethnic.control))
      system('ancestrymap -p param.reverse')

    system(paste('ancestrymap -p ', param.file.name, ' > ', out.prefix, 'out.file', sep=''))

    # remove anything un-wanted
    if(!is.null(ethnic.control))
      system('rm case.cont.reversed ')

    x <- read.amap(paste(out.prefix, 'out.file', sep=''), remove.files)

    return(x)
  }else{
    return()
  }
}


.preset <- function(parameter.name, value, preset)
{
  value <- as.list(value)

  if(!is.null(preset))
  {
    if(!(preset %in% c(-1:2)))
    {
      warning(paste("Invalid preset:", preset))
    }
  }

  # default
  if(is.null(preset))
    return(data.frame(parameter.name = c('risk', 'indivname', 'snpname', 'genotypename',
                        'badsnpname', 'tlreest', 'seed', 'splittau', 'fancyxtheta',
                        'output', 'trashdir', 'checkit', 'details', 'numburn', 'numiters',
                        'emiter', 'dotoysim', 'cleanit', 'reestiter', 'thetafilename',
                        'lambdafilename', 'freqfilename', 'indoutfilename',
                        'snpoutfilename', 'ethnicfilename', 'noxdata', 'fakespacing',
                        'thxpars', 'thpars', 'lampars', 'lamxpars', 'markersim',
                        'simnumindivs', 'risksim', 'tauscal', 'wrisk', 'lrisk',
                        'controlrisk', 'risk2', 'taulsdev', 'taulmean', 'allmale',
                        'allcases', 'usecontrols', 'pbfmodern', 'pbxname',
                        'genotoyoutfilename', 'indtoyoutfilename', 'unknowngender',
                        'ethnic.control',
                        parameter.name[51:max(51, length(parameter.name))]),
                      value = c(value$risk, value$indivname, value$snpname,
                        value$genotypename, value$badsnpname, value$tlreest, value$seed,
                        value$splittau, value$fancyxtheta, value$output, value$trashdir,
                        value$checkit, value$details, value$numburn, value$numiters,
                        value$emiter, value$dotoysim, value$cleanit, value$reestiter,
                        value$thetafilename, value$lambdafilename, value$freqfilename,
                        value$indoutfilename, value$snpoutfilename, value$ethnicfilename,
                        value$noxdata, value$fakespacing, value$thxpars, value$thpars,
                        value$lampars, value$lamxpars, value$markersim, value$simnumindivs,
                        value$risksim, value$tauscal, value$wrisk, value$lrisk,
                        value$controlrisk, value$risk2, value$taulsdev, value$taulmean,
                        value$allmale, value$allcases, value$usecontrols, value$pbfmodern,
                        value$pbxname, value$genotoyoutfilename, value$indtoyoutfilename,
                        value$unknowngender, value$ethnic.control,
                        value[51:max(51, length(value))])))


  # input file check
  if(preset == 0)
    return(data.frame(parameter.name = c('risk', 'indivname', 'snpname', 'genotypename',
                        'tlreest', 'seed', 'splittau', 'fancyxtheta', 'checkit',
                        'details', 'numburn', 'numiters', 'emiter', 'dotoysim',
                        'cleanit', 'reestiter', 'indoutfilename', 'snpoutfilename',
                        parameter.name[51:max(51, length(parameter.name))]),
                      value = c(value$risk, value$indivname, value$snpname, value$genotypename,
                        TRUE, value$seed, TRUE, value$fancyxtheta, TRUE, TRUE, 0, 0, 10, FALSE,
                        TRUE, 5, value$indoutfilename, value$snpoutfilename,
                        value[51:max(51, length(value))])))
  # data check
  if(preset == 1)
    return(data.frame(parameter.name = c('risk', 'indivname', 'snpname', 'genotypename',
                        'tlreest', 'seed', 'splittau', 'fancyxtheta', 'checkit', 'details',
                        'numburn', 'numiters', 'emiter', 'dotoysim', 'cleanit',
                        'reestiter', 'ethnicfilename',
                        parameter.name[51:max(51, length(parameter.name))]),
                      value = c(value$risk, value$indivname, value$snpname, value$genotypename,
                        TRUE, value$seed, TRUE, value$fancyxtheta, TRUE, TRUE, 5, 5, 10,
                        FALSE, TRUE, 5, value$ethnicfilename, value[51:max(51, length(value))])))
  # final run
  if(preset == 2)
    return(data.frame(parameter.name = c('risk', 'indivname', 'snpname', 'genotypename',
                        'badsnpname', 'tlreest', 'seed', 'splittau', 'fancyxtheta',
                        'output', 'trashdir', 'checkit', 'details', 'numburn', 'numiters',
                        'emiter', 'dotoysim', 'cleanit', 'reestiter', 'thetafilename',
                        'lambdafilename', 'freqfilename', 'indoutfilename',
                        'snpoutfilename', 'ethnicfilename',
                        parameter.name[51:max(51, length(parameter.name))]),
                      value = c(value$risk, value$indivname, value$snpname, value$genotypename,
                        value$badsnpname, TRUE, value$seed, TRUE, value$fancyxtheta,
                        value$output, value$trashdir, FALSE, FALSE, 100, 200, 30, FALSE,
                        TRUE, 5, value$thetafilename, value$lambdafilename, value$freqfilename,
                        value$indoutfilename, value$snpoutfilename, value$ethnicfilename,
                        value[51:max(51, length(value))])))
  # simulation
  if(preset == -1)
    return(data.frame(parameter.name = c('risk', 'indivname', 'snpname', 'genotypename',
                        'badsnpname', 'tlreest', 'seed', 'splittau', 'fancyxtheta',
                        'output', 'trashdir', 'checkit', 'details', 'numburn', 'numiters',
                        'simnumindivs', 'casecontrol', 'markersim', 'risksim', 'emiter',
                        'dotoysim', 'cleanit', 'reestiter', 'thetafilename',
                        'lambdafilename', 'freqfilename', 'indoutfilename',
                        'snpoutfilename', 'ethnicfilename'),
                      value = c(value$risk, value$indivname, value$snpname, value$genotypename,
                        value$badsnpname, TRUE, value$seed, TRUE, value$fancyxtheta,
                        value$output, value$trashdir, FALSE, FALSE, 10, 10, 1000, '500 200',
                        1200, value$risksim, 10, TRUE, TRUE, 5, value$thetafilename,
                        value$lambdafilename, value$freqfilename, value$indoutfilename,
                        value$snpoutfilename, value$ethnicfilename)))
}


### for debugging purposes only
#param.file.name = 'param'
#preset = NULL
#risk = 2.0
#indivname = 'indiv.dat'
#snpname = 'snps.dat'
#genotypename = 'geno.dat'
#badsnpname = NA
#tlreest = TRUE
#seed = NA
#splittau = TRUE
#fancyxtheta = TRUE
#output = NA
#trashdir = NA
#checkit = FALSE
#details = FALSE
#numburn = 1
#numiters = 5
#emiter = 10
#dotoysim = FALSE
#cleanit = TRUE
#reestiter = 5
#thetafilename = NA
#lambdafilename = NA
#freqfilename = NA
#indoutfilename = NA
#snpoutfilename = NA
#ethnicfilename = 'ethnic.out'
#noxdata = FALSE
#fakespacing = NA
#thxpars = NA
#thpars = NA
#lampars = NA
#lamxpars = NA
#markersim = NA
#simnumindivs = NA
#risksim = NA
#tauscal = NA
#wrisk = NA
#lrisk = NA
#controlrisk = NA
#risk2 = NA
#taulsdev = NA
#taulmean = NA
#allmale = NA
#allcases = NA
#usecontrols = NA
#pbfmodern = NA
#pbxname = NA
#genotoyoutfilename = NA
#indtoyoutfilename = NA
#unknowngender = NA
#ethnic.control = 'ethnic.control.out'
#parameter.name = NULL
#value = NULL
#run = TRUE
#in.prefix = NULL
#out.prefix = NULL
#remove.files = TRUE
