# summary.amap.R (ancestrymap)
# Summary method for an amap object
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick


summary.amap <- function(object, ...)
{
  #attach(object)

  cat("\n\n\nANCESTRYMAP Summary\n", ...)
  cat(paste("File: ", object$file.name, '\n', sep=''), ...)

  # Checkit section
  if(object$param$checkit)
  {
    cat("\ncheckit:\n", ...)
    cat("--------\n", ...)

    # Dups
    if(!is.null(object$checks$dups))
    {
      cat("\nDups\n", ...)
      print(object$checks$dups)
    }

    # genocheck
    if(!is.null(object$checks$counts$good.geno))
    {
      cat("\nGenotype Check\n", ...)
      print(data.frame(good.geno=object$checks$counts$good.geno,
                       bad.geno=object$checks$counts$bad.geno))
      if(object$checks$counts$bad.geno > 0)
        cat(paste("-- See ", object$file.name, " for more details", sep=''), ...)
    }

    # het.x check
    if(!is.null(object$checks$het.x))
    {
      x <- subset(object$checks$het.x, object$checks$het.x$num.het > 0)
      if(length(x$snp.id) > 0)
      {
        cat("\nMale X Heterozygousity\n", ...)
        print(x)
      }
    }

    # physcheck
    if(!is.null(object$checks$phys.check))
    {
      cat("\nPhysical Map Check\n", ...)
      print(object$checks$phys.check)
    }

    # hw check
    if(!is.null(object$checks$hw))
    {
      cat("\nTop 10 Hardy-Weibnurg Scores\n", ...)
      print(object$checks$hw[order(-object$checks$hw$hw.score),][1:10,])
    }

    # mapcheck
    if(!is.null(object$markers$ancestry.diff))
    {
      cat("\nTop 10 Mapcheck Scores\n", ...)
      x <- subset(object$markers, select=c('snp.id', 'snp.index', 'ancestry.diff'))
      print(x[order(-x$ancestry.diff),][1:10,])
    }

    # freqcheck
    if(!is.null(object$markers$score.all))
    {
      cat("\nFrequency Check -- All\n", ...)
      x <- subset(object$markers, select=c('snp.id', 'chr.num', 'score.all', 'score.controls',
                                           'amean', 'bmean', 'cal.af.fq', 'cal.eur.fq'))
      names(x) <- c('snp.id', 'chr', 'score.all', 'score.cont', 'input.A', 'input.B',
                    'calc.A', 'calc.B')
      print(x[order(-x$score.all),][1:10,])

      cat("\nFrequency Check -- Controls\n", ...)
      print(x[order(-x$score.cont),][1:10,])
    }

    # check.indiv
    if(!is.null(object$checks$check.indiv))
    {
      x <- object$checks$check.indiv

      cat("\nCheck Indiv\n", ...)
      print(x[order(-x$score),][1:10,]) # highest 10
      cat(".\n.\n.\n", ...)
      print(x[order(x$score),][1:10,]) # lowest 10
    }

    # genetic distance of chromosomes
    if(!is.null(object$checks$gen.dist))
    {
      cat("\nGenetic Distance of Chromosomes\n", ...)
      print(object$checks$gen.dist)
    }

    cat("\n-- end checkit --\n", ...)
  }
  # end checkit

  # Model Statistics
  if(!is.null(object$genome.score))
  {
    cat("\nModel Statistics\n", ...)
    print(object$genome.score)
  }

  # Chromosome statistics
  if(!is.null(object$stats$chr.stats))
  {
    cat("\nChromosome Statistics\n", ...)
    print(object$stats$chr.stats)
  }

  # Genome statistics
  if(!is.null(object$stats$best.scores))
  {
    cat("\nGenome Scores\n", ...)
    write.table(object$outf$stats$best.scores, sep='\t', row.names=FALSE)
  }

#  detach()

}
