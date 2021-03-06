% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xtest.R
\name{xtest}
\alias{xtest}
\title{Performs an exact test with full enumeration for Hardy-Weinberg proportions.}
\usage{
xtest(c, statName = "LLR", histobins = 0, histobounds = c(0, 0),
  showCurve = T, safeSecs = 100, detail = 2)
}
\arguments{
\item{c}{A matrix containing the genotype counts. It should be a square matrix, but only the lower-left half is used.}

\item{statName}{can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} depending on which one is to be ploted. Note that P values for all four are computed regardless of which one is specified with this parameter.}

\item{histobins}{If 0 no histogram is plotted. If 1 or \code{TRUE} a histogram with 500 bins is plotted. If set to a number greater than 1, a histogram with \code{histobins} is plotted.}

\item{histobounds}{A vector containing the left and right boundaries for the histogram's x axis. If you leave this as the default, \code{c(0,0)}, then \code{xtest} will compute reasonable bounds to include most of the distribution.}

\item{showCurve}{whether to show a blue curve indicating the asymptotic (chi squared) distribution. This only works for \code{LLR} and \code{Chisq}}

\item{safeSecs}{After this many seconds the calculation will be aborted. This is a safety valve to prevent attempts to compute impossibly large sets of tables.}

\item{detail}{Determines how much detail is printed. If set to 0, nothing is printed (useful if you use \code{xtest} programmatically.).}
}
\value{
\code{xtest} returns a list components
\item{$ Pvalues}{The four computed P values corresponding to the test statistics: \code{LLR}, \code{Prob}, \code{U} and \code{Chisq} in that order.}
\item{$ observed}{The four observed statistics in the same order as above}
\item{$ tableCount}{The number of tables examined during the calculation}
\item{$ ntrials}{placeholder}
\item{$ genotypes}{The input matrix of genotype counts}
\item{$ alleles}{The allele counts \eqn{m} corresponding to the input genotype counts}
\item{$ statID}{Which test statistic was used if a histogram was plotted}
\item{$ histobins}{If greater than zero, the number of bins to use for the histogram}
\item{$ histobounds}{The lower and upper limits of the test statistic in the histogram}
\item{$ histoData}{Vector of \eqn{histobins} values for the histogram}
\item{$ showCurve}{Whether the asymptotic curve should be plotted with the histogram}
}
\description{
Given a set of genotype counts, \code{xtest} examines all possible outcomes with the same set of allele counts. For each table, it computes four test statistics and compares them with the observed values. It returns the total probability of all tables with test statistics as \dQuote{extreme} or more so than the observed. It can also plot a histogram of one of the statitistics if \code{histobins} is greater than zero. More about these four test statistics and other information can be found in the vignette. This function will not normally be called directly. Instead, \code{\link{hwx.test}} calls either \code{xtest} or \code{\link{mtest}} depending on which method is to be used.
}
\references{
The methods are described by \href{http://dx.doi.org/10.1534/genetics.109.108977}{Engels, 2009. \bold{Genetics} 183:1431}.
}
\seealso{
\code{\link{hwx.test}}
}

