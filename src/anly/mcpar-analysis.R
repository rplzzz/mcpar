## mcpar-analysis.R: A set of (hopefully) useful functions for
## analyzing the results of Monte Carlo calculations.

library(reshape2)
library(ggplot2)

read.mc.data <- function(filename, varnames=NULL)
{
    ## This is just a thin wrapper around the read.table function.  If
    ## names are given, we assign them; if not, then we just tag the
    ## last column as the log likelihood value.
    data <- read.table(filename)

    if(!is.null(varnames))
        names(data) <- varnames
    else
        names(data)[ncol(data)] <- "LL"
    data
}


mcparam.density <- function(mc.data)
{
    ## Create a density plot for all of the MC variables (including likelihood).
    data.m <- melt(mc.data)

    ggplot(data.m, aes(x=value)) + facet_wrap(~variable, scales='free') + geom_density(fill='grey')
}

mcparam.sample <- function(mc.data, nsamp=100, func=NULL)
{
    ## Sample the MC results using bootstrap sampling.  Optionally,
    ## apply a function to the sampled values.
    ##
    ## Return value:  data frame of samples, if func=NULL; otherwise,
    ##       a list containing the data frame in the first element and
    ##       the function values in the second.
    mcsamp <- mc.data[sample.int(nrow(mc.data), size=nsamp, replace=TRUE),]

    if(!is.null(func)) {
        fvals <- apply(mcsamp, 1, func)
        list(mcsamp, fvals)
    }
    else {
        mcsamp
    }
}

mcparam.clip.tails <- function(mc.data, qlo=0.01, qhi=0.99)
{
    ## Return a vector: True for rows that are in the main body of the
    ##                  distribution; False for rows that are in the
    ##                  tails.
    ##
    ## qlo:  quantile for the lower tail boundary
    ## qhi:  quantile for the upper tail boundary

    quants <- apply(mc.data, 2, function(x) {quantile(x,probs=c(qlo,qhi))})
    ## We want to make an exception for the log-likelihood.  Don't
    ## clip its upper end.  Log-likelihood is in the last column, and
    ## its theoretical maximum is zero.
    quants[2,ncol(quants)] <- 0.0 

    keep <- apply(mc.data, 1, function(x){all(x>quants[1,] & x<quants[2,])})

    mc.data[keep,]
}

