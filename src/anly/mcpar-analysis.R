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

mcparam.ML <- function(mc.data)
{
    ## Return the maximum likelihood parameters as a vector
    nparam <- if('iter' %in% names(mc.data))
                  ncol(mc.data) - 2          # assumes iter column is at the end
              else
                  ncol(mc.data) - 1
    v <- mc.data[which.max(mc.data$LL),] %>% as.matrix %>% as.vector
    v[1:nparam]
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

mcparam.itercount <- function(niter, nproc, npset)
{
    ## Return a vector giving the iteration count for a dataset output
    ## by the monte carlo calc.  This is less straightforward than
    ## just numbering the rows of the dataset sequentially for two
    ## reasons.  First, because of the parallel markov chains, there
    ## are many data points at each iteration.  Second, the iterations
    ## are not output sequentially; they are gathered together into
    ## batches and dumped periodically.
    ##
    ## TODO: This is really brittle.  It would be better just to
    ## record the iteration number in the output.
    ##
    ## niter:  number of iterations in the mcpar loop
    ## nproc:  number of processors allocated to the calculation
    ## npset:  number of parameter sets per processor
    ##
    ## Return value: vector of iteration sequence numbers

    ntot <- niter*nproc*npset
    nchain <- nproc*npset
    outstep <- if(niter > 50)
                   niter/10
               else
                   5
    nbatch <- niter %/% outstep           # This will probably cause a problem if niter is not divisible by outstep.

    out.batch <- seq(0,ntot-1) %/% (nbatch*nchain) # output batch number for each output slot (start count at 0)

    ## A batch is a series of nproc blocks, each originating from one
    ## processor.  The structure within each block is identical, so we
    ## can build up the structure of a block and repeat it.  The
    ## batches themselves will be repeated below, so it isn't
    ## necessary to construct the batch structure; we can just repeat
    ## the block structure as required.
    nblock <- outstep*npset
    seq.block <- 1 + seq(0,nblock-1) %/% npset

    ## Now the iteration number is the batch sequence number plus outstep*(batch #)
    out.batch * outstep + seq.block
}
