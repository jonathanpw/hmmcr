# This script file produces trace plots and histograms of the mcmc output files

trace_plot = function( chain, n_post, steps, burnin){

	# Matrix row indices for the posterior sample to use for GFF computation
	index_post = (steps - burnin - n_post + 1):(steps - burnin)

	# ----------------------------------------------------------------------------
	# Create mcmc trace plots and histograms
	# ----------------------------------------------------------------------------

	# Plot and save the mcmc trace plots and histograms.
	par_mean = par_median = upper = lower = rep( NA, ncol(chain))
	pdf('trace_plot.pdf')
	par(mfrow=c(3, 2))
	for(r in c(3:ncol(chain),1:2)){

		plot( NULL, xlab=NA, ylab=NA, xlim=c(1,length(index_post)), 
		      ylim=range(chain[,r]) )
		lines( chain[,r], type='l', col='black')

		par_mean[r] = round( mean(chain[,r]), 4)
		par_median[r] = round( median(chain[,r]), 4)
		upper[r] = round( quantile( chain[,r], prob=.975), 4)
		lower[r] = round( quantile( chain[,r], prob=.025), 4)

		hist( chain[,r], ylab=NA, main=NA, freq=F,
			    breaks=sqrt(nrow(chain)),
		      xlab=paste0('Mean = ',toString(par_mean[r]),
				  ' Median = ',toString(par_median[r])))
		abline( v=upper[r], col='red', lwd=2, lty=2)
		abline( v=lower[r], col='purple', lwd=2, lty=2)
	}
	dev.off()

	# ----------------------------------------------------------------------------
}






