score_plots <- function(fitted_outrig, plot_together=TRUE){
  #' Plots the fitted score functions learnt in an \code{outrigger} fit.
  #'
  #' @param fitted_outrig a fitted object from \code{outrigger}
  #' @param plot_together a logical denoting whether all score function estimators (across all covariate bins) should be plotted toghether. Default is TRUE.
  #'
  #' @importFrom np npregbw
  #' @importFrom RColorBrewer brewer.pal
  #' @importFrom stats quantile
  #' @importFrom graphics hist legend lines
  #' @export
  #'

  hist(fitted_outrig$score_plot_metadata$all_pilot_resids,
       main="Histogram of all residuals in lambda*h ball about x",
       xlab="epsilon"
       )

  #xmin <- min(fitted_outrig$score_plot_metadata$all_pilot_resids)
  #xmax <- max(fitted_outrig$score_plot_metadata$all_pilot_resids)
  xmin <- quantile(fitted_outrig$score_plot_metadata$all_pilot_resids,0.01)
  xmax <- quantile(fitted_outrig$score_plot_metadata$all_pilot_resids,0.99)
  xvals <- seq(xmin,xmax,length.out=100)

  M <- length(fitted_outrig$score_plot_metadata$all_scores[[1]]$score)
  K <- length(fitted_outrig$score_plot_metadata$all_scores)

  coul <- brewer.pal(8, "Set1")
  for (k in seq_len(K)) {

    yvals <- matrix(NA, nrow=M, ncol=length(xvals))

    for (m in seq_len(M)) {
      yvals[m,] <- sapply(xvals, fitted_outrig$score_plot_metadata$all_scores[[k]]$score[[m]]$rho)
    }
    plot_ylimvals <- c(min(yvals), max(yvals))

    if (plot_together) {
      for (m in seq_len(M)) {
        if (m==1) {
          plot(xvals, yvals[m,], type="l", lwd=2, col=coul[1], ylim=plot_ylimvals,
               main=paste0("Score estimator. Fold number: ",k),
               xlab="error",
               ylab="score(error|x)"
               )
        } else {
          lines(xvals, yvals[m,], lwd=2, col=coul[m])
        }
        legend("topright",
               legend=c(1:M),
               col=c(coul[1:M]),
               title="Covariate bin no.",
               cex=0.8,
               lwd=2)
      }
    } else {
      for (m in seq_len(M)) {
        plot(xvals, yvals[m,], type="l", lwd=2, col=coul[m],
             main=paste0("Fold number: ",k, ", Covariate bin number: ",m),
             xlab="error",
             ylab="score(error|x)"
             )
      }
    }
  }
}
