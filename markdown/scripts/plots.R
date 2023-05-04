library(ggplot2)


plot_mtx <- function(Dx, xlab.title="Region", ylab.title="Region") {
  data <- melt(Dx)
  ggplot(data, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(name="ICC",
                         colours=(c("white", "#d45757")),
                         limits=c(min(Dx), max(Dx))) +
    xlab(xlab.title) +
    ylab(ylab.title) +
    theme_bw() + theme(aspect.ratio=1) + scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
}
