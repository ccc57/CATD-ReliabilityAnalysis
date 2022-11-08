library(ggplot2)
library(MESS)


powerfig <- data.frame(seq(0.01,1,.01),seq(0.01,1,.01))
colnames(powerfig) <- c('ICC','es')
adj_es <- matrix(nrow=length(powerfig$ICC),ncol=length(powerfig$es))
  
for(i in 1:length(powerfig$ICC)){
  adj_es[i,] <- powerfig$es * sqrt(powerfig$ICC[i])
}

ss <- matrix(nrow = length(powerfig$ICC),ncol = length(powerfig$es))

runtest <- function(d) {
  tryCatch(
    expr = {
      result<-power_t_test(delta=d, sd=1, sig.level =0.05, power=0.8, ratio=1, sd.ratio=1, type="two.sample", alternative="two.sided", df.method="classical")
      return(ceiling(result$n))
    },
    error = function(e){
      message('an error occurred')
      return(NA)
    }
  )
}

for(i in 1:nrow(ss)){
  for(j in 1:ncol(ss)){
    ss[i,j] <- runtest(adj_es[i,j])
  }
}
filled.contour(ss*2, ylim = c(0,1),levels = 1.15^(1:100), plot.title = title(main='Sample Size Required for 2-sample t-Test at 80% Power'), key.title= title('Sample Size'))

filled.contour(log(ss*2), ylim = c(0,1),xlab='ICC',ylab='True Effect Size', 
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(log(ss*2), add = TRUE)
               })

filled.contour(volcano, plot.axes = {
  axis(1)
  axis(2)
  contour(volcano, add = TRUE, lwd = 2)
}
)
