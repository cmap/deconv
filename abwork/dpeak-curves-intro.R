set.seed(1)
par(mar = c(0,0,1,1))
n <- 1e4
n.k <- 4
means <- sample(1:n.k)
sds <- rep(0.5, n.k)
plot(NA,NA, xlim = c(1 - max(sds) * 4, n.k + max(sds) * 4), ylim = c(0, n.k + 1), ann = F, axes = F)
abline(h = 1:n.k, lty = 2, col = gray(.90))
out <- sapply (1:n.k, function(k) {
  y <- rnorm(n, mean = means[k], sd = sds[k])
  y.pdens <- density(y)
  lines(y.pdens$x, y.pdens$y + k, lwd = 2)
  return(y)
})


x.scale <- range(out)
par(mfrow=c(3,3))
for( i in 1:4) {
  for (j in 1:4) {
    if (j > i) {
      x <- out[,c(i, j)]
      hist(x, breaks = seq(x.scale[1], x.scale[2]+1, by = 0.5), 
        , col = rgb(247,148,30, maxColorValue = 250), border = "white", ann =F, axes = F)
      #abline(v = apply(x, 2, mean), lwd = 2)
    }
  }
}

curve(dnorm(x, sd = 0.5) + dnorm(x, mean = 2, sd = 0.5), from = -5, to = 5, lwd = 5)
curve(dnorm(x, sd = 0.5), from = -5, to = 5, col = 2, add = T, lty = 2)
curve(dnorm(x, mean = 2, sd = 0.5), from = -5, to = 5, col = 2, add = T, lty = 2)
