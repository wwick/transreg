start <- Sys.time()

file.path <- "~\\projects\\transcriptomics\\data\\"
file.name <- "expressionMatrix.csv"

q = 0.970 # 95% confidence



# Dixon's Q Test
Outlier <- function(x) {
  x <- sort(x)
  q.max <- (x[3] - x[2]) / (x[3] - x[1])
  q.min <- (x[2] - x[1]) / (x[3] - x[1])
  if (!is.na(q.max) && q.max > q.min && q.max > q) {
    x[3] <- NA
    return(x)
  }
  if (!is.na(q.min) && q.min > q) {
    x[1] <- NA
    return(x)
  }
  return(x)
}

Graph <- function(fit) {
  title(xlab = 'mRNA (TPM)', ylab = 'RPF (TPM)')
  abline(fit, col = 'blue')

  r2 <- round(summary(fit)$r.squared, 2)
  text(3e-1, 1e5, paste('r2 =', r2))

  m <- round(coef(fit)[2], 3)
  b <- round(coef(fit)[1], 3)
  text(3e-1, 1e4, paste('y = ', m, 'x + ', b, sep = ''))

  f <- summary(fit)$fstatistic
  p <- signif(pf(f[1], f[2], f[3], lower.tail = FALSE), 3)
  text(3e-1, 1e3, paste('p =', p))
}

file <- paste(file.path, file.name, sep = '')
pre.data <- read.csv(file = file, row.names = 1)

for (i in 1:8) {
  apply(cbind(pre.data[,i * 3 - 2], pre.data[,i * 3 - 1], pre.data[,i * 3]),
    1, Outlier)
}

total.data <- pre.data
for (i in 1:24) {
  total.data <- total.data[total.data[,i] != 0,]
}
ribo.data <- total.data[!grepl('VNG', rownames(total.data)),]

par(mfrow = c(2, 4))
xlim <- c(1e-4, 1e6)
ylim <- xlim

for (i in 1:4) {
  x <- (c(total.data[,i * 3 - 2], total.data[,i * 3 - 1],
     total.data[,i * 3]))
  y <- (c(total.data[,i * 3 - 2 + 12], total.data[,i * 3 - 1 + 12],
    total.data[,i * 3 + 12]))
  fit <- lm(log10(y) ~ log10(x))
  plot(x, y, log = 'xy', xlab = '', ylab = '', xlim = xlim, ylim = ylim,
    main = paste('Total Time Point', i))
  Graph(fit)
}

for (i in 1:4) {
  x <- c(ribo.data[,i * 3 - 2], ribo.data[,i * 3 - 1], ribo.data[,i * 3])
  y <- c(ribo.data[,i * 3 - 2 + 12], ribo.data[,i * 3 - 1 + 12], ribo.data[,i * 3 + 12])
  fit <- lm(log10(y) ~ log10(x))
  plot(x, y, log = 'xy', xlab = '', ylab = '', xlim = xlim, ylim = ylim,
    main = paste('Ribosomal Time Point', i))
  Graph(fit)
}

end <- Sys.time()
print(end - start)
