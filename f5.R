start <- Sys.time()

file.path <- "~\\projects\\transcriptomics\\data\\"
file.name <- "expressionMatrix.csv"

file <- paste(file.path, file.name, sep = '')
pre.data <- read.csv(file = file, row.names = 1)

q = 0.970 # 95% condidence

Col <- function(x) {
  apply(cbind(pre.data[,x * 3 - 2], pre.data[,x * 3 - 1], pre.data[,x * 3]),
    1, Outlier)
}

# Dixon's Q Test
Outlier <- function(x) {
  x <- sort(x)
  q.max <- (x[3] - x[2]) / (x[2] - x[1])
  q.min <- (x[2] - x[1]) / (x[3] - x[2])
  if (!is.na(q.max) && q.max > q.min && q.max > q) {
    x[3] <- NA
    return(mean(x, na.rm = TRUE))
  }
  if (!is.na(q.min) && q.min > q) {
    x[1] <- NA
    return(mean(x, na.rm = TRUE))
  }
  return(mean(x))
}

post.data <- Col(1)
for (i in 2:8) {
  post.data <- cbind(post.data, Col(i))
}
rownames(post.data) <- rownames(pre.data)
for (i in 1:8) {
  post.data <- post.data[post.data[,i] != 0,]
}

par(mfrow = c(2, 2))

for (i in 1:4) {
  fit <- lm(log10(post.data[,i + 4]) ~ log10(post.data[,i]))
  plot(post.data[,i], post.data[,i + 4], log = 'xy', xlab = 'mRNA',
    ylab = 'RPF', main = paste('Time Point', i))
  abline(fit, col = 'red')
  r2 <- round(summary(fit)$r.squared, 2)
  text(1e-2, 1e4, r2)
  print(expression(r^'2'))
}

ribo.data <- post.data[]

end <- Sys.time()
print(end - start)
