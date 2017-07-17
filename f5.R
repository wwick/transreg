file.path <- "~/Documents/projects/transcriptomics/data/"
file.name <- "expressionMatrix.csv"

file <- paste(file.path, file.name, sep = '')
pre.data <- read.csv(file = file, row.names = 1)

q = 0.970 # 95% condidence

Cols <- function(x) {
  apply(cbind(pre.data[,x * 3 - 2], pre.data[,x * 3 - 1], pre.data[,x * 3]),
    1, Outlier)
}

# Dixon's Q Test
Outlier <- function(x) {
  x <- sort(x)
  q.max = (x[3] - x[2]) / (x[2] - x[1])
  q.min = (x[2] - x[1]) / (x[3] - x[2])
  if (!is.na(q.max) && q.max > q.min && q.max > q) {
    x[3] <- NA
    return(mean(x, na.rm = TRUE))
  }
  if (!is.na(q.min) && q.min > q) {
    x[1] <- NA
    return(mean(x, na.rm = TRUE))
  }
  return(mean(x, na.rm = TRUE))
}

start <- Sys.time()

mrna <- cbind(Cols(1), Cols(2), Cols(3), Cols(4))
rpf <- cbind(Cols(5), Cols(6), Cols(7), Cols(8))

rownames(mrna) <- rownames(pre.data)
rownames(rpf) <- rownames(pre.data)
colnames <- c('Time Point 1', 'Time Point 2', 'Time Point 3', 'Time Point 4')
colnames(mrna) <- colnames
colnames(rpf) <- colnames

end <- Sys.time()
print(end - start)

#par(mfrow = c(2,2))

for (i in 1:1) {
  fit = lm(log(mrna[,i]) ~ log(rpf[,i]), na.action = na.omit)
  plot(mrna[,i], rpf[,i], log = 'xy', xlab = '', ylab = '')
  abline(fit, col = 'red')
}
