start <- Sys.time()

file.path <- "~\\projects\\transcriptomics\\data\\"
file.name <- "expressionMatrix.csv"

q = 0.970 # 95% confidence
require('outliers')

Col <- function(x) {
  apply(cbind(pre.data[,x * 3 - 2], pre.data[,x * 3 - 1], pre.data[,x * 3]),
    1, Outlier)
}

# Dixon's Q Test
Outlier <- function(x) {
  x <- sort(x)
  q.max <- (x[3] - x[2]) / (x[3] - x[1])
  q.min <- (x[2] - x[1]) / (x[3] - x[1])
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

total.data <- list()
ribo.data <- list()
for (i in 1:4) {
  total.data[[i]] <- cbind(Col(i), Col(i + 4))
  rownames(total.data[[i]]) <- rownames(pre.data)
  total.data[[i]] <- total.data[[i]][total.data[[i]][,1] != 0,]
  total.data[[i]] <- total.data[[i]][total.data[[i]][,2] != 0,]
  ribo.data[[i]] <- total.data[[i]][!grepl('VNG', rownames(total.data[[i]])),]
}

par(mfrow = c(2, 4))
# xlim <- c(1e-4, 1e6)
# ylim <- xlim

Pairs <- function(x) {
  return (x[2] / x[1])
}

outs <- list()
for (i in 1:4) {
  x <- apply(total.data[[i]], 1, Pairs)
  boxplot(x, horizontal = TRUE, main = paste("Time Point", i), xlab = xlab)
  outliersValue <- boxplot.stats(x)$out
  outs[[i]] <- rownames(total.data[[i]])[x %in% outliersValue]
  y[i] <- mean(x)
}

names <- rownames(pre.data)
outs.mat <- matrix(nrow = length(names), ncol = 4)
rownames(outs.mat) <- names
colnames(outs.mat) <- c('TP 1', 'TP 2', 'TP 3', 'TP 4')
for (i in 1:dim(outs.mat)[1]) {
  for (j in 1:4) {
    outs.mat[i,j] <- rownames(outs.mat)[i] %in% outs[[j]]
  }
}
# outs.mat <- outs.mat[outs.mat[,1] == TRUE | outs.mat[,2] == TRUE |
#   outs.mat[,3] == TRUE | outs.mat[,4] == TRUE,]
# print(outs.mat)

# print(y)
# test <- dixon.test(y)
# print(test$statistic)
# print(test$p.value)

ribo.outs <- list()
for (i in 1:4) {
  x <- apply(ribo.data[[i]], 1, Pairs)
  boxplot(x, horizontal = TRUE, main = paste("Time Point", i), xlab = xlab)
  outliersValue <- boxplot.stats(x)$out
  #x <- x[!x %in% outliersValue]
  ribo.outs[[i]] <- rownames(total.data[[i]])[x %in% outliersValue]
  # y[i] <- mean(x)
}
ribo.names <- names[!grepl('VNG', names)]

ribo.mat <- matrix(nrow = length(ribo.names), ncol = 4)
rownames(ribo.mat) <- ribo.names
colnames(ribo.mat) <- c('TP 1', 'TP 2', 'TP 3', 'TP 4')
for (i in 1:dim(ribo.mat)[1]) {
  for (j in 1:4) {
    ribo.mat[i,j] <- rownames(ribo.mat)[i] %in% ribo.outs[[j]]
  }
}
# ribo.mat <- ribo.mat[ribo.mat[,1] == TRUE | ribo.mat[,2] == TRUE |
#   ribo.mat[,3] == TRUE | ribo.mat[,4] == TRUE,]
# print(ribo.mat)

end <- Sys.time()
print(end - start)
