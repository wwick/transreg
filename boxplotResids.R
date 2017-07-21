start <- Sys.time()

file.path <- "~\\projects\\transcriptomics\\data\\"
file.name <- "expressionMatrix.csv"

q = 0.970 # 95% confidence

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

par(mfrow = c(2, 2))

for (i in 1:4) {
  x <- total.data[[i]][,1]
  y <- total.data[[i]][,2]
  fit <- lm(log10(y) ~ log10(x))
  boxplot(resid(fit), main = paste('Total Time Point', i), horizontal = TRUE, xlab = 'Residual')

}

# for (i in 1:4) {
#   x <- ribo.data[[i]][,1]
#   y <- ribo.data[[i]][,2]
#   fit <- lm(log10(y) ~ log10(x))
#   boxplot(resid(fit), main = paste('Ribosomal Time Point', i))
#
# }

end <- Sys.time()
print(end - start)
