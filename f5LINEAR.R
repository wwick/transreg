start <- Sys.time() 

#file.path <- "~\\projects\\transcriptomics\\data\\"
file.path <- "~/transreg-master/"
file.name <- "expressionMatrix.csv"

file <- paste(file.path, file.name, sep = '')
pre.data <- read.csv(file = file, row.names = 1)

q = 0.970 # 95% confidence

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

total.data <- list()
ribo.data <- list()
for (i in 1:4) {
  total.data[[i]] <- cbind(Col(i * 2 - 1), Col(i * 2))
  rownames(total.data[[i]]) <- rownames(pre.data)
  total.data[[i]] <- total.data[[i]][total.data[[i]][,1] != 0,]
  total.data[[i]] <- total.data[[i]][total.data[[i]][,2] != 0,]
  ribo.data[[i]] <- total.data[[i]][!grepl('VNG', rownames(total.data[[i]])),]
}



#PLOTS BEGIN HERE

#mRNA-Data
par(pty="s",mfrow=c(2,4))
for (i in 1:4) {
  x <- total.data[[i]][,1]
  y <- total.data[[i]][,2]
  fit <- lm(y ~ x) #lm(log10(y) ~ log10(x))
  #print(summary(fit))

  plot(x, y, xlim = c(0,1000), ylim = c(0,1000),
       xlab = 'mRNA', ylab = 'RPF',
       main = paste('Total Time Point', i))
  
  m <- signif(coef(fit)[2],3)
  b <- signif(coef(fit)[1],3)
  rsq <- summary(lm(x~y)$r.squared)
  r2 <- round(summary(fit)$adj.r.squared, 2)
  text(300,750, paste('y=', m, 'x + ', b)) 
  text(300,900, paste('r2= ', r2))
  print(r2)
  }



#Ribo-Data

for (i in 1:4) {
  
  x <- ribo.data[[i]][,1]
  y <- ribo.data[[i]][,2]
  fit <- lm(y ~ x)
  #print(summary(fit))
  
  #check axes
  plot(x, y, xlim = c(0,1000), ylim = c(0,1000),
       xlab = 'mRNA', ylab = 'RPF',
       main = paste('Ribosomal Time Point', i))
  r2 <- round(summary(fit)$r.squared, 2)
  m <- signif(coef(fit)[2],3)
  b <- signif(coef(fit)[1],3)
  text(300,750, paste('y=', m, 'x + ', b)) # 1e2, 1e3, r2)
  text(300,900, paste('r2= ', r2) )
  print(r2)
  }

end <- Sys.time()
print(end - start)
