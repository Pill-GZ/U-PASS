library("MASS")
library("plotly")
library("dplyr")
library("shinyWidgets")

#### adjustments made to x and y axis for plotting ####

x.adj <- function(x) {ifelse(x < 0.5, log(x) + log(2), -log(1-x) - log(2))}
x.adj.inv <- function(x) {ifelse(x < 0, exp(x-log(2)), 1-exp(-x-log(2)))}
y.adj <- function(y) {y - 0.9}
y.adj.plotly <- function(y) {log(y-.9) - log(0.1)}
y.adj.plotly.inv <- function(y) {exp(y + log(0.1)) + 0.9}

#### convert RAF (risk allele frequency in control group) into marginal frequency ####

RAF_marginal <- function(f, R, phi) {
  # f can be a vector
  # R and phi are scalars
  # check that R > 1
  # check that 0 < phi < 1
  f*R*phi / (f*(R-1)+1) + f*(1-phi)
}

#### convert marginal frequency into RAF (risk allele frequency in control group) ####

RAF_control <- function(theta, R, phi) {
  # theta can be a vector
  # R and phi are scalars
  # check that R > 1
  # check that 0 < theta < 1
  ( phi*R - phi - R*theta + theta + 1 - 
      sqrt( (-phi*R + phi + R*theta - theta -1)^2 - 4*theta*(phi*R-phi-R+1) ) ) /
    (2*(phi-1)*(R-1))
}

#### analytical solution of the signal size ####

signal.size.ana.sol <- function(theta, phi, R){
  if (length(theta) > 1 && length(phi) > 1) {
    warning("only one of the marginals should be a vector!")
    return(NA)
  } else if (R==1) {
    lam <- rep(0, max(length(theta), length(phi)))
  } else {
    A <- theta*phi*(1-theta)*(1-phi)
    B <- (theta*phi+(1-theta)*(1-phi)) 
    C <- (theta*(1-phi)+(1-theta)*phi)
    w2 <- 1/(2*(R-1))^2 * 
      (B+C*R - sqrt((B+C*R)^2 - 4*A*(R-1)^2))^2 * 
      (1/(theta*phi)+1/((1-theta)*phi)+1/(theta*(1-phi))+1/((1-theta)*(1-phi)))
  }
  return(w2)
}

## analytical solution of signal size as a function of conditional (f)

signal.size.ana.sol.f <- function(f, phi, R){
  if (length(f) > 1 && length(phi) > 1) {
    warning("Only one of the marginals should be a vector!")
    return(NA)
  } else {
    theta <- (f*R*phi)/(f*R+1-f) + f*(1-phi)
    if (length(theta) != length(phi)) {
      warning("Length of the marginals are unequal!")
      return(NA)
    } else if (R==1) {
      w2 <- rep(0, max(length(theta), length(phi)))
    } else {
      A <- theta*phi*(1-theta)*(1-phi)
      B <- (theta*phi+(1-theta)*(1-phi)) 
      C <- (theta*(1-phi)+(1-theta)*phi)
      w2 <- 1/(2*(R-1))^2 * 
        (B+C*R - sqrt((B+C*R)^2 - 4*A*(R-1)^2))^2 * 
        (1/(theta*phi)+1/((1-theta)*phi)+1/(theta*(1-phi))+1/((1-theta)*(1-phi)))
    }
  }
  return(w2)
}


#### solve for (f, R), or (theta, R) combinations that achieve specified signal size at given phi ####

OR.finder <- function(phi = 1/2, signal.size) {
  theta.min <- (signal.size*phi)/((1-phi)+signal.size*phi)
  theta.mid <- phi
  theta.max <- (phi)/(signal.size*(1-phi)+phi)
  OR.root.finder <- function(OR, theta, phi, signal.size) {
    signal.size.ana.sol(theta, phi, OR) - signal.size
  }
  # left branch
  theta.vec.left <- exp(log(theta.min) + 1:100/100 * (log(theta.mid) - log(theta.min)))
  OR.vec.left <- vector(mode = "numeric", length = 100)
  for (i in 1:length(theta.vec.left)) {
    solution <- uniroot(f = OR.root.finder, interval = c(1.01, 1e6),
                        theta = theta.vec.left[i], phi = phi, signal.size)
    OR.vec.left[i] <- solution$root
  }
  # right branch
  theta.vec.right <- rev(1 - exp(log(1-theta.max) + 1:100/100 * (log(1-theta.mid) - log(1-theta.max))))
  OR.vec.right <- vector(mode = "numeric", length = 100)
  for (i in 1:length(theta.vec.right)) {
    solution <- uniroot(f = OR.root.finder, interval = c(1.01, 1e6),
                        theta = theta.vec.right[i], phi = phi, signal.size)
    OR.vec.right[i] <- solution$root
  }
  # combined
  theta <- c(theta.min, theta.vec.left, theta.vec.right, theta.max)
  R <- c(1e8, OR.vec.left, OR.vec.right, 1e8)
  f <- RAF_control(theta = theta, R = R, phi = phi)
  return(list(theta = theta, R = R, f = f))
}

#### required columns for user upload data ####

list_accessed_columns <- c("INITIAL.SAMPLE.SIZE", "X95..CI..TEXT.", "P.VALUE", "MAPPED_GENE", 
                           "PUBMEDID", "REPORTED.GENE.S.", "DATE.ADDED.TO.CATALOG", "STUDY", 
                           "JOURNAL", "FIRST.AUTHOR", "DISEASE.TRAIT", "REPLICATION.SAMPLE.SIZE")
list_required_columns <- c("OR.or.BETA", "RISK.ALLELE.FREQUENCY")

#### read datasets downloaded from EBI ####

# function to read and pre-process data
read_data <- function(dataset_directory_and_name) {
  tryCatch(
    {
      dataset <- read.delim(dataset_directory_and_name, 
                            na.strings = c("NA", "NR"), stringsAsFactors = F)
      missing_required_columns <- list_required_columns[!(list_required_columns %in% colnames(dataset))]
      missing_accessed_columns <- list_accessed_columns[!(list_accessed_columns %in% colnames(dataset))]
      if (length(missing_required_columns) > 0) {
        stop(paste("Required column", missing_required_columns, "is missing! "))
      }
      if(length(missing_accessed_columns) > 0) {
        stop(paste("Column", missing_accessed_columns, "is missing! "))
      }
      # filter out non-numerical stuff in RAF
      dataset$RISK.ALLELE.FREQUENCY <- as.numeric(gsub("\\(.*\\)", "", dataset$RISK.ALLELE.FREQUENCY))
      # separate beta's from OR's
      beta.indicator <- grepl(pattern = "increase|decrease", x = dataset$X95..CI..TEXT.)
      dataset$OR <- ifelse(beta.indicator, exp(dataset$OR.or.BETA), dataset$OR.or.BETA)
      dataset$index <- 1:nrow(dataset)
      dataset
    }, 
    # print error in a message box if a parsing error occurs
    error = function(e) {
      showModal(modalDialog(title = "Unable to read dataset!", toString(e)))
      return(NULL)
    }
  ) # end of try-catch data read
}

#### determine intersection for power calculation ####

determine.intersection <- function(x, y, target) {
  if (target > max(y)) {
    return(FALSE)
  } else {
    j <- min(which(y > target))
    i <- j - 1
    return(ceiling(exp(log(x[i]) + (target - y[i])/(y[j]-y[i])*(log(x[j]) - log(x[i])))))
  }
}

#### function to calculate rare variant region based on sample sizes ####

rare.variant.curves.calculator <- function(n1, n2, RV.threshold.left, RV.threshold.right,
                                           f.lim = c(1e-4, 1-5e-3), R.lim = c(1,110)) {
  n <- n1 + n2
  phi <- n1 / n
  RV.resolution <- 100
  # left branch
  R.min <- 1
  R.max <- R.lim[2]
  f.rare.variant.root.finder.left <- function(f.RAF, OR, phi, n, rare.variant.threshold) {
    f.RAF * n * (OR*phi/(f.RAF*OR+(1-f.RAF)) + (1-phi)) - rare.variant.threshold
  }
  OR.vec.left <- exp(log(R.min) + 0:RV.resolution/RV.resolution * (log(R.max) - log(R.min)))
  f.vec.left <- vector(mode = "numeric", length = RV.resolution+1)
  for (i in 1:length(OR.vec.left)) {
    if (f.rare.variant.root.finder.left(f.RAF = f.lim[1]/2, OR = OR.vec.left[i], 
                                        phi, n, RV.threshold.left) > 0) {
      f.vec.left[i] <- NA
    } else {
      solution <- uniroot(f = f.rare.variant.root.finder.left, interval = c(f.lim[1]/2, 0.90),
                          OR = OR.vec.left[i], phi, n, RV.threshold.left,
                          tol = .Machine$double.eps^0.5)
      f.vec.left[i] <- solution$root
    }
  }
  # right branch
  R.min <- 1
  R.max <- R.lim[2]
  f.rare.variant.root.finder.right <- function(f.RAF, OR, phi, n, rare.variant.threshold) {
    (1-f.RAF) * n * (OR*phi/(f.RAF*OR+(1-f.RAF)) + (1-phi)) - rare.variant.threshold
  }
  OR.vec.right <- exp(log(R.min) + 0:RV.resolution/RV.resolution * (log(R.max) - log(R.min)))
  f.vec.right <- vector(mode = "numeric", length = RV.resolution+1)
  for (i in 1:length(OR.vec.right)) {
    if (f.rare.variant.root.finder.right(f.RAF = f.lim[2], OR = OR.vec.right[i], 
                                         phi, n, RV.threshold.right) > 0) {
      f.vec.right[i] <- NA
    } else {
      solution <- uniroot(f = f.rare.variant.root.finder.right, interval = c(0.1, f.lim[2]),
                          OR = OR.vec.right[i], phi, n, RV.threshold.right,
                          tol = .Machine$double.eps^0.5)
      f.vec.right[i] <- solution$root
    }
  }
  # combined
  return(list(f.vec.left = f.vec.left, OR.vec.left = OR.vec.left, 
              f.vec.right = f.vec.right, OR.vec.right = OR.vec.right))
}

#### function to calculate minimum number of counts needed to calibrate Fisher's exact test ####

minimum.RV.Fishers.test <- function(n1, n2, p.val.threhold) {
  # left-hand side
  O21 <- 0
  O22 <- n2
  O11 <- 1
  O12 <- n1- O11
  
  test.res <- fisher.test(matrix(c(O11, O21, O12, O22), 2), alternative = "greater")
  while (test.res$p.value > p.val.threhold) {
    O11 <- O11 + 1
    O12 <- n1- O11
    test.res <- fisher.test(matrix(c(O11, O21, O12, O22), 2), alternative = "greater")
  }
  left.threshold <- O11
  
  # right-hand side
  O22 <- 1
  O21 <- n2 - O22
  O12 <- 0
  O11 <- n1
  
  test.res <- fisher.test(matrix(c(O11, O21, O12, O22), 2), alternative = "greater")
  while (test.res$p.value > p.val.threhold) {
    O22 <- O22 + 1
    O21 <- n2- O22
    test.res <- fisher.test(matrix(c(O11, O21, O12, O22), 2), alternative = "greater")
  }
  right.threshold <- O22
  
  return(list(left = left.threshold, right = right.threshold))
}

#### set plot limits and resolution ####

f.lim <- c(1e-4, 1-5e-3)
R.lim <- c(1,110)
resolution <- 100
