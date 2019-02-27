library("MASS")
library("plotly")
library("dplyr")
if (!require("shinyWidgets")) {
  install.packages("shinyWidgets")
  library("shinyWidgets")
}

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
  
  theta <- c(theta.min, theta.vec.left, theta.vec.right, theta.max)
  R <- c(1e8, OR.vec.left, OR.vec.right, 1e8)
  f <- RAF_control(theta = theta, R = R, phi = phi)
  return(list(theta = theta, R = R, f = f))
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
