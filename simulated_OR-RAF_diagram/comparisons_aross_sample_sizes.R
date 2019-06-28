
library("parallel")
library("abind")

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

#### Simulate power of statistical tests ####

n.vec <- c(as.vector(outer(c(1,1.5,2:9), 10^(2:5))), 1e6)
alpha.vec <- c(5e-5, 5e-8)
phi.vec <- c(0.05, 0.15, 0.25, 0.35, 0.5, 0.85)
R.vec <- c(1.05, 1.1, 1.2, 1.5, 3.0)
f.vec <- c(0.0001, 0.001, 0.01, 0.1, 0.5, 0.9)

RAF <- list('1e-04' = NULL, '0.001' = NULL, '0.01' = NULL, 
            '0.1' = NULL, '0.5' = NULL, '0.9' = NULL)
OR <- list('1.05' = RAF, '1.1' = RAF, '1.2' = RAF, '1.5' = RAF, '3.0' = RAF)
simulate_mat <- list('0.05' = OR, '0.15' = OR, '0.25' = OR, 
                     '0.35' = OR, '0.5' = OR, '0.85' = OR)

#### one-replication of each test ####

# Fisher's Exact test
test.FE <- function(n1, n2, p1, p2) {
  y1 <- rbinom(n = 1, size = n1, prob = p1)
  y2 <- rbinom(n = 1, size = n2, prob = p2)
  data.tabulated <- matrix(c(y1, y2, n1 - y1, n2 - y2), 2)
  fisher.test(data.tabulated)$p.value
}

# Welch's t-test
test.WT <- function(n1, n2, p1, p2) {
  y1 <- rbinom(n = 1, size = n1, prob = p1)
  y2 <- rbinom(n = 1, size = n2, prob = p2)
  if (y1 + y2 < 5) {
    test.res.WT <- 1
  } else {
    data.unrolled.x <- c(rep(1, y1), rep(0, n1 - y1))
    data.unrolled.y <- c(rep(1, y2), rep(0, n2 - y2))
    test.res.WT <- t.test(x = data.unrolled.x, y = data.unrolled.y, var.equal = F)$p.value
  }
  test.res.WT
}

# Student t-test
test.ST <- function(n1, n2, p1, p2) {
  y1 <- rbinom(n = 1, size = n1, prob = p1)
  y2 <- rbinom(n = 1, size = n2, prob = p2)
  if (y1 + y2 == 0) {
    test.res.ST <- 1
  } else {
    data.unrolled.x <- c(rep(1, y1), rep(0, n1 - y1))
    data.unrolled.y <- c(rep(1, y2), rep(0, n2 - y2))
    test.res.ST <- t.test(x = data.unrolled.x, y = data.unrolled.y, var.equal = T)$p.value 
  }
  test.res.ST
}

# Chi-square test
test.CS <- function(n1, n2, p1, p2) {
  y1 <- rbinom(n = 1, size = n1, prob = p1)
  y2 <- rbinom(n = 1, size = n2, prob = p2)
  if (y1 + y2 == 0) {
    test.res.CS <- 1
  } else {
    data.tabulated <- matrix(c(y1, y2, n1 - y1, n2 - y2), 2)
    test.res.CS <- chisq.test(data.tabulated, correct = F)$p.value 
  }
  test.res.CS
}

# Likelihood ratio test
test.LR <- function(n1, n2, p1, p2) {
  y1 <- rbinom(n = 1, size = n1, prob = p1)
  y2 <- rbinom(n = 1, size = n2, prob = p2)
  data.tabulated <- matrix(c(y1, y2, n1 - y1, n2 - y2), 2)
  LR.model <- glm(formula = data.tabulated ~ as.factor(c(1, 0)), family = binomial())
  LR.model0 <- glm(formula = data.tabulated ~ 1, family = binomial())
  anova(LR.model0, LR.model, test ="Chisq")$`Pr(>Chi)`[2]
}

#### Simulation starts here ####

# !!! change the number of cores in production !!! #
num.cores <- 2

alpha <- alpha.vec[2]
cutoff <- qchisq(p = 1 - alpha, df = 1, ncp = 0, lower.tail = T)

for (phi in 0.5) { # phi.vec) {
  for (R.idx in 3) { #1:length(R.vec)) {
    R <- R.vec[R.idx]
    for (f.idx in c(4, 5)) { # 1:length(f.vec)) {
      f <- f.vec[f.idx]
      print(paste0("phi = ", phi, ", R = ", R, ", f = ", f))
      p1 <- f*R / (f*R + 1 - f)
      p2 <- f
      temp.mat <- matrix(nrow = 6, ncol = length(n.vec))
      for (n.idx in 1:length(n.vec)) {
        n <- n.vec[n.idx] * 2
        n1 <- floor(phi * n)
        n2 <- n - n1
        res <- mclapply(X = 1:num.cores, function(...) {
          replicate(5, {
            # Fisher's Exact test
            test.res.FE <- test.FE(n1, n2, p1, p2)
            # Welch's t-test
            test.res.WT <- test.WT(n1, n2, p1, p2)
            # Student t-test
            test.res.ST <- test.ST(n1, n2, p1, p2)
            # Chi-square test
            test.res.CS <- test.CS(n1, n2, p1, p2)
            # Likelihood ratio test
            test.res.LR <- test.LR(n1, n2, p1, p2)
            c(test.res.FE, test.res.WT, test.res.ST, test.res.CS, test.res.LR) < alpha
          })
        }, mc.cores = num.cores)
        temp.mat[1:5, n.idx] <- apply(abind(res, along = 3), c(1), mean)
      } # end of n
      # calculate theoretical power. First the signal size per sample
      signal.size.per.sample <- signal.size.ana.sol.f(f = f, phi = phi, R = R)
      # then calculate power as a function of total sample size (assuming one allele pair per subject)
      variable.power.theoretical <- pchisq(q = cutoff, df = 1 , lower.tail = F, ncp = signal.size.per.sample * 2 * n.vec)
      temp.mat[6,] <- variable.power.theoretical
      simulate_mat[[as.character(phi)]][[as.character(R)]][[as.character(f)]] <- temp.mat 
    } # end of RAF
  } # end of OR
} # end of phi

#### save simulation results ####

# save(simulate_mat, file = "~/comparing_tests.Rdata")

load(file = "~/Research_office/MyWork/U-PASS/Shiny/power/simulated_OR-RAF_diagram/comparing_tests.Rdata")
# save(simulate_mat, file = "./simulate_diagram_mat_p-val_5e-5.Rdata")

#### plots ####


for (phi in phi.vec) {
  for (R in R.vec) {
    for (f in f.vec) {
      png(filename = paste0(gsub(x = paste0("phi=", format(phi, nsmall = 2), 
                                            "_R=", format(R, nsmall = 2), "_f=", f),
                                 pattern = "\\.", replacement = ""), ".png"),
          width = 300, height = 200)
      par(mar = c(2,2.5,0,0))
      plot(x = c(1, length(n.vec)), y = c(0, 1), type = 'n',
           ylim = c(0,1), axes = F, ylab = "")
      axis(side = 1, at = seq(1, length(n.vec), 10), 
           labels = sapply(2:6, FUN = function(x) {as.expression(bquote(10^.(x)))}))
      mtext(text = "n =", side = 1, line = 0.9, at = -3)
      axis(side = 2, at = seq(0, 1, 0.2), 
           labels = format(seq(0, 1, 0.2), digits = 1), las = 1)
      if (simulate_mat[[as.character(phi)]][[as.character(R)]][[as.character(f)]][6,21] < 0.5) {
        text.x <- 1; text.pos <- 4
      } else {
        text.x <- 41; text.pos <- 2
      }
      text(x = text.x, y = 0.8, labels = bquote(phi~'='~.(format(phi, nsmall = 1))), 
           pos =  text.pos, cex = 2, col = grey(0.7))
      text(x = text.x, y = 0.6, labels = bquote(R~'='~.(format(R, nsmall = 1))), 
           pos = text.pos, cex = 2, col = grey(0.7))
      text(x = text.x, y = 0.4, labels = bquote(f~'='~.(format(f, scientific = F))), 
           pos = text.pos, cex = 2, col = grey(0.7))
      matplot(t(simulate_mat[[as.character(phi)]][[as.character(R)]][[as.character(f)]][c(1,3:6),]), 
              ylim = c(0,1), pch = c(1,3:6), col =  c(1,3:6),
              type = 'b', cex = 0.8, axes = F, ylab = "", add = T)
      dev.off()
    }
  }
}
