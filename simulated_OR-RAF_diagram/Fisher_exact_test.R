library("parallel")

#### adjustments made to x and y axis for plotting ####

x.adj <- function(x) {ifelse(x < 0.5, log(x) + log(2), -log(1-x) - log(2))}
x.adj.inv <- function(x) {ifelse(x < 0, exp(x-log(2)), 1-exp(-x-log(2)))}
y.adj <- function(y) {y - 0.9}
y.adj.plotly <- function(y) {log(y-.9) - log(0.1)}
y.adj.plotly.inv <- function(y) {exp(y + log(0.1)) + 0.9}

#### grid points for simulation ####

f.lim <- c(1e-4, 1-5e-3)
R.lim <- c(1,110)
resolution <- 100
f.vec <- x.adj.inv(seq(from = x.adj(f.lim)[1], 
                       to = x.adj(f.lim)[2], 
                       length.out = resolution))
R.vec <- y.adj.plotly.inv(seq(from = y.adj.plotly(R.lim)[1], 
                              to = y.adj.plotly(R.lim)[2], 
                              length.out = resolution))

#### simulate the OR-RAF diagram - Fisher's exact test - n2000 phi0.15 ####

n.vec <- c(2e2, 2e3, 2e4, 2e5)
phi.vec <- c(0.15, 0.5, 0.85)
alpha <- 5e-5

temp <- list('0.15' = NULL, '0.5' = NULL, '0.85' = NULL)
simulate_mat <- list('200' = temp, '2000' = temp, '20000' = temp, '2e+05' = temp)

for (n in n.vec) {
  for (phi in phi.vec) {
    temp.mat <- matrix(nrow = length(R.vec), ncol = length(f.vec))
    n1 <- phi * n
    n2 <- n - n1
    
    for (i in 1:length(R.vec)) {
      print(i)
      R <- R.vec[i]
      for (j in 1:length(f.vec)) {
        f <- f.vec[j]
        p1 <- f*R / (f*R + 1 - f)
        p2 <- f
        res <- unlist(mclapply(X = 1:20, function(...) {
          replicate(25, {
            y1 <- rbinom(n = 1, size = n1, prob = p1)
            y2 <- rbinom(n = 1, size = n2, prob = p2)
            test.res <- fisher.test(matrix(c(y1, y2, n1 - y1, n2 - y2), 2))
            test.res$p.value < alpha
          })
        }, mc.cores = 20))
        
        temp.mat[i,j] <- mean(res)
      }
    }
    
    simulate_mat[[as.character(n)]][[as.character(phi)]] <- temp.mat 
  }
}


#### save simulation results ####

save(simulate_mat, file = "~/simulate_diagram_mat.Rdata")
# load(file = "~/simulate_diagram_mat.Rdata")
# save(simulate_mat, file = "./simulate_diagram_mat_p-val_5e-5.Rdata")

