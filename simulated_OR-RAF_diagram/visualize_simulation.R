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



#### calculate RV curves ####

rare.variant.curves <- function(n1, n2, p.val.threhold) {
  RV.thresholds <- minimum.RV.Fishers.test(n1, n2, p.val.threhold)
  rare.variant.curves.calculator(n1, n2, RV.thresholds$left, RV.thresholds$right, f.lim, R.lim)
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

#### generate theoreical power matrix for heatmaps ####

signal.mat <- function(n1, n2, phi) {
 temp.mat <- matrix(nrow = length(R.vec), ncol = length(f.vec))
  for (i in 1:length(R.vec)) {
    for (j in 1:length(f.vec)) {
      temp.mat[i,j] <- signal.size.ana.sol.f(f.vec[j], phi, R.vec[i]) * (n1 + n2)
    }
  }
  temp.mat
}

power.mat <- function(n1, n2, phi, p.val.threhold) {
  cutoff <- qchisq(p = 1 - p.val.threhold, df = 1, ncp = 0, lower.tail = T)
  pchisq(q = cutoff, df = 1, ncp = signal.mat(n1, n2, phi), lower.tail = F)
}

#### load results from simulations ####

#load(file = "~/Research_office/MyWork/U-PASS/Shiny/power/simulated_OR-RAF_diagram/simulate_diagram_mat_p5e-5.Rdata")
load(file = "~/Research_office/MyWork/U-PASS/Shiny/power/simulated_OR-RAF_diagram/simulate_diagram_mat_p5e-8.Rdata")

#### visualize simulation results in a OR-RAF diagram ####

library("plotly")

n <- 200000
phi <- 0.85
n1 <- n * phi; n2 <- n * (1-phi)
p.val.threhold <- 5e-8
plot.option <- 'difference'

if (plot.option == 'empirical') {
  plot.mat <- simulate_mat[[as.character(n)]][[as.character(phi)]]
} else if (plot.option == 'theoretical') {
  plot.mat <- power.mat(n1 = n1, n2 = n2, phi = phi, p.val.threhold = p.val.threhold)
} else if (plot.option == 'difference') {
  plot.mat  <- abs(power.mat(n1 = n1, n2 = n2, phi = phi, p.val.threhold = p.val.threhold) -
    simulate_mat[[as.character(n)]][[as.character(phi)]])
}

p <- plot_ly(x = x.adj(f.vec), y = y.adj.plotly(R.vec), 
             z = plot.mat, zmin = 0, zmax = 1,
             colors = "Greys", type = "contour",
             reversescale = ifelse(plot.option == 'difference', F, T),
             hoverinfo = "none", width = 750, height = 700,
             contours = list(showlabels = TRUE)) %>%
  # config(displayModeBar = F) %>% 
  layout(xaxis = list(range = x.adj(f.lim), # fixedrange=TRUE,
                      tickvals = x.adj(c(0.5, 10^(-1:-4), 1-10^(-1:-2))),
                      ticktext = c(0.5, 0.1, (10^{-2}), (10^{-3}), (10^{-4}),0.9, 0.99),
                      zeroline = FALSE, tickfont = list(size = 20),
                      showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
         yaxis = list(range = y.adj.plotly(R.lim), # fixedrange=TRUE,
                      tickvals = y.adj.plotly(c(1, 2, 5, 10, 20, 50, 100)),
                      ticktext = c(1, 2, 5, 10, 20, 50, 100), 
                      zeroline = FALSE, tickfont = list(size = 20),
                      showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
         showlegend = FALSE, 
         margin = list(l = 50, r = 20, b = 80, t = 50, pad = 4)) %>%
  add_annotations(yref = "paper", xref = "paper", y = 1.08, x = -0.07, 
                  text = "odds ratio", showarrow = F, font=list(size = 25)) %>%
  add_annotations(yref = "paper", xref = "paper", y = -0.12, x = 1, 
                  text = "risk allele frequency (in control group)", 
                  showarrow = F, font = list(size = 25)) %>%
  add_annotations(yref = "paper", xref = "paper", y = 1.04,  
                  text = ifelse(plot.option == 'difference', "diff", "power"), 
                  x = ifelse(plot.option == 'difference', 1.1, 1.15),
                  showarrow = F, font = list(size = 25)) %>%
  add_annotations(yref = "paper", xref = "paper", y = 0.8, x = 0.5, 
                  text = "<b>Cases/Controls</b>", showarrow = F, 
                  font=list(size = 40, color = toRGB("grey70"))) %>%
  add_annotations(yref = "paper", xref = "paper", y = 0.7, x = 0.5, 
                  text = paste0("<b>", format(n1/2, scientific = FALSE), 
                                "/", format(n2/2, scientific = FALSE), "</b>"), showarrow = F, 
                  font=list(size = 40, color = toRGB("grey70")))

RV.curves <- rare.variant.curves(n1, n2, p.val.threhold)

# if (plot.option == 'theoretical') {
#   p <- p %>% 
#     add_trace(x = x.adj(RV.curves$f.vec.left),
#               y = y.adj.plotly(RV.curves$OR.vec.left),
#               line = list(dash = 'dash', width = 2), type = "scatter", mode = "lines",
#               color = I("red"), hoverinfo = 'text',
#               text = paste('rare variant zone to the left',
#                            '\n', 'risk allele count < ',
#                            minimum.RV.Fishers.test(n1, n2, p.val.threhold)$left),
#               inherit = F) %>%
#     add_trace(x = x.adj(RV.curves$f.vec.right),
#               y = y.adj.plotly(RV.curves$OR.vec.right),
#               line = list(dash = 'dash', width = 2), type = "scatter", mode = "lines",
#               color = I("red"), hoverinfo = 'text',
#               text = paste('rare variant zone to the right',
#                            '\n', 'non-risk allele count < ', 
#                            minimum.RV.Fishers.test(n1, n2, p.val.threhold)$right),
#               inherit = F)
# }

p

#### calculate deviations between theory and simulation ####


n.vec <- 2*10^(2:5)
phi.vec <- c(0.05, 0.15, 0.25, 0.5, 0.85)
MAD.mat <- matrix(nrow = length(n.vec), ncol = length(phi.vec))
rownames(MAD.mat) <- n.vec
colnames(MAD.mat) <- phi.vec
MAD.sd.mat <- matrix(nrow = length(n.vec), ncol = length(phi.vec))
rownames(MAD.sd.mat) <- n.vec
colnames(MAD.sd.mat) <- phi.vec

p.val.threhold <- 5e-8
for (n in n.vec) {
  for (phi in phi.vec) {
    n1 <- n * phi; n2 <- n * (1-phi)
    MAD <- abs(power.mat(n1 = n1, n2 = n2, phi = phi, p.val.threhold = p.val.threhold) -
                 simulate_mat[[as.character(n)]][[as.character(phi)]])
    MAD.mat[as.character(n), as.character(phi)]  <- mean(MAD)
    MAD.sd.mat[as.character(n), as.character(phi)]  <- sd(MAD)
  }
}


