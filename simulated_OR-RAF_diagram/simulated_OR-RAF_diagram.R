
library("dplyr")
library("plotly")

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

#### simulate the OR-RAF diagram ####

n <- 2000
phi <- 0.85

temp.mat <- matrix(nrow = length(R.vec), ncol = length(f.vec))
for (i in 1:length(R.vec)) {
  print(i)
  R <- R.vec[i]
  for (j in 1:length(f.vec)) {
    f <- f.vec[j]
    p <- c(f*R*phi/(f*R+1-f), (1-f)*phi/(f*R+1-f), f*(1-phi), (1-f)*(1-phi))
    res <- replicate(1000, {
      x <- rmultinom(1, n, p)
      if (sum(x[c(1,3)]) > 2 && sum(x[c(2,4)]) > 2) {
        t.statistics <- t.test(x = c(rep(1, x[1]), rep(0, x[3])), 
                               y = c(rep(1, x[2]), rep(0, x[4])))
        t.statistics$p.value < 0.05/1000
      } else {
        FALSE
      }
    })
    temp.mat[i,j] <- mean(res)
  }
}

#### display the simulated OR-RAF diagram ####

# load("~/Research_office/MyWork/Chi-squared_model/Shiny/power/simulated_OR-RAF_diagram/drive-download-20190319T013843Z-001/temp.mat.1-33.Rdata")
# mat.part1 <- temp.mat
# load("~/Research_office/MyWork/Chi-squared_model/Shiny/power/simulated_OR-RAF_diagram/drive-download-20190319T013843Z-001/temp.mat.34-66.Rdata")
# mat.part2 <- temp.mat
# load("~/Research_office/MyWork/Chi-squared_model/Shiny/power/simulated_OR-RAF_diagram/drive-download-20190319T013843Z-001/temp.mat.67-100.Rdata")
# mat.part3 <- temp.mat
# 
# simulate_mat_n20000_phi05 <- matrix(nrow = 100, ncol = 100)
# simulate_mat_n20000_phi05[1:33,] <- mat.part1[1:33,]
# simulate_mat_n20000_phi05[34:66,] <- mat.part2[34:66,]
# simulate_mat_n20000_phi05[67:100,] <- mat.part3[67:100,]

simulate_mat_n2000_phi085 <-temp.mat 

# save(simulate_mat_n20000_phi05, simulate_mat_n2000_phi085, simulate_mat_n2000_phi015, file = "~/Research_office/MyWork/Chi-squared_model/Shiny/power/simulated_OR-RAF_diagram/simulate_diagram_mat.Rdata")
# load("~/Research_office/MyWork/Chi-squared_model/Shiny/power/simulated_OR-RAF_diagram/simulate_mat_n20000_phi05.Rdata")
load("~/Research_office/MyWork/Chi-squared_model/Shiny/power/simulated_OR-RAF_diagram/simulate_diagram_mat.Rdata")

image(t(simulate_mat_n20000_phi05))
image(t(temp.mat))

plot.mat <- simulate_mat_n20000_phi05
plot.mat <- temp.mat
# plot.mat[,1:5] <- temp.mat[,1:5]

plot_ly(x = x.adj(f.vec), y = y.adj.plotly(R.vec), 
        z = plot.mat, zmin = 0, zmax = 1,
        colors = "Greys", type = "contour",
        reversescale = T, hoverinfo = "none",
        # width = 750, height = 700,
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
         margin = list(l = 50, r = 20, b = 80, t = 50, pad = 4))


#### simulate Type I error ####

temp.vec.level <- vector(mode = "numeric", length = length(f.vec))
R <- 1
n <- 20000
phi <- 0.1
for (j in 1:length(f.vec)) {
  print(j)
  f <- f.vec[j]
  p <- c(f*R*phi/(f*R+1-f), (1-f)*phi/(f*R+1-f), f*(1-phi), (1-f)*(1-phi))
  res <- replicate(1000, {
    x <- rmultinom(1, n, p)
    if (sum(x[c(1,3)]) > 2 && sum(x[c(2,4)]) > 2) {
      t.statistics <- t.test(x = c(rep(1, x[1]), rep(0, x[3])), 
                             y = c(rep(1, x[2]), rep(0, x[4])), 
                             var.equal = T)
      t.statistics$p.value < 0.005
    } else {
      FALSE
    }
  })
  temp.vec.level[j] <- mean(res)
}

temp.vec.level5e3 <- temp.vec.level

plot(x = x.adj(f.vec), y = temp.vec.level, type = 'b', log = 'y', xaxt = 'n')
axis(side = 1, at = x.adj(c(0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99)), 
     labels = c(0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99))
abline(v = x.adj(5/n/10), lty = 2)
