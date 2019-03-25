n1 = 100000
n2 = 100000

RV.curves <- rare.variant.curves(n1, n2, rare.variant.threshold = 30,
                                 f.lim = c(1e-4, 1-5e-3), R.lim = c(1,110))

# f.rare.variant.root.finder.left(f.RAF = 1e-4, OR = 1, phi = 0.5, n = 2000, rare.variant.threshold = 0)

par(mar = c(4,3,3,1))
options(scipen=0)
y.plot.max <- 100
plot(x = x.adj(RV.curves$f.vec.left), y = y.adj(RV.curves$OR.vec.left),
     log = 'y', type = 'l', las = 1, xaxt = 'n', yaxt = 'n',
     xlim = x.adj(c(1e-4, 1-5e-3)), ylim = y.adj(c(1, y.plot.max)),
     xaxs = "i", yaxs = "i", xlab = "", ylab = "")
lines(x = x.adj(RV.curves$f.vec.right), y = y.adj(RV.curves$OR.vec.right))
mtext(text = "odds ratio", side = 3, at = x.adj(5e-5),
      las = 1, adj = 0, line = 1, cex = 1.5)
axis(side = 1, at = x.adj(c(0.5, 10^(-1:-4), 1-10^(-1:-2))), cex.lab = 1.5,
     labels = c(0.5, 0.1, expression(10^{-2}), expression(10^{-3}), expression(10^{-4}),
                0.9, 0.99), cex.axis = 1.5)
mtext(text = "risk allele frequency (in control group)", side = 1, at = x.adj(0.93),
      las = 1, line = 2.5, cex = 1.5)
axis(side = 2, at = y.adj(c(1, 2, 5, 10, 20, 50, 100)),
     labels = c(1, 2, 5, 10, 20, 50, 100), las = 1, cex.axis = 1.5)

text(x = x.adj(0.1), y = 20, labels = "Cases/Controls", cex = 5, pos = 1, col = "grey80")
text(x = x.adj(0.1), y = 10, labels = paste0(n1, "/", n2), cex = 5, pos = 1, col = "grey80")

