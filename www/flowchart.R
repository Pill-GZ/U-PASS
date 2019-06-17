# setEPS()
# postscript("www/flowchart.eps", width = 8, height = 3)
# pdf("www/flowchart.pdf", width = 8, height = 3)
png("www/flowchart.png", width = 8, height = 3, units = "in", res = 100)
# par(mar = c(2,2,1,1))
par(mar = c(0, 0, 0, 0))
plot(c(0,1), c(0,1), type = 'n', axes = F)
first_col <- 0.2
second_col <- 0.6
third_col <- 0.9
grey_background <- grey(0.95)

text(x = 0.425, y = 0.65, labels = bquote(. %=>% .), cex = 1.5)
text(x = 0.425, y = 0.18, labels = bquote(. %<=>% .), cex = 1.5)
text(x = 0.775, y = 0.65, labels = bquote(. %=>% .), cex = 1.5)
text(x = 0.775, y = 0.18, labels = bquote(. %=>% .), cex = 1.5)

## Disease model parameters
polygon(x = c(0, 0.4, 0.4, 0), y = c(0.38, 0.38, 1, 1), col = grey_background)
text(x = first_col, y = 0.90, "Disease model\n(dominant, additive etc.)")
text(x = first_col, y = 0.75, "Risk allele frequency in pop.")
text(x = first_col, y = 0.60, "Disease prevelance in pop.")
text(x = first_col, y = 0.45, "Genotype relative risk (GRR)")

## Core parameters 1
polygon(x = c(0.45, 0.75, 0.75, 0.45), y = c(0.5, 0.5, 0.82, 0.82), col = grey_background)
text(x = second_col, y = 0.72, "Risk allele frequency\nin control group (f)")
text(x = second_col, y = 0.55, bquote(paste("Odds ratio (", R, ")")))

## Sample sizes
polygon(x = c(0, 0.4, 0.4, 0), y = c(0.05, 0.05, 0.3, 0.3), col = grey_background)
text(x = first_col, y = 0.23, bquote(paste("Number of cases (", n[1], ")")))
text(x = first_col, y = 0.12, bquote(paste("Number of controls (", n[2], ")")))

## Core parameters 2
polygon(x = c(0.45, 0.75, 0.75, 0.45), y = c(0.05, 0.05, 0.3, 0.3), col = grey_background)
text(x = second_col, y = 0.23, bquote(paste("Fraction of cases (", phi, ")")))
text(x = second_col, y = 0.12, bquote(paste("Number of subjects (", n, ")")))

## Power analysis
polygon(x = c(0.80, 1, 1, 0.80), y = c(0.05, 0.05, 0.82, 0.82), col = grey_background)
text(x = third_col, y = 0.425, "Statistical\npower")

dev.off()

# if (!require("DiagrammeR")) {
#   install.packages("DiagrammeR")
#   library("DiagrammeR")
# }
# 
# DiagrammeR::grViz("www/flowchart.gv")
