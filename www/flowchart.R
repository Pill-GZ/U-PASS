# setEPS()
# postscript("www/flowchart.eps", width = 8, height = 3)
# pdf("www/flowchart.pdf", width = 8, height = 3)
png("www/flowchart.png", width = 8, height = 2.5, units = "in", res = 100)
# par(mar = c(2,2,1,1))
par(mar = c(0.1, 0, 0.1, 0))
plot(c(0, 1.15), c(0, 1.02), type = 'n', axes = F, yaxs="i")
first_col <- 0.2
second_col <- 0.6
third_col <- 0.9
grey_background <- grey(0.95)
grey_light <- grey(0.7)
grey_medium <- grey(0.4)
grey_dark <- grey(0.2)

text(x = 0.425, y = 0.60, labels = bquote(. %=>% .), cex = 1.5, col = grey_dark)
text(x = 0.425, y = 0.13, labels = bquote(. %<=>% .), cex = 1.5, col = grey_light)
text(x = 0.775, y = 0.60, labels = bquote(. %=>% .), cex = 1.5, col = grey_light)
text(x = 0.775, y = 0.13, labels = bquote(. %=>% .), cex = 1.5, col = grey_light)

## Disease model parameters
polygon(x = c(0, 0.40, 0.40, 0), y = c(0.33, 0.33, 0.95, 0.95), 
        col = "white", border = grey_light, lwd = 2)
text(x = first_col, y = 0.85, "Disease model\n(dominant, additive etc.)", col = grey_dark)
text(x = first_col, y = 0.70, "Risk allele frequency in pop.", col = grey_dark)
text(x = first_col, y = 0.55, "Disease prevelance in pop.", col = grey_dark)
text(x = first_col, y = 0.40, "Genotype relative risk (GRR)", col = grey_dark)

## Core parameters 1
polygon(x = c(0.45, 0.75, 0.75, 0.45), y = c(0.45, 0.45, 0.77, 0.77), 
        col = "white", border = grey_light, lwd = 2)
text(x = second_col, y = 0.67, "Risk allele frequency\nin control group (f)", col = grey_dark)
text(x = second_col, y = 0.50, bquote(paste("Odds ratio (", R, ")")), col = grey_dark)

## Sample sizes
polygon(x = c(0, 0.4, 0.4, 0), y = c(0.00, 0.00, 0.25, 0.25), col = grey_background, border = F)
text(x = first_col, y = 0.18, bquote(paste("Number of cases (", n[1], ")")), col = grey_medium)
text(x = first_col, y = 0.07, bquote(paste("Number of controls (", n[2], ")")), col = grey_medium)

## Core parameters 2
polygon(x = c(0.45, 0.75, 0.75, 0.45), y = c(0.00, 0.00, 0.25, 0.25), col = grey_background, border = F)
text(x = second_col, y = 0.18, bquote(paste("Fraction of cases (", phi, ")")), col = grey_medium)
text(x = second_col, y = 0.07, bquote(paste("Number of subjects (", n, ")")), col = grey_medium)

## Power analysis
polygon(x = c(0.80, 1, 1, 0.80), y = c(0.00, 0.00, 0.77, 0.77), col = grey_background, border = F)
text(x = third_col, y = 0.375, "Statistical\npower", col = grey_medium)

dev.off()

# if (!require("DiagrammeR")) {
#   install.packages("DiagrammeR")
#   library("DiagrammeR")
# }
# 
# DiagrammeR::grViz("www/flowchart.gv")
