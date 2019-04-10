
# check and install all required packages
required.packages <- c("shiny", "MASS", "plotly", "dplyr", "shinyWidgets", "knitr", "rintrojs")
missing.packages <- required.packages[!required.packages %in% rownames(installed.packages())]

for (missing.package in missing.packages) {
  install.packages(missing.package)
}

# we strongly recommend using our version of the load spinner
if (!"shinycssloaders" %in% rownames(installed.packages())) {
  NEEDS_LOAD_SPINNER <- TRUE
} else if (packageVersion("shinycssloaders") != package_version('0.2.1')) {
  NEEDS_LOAD_SPINNER <- TRUE
} else {
  NEEDS_LOAD_SPINNER <- FALSE
}

if (NEEDS_LOAD_SPINNER) {
  if (!"devtools" %in% rownames(installed.packages())) {
    install.packages("devtools")
  }
  devtools::install_github("Pill-GZ/shinycssloaders", upgrade = "never")
}
