U-PASS
======

U-PASS: a **U**nified **P**ower Analysis of **As**sociation **S**tudies.

GWAS power calculator / generic power calculations for large scale association tests.

U-PASS runs as an <a href="https://shiny.rstudio.com/" target="_blank">R Shiny application</a>. It is a free, open source software under the <a href="https://badges.mit-license.org/" target="_blank">MIT license</a>.

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

-   [Installation Guide](#installation-guide)
    -   [Download](#download)
    -   [Install dependencies](#install-dependencies)
    -   [Start / terminate the application](#start-terminate-the-application)
-   [User Guide](#user-guide)
    -   [The OR-RAF diagram](#the-or-raf-diagram)
    -   [Interactively explore reported findings in the NHGRI-EBI Catalog](#interactively-explore-reported-findings-in-the-nhgri-ebi-catalog)
    -   [Review and forensics of reported findings](#review-and-forensics-of-reported-findings)
    -   [Find optimal study designs](#find-optimal-study-designs)
-   [References](#references)

Installation Guide
==================

The application is hosted on <span target="_blank"><https://power.stat.lsa.umich.edu/u-pass/></span>.

You can run the application locally by following this three-steps installation guide.

Download
--------

You can download the project by running in your terminal:

``` shell
git clone https://github.com/Pill-GZ/U-PASS.git
```

or download as a zip file from its <a href="https://github.com/Pill-GZ/U-PASS" target="_blank">GitHub page</a>: <span target="_blank"><https://github.com/Pill-GZ/U-PASS></span>.

Install dependencies
--------------------

We have collected the required R packages inside `install_required_packages.R`.

You can install these packages by navigating to the project folder, and running in your terminal:

``` shell
Rscript install_required_packages.R 
```

or from inside R (RStudio):

``` r
source("install_required_packages.R")
```

Start / terminate the application
---------------------------------

You can start the application by running in your terminal:

``` shell
Rscript -e 'library(methods); shiny::runApp("./", launch.browser=TRUE)'
```

or from inside R (RStudio):

``` r
shiny::runApp()
```

The application can be terminated by simply closing the browser (or browser tab).

Alternatively, the application can be terminated by pressing `Ctrl + C` in the terminal, or by pressing the red stop button inside Rstudio.

User Guide
==========

The OR-RAF diagram
------------------

U-PASS calculates statistical power based the core parameters common to models of qualitative traits:

-   Sample sizes, i.e., the number of Cases, *n*<sub>1</sub>, and Controls, *n*<sub>2</sub>,
-   Conditional distribution of risk variant among Controls, i.e., risk allele frequency (RAF) in the Control group.
-   Odds ratio (OR) of having the defined trait among the genetic variants.

Users need only prescribe the sample sizes, by one of two ways provided in the first box, i.e., total sample size + fraction of Cases, or number of Cases + number of Controls.

Statistical power of familiar association tests, including ncluding the likelihood ratio test, chi-square test, Welch's t-test, and LR test for logistic regressions, have the same asymptotic power curves (see the <a href="https://power.stat.lsa.umich.edu/u-pass/U-PASS_documentation.html#a_test-independent_power_analysis" target="_blank">documentation</a> for details). This common power limit is calculated as a function of RAF and OR, and visualized as a heatmap in the OR-RAF diagram.

<center>
<img src="www/user_guide_files/basics.gif" style="width:90.0%" />
</center>
<br />

Interactively explore reported findings in the NHGRI-EBI Catalog
----------------------------------------------------------------

We provide options for users to load and overlay findings reported in the <a href="https://www.ebi.ac.uk/gwas/home" target="_blank">NHGRI-EBI GWAS Catalog</a>, or upload data from other sources compliant with the Catalog's data format.

A quick reference for the diagram with data overlay:

-   Circles: reported associations
    -   red: user selected loci
    -   orange: findings reported in the same study as the user selected loci
    -   blue: findings reported in studies other than the one selected
-   Greyscale heatmap: OR-RAF power diagram of association tests
-   red dashed lines: rare-variant threshold. We recommend specifying the threshold by the <a href="https://power.stat.lsa.umich.edu/u-pass/U-PASS_documentation.html#finite-sample_corrections" target="_blank">minimum calibration numbers</a>.
    -   left (if present): the minimum risk variant count needed for the asymptotic approximations to apply.
    -   right (if present): the minimum non-risk variant count needed for the asymptotic approximations to apply.

The initial sample sizes are dynamically adjusted, and automatically determined from texts of the article reporting the user selected loci.

Information of the selected loci and the study is also dynamically displayed below the diagram.

<center>
<img src="www/user_guide_files/interactive.gif" style="width:60.0%" />
</center>
<br />

Review and forensics of reported findings
-----------------------------------------

The unified power analysis allows us to examine results from different studies employing different models and applicable tests, in the same diagram, with the same power limits. It allows for a systematic review of reported findings for their statistical validity.

In particular, a reported association predicted to have low power given the study's sample size -- lying in the dark regions of the OR-RAF diagram -- while not impossible, invites further scrutiny. It should be noted that a reported association predicted to have high power is not automatically accurate, as survival bias induced by multiple testing may inflate the reported OR and RAF estimates.

Studies where reported associations show misalignment with the predicted powered curves may be further investigated for potential problems in the data curation process. The following figure shows one such study, where gross misalignment was identified.

<center>
<img src="www/user_guide_files/forensics.png" alt="optional caption text" style="width:60.0%" />
</center>
We reached out to the authors of the study (Domínguez-Cruz et al. 2018), who confirmed that the RAF reported in the Catalog were based on all subjects in the study, as opposed to only the control group, while the Catalog requires that <a href="https://www.ebi.ac.uk/gwas/docs/fileheaders" target="_blank">RAF be reported in the control group only</a>.As a consequence, the RAFs are systematically overestimated, shifting the reported findings to the right in the diagram.

In general, we expect this aspect of our software to be useful for discovering problems with data entry and catalog curation process, as well as for assessing the reproducibility and robustness of reported findings.

Find optimal study designs
--------------------------

We provide three ways to perform power analysis, depending on the contraint of the study design.

-   If the contraint is the total budget, i.e., total number of subjects recruited,
    -   power is calculated as a function of the fraction of Cases.
-   If the contraint is the number of Cases,
    -   power is calculated as a function of the number of Controls.
-   If the contraint is the fraction of Cases,
    -   power is calculated as a function of the total number of total subjects.

<center>
<img src="www/user_guide_files/design.gif" alt="optional caption text" style="width:100.0%" />
</center>
<br />

The power analysis tool uses the targeted RAF and OR directly to calculated optimal study designs. This allows us to bypass the disease models which define these two quantities implicitly. Indeed, when designing a study to replicate a reported finding, the core quantities RAF and OR are alaways available in GWAS catalogs, while disease models are often not reported in the literature. See <a href="https://power.stat.lsa.umich.edu/u-pass/U-PASS_documentation.html#a_model-invariant_parametrization" target="_blank">more arguments in the documentation</a> on why we prefer to specify these two quantities directly.

Type I error control criteria may be specified in terms of <a href="https://en.wikipedia.org/wiki/Family-wise_error_rate" target="_blank">family-wise error rate</a> or <a href="https://en.wikipedia.org/wiki/Type_I_and_type_II_errors" target="_blank">type I error</a>.

Target non-discovery rate may be specified in terms of <a href="https://en.wikipedia.org/wiki/Type_I_and_type_II_errors" target="_blank">power / type II error</a>, or the more stringent family-wise non-discovery rate, i.e., the probablity of not detecting any one of the loci with equal or stronger signal.

References
==========

Domínguez-Cruz, Miriam Givisay, María de Lourdes Muñoz, Armando Totomoch-Serra, María Guadalupe García-Escalante, Juan Burgueño, Nina Valadez-González, Doris Pinto-Escalantes, and Álvaro Díaz-Badillo. 2018. “Pilot Genome-Wide Association Study Identifying Novel Risk Loci for Type 2 Diabetes in a Maya Population.” *Gene* 677. Elsevier: 324–31.
