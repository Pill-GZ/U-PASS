# explicitly source the global.R file
source("./global.R")

#### start of server ####

server <- function(input, output, session) {
#### first tab: OR-RAF diagram #### 
  
  #### quick start guided tours with IntroJS ####
  
  # OR-RAF tab
  observeEvent(input$help_ORRAF, {
    rintrojs::introjs(session, options = list(
      steps = data.frame(element = c("#sample_size", "#false_discovery_control", "#plot_options", 
                                     "#data_options", "#OR-RAF_diagram", "#OR-RAF_diagram", "#gene_info", NA),
                         intro = c("Specify <b>sample sizes</b> here.",
                                   "Specify <b>false discovery control criteria</b> (FWER / Type I error rate) here.",
                                   "Select <b>plot options</b> here.",
                                   "Choose to <b>overlay reported findings</b> from the NHGRI-EBI GWAS Catalog, or upload your own data!",
                                   "Statistical power for association tests is visualized in the <b>OR-RAF diagram</b>. 
                                     Reported findings are also overlaid here.<br><br>
                                     It's <b>fully interactive</b>. Click on a reported locus to display detailed information.
                                     Sample sizes also automatically adapt to the study reporting the selected locus.",
                                   "We recommend specifying the <b>rare-variant regions</b> according to the minimum number of counts 
                                     needed for Fisher's exact test to be correctly calibrated.
                                     In this case, power calculations based on asymptotics are not applicable in the rare-variant region (outisde the red line(s));
                                     single-SNP-based association tests have 0 power, and are <b>not recommended</b>.<br><br>
                                     If a reported locus lies in the <b>low power region</b> (dark regions of the heatmap),
                                     the claim of statistical significance are dubious.",
                                   "When you select a reported locus in the OR-RAF diagram, detailed information is displayed here below the diagram.",
                                   "The power analysis is <b>model-free</b> and <b>test-independent</b>. 
                                     This means that you do not need to specify a disease model, or the test of association used.<br><br>
                                     Find out why in the Details tab."
                         ))
    ))
  }) # end of intro for OR-RAF tab
  
  # OR-RAF tab
  observeEvent(input$help_power_analysis, {
    # Auto fill in Steps 1-4 for the IntroJS demo
    updateSelectInput(session, "step1_target_OR_RAF", selected = "Allele frequency and odds ratio")
    updateSelectInput(session, "step2_fixed_quantity", selected = "Budget / total number of subjects")
    updateSelectInput(session, "step3_type_I_error_criteria", selected = "Family-wise error rate (FWER)")
    updateSelectInput(session, "step4_type_II_error_criteria", selected = "Type II error / non-discovery proportion (NDP)")
    # start IntroJS demo
    rintrojs::introjs(session, options = list(
      steps = data.frame(element = c("#step1", "#step2", "#step3", 
                                     "#step4", "#power_analysis_results", NA),
                         intro = c("Specify in here the <b>RAF</b> and <b>OR</b> of the loci you wish to target.",
                                   "Specify <b>constraint on sample sizes</b> here.<br><br>
                                     It could be total budget (i.e., total number of subjects), number of Cases, or fraction of Cases among recuited subjects.",
                                   "Specify <b>false discovery control criteria</b> (FWER / Type I error rate) here.",
                                   "Specify <b>non-discovery control criteria</b> (Type II error rate / FWNDR) here.",
                                   "Results from the power calculation is displayed here.<br><br>
                                     <ul>
                                       <li>If the contraint is <b>total budget</b>, power is shown as a function of the fraction of Cases.</li>
                                       <li>If the contraint is <b>number of Cases</b>, power is shown as a function of the number of Controls.</li>
                                       <li>If the contraint is <b>fraction of Cases</b>, power is shown as a function of the total number of subjects.</li>
                                     </ul>",
                                   "The power analysis is <b>model-free</b> and <b>test-independent</b>. 
                                     This means that you do not need to specify a disease model, or the test of association used.<br><br>
                                     Find out why in the Details tab."
                         ))
    ))
  }) # end of intro for design-my-study tab
  
  
  
  #### calculate number of cases (n1) and controls (n2), and the fraction of cases (phi) ####
  
  n1 <- reactive({
    if (input$sample_size_specification == 'Number of subjects + fraction of Cases') {
      validate(need({is.integer(input$n); input$n > 0}, "Number of subjects must be a positive integer"))
      floor(input$n * input$phi)
    } else {
      validate(need({is.integer(input$n1); input$n1 > 0}, "Number of cases must be a positive integer"))
      floor(input$n1)
    }
  })
  n2 <- reactive({
    if (input$sample_size_specification == 'Number of subjects + fraction of Cases') {
      validate(need({is.integer(input$n); input$n > 0}, "Number of subjects must be a positive integer"))
      input$n - floor(input$n * input$phi)
    } else {
      validate(need({is.integer(input$n1); input$n1 > 0}, "Number of cases must be a positive integer"))
      floor(input$n2)
    }
  })
  phi <- reactive({
    if (input$sample_size_specification == 'Number of subjects + fraction of Cases') {
      input$phi
    } else {
      validate(
        need({is.integer(input$n1); input$n1 > 0}, "Number of cases must be a positive integer"),
        need({is.integer(input$n2); input$n2 > 0}, "Number of controls must be a positive integer"))
      input$n1 / (input$n1 + input$n2)
    }
  })
  
  #### calculate rare variant region ####
  
  rare_variant_threshold_count <- reactive({
    if (input$rare_variant_zone_specification == "Minimum counts needed to calibrate Fisher's exact test") {
      minimum.RV.Fishers.test(n1 = 2 * n1(), n2 = 2 * n2(), p.val.threhold = p.val.cutoff())
    } else if (input$rare_variant_zone_specification == 'Absolute variant count in study') {
      threshold <- input$rare_variant_threshold_count
      list(left = threshold, right = threshold)
    } else if (input$rare_variant_zone_specification == 'Fraction of total subjects in study') {
      fraction <- as.numeric(gsub("%", "", input$rare_variant_threshold_fraction))/100
      threshold <- floor(fraction * (n1() + n2()))
      list(left = threshold, right = threshold)
    }
  })
  
  rare.variant.curves <- reactive({
    RV.thresholds <- rare_variant_threshold_count()
    rare.variant.curves.calculator(n1 = 2 * n1(), n2 = 2 * n2(), RV.thresholds$left, RV.thresholds$right, f.lim, R.lim)
  })

  #### calculate cutoff threshold for false discovery control ####
  
  p.val.cutoff <- reactive({
    if (input$type_I_error_criteria == 'Type I error') {
      input$alpha
    } else if (input$type_I_error_criteria == 'Family-wise error rate (FWER)') {
      validate(need({is.integer(input$p_FWER); input$p_FWER > 0}, "Effective dimension must be a positive integer"))
      input$alpha_FWER / input$p_FWER
    } 
  })
  
  output$p_val_cutoff <- renderText({ paste("<b>P-value threshold:</b>",p.val.cutoff()) })
  
  cutoff <- reactive({
    qchisq(p = 1 - p.val.cutoff(), df = 1, ncp = 0, lower.tail = T)
  })
  
  #### calculate power for the sample size inputs ####
  
  signal.size.vec <- as.vector(outer(c(1,2,5), 10^(-5:-2)))
  power.vec <- reactive({
    pchisq(q = cutoff(), df = 1, ncp = signal.size.vec * 2 * (n1() + n2()) , lower.tail = F)
  })
  output$power_vec <- renderText({ power.vec() })
  
  #### generate power matrix for heatmaps ####
  
  f.vec <- x.adj.inv(seq(from = x.adj(f.lim)[1], 
                         to = x.adj(f.lim)[2], 
                         length.out = resolution))
  R.vec <- y.adj.plotly.inv(seq(from = y.adj.plotly(R.lim)[1], 
                                to = y.adj.plotly(R.lim)[2], 
                                length.out = resolution))
  signal.mat <- reactive({
    temp.mat <- matrix(nrow = length(R.vec), ncol = length(f.vec))
    for (i in 1:length(R.vec)) {
      for (j in 1:length(f.vec)) {
        temp.mat[i,j] <- signal.size.ana.sol.f(f.vec[j], phi(), R.vec[i]) * 2 * (n1() + n2())
      }
    }
    temp.mat
  })
  
  power.mat <- reactive({
    pchisq(q = cutoff(), df = 1, ncp = signal.mat(), lower.tail = F)
  })
  
  #### OR-RAF plot in base R  (deprecated, but useful reference) ####
  
  output$OR.RAF.plot <- renderPlot({

    #OR <- 1.3^(-50:50)
    solution.vec <- OR.finder(phi = phi(), signal.size = 0.01)
    #hist(rnorm(input$obs), col = 'darkgray', border = 'white')
    par(mar = c(4,3,3,1))
    options(scipen=0)
    y.plot.max <- 100
    plot(x = x.adj(solution.vec$f), y = y.adj(solution.vec$R),
         log = 'y', type = 'n', las = 1, xaxt = 'n', yaxt = 'n',
         xlim = x.adj(c(1e-4, 1-5e-3)), ylim = y.adj(c(1, y.plot.max)),
         xaxs = "i", yaxs = "i", xlab = "", ylab = "")
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
    text(x = x.adj(0.1), y = 10, labels = paste0(n1(), "/", n2()), cex = 5, pos = 1, col = "grey80")

    for (signal.size in signal.size.vec) {
      solution.vec <- OR.finder(phi = phi(), signal.size)
      lines(x = x.adj(solution.vec$f), y = y.adj(solution.vec$R), lty = 2)
      if (signal.size %in% c(5e-4)) {lines(x = x.adj(solution.vec$f), y = y.adj(solution.vec$R), lty = 1, col = 2)}
      if (signal.size >= 1e-4){
        if (signal.size == 1e-1) {
          text(x = x.adj(solution.vec$p[1]*1.03), y = y.plot.max*0.7, pos = 4,
               labels = bquote(paste(lambda/n, " = ", 10^{.(log10(signal.size))})))
        } else if (signal.size %in% 10^(-4:-2)) {
          text(x = x.adj(solution.vec$p[1]*0.9), y = y.plot.max*0.7, pos = 4,
               labels = bquote(10^{.(log10(signal.size))}))
        }
      }
    }
  }, height = 750, width = 700)
  
  #### base OR-RAF plot in plotly ####
  
  OR_RAF_baseplot <- reactive({
    withProgress(message = 'Rendering OR-RAF phase diagram', detail = 'Rendering power heatmap',
                 value = 0, {
                   p <- plot_ly(x = x.adj(f.vec), y = y.adj.plotly(R.vec), 
                                z = power.mat(), zmin = 0, zmax = 1,
                                colors = "Greys", type = "contour",
                                reversescale = T, hoverinfo = "none",
                                width = 750, height = 700,
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
                     add_annotations(yref = "paper", xref = "paper", y = 1.04, x = 1.15, 
                                     text = "power", 
                                     showarrow = F, font = list(size = 25)) %>%
                     add_annotations(yref = "paper", xref = "paper", y = 0.8, x = 0.5, 
                                     text = "<b>Cases/Controls</b>", showarrow = F, 
                                     font=list(size = 40, color = toRGB("grey70"))) %>%
                     add_annotations(yref = "paper", xref = "paper", y = 0.7, x = 0.5, 
                                     text = paste0("<b>", format(n1(), scientific = FALSE), 
                                                   "/", format(n2(), scientific = FALSE), "</b>"), showarrow = F, 
                                     font=list(size = 40, color = toRGB("grey70")))
                   setProgress(value = 1)
                 }) # end of progress bar for base plot
    p
  }) # end of OR-RAF local plot
  
  #### overlay equi-power curves and rare variant zones on OR-RAF diagram ####
  
  OR_RAF_w_power_RV_curves <- reactive({
    withProgress(message = 'Rendering phase diagram with power and RV curves',
                 detail = 'Rendering base OR-RAF diagram', value = 0, {
                   # generate base OR-RAF diagram
                   p <- OR_RAF_baseplot()
                   # overlaying equi-power curves
                   if (input$overlay_equipower_curves == T) {
                     setProgress(value = 0.8, detail = "Overlaying equi-power curves")
                     for (i in 1:length(signal.size.vec)) {
                       signal.size <- signal.size.vec[i]
                       solution.vec <- OR.finder(phi = phi(), signal.size)
                       p <- add_trace(p, x = x.adj(solution.vec$f), y = y.adj.plotly(solution.vec$R), 
                                      line = list(dash = 'dash', width = 1), type = "scatter", mode = "lines", 
                                      color = I("gray50"), hoverinfo = 'text', 
                                      text = paste('signal size / sample =', signal.size, 
                                                   '\n', 'power = ', round(power.vec()[i], digits = 3)),
                                      inherit = F) 
                     }
                   }
                   # overlaying rare-variant zones
                   if (input$overlay_rare_variant_zones == T) {
                     setProgress(value = 0.9, detail = "Overlaying rare-variant region")
                     RV.thresholds <- rare_variant_threshold_count()
                     p <- p %>% 
                       add_trace(x = x.adj(rare.variant.curves()$f.vec.left),
                                 y = y.adj.plotly(rare.variant.curves()$OR.vec.left),
                                 line = list(dash = 'dash', width = 2), type = "scatter", mode = "lines",
                                 color = I("red"), hoverinfo = 'text',
                                 text = paste('rare variant zone to the left',
                                              '\n', 'risk allele count < ', RV.thresholds$left),
                                 inherit = F) %>%
                       add_trace(x = x.adj(rare.variant.curves()$f.vec.right),
                                 y = y.adj.plotly(rare.variant.curves()$OR.vec.right),
                                 line = list(dash = 'dash', width = 2), type = "scatter", mode = "lines",
                                 color = I("red"), hoverinfo = 'text',
                                 text = paste('rare variant zone to the right',
                                              '\n', 'non-risk allele count < ', RV.thresholds$right),
                                 inherit = F)
                   }
                   setProgress(value = 1)
    }) # end of progress bar
    return(p)
  })
  
  #### select and load dataset ####
  # plot redered and passed to output
  
  list_of_datasets <- c("Breast carcinoma" = "./data/breast_cancer.tsv",
                        "Coronary heart disease" = "./data/coronary_heart_disease.tsv", 
                        "Type II diabetes mellitus" = "./data/diabetes.tsv")
  
  dataset <- reactive({
    if (input$overlay_example_dataset == T) {
      if (input$choose_dataset %in% names(list_of_datasets)) {
        setProgress(value = 0.5, detail = 'Loading NHGRI-EBI GWAS Catalog')
        # print(list_of_datasets[input$choose_dataset])
        read_data(list_of_datasets[input$choose_dataset])
      } else {
        
        if (!is.null(input$my_data_upload)) {
          req(input$my_data_upload)
          return(read_data(input$my_data_upload$datapath))
        } else {
          return(NULL)
        }
      } # end of user upload
    } # end of data read
  })
  
  #### determine set of selected loci by capture click event ####
  
  selected_loci_idx <- reactive({
    click_data <- event_data("plotly_click")
    if ( length(click_data) == 0 | is.null(dataset()) ) {
      integer(0)
    } else {
      # there is no index or key values in plotly (although there is one in ggplot)
      # so we have to search for the selected gene in the data
      search_x <- x.adj.inv(click_data$x)
      search_y <- y.adj.plotly.inv(click_data$y)
      search_x_idx <- which(abs(dataset()$RISK.ALLELE.FREQUENCY-search_x) == min(abs(dataset()$RISK.ALLELE.FREQUENCY-search_x), na.rm =T))
      search_y_idx <- which(abs(dataset()$OR-search_y) == min(abs(dataset()$OR-search_y), na.rm =T))
      intersect(search_x_idx, search_y_idx)
    }
  })
  
  #### determine paper of selected loci by capture click event ####
  
  selected_paper_idx <- reactive({
    if (length(selected_loci_idx()) > 0) {
      which(dataset()$PUBMEDID == dataset()[selected_loci_idx()[1], "PUBMEDID"])
    } else (
      integer(0)
    )
  })
  
  #### display details of selected loci by capture click event ####
  
  output$selection <- renderText({
    if (length(selected_loci_idx()) > 0) {
      search_idx <- selected_loci_idx()
      paste("<b>Reported gene(s):</b> ", dataset()[search_idx, "REPORTED.GENE.S."], "<br>",
            "<b>Mapped gene:</b> ", dataset()[search_idx, "MAPPED_GENE"], "<br>",
            "<b>Estimated odds ratio:</b> ", dataset()[search_idx, "OR"], "<br>",
            "<b>Estimated allele frequency:</b> ", dataset()[search_idx, "RISK.ALLELE.FREQUENCY"], "<br>",
            "<b>p-value:</b> ", dataset()[search_idx, "P.VALUE"], "<br>",
            "<b>Date added to catalog:</b> ", dataset()[search_idx, "DATE.ADDED.TO.CATALOG"], "<br>",
            "<b>Study:</b> ", dataset()[search_idx, "STUDY"], "<br>",
            "<b>Journal:</b> ", dataset()[search_idx, "JOURNAL"], "<br>",
            "<b>First author:</b> ", dataset()[search_idx, "FIRST.AUTHOR"], "<br>",
            "<b>Disease trait:</b> ", dataset()[search_idx, "DISEASE.TRAIT"], "<br>",
            "<b>Initial sample size:</b> ", dataset()[search_idx, "INITIAL.SAMPLE.SIZE"], "<br>",
            "<b>Replication sample size:</b> ", dataset()[search_idx, "REPLICATION.SAMPLE.SIZE"], "<br><br>")
    } else if (is.null(dataset())) {
      "Choose a dataset in the NHGRI-EBI GWAS Catalog format."
    } else {
      "Click on data points to display details"
    }
  })
  
  #### update base plot sample sizes and false discovery control according to selected loci ####
  
  observe({
    if (input$adaptive_sample_size == T && length(selected_loci_idx()) > 0) {
      # determine sample sizes
      initial_sample_size <- strsplit(x = as.character(dataset()[selected_loci_idx()[1], "INITIAL.SAMPLE.SIZE"]), 
                                      split = ", ")
      initial_sample_size <- sapply(initial_sample_size, gsub, pattern = ",", replacement = "")
      cases_idx <- grepl(pattern = "cases", x = initial_sample_size)
      controls_idx <- grepl(pattern = "controls", x = initial_sample_size)
      n_cases <- sum(as.numeric(gsub(".*?([0-9]+).*", "\\1", initial_sample_size[cases_idx])))
      n_controls <- sum(as.numeric(gsub(".*?([0-9]+).*", "\\1", initial_sample_size[controls_idx])))
      # make changes to the inputs
      if (n_cases > 0 && n_controls > 0) {
        # Change the number of cases
        updateNumericInput(session, "n1", value = n_cases)
        # Change the number of controls
        updateNumericInput(session, "n2", value = n_controls)
        # Change selection for input$sample_size_specification
        updateSelectInput(session, "sample_size_specification",
                          selected = "Number of Cases + number of Controls")
      } else {
        showModal(modalDialog(
          title = "Unable to automatically determine initial sample sizes",
          paste("No information about number of", ifelse(n_cases > 0, "controls!", "cases!"))
        ))
      }
      # determine p-value cutoff point
      max_p_val <- max(dataset()[selected_paper_idx(), "P.VALUE"], na.rm = T)
      #print(max_p_val)
      if (max_p_val > -Inf) {
        dim_guess <- min(0.1^(ceiling(log10(max_p_val)) + 1), 1e4)
        FWER_guess <- max_p_val * dim_guess
        FWER_guess <- ifelse(FWER_guess <= 0.05, 0.05, 0.1)
        # Chage the FWER in UI
        shinyWidgets::updateSliderTextInput(session, "alpha_FWER", selected = FWER_guess)
        # Change the dimension of FWER control in UI
        updateNumericInput(session, "p_FWER", value = dim_guess)
        # Change false discovery control to FWER 
        updateSelectInput(session, "type_I_error_criteria",
                          selected = "Family-wise error rate (FWER)")
      } else {
        showModal(modalDialog(
          title = "Unable to guess p-value cutoff point"
        ))
      }
    }
  })
  
  #### overlay data points and render plotly ####
  
  output$OR.RAF.plotly <- renderPlotly({
    # overlay data points
    if ( input$overlay_example_dataset == T && !is.null(dataset()) ) {
      withProgress(message = 'Rendering phase diagram with points overlay', 
                   detail = 'Locating selected loci/article',
                   value = 0, {
                     with(data = dataset(), {
                       # separate layers for selected points (if any), points in the same paper (if any),
                       # indexed by binary (T/F) vectors
                       selected_loci_TF <- index %in% selected_loci_idx()
                       selected_paper_TF <- index %in% setdiff(selected_paper_idx(), selected_loci_idx())
                       not_selected_TF <- !selected_paper_TF & !selected_loci_TF
                       
                       # generate base OR-RAF diagram with equi-power and rare-variant curves
                       setProgress(value = 0.2, detail = "Rendering base plot")
                       p <- OR_RAF_w_power_RV_curves() 
                       
                       # Overlaying loci not selected / not in article
                       setProgress(value = 0.7, detail = "Overlaying loci not selected / not in article")
                       p <- p %>% add_markers(x = x.adj(RISK.ALLELE.FREQUENCY[not_selected_TF]), 
                                              y = y.adj.plotly(OR[not_selected_TF]), 
                                              marker = list(color = 'rgb(255, 255, 255)', size = 10, opacity = 0.5,
                                                            line = list(color = 'rgb(20, 100, 238)', width = 3)),
                                              hoverinfo = "text",
                                              text = paste("RAF: ", RISK.ALLELE.FREQUENCY[not_selected_TF],
                                                           "OR: ", round(OR[not_selected_TF], digits = 3),
                                                           '<br>MAPPED GENE:', MAPPED_GENE[not_selected_TF],
                                                           "<br>PUBMEDID: ", PUBMEDID[not_selected_TF]),
                                              type = "scatter", mode = 'markers', inherit = F) 
                       
                       # Overlaying loci in article
                       if (sum(selected_paper_TF) > 0) {
                         setProgress(value = 0.9, detail = "Overlaying loci selected / in article")
                         p <- p %>% add_markers(x = x.adj(RISK.ALLELE.FREQUENCY[selected_paper_TF]),
                                                y = y.adj.plotly(OR[selected_paper_TF]),
                                                marker = list(color = 'rgb(255, 255, 255)', size = 10, opacity = 0.6,
                                                              line = list(color = 'rgb(255, 165, 0)', width = 3)),
                                                hoverinfo = "text",
                                                text = paste("RAF: ", RISK.ALLELE.FREQUENCY[selected_paper_TF],
                                                             "OR: ", round(OR[selected_paper_TF], digits = 3),
                                                             '<br>MAPPED GENE:', MAPPED_GENE[selected_paper_TF],
                                                             "<br>PUBMEDID: ", PUBMEDID[selected_paper_TF]),
                                                type = "scatter", mode = 'markers', inherit = F) 
                       }
                       # Overlaying loci selected 
                       if (sum(selected_loci_TF) > 0) {
                         p <- p %>% add_markers(x = x.adj(RISK.ALLELE.FREQUENCY[selected_loci_TF]),
                                                y = y.adj.plotly(OR[selected_loci_TF]),
                                                marker = list(color = 'rgb(255, 165, 0)', size = 15, opacity = 0.9,
                                                              line = list(color = 'rgb(255, 0, 0)', width = 3)),
                                                hoverinfo = "text",
                                                text = paste("RAF: ", RISK.ALLELE.FREQUENCY[selected_loci_TF],
                                                             "OR: ", round(OR[selected_loci_TF], digits = 3),
                                                             '<br>MAPPED GENE:', MAPPED_GENE[selected_loci_TF],
                                                             "<br>PUBMEDID: ", PUBMEDID[selected_loci_TF]),
                                                type = "scatter", mode = 'markers', inherit = F)
                       }
                       
                       # suppress warnings  
                       # storeWarn<- getOption("warn")
                       # options(warn = -1) 
                       p <- p %>% layout( source = "dataset" )
                       # restore warnings, delayed so plot is completed
                       # shinyjs::delay(expr =({ 
                       #   options(warn = storeWarn) 
                       # }), ms = 100)
                       
                       # set progress bar to 1
                       setProgress(value = 1)
                       
                       p
                     }) # end of data overlay
                   }) # end of progress bar
    } else { 
      # don't overlay data points
      # suppress warnings  
      storeWarn<- getOption("warn")
      options(warn = -1)
      p <- OR_RAF_w_power_RV_curves()
      # restore warnings, delayed so plot is completed
      shinyjs::delay(expr =({
        options(warn = storeWarn)
      }), ms = 100)
      p
    }
  }) # end of OR-RAF plotly output
  
  
  
#### second tab: design my study ###################### #### 
  #### checking if design is complete ####
  
  waiting_for_design <- reactive({
    if (input$step1_target_OR_RAF != "Select a target" &&
        input$step2_fixed_quantity != "Select a constraint" &&
        input$step3_type_I_error_criteria != "Select a criterion" &&
        input$step4_type_II_error_criteria != "Select a criterion") {
      FALSE
    } else {
      TRUE
    }
  })
  output$waiting_for_design <- reactive({waiting_for_design()})
  outputOptions(output, 'waiting_for_design', suspendWhenHidden=FALSE)
  
  #### determine sample size specification ####
  fixed_n_design <- reactive({
    (!waiting_for_design()) && (input$step2_fixed_quantity == "Budget / total number of subjects")
  })
  output$fixed_n_design <- reactive({fixed_n_design()})
  outputOptions(output, 'fixed_n_design', suspendWhenHidden=FALSE)
  
  fixed_n1_design <- reactive({
    (!waiting_for_design()) && (input$step2_fixed_quantity == "Number of Cases")
  })
  output$fixed_n1_design <- reactive({fixed_n1_design()})
  outputOptions(output, 'fixed_n1_design', suspendWhenHidden=FALSE)
  
  fixed_phi_design <- reactive({
    (!waiting_for_design()) && (input$step2_fixed_quantity == "Fraction of Cases")
  })
  output$fixed_phi_design <- reactive({fixed_phi_design()})
  outputOptions(output, 'fixed_phi_design', suspendWhenHidden=FALSE)
  
  #### calculate rejection region for the design ####
  
  #### calculate cutoff threshold for false discovery control ####
  
  design.p.val.cutoff <- reactive({
    if (input$step3_type_I_error_criteria == 'Type I error') {
      input$design_alpha
    } else if (input$type_I_error_criteria == 'Family-wise error rate (FWER)') {
      validate(need({is.integer(input$p_FWER); input$p_FWER > 0}, "Effective dimension must be a positive integer"))
      input$design_alpha_FWER / input$design_p_FWER
    } 
  })
  
  output$design_p_val_cutoff <- renderText({ paste("<b>P-value threshold:</b>", design.p.val.cutoff()) })
  
  design.cutoff <- reactive({
    qchisq(p = 1 - design.p.val.cutoff(), df = 1, ncp = 0, lower.tail = T)
  })
  
  # design.cutoff <- reactive({
  #   if (input$step3_type_I_error_criteria == 'Type I error') {
  #     qchisq(p = 1 - input$design_alpha, df = 1, ncp = 0, lower.tail = T)
  #   } else if (input$step3_type_I_error_criteria == 'Family-wise error rate (FWER)') {
  #     validate(need({is.integer(input$design_p_FWER); input$design_p_FWER > 0}, "Effective dimension must be a positive integer"))
  #     qchisq(p = 1 - input$design_alpha_FWER / input$design_p_FWER, df = 1, ncp = 0, lower.tail = T)
  #   } 
  # })
  
  #### calculate target power of the study ####
  
  design.power <- reactive({
    if (input$step4_type_II_error_criteria == 'Type II error / non-discovery proportion (NDP)') {
      1 - input$design_power
    } else if (input$step4_type_II_error_criteria == 'Family-wise non-discovery rate (FWNDR)') {
      validate(need({is.integer(input$design_sparsity); input$design_sparsity > 0}, "Sparsity must be a positive integer"))
      1 - input$design_FWNDR / input$design_sparsity
    } 
  })
  
  #### fixed budget: power as a function of fraction variable.phi ####
  
  # input$fixed_budget
  output$optimal_design_fixed_n <- renderPlotly({
    if (fixed_n_design()) {
      validate(need({is.integer(input$fixed_budget); input$fixed_budget > 0}, "Number of subjects must be a positive integer"))
      fixed.n <- input$fixed_budget
      variable.phi <- 1:99/100
      # calculate signal size of the design
      design.signal.size.per.sample <- {
        if (input$step1_target_OR_RAF == 'Allele frequency and odds ratio') {
          validate(
            need({is.numeric(input$target_RAF) & input$target_RAF > 0 & input$target_RAF < 1}, 
                 "Target risk allele frequency must be between 0 and 1"),
            need({is.numeric(input$target_OR); input$target_OR > 1}, 
                 "Target odds ratio must be greater than 1")
          )
          signal.size.ana.sol.f(f = input$target_RAF, phi = variable.phi, R = input$target_OR)
        } else if (input$step1_target_OR_RAF == 'Signal size per sample (advanced user)') {
          input$target_w2
        }
      }
      # calculate power as a function of Controls (assuming one allele pair per subject)
      variable.power <- pchisq(q = design.cutoff(), df = 1 , lower.tail = F, 
                               ncp = design.signal.size.per.sample * 2 * fixed.n)
      # calculate optimal fraction of Cases
      optimal.fraction.of.cases <- variable.phi[which.max(design.signal.size.per.sample)]
      optimal.fraction.of.cases.prompt <- ifelse(max(variable.power) > design.power(), 
                                                 yes = paste("<b>Optimal Cases / Controls:</b>", 
                                                             ceiling(fixed.n * optimal.fraction.of.cases), "/",
                                                             fixed.n  - ceiling(fixed.n * optimal.fraction.of.cases)),
                                                 no = "<b>cannot be achieved</b>")
      
      plot_ly(x = variable.phi, y = variable.power, 
              type = 'scatter', mode = 'lines+markers',
              line = list(width = 4), marker = list(size = 8), 
              hoverinfo = 'text', 
              text = paste('Fraction of Cases =', format(variable.phi, scientific = F), 
                           '\n', 'Power =', round(variable.power, digits = 3)),
              width = 700, height = 700) %>%
        layout(xaxis = list(fixedrange=FALSE, range = c(-0.01, 1.01), 
                            tickvals = 0:10/10, 
                            ticktext = 0:10/10,
                            zeroline = FALSE, tickfont = list(size = 20),
                            showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
               yaxis = list(range = c(-0.01,1.01), fixedrange=TRUE,
                            tickvals = 0:10/10, ticktext = 0:10/10, 
                            zeroline = FALSE, tickfont = list(size = 20),
                            showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
               shapes = list(list(type = "line", line = list(color = toRGB("grey70"), dash = "dash"),
                                  xref = "x", yref = "y", 
                                  x0 = optimal.fraction.of.cases, x1 = optimal.fraction.of.cases, 
                                  y0 = 0, y1 = 1),
                             # list(type = "line", line = list(color = toRGB("grey20"), dash = "dash"),
                             #      xref = "x", yref = "y", 
                             #      x0 = 1/(1+sqrt(input$target_OR)), x1 = 1/(1+sqrt(input$target_OR)), 
                             #      y0 = 0, y1 = 1),
                             list(type = "line", line = list(color = toRGB("grey70"), dash = "dash"),
                                  xref = "x", yref = "y", 
                                  x0 = 0, x1 = 1, 
                                  y0 = design.power(), y1 = design.power())),
               showlegend = FALSE,
               margin = list(l = 50, r = 20, b = 80, t = 50, pad = 4)) %>%
        add_annotations(yref = "paper", xref = "paper", y = 1.08, x = -0.07, 
                        text = "power", showarrow = F, font=list(size = 25)) %>%
        add_annotations(yref = "paper", xref = "paper", y = -0.12, x = 1, 
                        text = "fraction of subjects in Case group", 
                        showarrow = F, font = list(size = 25))  %>%
        add_annotations(yref = "paper", xref = "paper", y = design.power() - 0.01, x = 0.05, 
                        text = paste("<b>Target power:</b>", format(design.power(), scientific = F)),
                        showarrow = F, font=list(size = 20, color = toRGB("grey60"))) %>%
        add_annotations(yref = "paper", xref = "paper", y = design.power() - 0.05, x = 0.05, 
                        text = optimal.fraction.of.cases.prompt,
                        showarrow = F, font=list(size = 20, color = toRGB("grey60")))
    }
  }) # end of fixed Cases plotly output
  
  
  #### fixed Cases: power as a function of number of Controls variable.n2 ####
  
  # input$fixed_cases
  output$optimal_design_fixed_n1 <- renderPlotly({
    if (fixed_n1_design()) {
      validate(need({is.integer(input$fixed_cases); input$fixed_cases > 0}, "Number of cases must be a positive integer"))
      variable.n2 <- as.vector(outer(c(1,1.5,2:9), 10^(1:6)))
      variable.n <- input$fixed_cases + variable.n2
      variable.phi <- input$fixed_cases / variable.n
      # calculate signal size of the design
      design.signal.size.per.sample <- {
        if (input$step1_target_OR_RAF == 'Allele frequency and odds ratio') {
          validate(
            need({is.numeric(input$target_RAF) & input$target_RAF > 0 & input$target_RAF < 1}, 
                 "Target risk allele frequency must be between 0 and 1"),
            need({is.numeric(input$target_OR); input$target_OR > 1}, 
                 "Target odds ratio must be greater than 1")
          )
          signal.size.ana.sol.f(f = input$target_RAF, phi = variable.phi, R = input$target_OR)
        } else if (input$step1_target_OR_RAF == 'Signal size per sample (advanced user)') {
          input$target_w2
        }
      }
      # calculate power as a function of Controls (assuming one allele pair per subject)
      variable.power <- pchisq(q = design.cutoff(), df = 1 , lower.tail = F, 
                               ncp = design.signal.size.per.sample * 2 * variable.n)
      # calculate required number of controls
      required.number.of.controls <- determine.intersection(variable.n2, variable.power, design.power())
      required.number.of.controls.prompt <- ifelse(required.number.of.controls, 
                                                   yes = paste("<b>Number of Controls needed:</b>", 
                                                               format(required.number.of.controls, scientific = F)),
                                                   no = "<b>cannot be achieved</b>")
      
      plot_ly(x = variable.n2, y = variable.power, 
              type = 'scatter', mode = 'lines+markers',
              line = list(width = 4), marker = list(size = 8), 
              hoverinfo = 'text', 
              text = paste('Number of Controls =', format(variable.n2, scientific = F), 
                           '\n', 'Power =', round(variable.power, digits = 3)),
              width = 700, height = 700) %>%
        layout(xaxis = list(fixedrange=FALSE,  type = "log", #range = log(range(variable.n2)), 
                            tickvals = c(as.vector(outer(c(1:9), 10^(1:6))), 1e7), 
                            ticktext = c(as.vector(rbind(sapply(10^(1:6), format, scientific = F), matrix("", 8, 6))), ""),
                            tickangle = 0,
                            zeroline = FALSE, tickfont = list(size = 20),
                            showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
               yaxis = list(range = c(-0.01,1.01), fixedrange=TRUE,
                            tickvals = 0:10/10, ticktext = 0:10/10, 
                            zeroline = FALSE, tickfont = list(size = 20),
                            showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
               shapes = list(list(type = "line", line = list(color = toRGB("grey70"), dash = "dash"),
                                  xref = "x", yref = "y", 
                                  x0 = min(variable.n2), x1 = 1e7, 
                                  y0 = design.power(), y1 = design.power())),
               showlegend = FALSE,
               margin = list(l = 50, r = 20, b = 80, t = 50, pad = 4)) %>%
        add_annotations(yref = "paper", xref = "paper", y = 1.08, x = -0.07, 
                        text = "power", showarrow = F, font=list(size = 25)) %>%
        add_annotations(yref = "paper", xref = "paper", y = -0.12, x = 1, 
                        text = "number of subjects in Control group", 
                        showarrow = F, font = list(size = 25))  %>%
        add_annotations(yref = "paper", xref = "paper", y = design.power() - 0.01, x = 0.05, 
                        text = paste("<b>Target power:</b>", format(design.power(), scientific = F)),
                        showarrow = F, font=list(size = 20, color = toRGB("grey60"))) %>%
        add_annotations(yref = "paper", xref = "paper", y = design.power() - 0.05, x = 0.05, 
                        text = required.number.of.controls.prompt,
                        showarrow = F, font=list(size = 20, color = toRGB("grey60")))
    }
  }) # end of fixed Cases plotly output
  
  #### fixed fraction of cases: power as a function of total subjects variable.n ####
  
  # input$fixed_phi
  output$optimal_design_fixed_phi <- renderPlotly({
    if (fixed_phi_design()) {
      # calculate signal size of the design
      design.signal.size.per.sample <- {
        if (input$step1_target_OR_RAF == 'Allele frequency and odds ratio') {
          validate(
            need({is.numeric(input$target_RAF) & input$target_RAF > 0 & input$target_RAF < 1}, 
                 "Target risk allele frequency must be between 0 and 1"),
            need({is.numeric(input$target_OR); input$target_OR > 1}, 
                 "Target odds ratio must be greater than 1")
          )
          signal.size.ana.sol.f(f = input$target_RAF, phi = input$fixed_phi, R = input$target_OR)
        } else if (input$step1_target_OR_RAF == 'Signal size per sample (advanced user)') {
          input$target_w2
        }
      }
      variable.n <- as.vector(outer(c(1,1.5,2:9), 10^(1:6)))
      # calculate power as a function of total sample size (assuming one allele pair per subject)
      variable.power <- pchisq(q = design.cutoff(), df = 1 , lower.tail = F, 
                               ncp = design.signal.size.per.sample * 2 * variable.n)
      
      # calculate required number of total samples
      required.number.of.total.samples <- determine.intersection(variable.n, variable.power, design.power())
      required.number.of.total.samples.prompt <- ifelse(required.number.of.total.samples, 
                                                   yes = paste("<b>Number of samples needed:</b>", 
                                                               format(required.number.of.total.samples, scientific = F),
                                                               "<br><b>Cases:</b>", 
                                                               ceiling(required.number.of.total.samples * input$fixed_phi), 
                                                               " <b>Controls:</b>",
                                                               required.number.of.total.samples - 
                                                                 ceiling(required.number.of.total.samples * input$fixed_phi)),
                                                   no = "<b>cannot be achieved</b>")
      
      plot_ly(x = variable.n, y = variable.power, 
              type = 'scatter', mode = 'lines+markers',
              line = list(width = 4), marker = list(size = 8), 
              hoverinfo = 'text', 
              text = paste('Number of samples =', format(variable.n, scientific = F), 
                           '\n', 'Power =', round(variable.power, digits = 3)),
              width = 700, height = 700) %>%
        layout(xaxis = list(fixedrange=FALSE,  type = "log", #range = log(range(variable.n)), 
                            tickvals = c(as.vector(outer(c(1:9), 10^(1:6))), 1e7), 
                            ticktext = c(as.vector(rbind(sapply(10^(1:6), format, scientific = F), matrix("", 8, 6))), ""),
                            tickangle = 0,
                            zeroline = FALSE, tickfont = list(size = 20),
                            showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
               yaxis = list(range = c(-0.01,1.01), fixedrange=TRUE,
                            tickvals = 0:10/10, ticktext = 0:10/10, 
                            zeroline = FALSE, tickfont = list(size = 20),
                            showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
               shapes = list(list(type = "line", line = list(color = toRGB("grey70"), dash = "dash"),
                                  xref = "x", yref = "y", 
                                  x0 = min(variable.n), x1 = 1e7, 
                                  y0 = design.power(), y1 = design.power())),
               showlegend = FALSE,
               margin = list(l = 50, r = 20, b = 80, t = 50, pad = 4)) %>%
        add_annotations(yref = "paper", xref = "paper", y = 1.08, x = -0.07, 
                        text = "power", showarrow = F, font=list(size = 25)) %>%
        add_annotations(yref = "paper", xref = "paper", y = -0.12, x = 1, 
                        text = "number of subjects total", 
                        showarrow = F, font = list(size = 25))  %>%
        add_annotations(yref = "paper", xref = "paper", y = design.power() - 0.01, x = 0.05, 
                        text = paste("<b>Target power:</b>", format(design.power(), scientific = F)),
                        showarrow = F, font=list(size = 20, color = toRGB("grey60"))) %>%
        add_annotations(yref = "paper", xref = "paper", y = design.power() - 0.05, x = 0.05, 
                        text = required.number.of.total.samples.prompt,
                        showarrow = F, font=list(size = 20, color = toRGB("grey60")))
    }
  }) # end of fixed phi plotly output
  
  
#### end of server ####
} # end of server