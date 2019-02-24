#### read datasets downloaded from EBI ####

# function to read and pre-process data
read_data <- function(dataset_directory_and_name) {
  dataset <- read.delim(dataset_directory_and_name, 
                        na.strings = c("NA", "NR"), stringsAsFactors = F)
  # filter out non-numerical stuff in RAF
  dataset$RISK.ALLELE.FREQUENCY <- as.numeric(gsub("\\(.*\\)", "", dataset$RISK.ALLELE.FREQUENCY))
  # separate beta's from OR's
  beta.indicator <- grepl(pattern = "increase|decrease", x = dataset$X95..CI..TEXT.)
  dataset$OR <- ifelse(beta.indicator, exp(dataset$OR.or.BETA), dataset$OR.or.BETA)
  dataset
}



#### start of server ####

server <- function(input, output, session) {
  
  #### calculate number of cases (n1) and controls (n2), and the fraction of cases (phi) ####
  
  n1 <- reactive({
    ifelse(input$sample_size_specification == 'Total number of subjects + fraction of controls',
           floor(input$n * input$phi), floor(input$n1))
  })
  n2 <- reactive({
    ifelse(input$sample_size_specification == 'Total number of subjects + fraction of controls',
           input$n - floor(input$n * input$phi), floor(input$n2))
  })
  phi <- reactive({
    ifelse(input$sample_size_specification == 'Total number of subjects + fraction of controls',
           input$phi, input$n1 / (input$n1 + input$n2))
  })
  
  #### calculate power for the sample size inputs ####
  
  cutoff <- reactive({
    if (input$type_I_error_criteria == 'Type I error') {
      qchisq(p = 1 - input$alpha, df = 1, ncp = 0, lower.tail = T)
    } else if (input$type_I_error_criteria == 'Family-wise error rate (FWER)') {
      qchisq(p = 1 - input$alpha.FWER / input$p.FWER, df = 1, ncp = 0, lower.tail = T)
    } 
  })
  
  signal.size.vec <- as.vector(outer(c(1,2,5), 10^(-5:-2)))
  power.vec <- reactive({
    pchisq(q = cutoff(), df = 1, ncp = signal.size.vec * (n1() + n2()) , lower.tail = F)
  })
  output$power.vec <- renderText({ power.vec() })
  
  #### generate power matrix for heatmaps ####
  
  f.lim <- c(1e-4, 1-5e-3)
  R.lim <- c(1,110)
  resolution <- 100
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
        temp.mat[i,j] <- signal.size.ana.sol.f(f.vec[j], phi(), R.vec[i]) * (n1() + n2())
      }
    }
    temp.mat
  })
  
  power.mat <- reactive({
    pchisq(q = cutoff(), df = 1, ncp = signal.mat(), lower.tail = F)
  })
  
  #### OR-RAF plot in base R  (deprecated) ####
  
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
  # not pass to output
  
  OR_RAF_baseplot <- reactive({
    p <- plot_ly(x = x.adj(f.vec), y = y.adj.plotly(R.vec), 
                 z = power.mat(), zmin = 0, zmax = 1,
                 colors = "Greys", type = "contour",
                 reversescale = T, hoverinfo = "none",
                 width = 750, height = 700, 
                 contours = list(showlabels = TRUE)) %>%
      # config(displayModeBar = F) %>% 
      layout(xaxis = list(range = x.adj(f.lim), fixedrange=TRUE,
                          tickvals = x.adj(c(0.5, 10^(-1:-4), 1-10^(-1:-2))),
                          ticktext = c(0.5, 0.1, (10^{-2}), (10^{-3}), (10^{-4}),0.9, 0.99),
                          zeroline = FALSE, tickfont = list(size = 20),
                          showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
             yaxis = list(range = y.adj.plotly(R.lim), fixedrange=TRUE,
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
      add_annotations(yref = "paper", xref = "paper", y = 1.04, x = 1.14, 
                      text = "power", 
                      showarrow = F, font = list(size = 25)) %>%
      add_annotations(yref = "paper", xref = "paper", y = 0.8, x = 0.5, 
                      text = "<b>Cases/Controls</b>", showarrow = F, 
                      font=list(size = 40, color = toRGB("grey70"))) %>%
      add_annotations(yref = "paper", xref = "paper", y = 0.7, x = 0.5, 
                      text = paste0("<b>", format(n1(), scientific = FALSE), 
                                    "/", format(n2(), scientific = FALSE), "</b>"), showarrow = F, 
                      font=list(size = 40, color = toRGB("grey70")))
    
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
    p
  }) # end of OR-RAF local plot
  
  
  #### select and load dataset ####
  # plot redered and passed to output
  
  list_of_datasets <- c("Breast carcinoma" = "./data/breast_cancer.tsv",
                        "Coronary heart disease" = "./data/coronary_heart_disease.tsv", 
                        "Type II diabetes mellitus" = "./data/diabetes.tsv")
  
  dataset <- reactive({
    if (input$overlay_example_dataset == T) {
      # print(input$choose_dataset)
      if (input$choose_dataset %in% names(list_of_datasets)) {
        # print(list_of_datasets[input$choose_dataset])
        read_data(list_of_datasets[input$choose_dataset])
      } else {
        # print(input$my_data_upload)
        if (!is.null(input$my_data_upload)) {
        req(input$my_data_upload)
        tryCatch(
          {
            read_data(input$my_data_upload$datapath)
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          }
        ) # end of try-catch data read
        } else {
          NULL
        }
      } # end of user upload
    } # end of data overlay
  })
  
  #output$debug <- renderText({ dataset()[1,1] })
  
  #### overlay data points and render plotly ####
  
  output$OR.RAF.plotly <- renderPlotly({
    if ( input$overlay_example_dataset == T && !is.null(dataset()) ) {
      with(data = dataset(), {
        OR_RAF_baseplot() %>% 
          add_trace(x = x.adj(RISK.ALLELE.FREQUENCY), 
                    y = y.adj.plotly(OR), 
                    marker = list(color = 'rgb(255, 255, 255)', size = 10, opacity = 0.7,
                                  line = list(color = 'rgb(20, 100, 238)', width = 3)), 
                    hoverinfo = "text",
                    text = paste("RAF: ", RISK.ALLELE.FREQUENCY,
                                 "OR: ", round(OR, digits = 3),
                                 '<br>MAPPED GENE:', MAPPED_GENE,
                                 "<br>PUBMEDID: ", PUBMEDID),
                    type = "scatter", mode = 'markers', inherit = F, source = "dataset")
      })
    } else {
      OR_RAF_baseplot()
    }
  }) # end of OR-RAF plotly output
  
  
  #### display details of selected loci by capture click event ####
  
  selected_data_idx <- reactive({
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
  
  output$selection <- renderText({
    if (length(selected_data_idx()) > 0) {
      search_idx <- selected_data_idx()
      paste("<b>Reported gene(s):</b> ", dataset()[search_idx, "REPORTED.GENE.S."], "<br>",
            "<b>Mapped gene:</b> ", dataset()[search_idx, "MAPPED_GENE"], "<br>",
            "<b>Estimated odds ratio:</b> ", dataset()[search_idx, "OR"], "<br>",
            "<b>Estimated allele frequency:</b> ", dataset()[search_idx, "RISK.ALLELE.FREQUENCY"], "<br>",
            "<b>Date added to catalog:</b> ", dataset()[search_idx, "DATE.ADDED.TO.CATALOG"], "<br>",
            "<b>Study:</b> ", dataset()[search_idx, "STUDY"], "<br>",
            "<b>Disease trait:</b> ", dataset()[search_idx, "DISEASE.TRAIT"], "<br>",
            "<b>Initial sample size:</b> ", dataset()[search_idx, "INITIAL.SAMPLE.SIZE"], "<br>",
            "<b>Replication sample size:</b> ", dataset()[search_idx, "REPLICATION.SAMPLE.SIZE"], "<br><br>")
    } else {
      "Click on data points to display details"
    }
  })
  
  
} # end of server