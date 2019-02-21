server <- function(input, output) {
  
  # plot in base R
  output$OR.RAF.plot <- renderPlot({
    
    #OR <- 1.3^(-50:50)
    solution.vec <- OR.finder(phi = input$phi, signal.size = 0.01)
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
    n1 <- floor(input$n * input$phi)
    n2 <- input$n - n1
    text(x = x.adj(0.1), y = 10, labels = paste0(n1, "/", n2), cex = 5, pos = 1, col = "grey80")
    
    for (signal.size in c(10^(-5:-1), 5*10^(-5:-2))) {
      solution.vec <- OR.finder(phi = input$phi, signal.size)
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
  }, height = 800, width = 800)
  
  # plot in plotly
  output$OR.RAF.plotly <- renderPlotly({
    n1 <- floor(input$n * input$phi)
    n2 <- input$n - n1
    m <- list(l = 50, r = 20, b = 80, t = 50, pad = 4)
    p <- plot_ly(width = 800, height = 800)  %>% 
      # config(displayModeBar = F) %>% 
      layout(xaxis = list(range = x.adj(c(1e-4,1-1e-3)),
                          tickvals = x.adj(c(0.5, 10^(-1:-4), 1-10^(-1:-2))),
                          ticktext = c(0.5, 0.1, (10^{-2}), (10^{-3}), (10^{-4}),0.9, 0.99),
                          zeroline = FALSE, tickfont = list(size = 20),
                          showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
             yaxis = list(range = y.adj.plotly(c(1,110)), 
                          tickvals = y.adj.plotly(c(1, 2, 5, 10, 20, 50, 100)),
                          ticktext = c(1, 2, 5, 10, 20, 50, 100), 
                          zeroline = FALSE, tickfont = list(size = 20),
                          showline = TRUE, linecolor = toRGB("black"), mirror = "ticks", linewidth = 3),
             showlegend = FALSE, margin = m) %>%
      add_annotations(yref = "paper", xref = "paper", y = 1.08, x = -0.05, 
                      text = "odds ratio", showarrow = F, font=list(size = 25)) %>%
      add_annotations(yref = "paper", xref = "paper", y = -0.1, x = 1, 
                      text = "risk allele frequency (in control group)", 
                      showarrow = F, font = list(size = 25)) %>%
      add_annotations(yref = "paper", xref = "paper", y = 0.6, x = 0.5, 
                      text = "<b>Cases/Controls</b>", showarrow = F, 
                      font=list(size = 40, color = toRGB("grey70"))) %>%
      add_annotations(yref = "paper", xref = "paper", y = 0.5, x = 0.5, 
                      text = paste0("<b>", n1, "/", n2, "</b>"), showarrow = F, 
                      font=list(size = 40, color = toRGB("grey70")))
      
    for (signal.size in c(10^(-5:-1), 5*10^(-5:-2))) {
      solution.vec <- OR.finder(phi = input$phi, signal.size)
      p <- add_lines(p, x = x.adj(solution.vec$f), y = y.adj.plotly(solution.vec$R),
                     line = list(dash = 'dash'), color = I("gray50"),
                     hoverinfo = 'text', text = paste('signal size =', signal.size)) 
    }
    p
  })
  
}