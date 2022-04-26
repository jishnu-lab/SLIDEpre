
# Define server logic required to draw a histogram
server <- function(input, output) {
  x_data <- reactive({
    ## check data files are valid
    x_file <- input$x_file
    ext <- tools::file_ext(x_file$datapath)
    req(x_file)
    validate(need(ext == "csv", "Please Upload .csv X File"))

    ## read in data
    read.csv(x_file$datapath,
             header = input$header,
             row.names = ifelse(input$row_names, 1, NULL))
  })

  y_data <- reactive({
    ## check data files are valid
    y_file <- input$y_file
    ext <- tools::file_ext(y_file$datapath)
    req(y_file)
    validate(need(ext == "csv", "Please Upload .csv Y File"))

    ## read in data
    read.csv(y_file$datapath,
             header = input$header,
             row.names = ifelse(input$row_names, 1, NULL))
  })

  ## run essential regression
  er_res <- eventReactive(input$run_button, {
    id <- showNotification("Running ER. . . ", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    plainER(y = y_data(),
            x = x_data(),
            sigma = cor(x_data()),
            delta = seq(input$delta_lbd, input$delta_ubd, by = 0.01),
            thresh_fdr = input$thresh_fdr,
            beta_est = input$beta_est,
            conf_int = input$conf_int,
            pred = input$pred,
            lambda = input$lambda,
            rep_cv = input$rep_cv,
            diagonal = input$diagonal,
            merge = input$merge,
            equal_var = input$equal_var,
            alpha_level = input$alpha,
            support = input$support,
            correction = input$correction,
            verbose = F)
  })

  output$pure_summ <- renderText({
    read_res <- readER(er_res())
    pure <- read_res$pure_vars
    c("Pure Variables: ", sort(pure, decreasing = F))
  })

  output$mix_summ <- renderText({
    read_res <- readER(er_res())
    mix <- read_res$mix_vars
    c("Mixed Variables: ", sort(mix, decreasing = F))
  })

  output$k <- renderText({
    c("Number of Clusters: ", er_res()$K)
  })

  output$betas <- renderText({
    betas <- er_res()$beta
    betas <- round(betas, 3)
    paste0(betas, "\n")
  })

  output$samp_corr <- renderPlot({
    makeHeatmap(cor(x_data()), "Sample Correlation Heatmap", T, T)
  })

  output$thresh_corr <- renderPlot({
    makeHeatmap(er_res()$thresh_sigma, "Thresholded Correlation Heatmap", T, T)
  })
}
