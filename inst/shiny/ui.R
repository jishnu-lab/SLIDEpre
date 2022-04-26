# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Essential Regression"),

  # Sidebar with a slider input for number of bins
  tabsetPanel(
    tabPanel("Plain Essential Regression",
             column(4,
                    h4("File Input:"),
                    # File input
                    fileInput("x_file", "Upload X Data", multiple = TRUE, accept = ".csv"),
                    fileInput("y_file", "Upload Y Data", multiple = TRUE, accept = ".csv"),

                    # Check box if file has header
                    checkboxInput("header", "Header?", TRUE),

                    # Check box if file has row names
                    checkboxInput("row_names", "Row Names?", TRUE)),
             column(4,
                    h4("Parameters:"),
                    # Argument selection
                    numericInput("delta_lbd", "Delta - Lower Bound:", value = 0.01, min = 0.001, max = 0.3),
                    numericInput("delta_ubd", "Delta - Upper Bound:", value = 0.1, min = 0.001, max = 0.3),
                    numericInput("lambda", "Lambda", value = 0.1, min = 0.01, max = 1),
                    numericInput("alpha", "Alpha Level", value = 0.05, min = 0.01, max = 0.1),
                    numericInput("thresh_fdr", "FDR Threshold:", value = 0.2, min = 0.01, max = 0.4),
                    numericInput("rep_cv", "Cross-Validation Replicates", value = 50, min = 1, max = 300)),
             column(4,
                    h4("ER Flags:"),
                    # ER flags
                    checkboxInput("merge", "Merge With Union?", FALSE),
                    checkboxInput("beta_est", "Least Squares Beta Estimation?", TRUE),
                    checkboxInput("conf_int", "Beta Estimation Confidence Intervals?", TRUE),
                    checkboxInput("pred", "Prediction?", TRUE),
                    checkboxInput("diagonal", "Diagonal?", FALSE),
                    checkboxInput("equal_var", "Data Has Equal Variance?", FALSE),
                    checkboxInput("support", "Use Support?", FALSE),
                    checkboxInput("correction", "Use Multiple Testing Correction?", TRUE)),

             # Run button
             actionButton("run_button", "Run!")
    ),
    tabPanel("Prior Essential Regression")
  ),

  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel("Cluster Summary",
               hr(),
               textOutput("k"),
               hr(),
               textOutput("pure_summ"),
               hr(),
               textOutput("mix_summ"),
               hr()),
      tabPanel("Beta Results",
               textOutput("betas")),
      tabPanel("Sample Correlation Heatmap",
               plotOutput("samp_corr")),
      tabPanel("Thresholded Correlation Heatmap",
               plotOutput("thresh_corr"))
    )
  )

)
