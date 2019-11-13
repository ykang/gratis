library(gratis)
library(tsfeatures)
library(Mcomp)
library(tscompdata)
features <- c(
  "ndiffs",
  "nsdiffs1",
  "acf_features",
  "pacf_features",
  "entropy",
  "nonlinearity",
  "hurst",
  "stability",
  "lumpiness",
  "unitroot_kpss",
  "unitroot_pp",
  "max_level_shift",
  "max_kl_shift",
  "max_var_shift",
  "stl_features",
  "heterogeneity"
)
featureNames <- suppressWarnings(
  unique(c(
    names(tsfeatures(M3[[1907]]$x, features = features)),
    names(tsfeatures(tscompdata::nn5[[1]], features = features, na.action = na.interp))
  ))
)
featureNames <- setdiff(featureNames, c("seasonal_period", "seasonal_period1", "seasonal_period2", "nperiods"))
ui <- fluidPage(
  tags$head(tags$style(HTML("
                            #text {
                            font-size: 25px;
                            }
                            "))),
  tags$head(tags$style(HTML("
                            #text2 {
                            font-size: 25px;
                            }
                            "))),

  titlePanel("Generation of time series with controllable features"),
  fluidRow(
    column(3,
      style = "font-size: 15pt",
      wellPanel(
        selectInput("seasonal_period",
          label = "Please select a seasonal pattern:",
          choices = c(
            "Yearly", "Quarterly", "Monthly",
            "Weekly", "Daily"
          ),
          selected = "Yearly"
        ),
        numericInput("Length",
          label = "Length:",
          min = 30,
          max = 150,
          value = 60
        ),
        numericInput("Number",
          label = "How many time series do you want?",
          min = 1,
          max = 100,
          value = 1
        ),
        h4("Please set values for the features you are interested in:"),
        eval(parse(text = paste("wellPanel(", paste(
          "numericInput(featureNames[", 1:9, "], label = featureNames[", 1:9, "], value = 0)",
          sep = "", collapse = ","
        ), ")")))
      ),
      actionButton("do", "Generate"),
      # Button
      downloadButton("downloadData", "Download")
    ),
    column(3,
      style = "font-size: 15pt",
      eval(parse(text = paste("wellPanel(", paste(
        "numericInput(featureNames[", 10:21, "],label = featureNames[", 10:21, "], value = 0)",
        sep = "", collapse = ","
      ), ")")))
    ),
    column(3,
      style = "font-size: 15pt",
      eval(parse(text = paste("wellPanel(", paste(
        "numericInput(featureNames[", 22:33, "],label = featureNames[", 22:33, "], value = 0)",
        sep = "", collapse = ","
      ), ")")))
    ),
    column(3,
      style = "font-size: 15pt",
      eval(parse(text = paste("wellPanel(", paste(
        "numericInput(featureNames[", 34:45, "],label = featureNames[", 34:45, "], value = 0)",
        sep = "", collapse = ","
      ), ")")))
    )
  ),
  fluidRow(
    column(12, mainPanel(
      h3("Please see below the generated time series."),
      plotOutput("plot"),
      verbatimTextOutput("text"),
      p("\n"),
      verbatimTextOutput("text2")
    ))
  )
)

server <- function(input, output) {
  evolved.ts <- eventReactive(input$do, {
    selected.features <- unlist(lapply(featureNames, function(x) {
      ifelse(input[[x]] != "0", x, NA)
    })) %>% na.omit()
    n.features <- length(selected.features)
    target <- rep(NA, n.features)
    for (i in 1:n.features) {
      target[i] <- as.numeric(input[[selected.features[i]]])
    }
    seasonal <-
      ifelse(input$seasonal_period == "Daily", 2, ifelse(input$seasonal_period == "Yearly", 0, 1))
    ga_min <-
      if (seasonal == 0) {
        c(rep(0, 10))
      } else if (seasonal == 1) {
        c(rep(0, 17))
      } else {
        c(rep(0, 35))
      }
    ga_max <-
      if (seasonal == 0) {
        c(rep(1, 10))
      } else if (seasonal == 1) {
        c(rep(1, 17))
      } else {
        c(rep(1, 35))
      }
    freqMat <- as.list(c(1, 4, 12, 52, "multiple"))
    names(freqMat) <- c("Yearly", "Quarterly", "Monthly", "Weekly", "Daily")
    freqMat[["Daily"]] <- c(7, 365)
    freqMat <- lapply(freqMat, as.numeric)
    freq <- freqMat[[input$seasonal_period]]
    evolved.ts <- c()
    withProgress(message = "Generating data", detail = "0%", {
      ptm <- proc.time()
      while (ifelse(is.null(dim(evolved.ts)), 0 < 1, dim(evolved.ts)[2] < input$Number)) {
        GA <- ga_ts(
          type = "real-valued", fitness = fitness_ts, features = features, seasonal = seasonal,
          input$Length, freq, target, 3, selected.features,
          n = input$Length,
          min = ga_min,
          max = ga_max,
          parallel = TRUE, popSize = 30, maxiter = 100,
          pmutation = 0.3, pcrossover = 0.8, maxFitness = -0.05,
          run = 30, keepBest = TRUE, monitor = GA::gaMonitor
        )
        evolved.ts.new <-
          unique(do.call(
            cbind,
            eval(parse(text = paste("list(", paste("GA@bestSol[[GA@iter - ", 0:(GA@run - 1), "]]", sep = "", collapse = ","), ")")))
          ), MARGIN = 2)
        evolved.ts <- cbind(evolved.ts, evolved.ts.new)
        n.ts.evolved <- dim(evolved.ts)[2]
        ## figures here are only using the last GA process
        # output$text <- renderPrint({
        #   sprintf('The fitness values are within [%.3f, %.3f]', max(GA@fitness) - 0.03, max(GA@fitness))
        # })
        ifelse(length(freq) == 1, optimal.ts <- ts(evolved.ts[, 1], frequency = freq),
          optimal.ts <- msts(evolved.ts[, 1], seasonal.periods = freq)
        )
        output$text2 <- renderPrint({
          as.data.frame(round(tsfeatures(optimal.ts,
            features = features
          )
          %>%
            dplyr::select(selected.features), 3))
        })
        incProgress(dim(evolved.ts.new)[2] / input$Number, detail = paste0(min(round(n.ts.evolved / input$Number, 2) * 100, 100), "%"))
      }
      print(proc.time() - ptm)
    })
    if (length(freq) == 1) {
      ts(evolved.ts[, 1:(input$Number)], frequency = freq)
    } else {
      msts(evolved.ts[, 1:(input$Number)], seasonal.periods = freq)
    }
  })
  output$plot <- renderPlot({
    autoplot(evolved.ts(), ylab = "Evolved Time Series") + ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 20))
  })
  output$downloadData <- downloadHandler(
    filename <- function() {
      paste("test", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(evolved.ts(), file, row.names = FALSE)
    }
  )
}
shinyApp(ui = ui, server = server)
