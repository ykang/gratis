library(shiny)
library(purrr)
library(rlang)

shinyServer(
  function(input, output, session) {
    # Constants
    avgperiods <- c(Yearly = 31557600, Quarterly = 7889400, Monthly = 2629800, Weekly = 604800, Daily = 86400, Hourly = 3600)
    freq_sec <- c(Year = 31557600, Week = 604800, Day = 86400, Hour = 3600)
    acf_features <- c("x_acf1", "diff1_acf1", "diff2_acf1", "x_acf10", "diff1_acf10", "diff2_acf10")
    pacf_features <- c("x_pacf5", "diff1x_pacf5", "diff2x_pacf5")
    stl_features <- c("trend", "spike", "linearity", "curvature", "e_acf1", "e_acf10")
    stl_seas_features <- c("seasonal_strength", "peak", "trough")
    shift_features <- c("max_level_shift", "time_level_shift",
                        "max_kl_shift", "time_kl_shift",
                        "max_var_shift", "time_var_shift")
    behaviour_features <- c("entropy", "nonlinearity", "hurst", "stability", "lumpiness",
                            "unitroot_kpss", "unitroot_pp")
    heterogeneity_features <- c("arch_acf", "garch_acf", "arch_r2", "garch_r2")

    all_features <- reactive({
      out <- c(acf_features, pacf_features,
        stl_features,
        shift_features,
        behaviour_features,
        heterogeneity_features
      )
      if(!is_empty(seasonal_freq())){
        out <- c(out,
                 "seas_acf1",
                 "seas_pacf",
                 outer(stl_seas_features, seq_along(seasonal_freq()), "paste0")
                 )
      }
      out
    })


    interval_seconds <- reactive({
      req(input$data_period)
      avgperiods[input$data_period]
    })

    seasonal_freq <- reactive({
      freq_sec[input$data_freq] / interval_seconds()
    })

    output$series_period <- renderUI({
      selectInput("data_period",
                  label = "Time series observation frequency:",
                  choices = names(avgperiods),
                  selected = avgperiods[1]
      )
    })

    output$seasonal_patterns <- renderUI({
      seconds <- freq_sec / interval_seconds()
      seasonal_periods <- names(seconds)[seconds > 1]
      if(!is_empty(seasonal_periods)){
        checkboxGroupInput("data_freq", label = "Seasonal period(s):",
                           choices = seasonal_periods,
                           selected = seasonal_periods,
                           inline = TRUE)
      }
    })

    output$feature_diff <- renderUI({
      do.call("tagList", c(
          list(numericInput("par_ndiff", "Number of differences:", value = 0, min = 0, max = 2)),
          map(names(seasonal_freq()),
              ~ numericInput(
                paste0("par_nsdiff_", .x),
                paste0("Number of seasonal differences [", .x, "]:"),
                value = 0, min = 0, max = 2)
          )
        )
      )
    })

    output$feature_acf <- renderUI({
      features <- acf_features
      if(!is_empty(seasonal_freq())){
        features <- c(features, "seas_acf1")
      }

      do.call("tagList",
              map(features,
                  ~ numericInput(
                    paste0("par_", .x),
                    paste0(.x),
                    value = 0, min = -1, max = 1, step = 0.01)
              )
      )
    })

    output$feature_pacf <- renderUI({
      features <- pacf_features
      if(!is_empty(seasonal_freq())){
        features <- c(features, "seas_pacf")
      }

      do.call("tagList",
              map(features,
                  ~ numericInput(
                    paste0("par_", .x),
                    paste0(.x),
                    value = 0, min = -1, max = 1, step = 0.01)
              )
      )
    })

    output$feature_stl <- renderUI({
      features <- stl_features

      io <- map(stl_features,
          ~ numericInput(
            paste0("par_", .x),
            .x,
            value = 0, step = 0.01)
      )

      if(!is_empty(seasonal_freq())){
        seas <- cross2(names(seasonal_freq()), stl_seas_features) %>%
          map(~ numericInput(
            paste0("par_", .x[[2]], match(.x, names(seasonal_freq()))),
            paste0(.x[[2]], " [", .x[[1]], "]"),
            value = 0, step = 0.01))

        io <- c(io, seas)
      }

      do.call("tagList", io)
    })

    output$feature_shift <- renderUI({
      timevars <- grepl("time", shift_features, fixed = TRUE)

      time_shift <- map(shift_features[timevars],
                  ~ numericInput(
                    paste0("par_", .x),
                    paste0(.x),
                    value = 0, min = 0, max = input$data_length, step = 1)
      )

      max_shift <- map(shift_features[!timevars],
                      ~ numericInput(
                        paste0("par_", .x),
                        paste0(.x),
                        value = 0, step = 0.01)
      )

      do.call("tagList",
              c(time_shift, max_shift)
      )
    })

    output$feature_behave <- renderUI({
      do.call("tagList",
              map(behaviour_features,
                  ~ numericInput(
                    paste0("par_", .x),
                    paste0(.x),
                    value = 0, step = 0.01)
              )
      )
    })

    output$feature_behave <- renderUI({
      do.call("tagList",
              map(behaviour_features,
                  ~ numericInput(
                    paste0("par_", .x),
                    paste0(.x),
                    value = 0, step = 0.01)
              )
      )
    })

    output$feature_heterogeneity <- renderUI({
      do.call("tagList",
              map(heterogeneity_features,
                  ~ numericInput(
                    paste0("par_", .x),
                    paste0(.x),
                    value = 0, min = -1, max = 1, step = 0.01)
              )
      )
    })

    observeEvent(input$btn_gen, {
      # features: list of relevant features
      # target: target features for fitness in GA
      # seasonal: number of seasonal components
      # n: length of series to generate
      # freq: frequencies of time series seasonality
      # nComp: Number of components in mixture models
      # selected.features: Features actually used

      ga_len <- c(10,17,35)[min(2, length(seasonal_freq())+1)]
      ga_min <- rep(0, ga_len)
      ga_max <- rep(1, ga_len)

      ga_ts(
        type = "real-valued", fitness = fitness_ts, features = features, seasonal = seasonal,
        input$data_length, # n for fitness_ts
        freq = seasonal_freq(), target = target, nComp = 3, selected.features = all_features(),
        min = ga_min,
        max = ga_max,
        parallel = TRUE, popSize = 30, maxiter = 100,
        pmutation = 0.3, pcrossover = 0.8, maxFitness = -0.05,
        run = 30, keepBest = TRUE, monitor = GA::gaMonitor
      )
    })

    output$export <- downloadHandler(
      filename = function() {
        paste0("tsgen-", Sys.time(), ".csv")
      },
      content = function(file){
        write_csv(generated_ts(), path = file)
      }
    )
  }
)
