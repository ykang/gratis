library(shiny)
library(purrr)
library(rlang)

shinyServer(
  function(input, output, session) {
    # Constantns
    avgperiods <- c(Yearly = 31557600, Quarterly = 7889400, Monthly = 2629800, Weekly = 604800, Daily = 86400, Hourly = 3600)
    freq_sec <- c(Year = 31557600, Week = 604800, Day = 86400, Hour = 3600)

    interval_seconds <- reactive({
      req(input$data_period)
      avgperiods[input$data_period]
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
