library(shiny)
library(purrr)
library(rlang)

shinyServer(
  function(input, output, session) {

    output$export <- downloadHandler(
      filename = function() {
        paste0("tsgen-", Sys.time(), ".csv")
      },
      content = function(file){
        write_csv(grade_data(), path = file)
      }
    )

    onStop(function(){
      cache <- isolate(reactiveValuesToList(v))
      save(cache, file = "cache.RData")
    })
  }
)
