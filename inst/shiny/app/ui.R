library(shiny)
library(shinydashboard)

shinyUI(
  dashboardPage(
    dashboardHeader(
      title = "tsgeneration",
      titleWidth = "200px",
      tags$li(class = "dropdown", a(href="https://github.com/ykang/tsgeneration", target="_blank", span(icon("github"), " GitHub")))
    ),
    dashboardSidebar(
      width = "200px",
      sidebarMenu(
        menuItem("Structure", tabName = "tab_struct", icon = icon("line-chart")),
        menuItem("Features", tabName = "tab_features", icon = icon("bar-chart")),
        menuItem("Visualise", tabName = "tab_vis", icon = icon("eye")),
        br(),
        tags$li(
          downloadLink("export",
                     style = "margin: 0;",
                     label = NULL,
                     class = "",
                     icon("download"),
                     span("Download series")
          )
        )
      )
    ),
    dashboardBody(
      tabItems(
        tabItem("tab_struct",
                box(
                  title = "Time series structure",
                  selectInput("data_frequency",
                              label = "Time series observation frequency:",
                              choices = c(
                                "Yearly", "Quarterly", "Monthly",
                                "Weekly", "Daily"
                              ),
                              selected = "Yearly"
                  ),
                  uiOutput("seasonal_patterns"),
                  numericInput("data_length",
                               label = "Series length:",
                               min = 30,
                               max = 150,
                               value = 60
                  ),
                  numericInput("data_ngen",
                               label = "Number of series generated:",
                               min = 1,
                               max = 100,
                               value = 1
                  )
                )

        ),

        tabItem("tab_features",
                uiOutput("out_features")
        ),

        tabItem("tab_vis",
                uiOutput("out_plot")
        )
      )
    )
  )
)
