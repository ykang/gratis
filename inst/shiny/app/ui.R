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
                     icon("floppy-o"),
                     span("Export Grades")
          )
        )
      )
    ),
    dashboardBody(
      tabItems(
        tabItem("tab_struct",
                "ts structure"

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
