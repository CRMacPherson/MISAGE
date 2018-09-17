custom.dashboardMenu = function (..., title = NULL, titleWidth = NULL, disable = FALSE, .list = NULL) 
  {
    items <- c(list(...), .list)
    lapply(items, shinydashboard:::tagAssert, type = "li", class = "dropdown")
    titleWidth <- validateCssUnit(titleWidth)
    custom_css <- NULL
    if (!is.null(titleWidth)) {
      custom_css <- tags$head(tags$style(HTML(gsub("_WIDTH_", 
                                                   titleWidth, fixed = TRUE, "\n      @media (min-width: 200px) {\n        .main-header > .navbar {\n          margin-left: _WIDTH_;\n        }\n        .main-header .logo {\n          width: _WIDTH_;\n        }\n      }\n    "))))
    }
    tags$header(class = "main-header", custom_css, style = if (disable) 
      "display: none;", span(class = "logo", title), tags$nav(class = "navbar navbar-static-top", 
                                                              role = "navigation", span(shiny::icon("bars"), style = "display:none;"),
                                                              div(class = "navbar-custom-menu", tags$ul(class = "nav navbar-nav", 
                                                                                                        items))))
  }