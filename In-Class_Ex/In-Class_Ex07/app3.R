pacman::p_load(shiny, sf, tmap, bslib, tidyverse,
               sfdep, shinydashboard, shinythemes)

hunan <- st_read(dsn = "data/geospatial", 
                layer = "Hunan")
data <- read_csv("data/aspatial/Hunan_2012.csv")
hunan_profile <- left_join(hunan, data,
                        by = c("County" = "COUNTY"))

#========================#
###### Shiny UI ######
#========================#  

ui <- navbarPage(
  title = "GLSA Application",
  fluid = TRUE,
  theme=shinytheme("flatly"),
  id = "navbarID",
  tabPanel("GeoVisualisation",
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = "variable",
                           label = "Mapping variable",
                           choices = list("Gross Domestic Product, GDP" = "GDP",
                                          "Gross Domestic Product Per Capita" = "GDPPC",
                                          "Gross Industry Output" = "GIO",
                                          "Output Value of Agriculture" = "OVA",
                                          "Output Value of Service" = "OVS"),
                           selected = "GDPPC"),
               selectInput(inputId = "classification",
                           label = "Classification method:",
                           choices = list("sd" = "sd", 
                                          "equal" = "equal", 
                                          "pretty" = "pretty", 
                                          "quantile" = "quantile", 
                                          "kmeans" = "kmeans", 
                                          "hclust" = "hclust", 
                                          "bclust" = "bclust", 
                                          "fisher" = "fisher", 
                                          "jenks" = "jenks"),
                           selected = "pretty"),
               sliderInput(inputId = "classes",
                           label = "Number of classes",
                           min = 5,
                           max = 10,
                           value = c(6)),
               selectInput(inputId = "colour",
                           label = "Colour scheme:",
                           choices = list("blues" = "Blues", 
                                          "reds" = "Reds", 
                                          "greens" = "Greens",
                                          "Yellow-Orange-Red" = "YlOrRd",
                                          "Yellow-Orange-Brown" = "YlOrBr",
                                          "Yellow-Green" = "YlGn",
                                          "Orange-Red" = "OrRd"),
                           selected = "YlOrRd"),
               sliderInput(inputId = "opacity",
                           label = "Level of transparency",
                           min = 0,
                           max = 1,
                           value = c(0.5))
               ),
             mainPanel(
               tmapOutput("mapPlot",
                          width = "100%", 
                          height = 580)
               )
             )
           ),
  navbarMenu("Global Measures",
             tabPanel("Moran's I"),
             tabPanel("Geary's c"),
             tabPanel("Getis-Ord Global G")
             ),
  navbarMenu("Local Measures",
             tabPanel("Local Moran",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput(inputId = "variable",
                                      label = "Mapping variable",
                                      choices = list("Gross Domestic Product, GDP" = "GDP",
                                                     "Gross Domestic Product Per Capita" = "GDPPC",
                                                     "Gross Industry Output" = "GIO",
                                                     "Output Value of Agriculture" = "OVA",
                                                     "Output Value of Service" = "OVS"),
                                      selected = "GDPPC"),
                          radioButtons(inputId = "Contiguity1",
                                       label = "Contiguity Method",
                                       choices = c("Queen" = TRUE, 
                                                   "Rook" = FALSE),
                                       selected = "TRUE",
                                       inline = TRUE),
                          selectInput("MoranWeights", "Spatial Weights Style",
                                      choices = c("W: Row standardised" = "W",
                                                  "B: Binary" = "B",
                                                  "C: Globally standardised" = "C",
                                                  "U: C / no of neighbours" = "U",
                                                  "minmax" = "minmax",
                                                  "S: Variance" = "S"),
                                      selected = "W"),
                          sliderInput(inputId = "MoranSims", 
                                      label = "Number of Simulations:", 
                                      min = 99, max = 499,
                                      value = 99, step = 100),
                          actionButton("MoranUpdate", "Update Plot"),
                          hr(),
                          radioButtons(inputId = "MoranConf",
                                       label = "Select Confidence level",
                                       choices = c("0.95" = 0.05, 
                                                   "0.99" = 0.01),
                                       selected = 0.05,
                                       inline = TRUE),
                          selectInput("LisaClass", "Select Lisa Classification",
                                      choices = c("mean" = "mean",
                                                  "median" = "median",
                                                  "pysal" = "pysal"),
                                      selected = "mean"),
                          selectInput("localmoranstats", "Select Local Moran's Stat:",
                                      choices = c("local moran(ii)" = "local moran(ii)",
                                                  "expectation(eii)" = "expectation(eii)",
                                                  "variance(var_ii)" = "variance(var_ii)",
                                                  "std deviation(z_ii)" = "std deviation(z_ii)",
                                                  "P-value" = "p_value"),
                                      selected = "local moran(ii)")
                        ),
                        mainPanel(
                          fluidRow(
                            column(6, tmapOutput("LocalMoranMap")),
                            column(6, tmapOutput("LISA"))
                          )
                        )
                      )
                      ),
             tabPanel("Local Gi")
             ),
  navbarMenu("Emerging Hot Spot Analysis")
)

#========================#
###### Shiny Server ######
#========================# 

server <- function(input, output){
    output$mapPlot <- renderTmap({
      tmap_options(check.and.fix = TRUE) +
        tm_shape(hunan_profile)+
        tm_fill(input$variable,
                n = input$classes,
                style = input$classification,
                palette = input$colour,
                alpha = input$opacity) +
        tm_borders(lwd = 0.1,  alpha = 1) +
        tm_view(set.zoom.limits = c(6.5, 8)
                )
    })
    
    #==========================================================
    # Local Measures of Spatial AutoCorrelation
    #==========================================================   
    
    localMIResults <- eventReactive(input$MoranUpdate,{
      
      if(nrow(hunan_profile) == 0) return(NULL)  # Exit if no data
      
      # Computing Contiguity Spatial Weights
      wm_q <- hunan_profile %>%
        mutate(nb = st_contiguity(geometry, 
                                  queen = !!input$Contiguity1),
               wt = st_weights(nb,
                               style = input$MoranWeights))

      # Computing Local Moran's I

      lisa <- wm_q %>%
        mutate(local_moran = local_moran(
          hunan_profile$GDPPC, nb, wt, 
          nsim = as.numeric(input$MoranSims)),
          .before = 5) %>%
        unnest(local_moran)

      lisa <- lisa %>%
        rename("local moran(ii)" = "ii", "expectation(eii)" = "eii",
               "variance(var_ii)" = "var_ii", "std deviation(z_ii)" = "z_ii",
               "p_value" = "p_ii")
      
      return(lisa)       
    })
    
    #==========================================================
    # Render output maps
    #==========================================================
    
    #Render local Moran I statistics
    output$LocalMoranMap <- renderTmap({
      df <- localMIResults()
      
      if(is.null(df) || nrow(df) == 0) return()  # Exit if no data
      
      # Map creation using tmap
      localMI_map <- tm_shape(df) +
        tm_fill(col = input$localmoranstats, 
                style = "pretty", 
                palette = "RdBu", 
                title = input$localmoranstats) +
        tm_borders() +
        tm_view(set.zoom.limits = c(6, 7))
      
      localMI_map 
    })

    #Render LISA map 
    output$LISA <- renderTmap({
      df <- localMIResults()
      if(is.null(df)) return()
      
      
      lisa_sig <- df  %>%
        filter(p_value < as.numeric(input$MoranConf))  
      
      lisamap <- tm_shape(df) +
        tm_polygons() +
        tm_borders() +
        
        tm_shape(lisa_sig) +
        tm_fill(col = input$LisaClass,  
                palette = "-RdBu",  
                title = (paste("Significance:", input$LisaClass))) +
        tm_borders(alpha = 0.4) +
        tm_view(set.zoom.limits = c(6, 7))

      lisamap 
    })
}

shinyApp (ui=ui, server=server)

