setwd("~/Documents/R/netbid_app/netbid_analysis_pipeline_lite/")
library(shiny)
library(shinycssloaders)
source("./drl_control_panel.R")
source("./utility_functions.R")
# utility functions ------------
hover_delay <- 50
# rest of app ---------------

source("./drl_helper_functions.R")
source("./drl_main_functions.R")

#dc <- readRDS("./DATA/dc_example_corrected.RDS") #instead of running the pipeline for now
dc <- readRDS("./DATA/dc_example_corrected_filtered.RDS") #instead of running the pipeline for now

categorical_colors <- readRDS(categorical_color_file)

side_panel <- function(){
  sidebarPanel(width = 2,
    actionButton(inputId = "layout_graph", label = "Run Network Layout"),
    numericInput(inputId = "hdbscan_n", label = "N for hdbscan", min = 5, max = 1000, value = 100, step = 1) ,
    checkboxInput(inputId = "force_hard_cluster", label = "Force Hard Cluster", value = FALSE),
    actionButton(inputId = "run_cluster",label = "Cluster Network"),
    radioButtons(inputId = "activity_type", label = "Select Activity Type", choiceNames = c("Transcription Factor","Signalling Factor"), choiceValues =c("TF","SIG")),
    selectInput(inputId = "numerator",label = "Numerator",choices = NULL),
    selectInput(inputId = "denominator",label = "Denominator",choices = NULL),
    fluidRow(
      column(width = 6,numericInput(inputId = "color_max", label = "color max", value = .005)),
    column(width = 6,numericInput(inputId = "color_min", label = "color min", value = -.005))
    )
  ) 
}

first_row <- function(){
  fluidRow(column(6, 
                  plotOutput("graph_layout_controller",
                             brush = brushOpts(
                               id = "layout_brush",
                               resetOnNew = TRUE
                             ),
                             hover = hoverOpts(id = "hover_main", delay = hover_delay)
                  ) 
  ),
  column(6, 
         plotOutput("graph_layout",dblclick = "graph_dblclick_2", 
                    brush = brushOpts(id = "layout_brush_2", resetOnNew = TRUE),
                    hover = hoverOpts(id = "hover_zoom",hover_delay))
  ),
		   column(6,
          tableOutput("nearest_gene"))
  
  )
  
}


second_row <- function(){
  
  fluidRow(
    column(6,
           plotOutput(
             "colored_graph_layout",dblclick = "graph_dblclick_3", 
             brush = brushOpts(id = "layout_brush_3", resetOnNew = TRUE),
			 hover = hoverOpts(id = "hover_activity",hover_delay)
           )),

    column(6, plotOutput("beeswarm_plot",
                         dblclick = "beeswarm_dblclick",
                         brush = brushOpts(
                           id = "beewarm_ybrush",
                           resetOnNew = TRUE,direction = "y"),
                         hover = hoverOpts(id = "hover_beeswarm", delay = hover_delay)
    )
    )
  )
}

third_row <- function(){
  fluidRow(
    column(12,
           tableOutput(outputId = "nearest_gene")
           )
    
    )
}



main_app <- function() {
  sidebarLayout(
    
    # Sidebar with a slider input
    side_panel(),
    
    # Show a plot of the generated distribution
    mainPanel(
      first_row(),
      second_row()
    )
  )
  
}


ui <- fluidPage(
  
  # Application title
  titlePanel("Functional Network Explorer"),
  tabsetPanel(
    tabPanel("Introduction",fluidRow()),
    tabPanel("Explorer",main_app())
  )
)


server <- function(input, output, session) {
  #can try to get rid of all the dcs that are getting passed around, but I mightbe limited
  # by the architecture of the original transformation pipeline
  
  dc_parameters <- reactive( # just an empty data container that has the parameters
    {
      dc_used <- dc    
      
      dc_used$data$input$sample_info <- dc_used$data$input$sample_info %>%
        tidyr::unite("sample_type",Met.Site,Source,sep = " ")
      dc_used
    }
  )
  
  tissue_types <- reactive(unique(dc_parameters()$data$input$sample_info$sample_type))
  
  observeEvent(dc_parameters(),{
    choices_used <- tissue_types()
    freezeReactiveValue(input, "numerator")
    updateSelectInput(inputId = "numerator", choices = choices_used) 
    freezeReactiveValue(input, "denominator")
    updateSelectInput(inputId = "denominator", choices = choices_used) 
  },priority = 1000000)
  
  
  #can structure this better so that the clusters/layout are dynamically updated in the same reactive 
  layout_dc <- eventReactive(eventExpr = input$layout_graph,
                             {
                               
                               out_dc <- dc_parameters()
                               out_dc$network <- out_dc$network %>%  activate(nodes) %>% mutate(Cluster = factor(NA))
                               out_dc
                             })
  
  layout <- reactive({
    layout_dc()$network
  })
  
  # RUNS CLUSTERING ON THE LAYOUT -----------------
  cluster_dc <- eventReactive(eventExpr = input$run_cluster,
                              {
                                clust_alg_used_hardcode <- "hdbscan_knn"
                                to_clust <- layout_dc()
                                to_clust$parameters$cluster_parameters$hdbscan_knn$force_hard_cluster <- input$force_hard_cluster
                                to_clust$parameters$cluster_parameters$hdbscan_knn$hdbscan_min_pts <- input$hdbscan_n
                                to_clust$parameters$cluster_algorithms <- clust_alg_used_hardcode
                                
                                dc_out <- cluster_network_igraph(
                                  dc_in = to_clust
                                )
                                dc_out})
  
  cluster <- reactive(cluster_dc()$network)
  # can try to have a single ggplot functin that just reacts to changing data; 
  # then can try to have a particular reactivevalue that updates, rather than separate reactives; but this might just be doing what I could o with reactives anyway
  data_used <- reactiveValues(network_used = NULL)
  
  observeEvent(layout() , {
    data_used$network_used <- layout() 
  },priority = 1000
  )
  
  observeEvent(cluster() , {
    data_used$network_used <- cluster() 
  },priority = 1000
  )
  
  
  layout_gg_unified <- eventReactive(data_used$network_used,{
    
    df_used <- data_used$network_used #making the eventreactive here
    nodes <- df_used %>% activate(nodes) %>% as_tibble()
    gg_unified <- nodes %>% ggplot(aes(x = x, y = y)) + 
      theme_prism() + geom_point(aes(color = Cluster, ), alpha = .333333, size = 1) + 
      scale_color_manual(values = categorical_colors, na.value = "black", na.translate  = TRUE) + guides(colour = guide_legend(override.aes = list(size=10, alpha = 1)))
    return(gg_unified)
    
  })
  
  
  #------------------------- HERE IS THE FIRST ACTUAL DRAWING of INPUT -----
  # draiwng the  cntoroller for raw layout -------
  #the code below deals with the zooming options for the plot
  ranges_used_layout <- reactiveValues(x = NULL, y = NULL)
  
  observe({ 
    #here is where we watch to see if they are zooming in with the brush
    brush <- input$layout_brush
    if (!is.null(brush)) {
      ranges_used_layout$x <- c(brush$xmin, brush$xmax)
      ranges_used_layout$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges_used_layout$x <- NULL
      ranges_used_layout$y <- NULL
    }
  })
  
  ranges_controller_second <- reactiveValues(x = NULL, y = NULL) #the first controller
  observeEvent(input$graph_dblclick_2,{
    brush <- input$layout_brush_2
    brush_previous <- input$layout_brush
    if (!is.null(brush)) {
      ranges_used_layout$x <- c(brush$xmin, brush$xmax)
      ranges_used_layout$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_used_layout$x <- c(brush_previous$xmin, brush_previous$xmax)
      ranges_used_layout$y <- c(brush_previous$ymin, brush_previous$ymax)
    }
  })
  
  observeEvent(input$graph_dblclick_3,{
    brush <- input$layout_brush_3
    brush_previous <- input$layout_brush
    if (!is.null(brush)) {
      ranges_used_layout$x <- c(brush$xmin, brush$xmax)
      ranges_used_layout$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_used_layout$x <- c(brush_previous$xmin, brush_previous$xmax)
      ranges_used_layout$y <- c(brush_previous$ymin, brush_previous$ymax)
    }
  })
  
  
  output$graph_layout_controller <- renderPlot({ #this you can drag around on to move the cursor
    layout_gg_unified()}) 
  
  output$graph_layout <- renderPlot({#this one is the main plot
    layout_gg_unified() +   coord_cartesian(xlim = ranges_used_layout$x, ylim = ranges_used_layout$y )})
  
  
  
  #SECOND ACTUAL DRAWING --------------------------
  # drawing the controller for the clustered layout  ------
  
  
  
  # reactive that returns a tidy activity dataframe
  source("./activity_chain.R") #where I have a bunch of activity calculations written
  
  
  activity_long <- reactive(
    {
      dc_used <- dc_parameters()
      activity_long_fun(dc_used)
    }
  )
  
  activity_wide <-  reactive(
    {
      tidy_activity <- activity_long()# dc_used$data$input$tidy_activity_long
      activity_wide_fun(tidy_activity)
    }
  )
  activity_wide_filtered <- reactive({
    
    tidy_activity_wide <-  activity_wide()  
    out <- activity_wide_filtered_fun(tidy_activity_wide = tidy_activity_wide,
                                      activity_type = input$activity_type,
                                      numerator = input$numerator,
                                      denominator = input$denominator
    )
    out
  })
  
  
  network_eventreactive <- eventReactive(data_used$network_used,{ #using the eventreactive for the ignore null property. 
    out <- data_used$network_used                                        #could have use req instead
    out
  })
  
  activity_wide_joined <-  reactive({
    network <- network_eventreactive()
    #browser()
    out <- activity_wide_joined_fun(activity_wide_filtered_in = activity_wide_filtered(),
                                    nodes_in = network)
    #browser()
    #here is where we filter on activity type 
    out
  })
  
  scale_turbo <- reactive({
    scale_color_gradientn(colors =  viridis::viridis_pal(option = "turbo", begin = 0, end = 1)(1000), limits = c(input$color_min,input$color_max), na.value = "#787878", oob = scales::oob_squish)
  })
  
  # reactive that returns plot with activity
  layout_gg_activity <- reactive({
    #here is where we filter on activity type 
    activity_table <-  activity_wide_joined() #chains back to an eventreactive
    
    gg_activity <- activity_table %>% ggplot(aes(x = x, y = y, color = differential_activity)) + geom_point(size = 1,alpha = .95) +
      theme_prism() +   coord_cartesian(xlim = ranges_used_layout$x, ylim = ranges_used_layout$y ) + scale_turbo()
    
    gg_activity
    
  }
  )
  
  output$colored_graph_layout <- renderPlot(layout_gg_activity())
  
  
  beeswarm_yrange <- reactiveValues(y = NULL) # the place where we define the zoom control for the activity plot
  observeEvent(input$beeswarm_dblclick, {
    
    brush <- input$beewarm_ybrush
    if (!is.null(brush)) {
      
      beeswarm_yrange$y <- c(brush$ymin, brush$ymax)
      
    } else {
      
      beeswarm_yrange$y <- NULL
    }
  })
  
  beeswarm_plot <- reactive({
    plot_activity <- activity_wide_joined()  #goes back to an event reactive 
    #browser()
    plot_activity <- plot_activity %>%  dplyr::filter((Cluster != "0") | (is.na(Cluster)))# %>%     mutate(Cluster = factor(Cluster))
    ##browser()
    
    plot_activity %>% ggplot(aes(x = Cluster, y = differential_activity)) +
      geom_boxplot(outlier.shape = NA, aes( group = factor(Cluster))) +
      ggbeeswarm::geom_quasirandom(size = 0.75, alpha = 0.9, aes(color = differential_activity)) + scale_turbo() +
      coord_cartesian(xlim = beeswarm_yrange$x, ylim = beeswarm_yrange$y, expand = FALSE) +
      geom_hline(yintercept = 0, linetype = "dotdash") +
      theme_prism() 
    
  })
  
  
  output$beeswarm_plot <- renderPlot(beeswarm_plot())
 table_nearest_gene <- reactiveValues(table = NULL) 
 output$nearest_gene <- renderTable(table_nearest_gene$table)
  observe({ 
    #here is where we watch to see if they are zooming in with the brush
    hover_main <- input$hover_main
    hover_zoom <- input$hover_zoom
    hover_activity <- input$hover_activity
    
    dataframe_main <- data_used$network_used
    dataframe_zoom <- data_used$network_used
    
    if (!is.null(hover_main) & !is.null(dataframe_main)) {
		dataframe_main <- dataframe_main %>% activate(nodes) %>% as_tibble()
      table_nearest_gene$table <- nearPoints(df = dataframe_main,coordinfo = hover_main,maxpoints = 1, threshold = 100)
    } 
  })
  
  
  observe({ 
    #here is where we watch to see if they are zooming in with the brush
    hover_zoom <- input$hover_zoom
    dataframe_zoom <- data_used$network_used
    
     if(!is.null(hover_zoom) & !is.null(dataframe_zoom)) {
		dataframe_zoom <- dataframe_zoom %>% activate(nodes) %>% as_tibble()
      table_nearest_gene$table <- nearPoints(df = dataframe_zoom,coordinfo = hover_zoom,maxpoints = 1, threshold = 100)
    }
  })
  
  observe({ 
    #here is where we watch to see if they are zooming in with the brush
    
    hover_activity <- input$hover_activity
  
    tissues_unwanted <- tissue_types() 
    tissues_unwanted <- tissues_unwanted[!(tissues_unwanted %in% c(input$numerator, input$denominator))]
    dataframe_activity <- activity_wide_joined()
    
    if(!is.null(hover_activity) & !is.null(dataframe_activity)) {
      table_nearest_gene$table <- nearPoints(df = dataframe_activity,coordinfo = hover_activity,maxpoints = 1, threshold = 100) %>% 
        dplyr::select(-tissues_unwanted)
    }
  })
  
  observe({ 
    #here is where we watch to see if they are zooming in with the brush
    
    hover_beeswarm <- input$hover_beeswarm
    
    tissues_unwanted <- tissue_types() 
    tissues_unwanted <- tissues_unwanted[!(tissues_unwanted %in% c(input$numerator, input$denominator))]
    dataframe_beeswarm <- activity_wide_joined()
    
    if(!is.null(hover_beeswarm) & !is.null(dataframe_beeswarm)) {
      table_nearest_gene$table <- nearPoints(df = dataframe_beeswarm,coordinfo = hover_beeswarm,maxpoints = 1, threshold = 100) %>% 
        dplyr::select(-tissues_unwanted)
    }
  })
  
  
}

shinyApp(ui, server,options = list(port = 7517))
