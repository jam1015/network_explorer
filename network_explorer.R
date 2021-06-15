#adding comment
setwd("~/Documents/R/netbid_app/netbid_analysis_pipeline_lite/")
library(shiny)
library(shinycssloaders)
source("./drl_control_panel.R")
source("./utility_functions.R")
# utility functions ------------

# rest of app ---------------

source("./drl_helper_functions.R")
source("./drl_main_functions.R")

dc <- readRDS("./DATA/dc_example_corrected.RDS") #instead of running the pipeline for now

categorical_colors <- readRDS(categorical_color_file)
ui <- fluidPage(
  
  # Application title
  titlePanel("Functional Network Explorer"),
  
  sidebarLayout(
    
    # Sidebar with a slider input
    sidebarPanel(
      sliderInput("filter_percentage", label = h3("Filter Percentage"), min = 0, 
                  max = 100, value = 100),
      actionButton(inputId = "layout_graph", label = "Run Network Layout"),
      numericInput(inputId = "hdbscan_n", label = "N for hdbscan", min = 5, max = 1000, value = 500, step = 1) ,
      checkboxInput(inputId = "force_hard_cluster", label = "Force Hard Cluster", value = FALSE),
      actionButton(inputId = "run_cluster",label = "Cluster Network"),
      radioButtons(inputId = "activity_type", label = "Select Activity Type", choiceNames = c("Transcription Factor","Signalling Factor"), choiceValues =c("TF","SIG")),
      selectInput(inputId = "numerator",label = "Numerator",choices = NULL),
      selectInput(inputId = "denominator",label = "Denominator",choices = NULL),
      numericInput(inputId = "color_max", label = "color max", value = .5),
      numericInput(inputId = "color_min", label = "color min", value = -.5)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(column(6, 
                      plotOutput("graph_layout_controller",
                                 brush = brushOpts(
                                   id = "layout_brush",
                                   resetOnNew = TRUE
                                 )
                      ) 
      ),
      column(6, 
             plotOutput("graph_layout")
      )
      
      ),
      fluidRow(
        column(6,
               plotOutput(
                 "colored_graph_layout"
               )),
        column(6, plotOutput("beeswarm_plot",
                             dblclick = "plot1_dblclick",
                             brush = brushOpts(
                               id = "plot1_brush",
                               resetOnNew = TRUE,direction = "y")
        )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
  dc_parameters <- reactive( # just an empty data container that has the parameters
    {
      dc_used <- dc    
      dc_used$data$input$sample_info <- dc_used$data$input$sample_info %>%
        tidyr::unite("sample_type",Met.Site,Source,sep = " ")
      dc_used
    }
  )
  
  dc_filtered <-  reactive(# to try to cut down on the data used
    {
      dc_used <- dc_parameters()
      sample_percentage <- input$filter_percentage
      nodes <- dc_used$network %>% activate(nodes) 
      n_points <- dim(as_tibble(nodes))[1]
      sample_n <- ceiling((sample_percentage/100)*n_points)
      indices <- sample.int(n = n_points,size = sample_n)
      dc_used$network <- dc_used$network %>% activate(nodes) %>% slice(indices)
      dc_used
    }    
  )
  
  dc_layout_reactive <- reactive({dc_filtered()}) #will eventually run the actual layout
  
  activity_plot_trigger <- reactive(list(input$run_cluster,input$layout_graph))
  
  dc_layout <- eventReactive(eventExpr = input$layout_graph,
                             {
                               dc_layout_reactive() 
                             }
  )
  # RUNS CLUSTERING ON THE LAYOUT -----------------
  dc_cluster <- eventReactive(eventExpr = input$run_cluster,
                              {
                                to_clust <- dc_layout_reactive()
                                to_clust$parameters$cluster_parameters$hdbscan_knn$force_hard_cluster <- input$force_hard_cluster
                                to_clust$parameters$cluster_parameters$hdbscan_knn$hdbscan_min_pts <- input$hdbscan_n
                                dc_out <- cluster_network_igraph(
                                  dc_in = to_clust
                                )
                                dc_out
                              } 
  )
  
  dc_layout_gg_raw <- reactive({ #Called raw because there is no clustering
    dc_used <- dc_layout() 
    nodes <- dc_used$network %>% activate(nodes) %>% as_tibble()
    gg_raw <- nodes %>% ggplot(aes(x = x, y = y)) + 
      theme_prism() + geom_point(color = "black", alpha = .25, size = .5) + 
      scale_color_manual(values = categorical_colors)# + theme(aspect.ratio = 1)
    gg_raw
#    dc_used$data$output$gg_layout <- gg_raw
#    dc_used
  }
  )
  #------------------------- HERE IS THE FIRST ACTUAL DRAWING of INPUT -----
  # draiwng the  cntoroller for raw layout -------
  plots_out <- reactiveValues(plot_used = NULL)
  observeEvent(input$layout_graph ,{ plots_out$plot_used <- dc_layout_gg_raw()#$data$output$gg_layout
  })
  output$graph_layout_controller <- renderPlot(plots_out$plot_used)
  output$graph_layout <- renderPlot(plots_out$plot_used +   coord_cartesian(xlim = ranges2$x, ylim = ranges2$y ))
  
  dc_layout_gg_forced <- reactive({ #a reactive that just makes a ggplot layout of the clustered output with forcing of hard clustering
    dc_used <- dc_cluster() 
    nodes <- dc_used$network %>% activate(nodes) %>% as_tibble()
    gg_forced <- nodes %>% ggplot(aes(x = x, y = y,color = factor(Cluster))) + 
      theme_prism() + geom_point(size = .5) + 
      scale_color_manual(values = categorical_colors)#+ theme(aspect.ratio = 1)
    gg_forced
  }
  )
  
  dc_layout_gg_natural <- reactive({ #a reactive that just makes a ggplot layout of the clustered output, but no forcing of hard clustering
    dc_used <- dc_cluster() 
    nodes <- dc_used$network %>% activate(nodes) %>% as_tibble()
    no_zero_cluster <- nodes %>% dplyr::filter(Cluster != 0)
    zero_cluster   <- nodes %>% dplyr::filter(Cluster == 0)
    gg_natural <- no_zero_cluster %>% ggplot(aes(x = x, y = y, color = factor(Cluster))) + geom_point(size = .5) + 
      scale_color_manual(values = categorical_colors) + 
      geom_point(data = zero_cluster, size = .5,alpha = .25) + theme_prism()
    gg_natural
  }
  )
  
  #SECOND ACTUAL DRAWING ---------------------------
  # drawing the controller for the clustered layout  ------
  
  observeEvent(input$run_cluster,{
    if(input$force_hard_cluster){
     plots_out$plot_used <- dc_layout_gg_forced()
    } else{
     plots_out$plot_used <-  dc_layout_gg_natural()
    }
  }
  )
  
  
  # reactive that returns data container with a tidy activity dataframe
  activity_dc_long <- reactive(
    {
      dc_used <- dc_layout()
      
      tidy_activity <- dc_used$data$input$activity %>% exprs() %>%
        as_tibble(rownames = "id") %>% 
        pivot_longer(cols = -id, names_to = "sample", values_to = "activity") %>%
        separate(col = id, sep = "_", into = c("id","type"))  
      
      samples_to_join <- dc_used$data$input$sample_info #selecting the sample info to join on
      
      tidy_activity <- tidy_activity %>%
        full_join(samples_to_join, by = "sample")  
      
      tidy_activity <-    dc$network %>% as_tibble() %>% full_join(tidy_activity, by = "id")
      
      #dc_used$data$input$tidy_activity_long <- tidy_activity
      #dc_used
    }
  )
  
  activity_dc_wide <- reactive(
    {
      dc_used <- activity_dc_long()
      
      tidy_activity <- dc_used$data$input$tidy_activity_long
      
      tidy_activity_wide <- tidy_activity %>% filter(!is.na(activity)) %>%  group_by(id,sample_type,type) %>%
        summarize(activity = mean(activity,na.rm = TRUE)) %>% 
        pivot_wider(names_from = sample_type, values_from = activity) 
      
      dc_used$data$input$tidy_activity_wide <- tidy_activity_wide
      dc_used
    }
  )
  
  activity_dc_wide_joined <- reactive({
    dc_used <-  activity_dc_wide()  
    #here is where we filter on activity type 
    activity_table <- dc_used$data$input$tidy_activity_wide %>% dplyr::filter(type %in% input$activity_type)
    activity_table <- activity_table %>% mutate(differential_activity = .data[[input$numerator]]-.data[[input$denominator]]) %>% dplyr::filter(!is.na(differential_activity))
    activity_table <- activity_table %>% dplyr::filter(!is.na(differential_activity))
    nodes <- dc_cluster()$network %>% activate(nodes) %>% as_tibble() 
    
    
    activity_table <-  activity_table %>% full_join(nodes, by = "id") %>% dplyr::filter(!is.na(differential_activity))
    dc_used$data$input$tidy_activity_wide_joined <- activity_table
    dc_used
  })
  
  scale_turbo <- reactive({
  scale_color_gradientn(colors =  viridis::viridis_pal(option = "turbo", begin = 0, end = 1)(1000), limits = c(input$color_min,input$color_max), na.value = "#787878")
    
  })
  # reactive that returns plot with activity
  dc_layout_gg_activity <- reactive({
    dc_used <-  activity_dc_wide_joined()  
    #here is where we filter on activity type 
    activity_table <-  dc_used$data$input$tidy_activity_wide_joined
    gg_activity <- activity_table %>% ggplot(aes(x = x, y = y, color = differential_activity)) + geom_point(size = 1,alpha = .95) + scale_turbo() + 
      theme_prism()
    dc_used$data$output$gg_layout <- gg_activity
    dc_used
  }
  )
  
  #    drawing the layout colored by activity----------- 
  observeEvent(activity_plot_trigger(),
               {output$colored_graph_layout <- renderPlot(
                 expr = { 
                   gg <- dc_layout_gg_activity()$data$output$gg_layout
                   x_limit <- ranges2$x
                   y_limit <- ranges2$y
                   gg +  coord_cartesian(xlim = x_limit, ylim = y_limit )
                 }
               )
               }
  )
  
  observeEvent(input$run_cluster,{
    output$beeswarm_plot <- renderPlot(
      expr = { 
        dc_used <- activity_dc_wide_joined()
        plot_activity <- dc_used$data$input$tidy_activity_wide_joined %>% dplyr::filter(Cluster != 0)# %>%     mutate(Cluster = factor(Cluster))
      
         plot_activity %>% ggplot(aes(x = factor(Cluster), y = differential_activity)) + geom_boxplot(outlier.shape = NA, aes( group = factor(Cluster))) + ggbeeswarm::geom_quasirandom(size = .25, alpha = .5,aes(color = differential_activity)) + scale_turbo() +
         coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) + theme_prism()

      }
    )
    
  })
  
  
  #drawing the numerator for the color control
  observeEvent(activity_dc_wide(),{
    freezeReactiveValue(input,"numerator")
    choices_used <- unique(activity_dc_wide()$data$input$sample_info$sample_type)
    updateSelectInput(inputId = "numerator", choices = choices_used) 
  }
  )
  #drawing the denominator for the color control
  observeEvent(activity_dc_wide(),{
    freezeReactiveValue(input,"denominator")
    choices_used <- unique(activity_dc_wide()$data$input$sample_info$sample_type)
    updateSelectInput(inputId = "denominator", choices = choices_used) 
  }
  )
  
  
  ranges <- reactiveValues(y = NULL)
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {

      ranges$y <- c(brush$ymin, brush$ymax)

    } else {

      ranges$y <- NULL
    }
  })

  
  #the code below deals with the zooming options for the plot
  ranges2 <- reactiveValues(x = NULL, y = NULL)
  observe({ #here is where we watch to see if they are zooming in with the brush
    brush <- input$layout_brush
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      ranges2$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
}

shinyApp(ui, server)
