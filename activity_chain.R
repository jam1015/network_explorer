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
		tidy_activity
	}
)

activity_dc_wide <- reactive(
	{
		#
		#    dc_used <- activity_dc_long()
		
		tidy_activity <- activity_dc_long()# dc_used$data$input$tidy_activity_long
		
		tidy_activity_wide <- tidy_activity %>% filter(!is.na(activity)) %>%  group_by(id,sample_type,type) %>%
			summarize(activity = mean(activity,na.rm = TRUE)) %>% 
			pivot_wider(names_from = sample_type, values_from = activity) 
		tidy_activity_wide
		#dc_used$data$input$tidy_activity_wide <- tidy_activity_wide
		#dc_used
	}
)

activity_dc_wide_joined <- reactive({
	
	tidy_activity_wide <-  activity_dc_wide()  
	#here is where we filter on activity type 
	
	activity_table <- tidy_activity_wide %>% dplyr::filter(type %in% input$activity_type)
	activity_table <- activity_table %>% mutate(differential_activity = .data[[input$numerator]]-.data[[input$denominator]]) %>% dplyr::filter(!is.na(differential_activity))
	
	activity_table <- activity_table %>% dplyr::filter(!is.na(differential_activity))
	nodes <- dc_layout()$network %>% activate(nodes) %>% as_tibble() 
	
	
	activity_table <-  activity_table %>% full_join(nodes, by = "id") %>% dplyr::filter(!is.na(differential_activity))
	activity_table
})