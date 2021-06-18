activity_dc_long_fun <- function(dc_used){
	
		tidy_activity <- dc_used$data$input$activity %>% exprs() %>%
			as_tibble(rownames = "id") %>% 
			pivot_longer(cols = -id, names_to = "sample", values_to = "activity") %>%
			separate(col = id, sep = "_", into = c("id","type"))  
		
		samples_to_join <- dc_used$data$input$sample_info #selecting the sample info to join on
		
		tidy_activity <- tidy_activity %>%
			full_join(samples_to_join, by = "sample")  
		
		tidy_activity <-    dc_used$network %>% as_tibble() %>% full_join(tidy_activity, by = "id")
		
		#dc_used$data$input$tidy_activity_long <- tidy_activity
		#dc_used
		tidy_activity
	
}


	activity_dc_wide_fun <- function(tidy_activity){
		
		tidy_activity_wide_out <- tidy_activity %>% filter(!is.na(activity)) %>%  group_by(id,sample_type,type) %>%
			summarize(activity = mean(activity,na.rm = TRUE)) %>% 
			pivot_wider(names_from = sample_type, values_from = activity) 
		
		tidy_activity_wide_out
		
	}
	
activity_dc_wide_filtered_fun <- function(tidy_activity_wide, activity_type, numerator, denominator){
	
	activity_table <- tidy_activity_wide %>% dplyr::filter(type %in% activity_type)
	activity_table <- activity_table %>% mutate(differential_activity = .data[[numerator]]-.data[[denominator]]) %>% dplyr::filter(!is.na(differential_activity))
	
	activity_table <- activity_table %>% dplyr::filter(!is.na(differential_activity))
	activity_table
}

activity_dc_wide_joined_fun <- function(activity_wide_filtered_in, nodes_in ){
	nodes <- nodes_in %>% activate(nodes) %>% as_tibble() 
	tidy_activity_wide_filtered <-  activity_wide_filtered_in %>%
		full_join(nodes, by = "id") %>%
		dplyr::filter(!is.na(differential_activity))
	tidy_activity_wide_filtered
	
}
