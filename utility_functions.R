#making a change just to have it in git
assign_parameters <- function(df_in = NA , env = parent.frame()){
	#checking if there is a parameters dataframe in global environment
	if(is_na(df_in) || is_null(df_in)) {
		if ("parameters" %in% ls(envir = .GlobalEnv)) {
			df_in <- get("parameters", envir = .GlobalEnv)
		}  else {
			stop("parameters df missing from global environment")
		}
	}
	
	passed_args <- as.list(env)[names(as.list(env)) ] #valuates the values in the calling environment and filters to include only the args to the function
	
	missing_vals <-  map_lgl(passed_args,is_missing) #seeing which are missing
	wanted_params <- names(missing_vals)[missing_vals]   #keeping the names of just the missing values
	
	##### searching through misssing values in parent environments
	values_collector <- vector("list", length(wanted_params))
	names(values_collector) <- wanted_params
	
	for(nam in names(values_collector)){
		n <- 1
		found <- FALSE
		reached_global <- FALSE
		
		while( !reached_global & !found ){
			
			
			candidate_value <-  mget(nam, envir = parent.frame(n),ifnotfound = "varnotfound", mode = "any",inherits = FALSE)
			if(length(candidate_value) == 1 && is.list(candidate_value) && candidate_value[[1]] == "varnotfound"){
				found_var_in_env <- FALSE
			} else {
				found_var_in_env <- TRUE
			}
			
			if(!map_lgl(.x = candidate_value,.f = is_null) & 
				 !map_lgl(.x = candidate_value,.f = is_missing) & 
				 found_var_in_env){
				
				values_collector[[nam]] <- candidate_value[[nam]]
				found <- TRUE  
				
				
			}  else if(identical(parent.frame(n),globalenv())) {
				reached_global <- TRUE
			} else {
				n <- n + 1 
			}
			
		}
		
	}
	
	to_assign_from_env <- values_collector[!map_lgl(.x = values_collector, .f = is_null)]
	to_assign_from_params_df <- names(values_collector)[!(names(values_collector) %in% names(to_assign_from_env))]
	
	
	df_filtered <- df_in %>%  dplyr::filter(parameter %in% to_assign_from_params_df) %>% dplyr::select(parameter,value) #unnests and filters the df to have only the parameters we want
	assign_var <- function(parameter,value,envi){  assign(x = parameter, value = value, envir = envi)  }  #function to the variables to the parent environment
	pwalk(.l =   df_filtered ,                                                          .f = assign_var, envi = env)  #doing the actual assignment of the variables to the parent environment
	pwalk(.l = list(parameter = names(to_assign_from_env), value = to_assign_from_env), .f = assign_var, envi = env)
	
	
	#df in is a dataframe with the name of the parameter, 
	#calling_fun <- deparse(sys.call(-1)) #gets the name of the calling fucntion
	#calling_fun <- calling_fun %>% str_split("\\(") %>% unlist() %>% `[`(1)   #gets name of  calling fun without the parentheses
	#wanted_params <- formalArgs(calling_fun)  #getting the names of the arguments of the calling function
	
}



#exponentiates base two and decrements by 1
exp_dec <- function(x) {
	(2^x)-1
}

#symmetric set difference
symdiff <- function( x, y) {
	setdiff( union(x, y), intersect(x, y))
}

# takes a vector and checks for duplicates, INCLUDING the first instance
all_dup <- function (value){
	duplicated(value) | duplicated(value, fromLast = TRUE)
} 

#joins and coalesces
coalesce_join <- function(x, y,  by = NULL, suffix = c(".x", ".y"),   join = dplyr::full_join, ...) {
	assign_parameters()
	
	joined <- join(x, y, by = by, suffix = suffix, ...)
	
	# names of desired output
	cols <- union(names(x), names(y))  #keeps original names w/o the suffixes added
	
	to_coalesce <- names(joined)[!names(joined) %in% cols]  #names not in cols because there is a suffix addd
	suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)] #the suffixes of to_coalesce (selecting from the input suffix vector)
	# remove suffixes and deduplicate
	to_coalesce <- unique(substr(
		to_coalesce, 
		1, 
		nchar(to_coalesce) - nchar(suffix_used)
	))  #explicitly discarding suffix of what we are coalescing and then, gets unique values 
	
	coalesced <- to_coalesce %>% purrr::map_dfc(
		~dplyr::coalesce(
			joined[[paste0(.x, suffix[1])]],   
			joined[[paste0(.x, suffix[2])]]
		)
	)    #pipes in unique values of to_coalesce, reattaches suffixes, coalesces, colects into dataframe 
	names(coalesced) <- to_coalesce  #can assign names because it happens in the same order
	
	dplyr::bind_cols(joined, coalesced)[cols] #binds coalesced columns (no suffixes) with new ones, selects all with now suffixes
}

#makies things title cae and does a replacement
make_title_case <- function(str){
	str <- sub(pattern = "_",replacement = " ",fixed = TRUE, x = str)
	str <- toTitleCase(str)
}

check_if_missing <- function(env,var){
	passed_args <- as.list(env)[names(as.list(env)) ] #valuates the values in the calling environment and filters to include only the args to the function
	
	if(var %in% names(passed_args)){
		missing_vals <-  map_lgl(passed_args,is_missing) #seeing which are missing
		missing_names <- names(missing_vals)[missing_vals]   #keeping the names of just the missing values
		return(var %in% missing_names)
	}  else{
		return(TRUE)
	}
	
}





standardize_cols <- function(x) { #make columns sum to 1
	x <- sweep(x,2,colSums(x),"/")
	x[is.na(x)] <- 0
	return(x)
}

standardize_rows <- function(x) {  #make rows sum to 1 
	x <- sweep(x,1,rowSums(x),"/")
	x[is.na(x)] <- 0
	return(x)
}


#centering mean and projecting onto unit hypersphere

spherize_cols <- function(x) {
	centroid <- rowMeans(x)
	x <- sweep(x,1,centroid,"-")
	radius <- sqrt(colSums(x^2))
	x <- sweep(x,2,radius,"/")
	return(x)
}



#centering mean and projecting onto unit hypersphere
spherize_rows <- function(x) {
	centroid <- colMeans(x)
	x <- sweep(x,2,centroid,"-")
	radius <- sqrt(rowSums(x^2))
	x <- sweep(x,1,radius,"/")
	return(x)
}



order_rows <- function(x) {
	ordered <- rowRanks(x = x, ties.method = "average", preserveShape = TRUE)
	rownames(ordered) <- rownames(x)
	colnames(ordered) <- colnames(x)
	return(ordered)
}

order_cols <- function(x) {
	ordered <- colRanks(x = x, ties.method = "average", preserveShape = TRUE)
	rownames(ordered) <- rownames(x)
	colnames(ordered) <- colnames(x)
	return(ordered)
}


# freduce_with_param <- function(value, function_list, ...)
# {
#   if (length(function_list) == 1L)
#     function_list[[1L]](value, ...)
#   else 
#     Recall(function_list[[1L]](value, ...), function_list[-1L])
# }



freduce_with_param <- function (value, function_list,...) {
	k <- length(function_list)
	
	for (i in 1:k) {
		
		value <- function_list[[i]](value,...)
	}
	
	return(value)
}


identity_param <- function(x,...) {
	#identity function that returns first argument
	return(x)
}

theme_prism  <- function(){
	#theme that makes graphs look like output from graphpad prism 
	theme_classic()+
		theme(  #aspect.ratio = 1,
			text = element_text(face ="bold",colour = "black") ,
			line = element_line(size = 1.5,lineend = "square",colour = "black"),
			axis.text = element_text(face ="bold",color = "black"),
			axis.ticks = element_line(lineend = "square", color = "black"),
			plot.title = element_text(hjust = 0.5),
			#panel.border=element_rect(size=2),
			#strip.placement = "outside",
			strip.background=element_rect(size=3) 
			
		) 
}

