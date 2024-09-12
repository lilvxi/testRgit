# Load necessary libraries
library(shiny)
library(dplyr)
library(survival)
library(lme4)
library(GLMMadaptive)
library(ggplot2)

library(bslib)
library(DT)

#####
cvp<-function(xx)
{
card(
card_header(paste0('Variable', xx)),
layout_columns(
  numericInput(paste0("var",xx,"_mean"), paste0("Mean"), value = 10),
  numericInput(paste0("var",xx,"_min"),  paste0("Min"), value = 5),
  numericInput(paste0("var",xx,"_max"),  paste0("Max"), value = 18)
),
class='bg-success & text-white'
)
}
  
# Function to generate variable values
generate_variable <- function(var_type, var_mean, var_min, var_max) {
      if (var_type == "Continuous") {
        var <- rnorm(1, mean = var_mean, sd = (var_max - var_min) / 6)
        return(pmin(pmax(var, var_min), var_max))
      } else {
        return(rbinom(1, 1, 0.5))
      }
    }


# Define UI for the app

ui<-page_navbar(

title="Stepped Wedge Design Simulation",

sidebar=sidebar(
width=350,
open=NA,

actionButton("runSim", "Run Simulation", icon = icon("play")),


accordion(
open=c('Model Selection','Beta Coefficients'),
accordion_panel(
'Model Selection',
radioButtons("model_choice",NULL, choices = c("glimmix", "glmmadaptive")),
),


accordion_panel(
'Beta Coefficients',
numericInput("beta", "Overall Beta:", value = 1.5),
selectInput("num_vars", "Number of Variables:", choices = c(1, 2, 3)),
uiOutput('bcv3')
),
accordion_panel(
'Sample Size and Cluster Information',
      numericInput("samplesize", "Total Sample Size:", value = 200, min = 1, step = 50),
      numericInput("steps", "Number of Clusters (Steps):", value = 4, min = 1),
      numericInput("sim_num", "Number of Simulations:", value = 3, min = 1),
),
accordion_panel(
'Interval and Switch Parameters',
      numericInput("intv_length", "Interval Length (days):", value = 30, min = 1),
      numericInput("switch_interval", "Switch Interval (days):", value = 180, min = 1),
),
accordion_panel(
'Time Parameters',
      dateRangeInput("study_time", "Study Time Range:", 
                     start = "2021-09-01", end = "2025-02-28"),
      dateInput("switch_begin", "Switch Begin Time:", value = "2021-09-01"),
      dateInput("recruit_end", "Recruitment End Time:", value = "2024-02-28"),
),
accordion_panel(
'Distribution Parameters',
      numericInput("scale", "Scale:", value = 0.003),
      numericInput("shape", "Shape:", value = 0.441),
),

accordion_panel(
'Standard Errors',
      numericInput("sigma_c", "Standard Error for Error (sigma_c):", value = 0.3, min = 0),
      # Conditional input for sigma_tau, only show if glmmadaptive is selected
      conditionalPanel(
        condition = "input.model_choice == 'glmmadaptive'",
        numericInput("sigma_tau", "Standard Error for Tau-bti (sigma_tau):", value = 0.15, min = 0)
      ),
),
),

#textOutput("estimatedTime"),
#tags$small("Note: The estimated time includes the total running time for simulation and model fitting.")
), ####end sidebar



nav_panel("SimulationData",
uiOutput('logDetail1'),
uiOutput('logDetail2'),
uiOutput('logDetail3'),
uiOutput('simuResu')
),


nav_panel("Density Plot", 
plotOutput("densityPlot")),

nav_panel("Mean Beta and Power",
         uiOutput("simulationResults")
),

nav_panel("Optional Histogram", 
         checkboxInput("showHistogram", "Show Histogram", value = FALSE),
         conditionalPanel(
           condition = "input.showHistogram == true",
           plotOutput("histogramPlot")
         ),
),


#footer=
#tags$div(
##downloadButton("downloadData", "Download Simulation Data"),
#tags$h6("Notes on Simulation Process"),
#tags$ul(style='font-size:75%',
#tags$li("The simulation process involves generating time-to-event data for a stepped wedge design study based on the specified input parameters. The resulting data are then used to fit a mixed-effects model using either the glimmix or glmmadaptive approach, as chosen by the user."),
#tags$li("The simulation time displayed includes both the data generation and the model fitting processes."),
#)
#)

 ##end footer
	
)


server <- function(input, output, session) {

output$bcv3<-renderUI({
input$num_vars->nv

lapply(1:nv,function(xx)
{
tagList(
selectInput(paste0("var",xx,"_type"), paste("Variable",xx,"Type:"), choices = c("Continuous", "Binary")),
numericInput(paste0("beta",xx),paste0("Beta for Variable ",xx," :"), value = 0.5),
uiOutput(paste0('bcvv',xx))
)
}
)->inter.bcv

tagList(inter.bcv)
})

################
output$bcvv1<-renderUI({
if(input$var1_type=="Continuous") cvp(1) else NULL
})

output$bcvv2<-renderUI({
if(input$num_vars>=2 & input$var2_type=="Continuous") cvp(2) else NULL
})

output$bcvv3<-renderUI({
if(input$num_vars>=3 & input$var3_type=="Continuous") cvp(3) else NULL
})


########
observeEvent(input$runSim, {
    set.seed(123)


    # Initialize an empty data frame to store results
    results <- data.frame()
    
    # Function to generate variable values
    generate_variable <- function(var_type, var_mean, var_min, var_max) {
      if (var_type == "Continuous") {
        var <- rnorm(1, mean = var_mean, sd = (var_max - var_min) / 6)
        return(pmin(pmax(var, var_min), var_max))
      } else {
        return(rbinom(1, 1, 0.5))
      }
    }
	
	# Progress tracking
	start.sim_time <- Sys.time()
	
	progress <- Progress$new(session, min = 1, max = input$steps*input$sim_num)
    progress$set(message = paste('Simulation starts at:',format(start.sim_time,"%H:%M:%S"),'...in progress...'), value = 1)
    
    # Simulation loop
    for (sim in 1:input$sim_num) {
      for (cluster in 1:input$steps) {
        error <- rnorm(1, mean = 0, sd = input$sigma_c)
        tau_bti <- if (input$model_choice == "glmmadaptive") rnorm(1, mean = 0, sd = input$sigma_tau) else 0
        
        for (id in 1:input$samplesize) {
          # Generate variables based on the number of variables selected by the user
          var1 <- if (input$num_vars >= 1) {
            generate_variable(input$var1_type, input$var1_mean, input$var1_min, input$var1_max)
          } else {
            NA
          }
          
          if (input$num_vars >= 2) {
            var2 <- generate_variable(input$var2_type, input$var2_mean, input$var2_min, input$var2_max)
          } else {
            var2 <- NULL
          }
          
          if (input$num_vars == 3) {
            var3 <- generate_variable(input$var3_type, input$var3_mean, input$var3_min, input$var3_max)
          } else {
            var3 <- NULL
          }
          
          # Define study dates
          studybegin <- as.Date(input$study_time[1])
          recuitend <- as.Date(input$recruit_end)
          studyend <- as.Date(input$study_time[2])
          
          # Calculate entry/switch/drop date and switch time
          entertemp <- runif(1, min = 0, max = as.numeric(difftime(recuitend, studybegin, units = "days")))
          entrydate <- studybegin + as.integer(entertemp)
          
          switchdate <- as.Date(input$switch_begin) + as.integer(cluster * input$switch_interval)
          switchtime <- as.numeric(pmax(switchdate - entrydate + 1, 0))
          
          c <- rexp(1, rate = exp(input$shape)) * as.numeric(difftime(studyend, studybegin, units = "days"))
          dropdate <- entrydate + as.integer(c)
          
          # Simulate survival time
          logu <- -log(runif(1))
          linpre <- sum(c(if (!is.null(var1)) var1 * input$beta1 else 0,
                          if (!is.null(var2)) var2 * input$beta2 else 0,
                          if (!is.null(var3)) var3 * input$beta3 else 0)) + error
          linpre2 <- linpre + tau_bti
          scalee <- input$scale * exp(linpre)
          scalee2 <- input$scale * exp(linpre2)
          
          nom <- logu - scalee * (switchtime^input$shape) + scalee2 * exp(input$beta) * (switchtime^input$shape)
          denom <- scalee2 * exp(input$beta)
          
          if (logu < scalee * (switchtime^input$shape)) {
            t <- (logu / scalee)^(1 / input$shape)
          } else {
            t <- (nom / denom)^(1 / input$shape)
          }
          
          disclosuredate <- entrydate + as.integer(t)
          
          # Determine censoring status and observed end date
          if (disclosuredate > pmin(dropdate, studyend, entrydate + 365)) {
            censor <- 0
            obenddate <- pmin(dropdate, studyend, entrydate + 365)
          } else {
            censor <- 1
            obenddate <- disclosuredate
          }
          
          # Calculate survival time
          disclosuretime <- as.numeric(difftime(obenddate, entrydate, units = "days")) + 1
          
          # Store the results
          results <- rbind(results, data.frame(
            sim, cluster, id, 
            var1 = if (!is.null(var1)) var1 else NA, 
            var2 = if (!is.null(var2)) var2 else NA, 
            var3 = if (!is.null(var3)) var3 else NA, 
            studybegin, entrydate, switchdate, dropdate, obenddate, censor, disclosuretime
          ))
        }
		progress$inc(1,message=paste('Simmulating',scales::percent(progress$getValue()/(input$steps*input$sim_num),0.01)))
      }
     
    }
	
	progress$close()
	
	
	#########
	end.sim_time<-Sys.time()
	output$logDetail1<-renderUI({
	tagList(
	tags$p(tags$em(tags$b('Simulation'),' starts at:'),start.sim_time,
	tags$em('ends at:'),end.sim_time,
	tags$em('takes:',sprintf('%4.3f',end.sim_time-start.sim_time),'s'))
	)
	})
	
    
    # Convert results to a data frame
    results <- as.data.frame(results)
	
	
	###
	start.wrangle.time<-Sys.time()
    progress2 <- Progress$new(session, min = 1, max = 10)
    progress2$set(message = paste('wrangling of simulation data starts at:',format(start.wrangle.time,"%H:%M:%S"),'...in progress...'), value = 1)
	
	
    # Create interval data
    intvdata <- results %>%
      mutate(entryintv = as.integer(difftime(entrydate, studybegin, units = "days") / input$intv_length) + 1,
             switchintv = as.integer(difftime(switchdate, studybegin, units = "days") / input$intv_length) + 1,
             surveintv = as.integer(difftime(obenddate, studybegin, units = "days") / input$intv_length) + 1,
             survtimediscrete = surveintv - entryintv + 1)
			 
	####density Plot
	output$densityPlot<-renderPlot({
	 ggplot(intvdata, aes(x =survtimediscrete, fill = factor(censor))) +
          geom_density(alpha = 0.5) +
          labs(title = "Density Plot of Month by Outcome",
               x = "Interval",
               y = "Density",
               fill = "Outcome") +
          theme_minimal()
    })
    
    # Event history data processing
    eventhistorydata <- intvdata %>% 
      rowwise() %>% 
      do({
        interval_data <- lapply(1:.$survtimediscrete, function(interval) {
          event <- ifelse(interval == .$survtimediscrete, .$censor, 0)
          trt <- ifelse(interval < (.$switchintv - .$entryintv + 1), 0, 1)
          data.frame(
            sim = .$sim, 
            cluster = .$cluster, 
            id = .$id, 
            var1 = if (!is.null(var1)) .$var1 else NA,
            var2 = if (!is.null(var2)) .$var2 else NA,
            var3 = if (!is.null(var3)) .$var3 else NA,
            disclosuretime = .$disclosuretime, 
            interval = interval, 
            event = event, 
            trt = trt
          )
        })
        bind_rows(interval_data)
      }) %>% 
      ungroup()
    
    eventhistorydata <- eventhistorydata %>% 
                        arrange(sim, cluster)
    
	
	end.wrangle.time<-Sys.time()
	progress2$close()
	
	####
	output$logDetail2<-renderUI({
	tagList(
	tags$p(tags$em(tags$b('Wrangling'),' starts at:'),start.wrangle.time,
	tags$em('ends at:'),end.wrangle.time,
	tags$em('takes:',sprintf('%4.3f',end.wrangle.time-start.wrangle.time),'s'))
	)
	})
	
	
	
	###
	start.model.time<-Sys.time();
	progress3 <- Progress$new(session, min = 1, max = 10)
    progress3$set(message = paste('Modelling starts at:',format(start.model.time,"%H:%M:%S"),'...in progress...'), value = 1)
	
	
    # Check for multicollinearity
    if (input$num_vars > 1) {
      selected_vars <- eventhistorydata %>% select(var1, var2, var3) %>% select_if(~ sum(!is.na(.)) > 0)
      corr_matrix <- cor(selected_vars, use = "complete.obs")
      
      if (any(abs(corr_matrix[lower.tri(corr_matrix)]) > 0.9)) {
        showModal(modalDialog(
          title = "Warning: Multicollinearity Detected",
          "High correlation between variables detected. Consider reducing the number of variables or modifying them.",
          easyClose = TRUE
        ))
      }
    }
	
    # Prepare model fitting
    variable_names <- c("var1", "var2", "var3")[1:input$num_vars]
    variable_names <- variable_names[!sapply(variable_names, function(x) all(is.na(eventhistorydata[[x]])))]
    
    # Model fitting based on user selection
    betas <- numeric(input$sim_num)
    pvalues <- numeric(input$sim_num)
    
    for (sim in 1:input$sim_num) {
      # Create formula dynamically
      formula <- as.formula(paste("event ~", 
                                  paste(c(variable_names, "trt", "interval"), collapse = " + "), 
                                  "+ (1 | cluster)"))
      
      if (input$model_choice == "glimmix") {
        glimmix_model <- tryCatch({
          glmer(formula, 
                data = eventhistorydata %>% filter(sim == !!sim), 
                family = binomial(link = "cloglog"), nAGQ = 1)
        }, error = function(e) e)
        
        if (inherits(glimmix_model, "error")) {
          showModal(modalDialog(
            title = "Model Fitting Error",
            paste("Simulation", sim, "encountered an error during fitting:", glimmix_model$message),
            easyClose = TRUE
          ))
          next
        }
        
        if (!is.null(glimmix_model@optinfo$conv$lme4$messages)) {
          showModal(modalDialog(
            title = "Model Fitting Warning",
            paste("Simulation", sim, "encountered a fitting issue:", glimmix_model@optinfo$conv$lme4$messages),
            easyClose = TRUE
          ))
          next
        }
        
        betas[sim] <- fixef(glimmix_model)["trt"]
        pvalues[sim] <- summary(glimmix_model)$coefficients["trt", "Pr(>|z|)"]
        
      } else if (input$model_choice == "glmmadaptive") {
        glmmadaptive_model <- mixed_model(
          fixed = formula, 
          random = ~ trt | cluster, 
          family = binomial(link = "cloglog"), 
          data = eventhistorydata %>% filter(sim == !!sim),
          control = list(nAGQ = 1)
        )
        
        betas[sim] <- fixef(glmmadaptive_model)["trt"]
        summary_model <- summary(glmmadaptive_model)
        coefficients_table <- summary_model$coef_table
        p_value_trt <- coefficients_table["trt", "p-value"] 
        
        pvalues[sim] <- p_value_trt
      }
    }
	
	####
	end.model.time<-Sys.time()
	progress3$close()
	
	####
	output$logDetail3<-renderUI({
	tagList(
	tags$p(tags$em(tags$b('Modelling'),' starts at:'),start.model.time,
	tags$em('ends at:'),end.model.time,
	tags$em('takes:',sprintf('%4.3f',end.model.time-start.model.time),'s'))
	)
	})
    
    # Calculate the mean of the beta coefficients and estimated power
    mean_beta <- mean(betas)
    power <- mean(pvalues < 0.05)
	end_time<-Sys.time()
    
    # Output the mean beta and estimated power
	output$simulationResults <- renderUI({
	  tagList(
	  tags$h6("Simulation Results:"),
	  tags$p("Mean of Beta Coefficients:",tags$em(mean_beta)),
	  tags$p("Estimated Power:",tags$em(power)),
	  #tags$p('beta:',paste(sprintf("%4.3f",betas),collapse=", ")),
	  #tags$p(tags$em(tags$b('P')),'values:',paste(sprintf("%4.3e",pvalues),collapse=", ")),
	  )
    })
    
    # Convert betas to a data frame for plotting
    betas_df <- data.frame(beta = betas)
    
    # Render histogram plot
    output$histogramPlot <- renderPlot({
      ggplot(betas_df, aes(x = beta)) +
        geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
        geom_vline(aes(xintercept = mean_beta), color = "red", linetype = "dashed", size = 1) +
        labs(
          title = paste("Histogram of Treatment Effect Coefficients (Beta =", input$beta, ")"),
          x = "Beta (Treatment Effect)",
          y = "Frequency"
        ) +
        theme_minimal()
    })
    
    # Make results available for download
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("simulation_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {

        write.csv(results, file, row.names = FALSE)
      }
    )
	
	################
    output$simuResu<-renderUI({
    tagList(
    layout_columns(col_widths=c(3),downloadButton("downloadData", "Download Simulation Data")),
    DTOutput('simdata')
    )
    })
    
    output$simdata<-renderDT({
	#nocols<-setdiff(c('var1','var2','var3'),paste0('var',1:input$num_vars))
	nocols<-''
	colss<-setdiff(names(results),nocols)
    datatable(results[,colss],options=list(dom='tp'))
    })
    



	
    
  })

############




}




# Run the application 
shinyApp(ui = ui, server = server)
