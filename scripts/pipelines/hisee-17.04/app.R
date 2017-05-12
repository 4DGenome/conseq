#==================================================================================================
# Created on: 2016-04-04
# Usage: in R: shinyApp(ui, server)
# Author: Javier Quilez (GitHub: jaquol)
# Goal: launches hisee-XX.XX --a visualisation tool for the Hi-C data generated with the hic-XX.XX pipeline
#==================================================================================================

CONSEQ <- '/users/GR/mb/jquilez/projects/conseq'



#==================================================================================================
# CONFIGURATION
#==================================================================================================

# Load R packages 
library("htmltools")
require('gtools')
library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library("RSQLite")

# Load functions

source("functions.r")

# Paths
DATA <- paste0(CONSEQ, '/data/hic/samples')
metadata <- paste0(CONSEQ, '/metadata/metadata.db'

# Samples
samples <- mixedsort(dir(DATA))

# Types of reads
reads_types <- c("Mapped" = "mapped",
					"Filtered" = "filtered",
					 "Dangling-end" = "dangling_end",
					 "Self-circle end" = "self_circle")

# Load metadata
db <- src_sqlite(metadata)
input_metadata <- tbl(db, "input_metadata")
quality_control_raw_reads <- tbl(db, 'quality_control_raw_reads')
hic <- tbl(db, "hic")
jobs <- tbl(db, "jobs")

# Filters

filters <- list(self_circle = 1,
                dangling_end = 2,
                error = 4,
                extra_dangling_end = 8,
                too_close_from_RES = 16,
                too_short = 32,
                too_large = 64,
                over_represented = 128,
                duplicated = 256,
                random_breaks = 512,
                trans = 1024)

# Color palettes

rwb_pal <- colorRampPalette(c("red", "white", "blue"))(100)
viridis_pal <- viridis(100)
reds_pal <- colorRampPalette(c("white", "red"))(100)

palettes <- list(viriris = viridis_pal,
                 reds = reds_pal,
                 red_white_blue = rwb_pal)

#==================================================================================================
# USER INTERFACE (UI)
#==================================================================================================

# Header
header <- dashboardHeader(title = "hisee-17.04 (4DGenome visualisation tool)")


# Side bar
sidebar <-  dashboardSidebar(
	sidebarMenu(
		menuItem("Sample information", tabName = "sample_info"),
		menuItem("Sample summary", tabName = "sample_summary"),
#		menuItem("Coverage", tabName = "coverage"),
		menuItem("Interaction matrices", tabName = "interaction_matrices"),
#		menuItem("Replicability", tabName = "replicability"),
#		menuItem("Metrics", tabName = "metrics"),
#		menuItem("Project progress", tabName = "project_progress"),
		menuItem("Performance", tabName = "performance")	
	)
)


# Body
body <- dashboardBody(

    includeCSS("www/styles.css"),
    
	tabItems(

		# Sample summary
		tabItem(tabName = "sample_summary",

			# Select sample
			fluidRow(column(3, selectInput("sample", label = h2("Select sample"), choices = samples, selected = "TE_S_T"))),
     			
			# last time sample was run on the pipeline
			fluidRow(column(12, h2("Last run on"))),
			fluidRow(column(12, box(tableOutput("last_run_on"), width = 12))),

			# Summary of input metadata
			fluidRow(column(12, h2("Input metadata"))),
			fluidRow(column(12, box(tableOutput("input_metadata"), width = 12))),

			# Pre-mapping quality control
			fluidRow(column(12, h2("Pre-mapping quality control"))),
			# FastQC reports of raw reads, and number of raw/processed reads
			fluidRow(column(12, box(title = "FastQC flags", tableOutput("fastqc_flags"), width = 12))),
			# Quality metrics and plots
			fluidRow(column(12, box(title = "Quality metrics and plots", width = 12,
						h4("Read1"), imageOutput("read1"),
						h4("Read2"), imageOutput("read2")))),

			# Post-mapping Statistics
			fluidRow(column(12, h2("Post-mapping statistics"))),
			fluidRow(column(12, box(title = "Mapping efficiency", width = 6, imageOutput("proportion_mapped_reads")),
								box(title = "Dangling-ends", width = 6, imageOutput("distribution_dangling_ends_lengths")))),
			fluidRow(column(12, box(title = "Decay of interactions with distance", width = 6, imageOutput("decay_interaction_counts_genomic_distance")),
								box(title = "Interaction matrix of mapped reads", width = 6, imageOutput("matrix_mapped")))),

			# Post-filtering statistics
			fluidRow(column(12, h2("Post-filtering statistics"))),
			fluidRow(column(12, box(title = "Summary of excluded reads", width = 8, tableOutput("excluded_reads")))),
			fluidRow(column(12, box(title = "Interaction matrix of filtered reads", width = 12, imageOutput("matrix_filtered", height = 750)))),

			# Overall sequencing coverage distribution
			fluidRow(column(9, box(title = "1-Mb coverage distribution", plotOutput("overall_coverage"), width = 12)),
					column(3, checkboxGroupInput("overall_coverage_read_type", label = h3("Select read type"), choices = reads_types, selected =reads_types[1]),
						h3("Display options"),
						sliderInput("overall_coverage_xlim", label = h5("X-axis range"), min = 0, max = 1000, value = 100),
						sliderInput("overall_coverage_ylim", label = h5("Y-axis range"), min = 0, max = 1500, value = 150))),

			# Scatter plot of coverage
			fluidRow(column(9, box(title = "Scatter plot of 1-Mb coverages", plotOutput("scatter_coverage"), width = 12)),
					column(3, checkboxGroupInput("scatter_coverage_read_type", label = h3("Selec read type"), choices = reads_types, selected = c(reads_types[1], reads_types[2])),
						h3("Display options"),
						sliderInput("scatter_coverage_xlim", label = h5("X-axis range"), min = 0, max = 1500, value = 150),
						sliderInput("scatter_coverage_ylim", label = h5("Y-axis range"), min = 0, max = 1500, value = 150))),

			# Chromosome plots of coverage
			fluidRow(column(9, box(title = "Chromosome plots of coverage", plotOutput("chromosome_coverage", height = 1250), width = 12, height = 1500)),
					column(3, h3("Display options"), sliderInput("chromosome_coverage_ylim", label = h5("Y-axis range"), min = 0, max = 500, value = 150)))

		),

		# Performance
		tabItem(tabName = "performance"),

        # raw contacts maps
		tabItem(tabName = "sample_info",

                verticalLayout(
                     
                    wellPanel(
                        
                        selectInput("meta_columns",
                                    label = "Select columns to show",
                                    choices = collect(input_metadata) %>% names,
                                    multiple = TRUE,
                                    selected = c("SAMPLE_ID", "SAMPLE_NAME", "SPECIES", "CELL_TYPE", "PRE_TREATMENT", "TREATMENT", "TREATMENT_TIME", "NOTES")),
                        
                        width = 12),
                    
                    DT::dataTableOutput("meta")
                    
                )
                
                
                ),
        
        # raw contacts maps
		tabItem(tabName = "interaction_matrices",

                # Input
                
                fluidRow(
                    
                    column(4,
                           
                           box(title = "Input info", width = NULL,
                               
                               # Select sample
                
                               fluidRow(column(6, selectInput("sample1",
                                                              label = "Select sample(s) [A]",
                                                              choices = samples,
                                                              multiple = TRUE)),
                                        column(6, selectInput("sample2",
                                                              label = "Select sample(s) to compare [B]",
                                                              choices = samples,
                                                              multiple = TRUE))),

                               # Select chromosome and resolution
                    
                               fluidRow(column(6, textInput("chrom", "Chromosome", "chr6")),
                                        column(6, numericInput("resolution", "Resolution (Kbp)",
                                                               500, 10, 5000, 10))),

                               # Buttons to show filters and fetch data
                    
                               fluidRow(column(4,
                                               align = "left",
                                               actionButton("show_filters", "Active filters")),
                                        column(4,
                                               align = "right",
                                               actionButton("fetch", "Fetch data"),
                                               offset = 4)),

                               # Filters

                               conditionalPanel(
                                   
                                   condition = "input.show_filters % 2 == 1",
                                   
                                   br(),
                                   
                                   fluidRow(column(6, checkboxGroupInput("filterex",
                                                                         label = "Excluding filters",
                                                                         choices = filters,
                                                                         selected = c(1, 2, 4, 8, 256, 512, 1024))),
                                            column(6, checkboxGroupInput("filterin",
                                                                         label = "Including filters",
                                                                         choices = filters,
                                                                         selected = 0))),
                               
                                   fluidRow(column(6, checkboxInput(inputId = "onlymulti",
                                                                    label = strong("Show only multicontacts"),
                                                                    value = FALSE))))
                               ),

                           box(title = "Download whole chromosome", width = NULL,

                               # Download whole chromosome
                    
                               #fluidRow(column(12, uiOutput("download_whole_chromosome"))),

                               fluidRow(column(4, align = "left", uiOutput("button_down_raw")),
                                        column(4, align = "center",  uiOutput("button_down_nor")),
                                        column(4, align = "right",  uiOutput("button_down_cor")))
                              ),

                           box(title = "Select region", width = NULL,

                               # Select region

                               #fluidRow(column(12, uiOutput("select_region"))),
                               
                               #fluidRow(column(12, uiOutput("range"))),

                               fluidRow(column(12, uiOutput("range2"))),

                               fluidRow(column(12, uiOutput("real_plot"))),

                               fluidRow(column(12, uiOutput("palette"))),

                               fluidRow(column(4, align = "right", uiOutput("update_plot"),
                                               offset = 8))
                              ),

                           box(title = "Download selected region", width = NULL,
                               
                               # fluidRow(column(12, uiOutput("download_region"))),
                               
                               fluidRow(column(4, align = "left",
                                               uiOutput("button_down_raw_region")),
                                        column(4, align = "center",
                                               uiOutput("button_down_nor_region")),
                                        column(4, align = "right",
                                               uiOutput("button_down_cor_region")))
                              )
                          ),

                    column(8,
                          
                           # Output
                           tabsetPanel(
                               
                               tabPanel("Raw",
                                        #textOutput("debug1"),
                                        #textOutput("debug3"),
                                        plotOutput("plot_raw")),
                               tabPanel("Normalized", plotOutput("plot_nor")),
                               tabPanel("Correlation", plotOutput("plot_cor")),
                               type = "pill")

                           )
                    )
                )
        )
    )

ui <- dashboardPage(header, sidebar, body)



#==================================================================================================
# SERVER
#==================================================================================================

# Plotting parameters
read_type_to_color <- list('mapped' = rgb(190/255, 190/255, 190/255, 0.5),
							'filtered' = rgb(0, 128/255, 1, 0.5),
							'dangling_end' = rgb(1, 102/255, 102/255, 0.5),
							'self_circle' = rgb(1, 102/255, 178/255, 0.5))

# Load existing images/plots
call_renderImage <- function(name, my_session = session, my_input = input){
	renderImage({
        width  <- my_session$clientData[[paste0("output_", name, "_width")]]
        height <- my_session$clientData[[paste0("output_", name, "_height")]]
		filename <- normalizePath(Sys.glob(file.path(DATA, my_input$sample, "*", "*", "*",
			paste0("*", name, "*.png"))))
		print(filename)
		list(src = filename, alt = "image number", width = width, height = height)
	}, deleteFile = FALSE)
}

server <- function(input, output, session, plot) {

	# Last run on
	output$last_run_on <- renderTable({
		df <- data.frame(t(data.frame(filter(hic, SAMPLE_ID == input$sample))))
		names(df) <- c("Value")
	df$field <- rownames(df)
	my_order <- c("SAMPLE_ID", "LAST_RUN_DATE")
	df[match(my_order, df$field), 1:2]
	})

	# Summary of input metadata
	output$input_metadata <- renderTable({
		df <- data.frame(t(data.frame(filter(input_metadata, SAMPLE_ID == input$sample))))
		names(df) <- c("Value")
	df$field <- rownames(df)
	my_order <- c("Timestamp", "SAMPLE_ID", "CELL_TYPE", "PRE_TREATMENT", "PRE_TREATMENT_TIME", "TREATMENT", "TREATMENT_TIME",
					"CONTROL", "EXPERIMENT_ID", "HIC", "SEQUENCING_PRE_LIBRARY", "SEQUENCING_LIBRARY", "SEQUENCING_CORE",
					"SEND_FOR_SEQUENCING_ON", "USER", "SAMPLE_NAME", "EXPERIMENT", "SPECIES", "DNA_CONCENTRATION", "SEQUENCING_INDEX",
					"ILLUMINA_MACHINE", "APPLICATION", "READ_LENGTH", "RESTRICTION_ENZYME", "SEQUENCING_TYPE", "NOTES")
	df[match(my_order, df$field), 1:2]
	})

	# FastQC reports of raw reads
	output$fastqc_flags <- renderTable({
		read1 <- data.frame(filter(select(quality_control_raw_reads, SAMPLE_ID, starts_with("READ1_")), SAMPLE_ID == input$sample))
		read2 <- data.frame(filter(select(quality_control_raw_reads, SAMPLE_ID, starts_with("READ2_")), SAMPLE_ID == input$sample))
		names(read1) <- gsub('READ1_', '', names(read1))
		names(read2) <- gsub('READ2_', '', names(read2))
		df <- data.frame(t(read1), t(read2))
		names(df) <- c("READ1", "READ2")
		df <- df[-c(1), ]
		df
	})

	# Quality metrics
	output$quality_metrics <- renderTable({
		df1 <- data.frame(hic %>%
			filter(SAMPLE_ID == input$sample) %>%
			select(starts_with("PERCENTAGE_")) %>%
			select(contains("READ1"))
			)
		df2 <- data.frame(hic %>%
			filter(SAMPLE_ID == input$sample) %>%
			select(starts_with("PERCENTAGE_")) %>%
			select(contains("READ2"))
			)
		names(df1) <- gsub('_READ1', '', names(df1))
		names(df2) <- gsub('_READ2', '', names(df2))
		df <- data.frame(t(df1), t(df2))
		names(df) <- c("READ1", "READ2")
		df		
	})

	# Load existing images/plots
	output$read1 <- call_renderImage("read1", session, input)
	output$read2 <- call_renderImage("read2", session, input)
	output$proportion_mapped_reads <- call_renderImage("proportion_mapped_reads", session, input)
	output$distribution_dangling_ends_lengths <- call_renderImage("distribution_dangling_ends_lengths", session, input)
	output$decay_interaction_counts_genomic_distance <- call_renderImage("decay_interaction_counts_genomic_distance", session, input)
	output$matrix_mapped <- call_renderImage("matrix_mapped", session, input)
	output$matrix_filtered <- call_renderImage("matrix_filtered", session, input)

	# Summary of exclude reads
	output$excluded_reads <- renderTable({
		df <- data.frame(t(data.frame(hic %>%
				filter(SAMPLE_ID == input$sample) %>%
				select(starts_with("EXCLUDED_"))
				)))
		names(df) <- 'values'
	n <- c()
	f <- c()
	for (i in df$values) {
		nc <- strsplit(i, ";")
		n <- c(n, as.integer((nc[[1]][1])))
		f <- c(f, as.numeric((nc[[1]][2])))
	}
	df$N_EXCLUDED_READS <- n
	df$FRACTION_EXCLUDED_READS <- f
	df$IS_APPLIED <- 'no'
	df["EXCLUDED_SELF_CIRCLE", "IS_APPLIED"] <- 'yes'
	df["EXCLUDED_DANGLING_END", "IS_APPLIED"] <- 'yes'
	df["EXCLUDED_EXTRA_DANGLING_END", "IS_APPLIED"] <- 'yes'
	df["EXCLUDED_ERROR", "IS_APPLIED"] <- 'yes'
	df["EXCLUDED_DUPLICATED", "IS_APPLIED"] <- 'yes'
	df["EXCLUDED_RANDOM_BREAKS", "IS_APPLIED"] <- 'yes'	
	df[, 2:4]
	})

	# Overall coverage distribution for different types of reads
	output$overall_coverage <- renderPlot({
		a = 0
		for (r in input$overall_coverage_read_type) {
			infile <- Sys.glob(file.path(DATA, input$sample, "*/*/genomic_coverages", paste0("*", r, "*.bed")))
			df <- read.delim(infile)
			bins <- seq(0, 2e3)
			my_colors = c()
			hist(df$read_counts/1e3,
					col = read_type_to_color[[r]],
					add = ifelse(a == 0, F, T),
					breaks = bins,
					border = F,
					xlim = c(0, input$overall_coverage_xlim),
					ylim = c(0, input$overall_coverage_ylim),
					main = input$sample,
					xlab = "No. reads x 1,000 / Mb",
					ylab = "Number of 1-Mb bins")
			my_colors = c(my_colors, read_type_to_color[[r]])
			a = a + 1
		}
		legend(input$overall_coverage_xlim * 0.8, input$overall_coverage_ylim * 0.8,
				names(read_type_to_color),
				fill = unlist(read_type_to_color, use.names=F),
				bty = "n")
	})

	# 2-way scatter plot of coverage (e.g. mapped vs filtered reads)
	output$scatter_coverage <- renderPlot({
		infile1 <- Sys.glob(file.path(DATA, input$sample, "*/*/genomic_coverages", paste0("*", input$scatter_coverage_read_type[[1]], "*.bed")))
		infile2 <- Sys.glob(file.path(DATA, input$sample, "*/*/genomic_coverages", paste0("*", input$scatter_coverage_read_type[[2]], "*.bed")))
		df1 <- read.delim(infile1)
		df2 <- read.delim(infile2)
		df1$read_counts <- df1$read_counts / 1e3
		df2$read_counts <- df2$read_counts / 1e3
		label1 <- input$scatter_coverage_read_type[[1]]
		label2 <- input$scatter_coverage_read_type[[2]]
		dat <- merge(df1, df2, by = names(df1)[1:3])
		names(dat) <- c("chrom", "start", "end", label1, label2)
		ggplot(dat, aes_string(x = label1, y = label2)) + 
				geom_point(shape = 1) +	
				geom_smooth(method = lm) +
				xlab(paste(label1, "(No. reads x1,000 / Mb)")) + xlim(0, input$scatter_coverage_xlim) +
				ylab(paste(label2, "(No. reads x1,000 / Mb)")) + ylim(0, input$scatter_coverage_ylim) +
				ggtitle(input$sample)
	})

	# Chromosome plots of coverage
	output$chromosome_coverage <- renderPlot({
		df_combined <- data.frame()
		for (r in reads_types) {
			infile <- Sys.glob(file.path(DATA, input$sample, "*/*/genomic_coverages", paste0("*", r, "*.bed")))
			df <- read.delim(infile)
			df$read_type <- r
			df_combined <- rbind(df_combined, df)
		}
		#df_combined$start_mb <- df_combined$start / 1e6
		df_combined$chrom <- factor(df_combined$chrom, levels = unique(as.vector(df_combined$chrom)))
		ggplot(df_combined, aes(x = start/1e6, y = read_counts/1e3, group = read_type, colour = read_type)) +
		geom_line() +
		facet_grid(chrom ~ .) +
		xlab("Position on chromosome (Mb)") +
		ylab("No. reads x1,000 / Mb") + ylim(0, input$chromosome_coverage_ylim) +
		ggtitle(input$sample)
	})

    # Reactive elements for contact maps

    chrom <- eventReactive(input$fetch, input$chrom)

    resolution <- eventReactive(input$fetch, input$resolution)
    
    samples_1 <- eventReactive(input$fetch, get_bam_path(input$sample1))

    samples_2 <- eventReactive(input$fetch, get_bam_path(input$sample2))

    filterin <- eventReactive(input$fetch, sum(as.numeric(input$filterin)))

    filterex <- eventReactive(input$fetch, sum(as.numeric(input$filterex)))

    onlymulti <- eventReactive(input$fetch, input$onlymulti)
    
    #region0 <- eventReactive(input$update_plot + 1, input$range * 1e3)

    region <- eventReactive(input$update_plot + 1,{

        if(is.null(input$range2)) return(NULL)

        strsplit(input$range2, "-") %>% unlist %>% as.numeric

        })

    m1 <- reactive(if(length(samples_2()) == 0) return(NULL) else return("[A]"))
    m2 <- reactive(if(length(samples_2()) == 0) return(NULL) else return("[B]"))

    transformation <- eventReactive(input$update_plot + 1, {
        
        if(input$realplot){
            
            return(logfinite)
            
        }else{
            
            return(function(x) log10(x + min(x[x>0], na.rm = TRUE)/2))
            
        }
        
    })

    the_palette <- eventReactive(input$update_plot + 1, palettes[[input$palette]])

    # fetch data
    
    dat_1 <- eventReactive(input$fetch, {

        selected_samples <- samples_1()
        
        progress <- shiny::Progress$new(max = 6 + length(selected_samples))
        on.exit(progress$close())
        progress$set(message = "Fetching data", value = 0)  
        
        wrap(selected_samples,
             chrom(),
             resolution() * 1e3,
             filterin(),
             filterex(),
             onlymulti(),
             progress)
        
    })

    dat_2 <- eventReactive(input$fetch, {

        selected_samples <- samples_2()

        if(length(selected_samples) == 0) return(NULL)
        
        progress <- shiny::Progress$new(max = 6 + length(selected_samples))
        on.exit(progress$close())
        progress$set(message = "Fetching data", value = 0)  
        
        wrap(selected_samples,
             chrom(),
             resolution() * 1e3,
             filterin(),
             filterex(),
             onlymulti(),
             progress)
        
    })

    # combine matrices
    
    dat_raw <- reactive({
        out <- dat_1()$raw
        if(!is.null(dat_2())) out <- combine_contacts(out, dat_2()$raw)
        out
    })

    dat_nor <- reactive({
        out <- dat_1()$nor
        if(!is.null(dat_2())) out <- combine_contacts(out, dat_2()$nor)
        out
    })

    dat_cor <- reactive({
        out <- dat_1()$cor
        if(!is.null(dat_2())) out <- combine_contacts(out, dat_2()$cor)
        out
    })
   
    # plots
    
    output$plot_raw <- renderPlot(plot_matrix(dat_raw(),
                                              region(),
                                              resolution() * 1e3,
                                              m1 = m1(), m2 = m2(),
                                              transformation = transformation(),
                                              color = the_palette(), sym = FALSE, trim = .01,
                                              unit_x_axis = 1e6,
                                              label_x_axis = "Genomic Position / Mbp"),
                                  width = 600, height = 600)
    
    output$plot_nor <- renderPlot(plot_matrix(dat_nor(),
                                              region(),
                                              resolution() * 1e3,
                                              m1 = m1(), m2 = m2(),
                                              transformation = transformation(),
                                              color = the_palette(), sym = FALSE, trim = .01,
                                              unit_x_axis = 1e6,
                                              label_x_axis = "Genomic Position / Mbp"),
                                  width = 600, height = 600)

    output$plot_cor <- renderPlot(plot_matrix(dat_cor(),
                                              region(),
                                              resolution() * 1e3,
                                              m1 = m1(), m2 = m2(),
                                              transformation = I,
                                              color = palettes$red_white_blue,
                                              sym = FALSE, trim = .01,
                                              unit_x_axis = 1e6,
                                              label_x_axis = "Genomic Position / Mbp"),
                                  width = 600, height = 600)
    
    # range of interest
    
    output$range <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        h2("Select region")
        
        sliderInput("range", "Region of interest (Kbp)",
                    min = 0, max = round(size / 1e3),
                    value = round(size / 1e3) / 2 + c(0, res * 10),
                    step = res, ticks = F)
        
    })

    output$range2 <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        h2("Select region")
        
        textInput("range2", "Region of interest (start-end)",
                    value = paste(round(size) / 2 + c(0, res * 1e3 * 10), collapse = "-"))
        
    })

    output$update_plot <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        actionButton("update_plot", "Update plot")
        
    })

    output$real_plot <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()
        checkboxInput("realplot", label = "Show holes", value = FALSE)
        
    })

    output$palette <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()
        selectInput("palette",
                    label = "Select color key",
                    choices = names(palettes),
                    multiple = FALSE)
        
    })

    output$button_down_raw <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        downloadButton('down_raw', 'Raw')
        
    })

    output$button_down_nor <- renderUI({
        
        size <- attr(dat_1(), "size")
        res <- resolution()
        
        downloadButton('down_nor', 'Normalized')
        
    })
    
    output$button_down_cor <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        downloadButton('down_cor', 'Correlation')
        
    })

    output$download_whole_chromosome <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        h2("Download whole chromosome")

    })

    output$download_region <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        h2("Download region")

    })

    output$select_region <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        h2(strong("Select region"))

    })

    output$button_down_raw_region <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        downloadButton('down_raw_region', 'Raw')
        
    })

    output$button_down_nor_region <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        downloadButton('down_nor_region', 'Normalized')
        
    })

    output$button_down_cor_region <- renderUI({

        size <- attr(dat_1(), "size")
        res <- resolution()

        downloadButton('down_cor_region', 'Correlation')
        
    })

    
    # metadata table

    output$meta <- DT::renderDataTable(collect(input_metadata)[, input$meta_columns] %>%
                                       filter(SAMPLE_ID %in% samples),
                                       rownames = FALSE,
                                       filter = "top")
    
    # download whole chromosome

    filename_raw <- reactive(make_download_filename(chrom(),
                                                    resolution()))
    content_raw <- function(filename) make_download_content(dat_raw(),
                                                            filename,
                                                            resolution())
    
    output$down_raw <- downloadHandler(filename = function() filename_raw(),
                                       content = function(filename) content_raw(filename))

    filename_nor <- reactive(make_download_filename(chrom(),
                                                    resolution()))
    content_nor <- function(filename) make_download_content(dat_nor(), filename, resolution())
    
    output$down_nor <- downloadHandler(filename = function() filename_nor(),
                                       content = function(filename) content_nor(filename))

    filename_cor <- reactive(make_download_filename(chrom(),
                                                    resolution()))
    content_cor <- function(filename) make_download_content(dat_cor(), filename, resolution())
    
    output$down_cor <- downloadHandler(filename = function() filename_cor(),
                                       content = function(filename) content_cor(filename))

    # download region

    filename_raw_region <- reactive(make_download_filename(chrom(),
                                                           resolution(),
                                                           region()))
    
    content_raw_region <- function(filename) make_download_content(dat_raw(),
                                                                   filename,
                                                                   resolution() * 1e3,
                                                                   region())
    
    output$down_raw_region <- downloadHandler(filename = function() filename_raw_region(),
                                              content = function(filename) content_raw_region(filename))

    filename_nor_region <- reactive(make_download_filename(chrom(),
                                                           resolution(),
                                                           region()))
    
    content_nor_region <- function(filename) make_download_content(dat_nor(),
                                                                   filename,
                                                                   resolution() * 1e3,
                                                                   region())
    
    output$down_nor_region <- downloadHandler(filename = function() filename_nor_region(),
                                              content = function(filename) content_nor_region(filename))

    filename_cor_region <- reactive(make_download_filename(chrom(),
                                                           resolution(),
                                                           region()))
    
    content_cor_region <- function(filename) make_download_content(dat_cor(),
                                                                   filename,
                                                                   resolution() * 1e3,
                                                                   region())
    
    output$down_cor_region <- downloadHandler(filename = function() filename_cor_region(),
                                              content = function(filename) content_cor_region(filename))

    
    # Output debug

    output$debug1 <- renderText({
        
        paste(samples_1(),
              chrom(),
              resolution() * 1e3,
              filterin(),
              filterex())
              
    })

    output$debug2 <- renderText({
        
        c(resolution(), region())
        
    })

    output$debug3 <- renderText({
        
        sum(dat_raw())
        
    })

    

}


shinyApp(ui, server)
