options(shiny.maxRequestSize = 20*1024^2)

library(shiny)
library(shinyjs)

#library(tidyverse)
library(purrr)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

library(markdown)
library(rmarkdown)

## cloupe df to segex df (k76 - dataset from segex)
cloupe2segex <- function(input_df, df_76k, not_in_76k) {

    ## from 8 cloupe columns to 6 SEGEX columns
    conversion <- function(df){
        df %>% mutate(ratio = 2^(.[[5]]),
                      fc = ifelse(ratio < 1, -1/ratio, ratio),
                      intensity_1 = .[[7]],
                      intensity_2 = .[[4]],
                      pvalue_1 = .[[6]]) %>% 
            select(probe_id, ratio, fc, intensity_1, intensity_2, pvalue_1) %>% 
            replace_na(list(ratio = 1, fc = 1, intensity_1 = 0, intensity_2 = 0, pvalue_1 = 1)) %>% 
            as.data.frame(.) %>% 
            format(., digits = 8) 
    }
    
    segex_df <- input_df %>% 
        left_join(df_76k, ., by = "FeatureName") %>% 
        conversion(.)

    
    segex_df_appendix <- input_df %>% 
        left_join(not_in_76k, ., by = "FeatureName") %>% 
        conversion(.)
    
    list(tosegex = segex_df, appendix = segex_df_appendix)
}

convert_cloupe_to_segex <- function(old, new, order, all_76k, not_in_76k){
    #message <- str_interp("${old} ===> ${new} conversion")
    #print(message)
    input <- read_csv(old, col_names = T, show_col_types = FALSE) 
    
    ## if order = 1 (Cond1/Cond2 situation) - we have to change order of columns in input file, to make it easier to import with segex
    ## by default segex import data as Cond2/Cond1, so we have to preserve the same order of columns
    ## loupe browser export data in different order (Cond1/Cond2 or Cond2/Cond1) with any attention to order
    if (order == 1) {
      input <- input %>% 
        ## change order of columns
        relocate(c(3,4,5),.after = c(6,7,8))
    }
    
    input <- cloupe2segex(input, all_76k, not_in_76k) 
    
    write_tsv(input$tosegex, new, col_names = T)
    write_tsv(input$appendix, paste0("APPENDIX_",new), col_names = T)
}

is_correct_leaderfile <- function(df_path) {
    ## file should contain at least 4 columns
    if (ncol(read_csv(df_path, col_names = T)) != 5) return(FALSE)
    TRUE
}

# 
is_correct_datafiles <- function(datafiles) {
    # datafiles = 2 columns - filename(1), path(2)

    checkone <- function(df_path) {
        if (ncol(read_csv(df_path, col_names = T, n_max = 10)) != 8) return(FALSE)
        TRUE
    }
    
    ## any problems with files?
    datafiles_vector_of_states <- sapply(datafiles[,2],checkone)
    datafiles_state <- all(datafiles_vector_of_states)
    
    ## if state is not TRUE (we have a problem in one or more files)
    if (!datafiles_state) {
        
        ## detect in which files and print them
        fnames <- which(datafiles_vector_of_states == FALSE)
        fnames <- datafiles[,1][fnames]
        fnames <- paste0(fnames, collapse = "\n")
        return(list(status = datafiles_state, names = fnames))
    }

    return(list(status = TRUE, names = NULL))
}

ui <- fluidPage(
    useShinyjs(),
    titlePanel("cloupe to SEGEX converter"),
    
    ## change color for error messages
    tags$head(
        tags$style(HTML("
      #error {
        color: #8B0000;
        white-space:pre-wrap;
      }
    "))
    ),
    
    sidebarLayout(
        sidebarPanel(width=6,
            
            includeMarkdown("leader.md"),
            
            fileInput("leader", label="Setup file upload:"),
            
            includeMarkdown("cloupe.md"),
            
            fileInput("csvs", label="Upload of Cloupe file(s) to be converted:", multiple = TRUE),
            
            actionButton(inputId = "go", label = "Convert"),
            
            #downloadButton(outputId = "downloadData", "Download combined pdf")
            conditionalPanel(
                "false", # always hide the download button
                downloadButton(outputId = "downloadData", "Download combined pdf")
            )
            
    ),

    mainPanel(
        verbatimTextOutput("error")
    ))
)

server <- function(input, output, session) {
    
    onSessionEnded(function() {
        list.dirs(recursive=FALSE) %>% 
            keep(function(x) grepl("^./tmp", x)) %>% 
            unlink(recursive = TRUE)
    })
    
    resulted_zip <- reactiveVal()

    observeEvent(input$go,{
        
        # show error if files in not presented
        output$error <- renderText({
            validate(
                need(input$leader != '', 'Please provide file with the description'),
                need(input$csvs != '', 'Please provide file(s) exported from loupe browser')
            )
        })
        
        ## field should not be empty
        shiny::req(input$csvs)
        shiny::req(input$leader)

        ## -----
        
        ## get list of filenames
        input_files <- input$csvs %>% 
            select(old_name = name, datapath)

        ## load leader filename
        leader_fname <- input$leader$datapath
        
        ## TODO add dataframe checking here
        ## add package validate
        ## https://cran.r-project.org/web/packages/validate/vignettes/cookbook.html

        if (!is_correct_leaderfile(leader_fname)) {
            output$error <- renderText({"Description file format is not correct, should be 5 columns"})
            return()
        }
        
        checking_datafiles <- is_correct_datafiles(input_files)
        if (!checking_datafiles$status) {
            output$error <- renderText({paste0("The format of the following files is not correct (8 columns are required):\n", 
                                               checking_datafiles$names)})
            return()
        }    

        main_dir <- getwd()
        
        ## read all platforms
        all_mm9mm10 <- read_csv("./data/mm10_to_segex.csv", col_names = T, show_col_types = FALSE)
        
        ## read all probe_id which are presented in 76K platform
        all_76k <- all_mm9mm10 %>% 
          filter(in_segex == T) %>% 
          select(probe_id, FeatureName = gene_short_mm10)
        
        ## read all probe_id which is NOT presented in 76k (excluding ERCC genes)
        not_in_76k <- all_mm9mm10 %>% 
          filter(in_segex == F & !grepl("ERCC", probe_id)) %>% 
          select(probe_id, FeatureName = gene_short_mm10) %>% 
          arrange(probe_id)
        
        ## read leader file
        leader_file <- read_csv(leader_fname, col_names = T) %>% 
            select(current = 1, cond2 = 2, cond1 = 3, feature = 4, order = 5) %>% # order = 0/1. 0 -> Cond2/Cond1 - normal, 1 -> Cond1/Cond2 - has to be flipped
            mutate(number = str_extract(current,"(^[:digit:]+)(?=_)")) %>% 
            rowwise() %>% 
            mutate(new_name = str_interp("${number}_scLoupe_${cond2}_vs_${cond1}_DiffExp_${feature}.tsv"),
                   old_name = str_interp("${current}.csv")) %>% 
            select(-number) %>% 
            select(old_name, new_name, order)

        ## just to check if all filenames in the leader file are correct
        list_in_file <- leader_file %>% select(old_name)
        
        list_intersection <- inner_join(leader_file, input_files, by = "old_name")
     
        if (nrow(list_in_file) != nrow(list_intersection)) {
            if(ncol(list_intersection) == 0) {
                output$error <- renderText({"There is nothing to process. The filenames in the file with description do not match the input filenames"})
                return()
            }
            else {
                output$error <- renderText({paste0("WARNING: partial processing. Only the following filenames were processed:\n",
                                                   paste0(list_intersection$old_name, collapse = "\n"))})
            }
        }
        
        leader_file <- inner_join(leader_file, input_files, by = "old_name") %>% 
            select(-old_name) %>% 
            select(old_name = datapath, new_name, order)
        
        ## create temporary directory
        zipdir <- tempfile("tmp", tmpdir = "./")
        if (!dir.exists(zipdir)) {
            dir.create(zipdir)  
        }
        setwd(zipdir)
        
        ## conversion
        withProgress(message = "Processing files", 
                     value = 0,
                     {
                       
                         n <- length(leader_file$new_name)

                         ## process all loupe files
                         pwalk(leader_file, function (old_name,new_name,order) {
                           incProgress(1/n)
                           ## create files for SEGEX deploy and appendix files
                           
                           convert_cloupe_to_segex(old_name, new_name, order, all_76k, not_in_76k)
                         })
                         
                         
                         ## combine all appendix files to one
                         ## add extra header for each filename
                         ## Header1\t\t\t\t\t\tHeader2\t\t\t\t\t\t
                         extraheader <- paste0(leader_file$new_name, "\t\t\t\t\t\t", collapse = "")
                         writeLines(extraheader, "COMBINED_APPENDICES.tsv")
                         
                         combined_file <- dplyr::bind_cols(lapply(paste0("APPENDIX_",leader_file$new_name), read_tsv, col_names = T),
                                                           .name_repair = "minimal") %>% 
                             write_tsv(file = "COMBINED_APPENDICES.tsv", col_names = T, append = T)

                         setProgress(detail = "Zipping files...")
                         ## zip files
                         ## files for segex and appendix files and combined files with appendices
                         fnames_for_archive <- c(leader_file$new_name, paste0("APPENDIX_",leader_file$new_name), "COMBINED_APPENDICES.tsv")
                         #zip("data.zip", leader_file$new_name)
                         zip("data.zip", fnames_for_archive)
                         #combined_name <- paste0(zipdir,"/test1.zip")
                     }
        )
        
        resulted_zip(paste0(zipdir,"/data.zip"))
        setwd(main_dir)
        
        runjs("$('#downloadData')[0].click();")

    })
    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("for_segex-", format(Sys.time(), "%a_%b_%d_%Y_%I-%M%p"), ".zip", sep="")
        },
        content = function(file) {
            file.copy(resulted_zip(), file)
        },
        contentType = "application/zip")
}

shinyApp(ui = ui, server = server)