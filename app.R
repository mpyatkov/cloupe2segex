library(shiny)
library(shinyjs)
library(tidyverse)
library(markdown)
library(rmarkdown)

## cloupe df to segex df (k76 - dataset from segex)
cloupe2segex <- function(input_df, k76_df) {
    # reducing name for SEGEX
    replace_name <- function(name) {
        name %>% 
            str_replace(.,"Gt\\(ROSA\\)26Sor","Gt_ROSA_26Sor") %>% 
            str_replace(.,"^lnc_","nc_") %>%
            str_replace(.,"^ncRNA_","nc_") %>% 
            str_replace(.,"intra-as","intra") %>% 
            str_replace(.,"_chr","_c") %>% 
            str_replace(., "@", "_") %>% 
            str_replace(., "\\(\\+\\)", "p") %>% 
            str_replace(., "\\(-\\)", "m") %>% 
            str_replace(., "random", "") %>% 
            str_replace(., "_JH[0-9]+","_JH") %>% 
            str_replace(.,"LOC100861978_chr4_JH_m", "LOC100861978_chr4m") %>% 
            str_replace(., "Fam205a2_chr4_GL456350_m", "Fam205a2_chr4_GL_m") %>% 
            ifelse(str_detect(., "^nc_"), str_to_lower(.), .) 
    }
    
    segex_df <- input_df %>% 
        mutate(probe_id = replace_name(FeatureID)) %>% 
        left_join(k76_df, ., by = "probe_id") %>% 
        mutate(ratio = 2^(.[[5]]),
               fc = ifelse(ratio < 1, -1/ratio, ratio),
               intensity_1 = .[[7]],
               intensity_2 = .[[4]],
               pvalue_1 = .[[6]]) %>% 
        select(probe_id, ratio, fc, intensity_1, intensity_2, pvalue_1) %>% 
        
        replace_na(list(ratio = 1, fc = 1, intensity_1 = 0, intensity_2 = 0, pvalue_1 = 1)) %>% 
        as.data.frame(.) %>% 
        format(., digits = 8) 
}

convert_cloupe_to_segex <- function(old,new,all_76k){
    #message <- str_interp("${old} ===> ${new} conversion")
    #print(message)
    input <- read_csv(old, col_names = T, show_col_types = FALSE) %>% 
        cloupe2segex(., all_76k) %>% 
        #slice(1:1) %>% 
        write_tsv(new, col_names = T)
}

ui <- fluidPage(
    useShinyjs(),
    titlePanel("cloupe to SEGEX converter"),
    
    sidebarLayout(
        sidebarPanel(
            
            includeMarkdown("leader.md"),
            
            fileInput("leader", label="Description"),
            
            includeMarkdown("cloupe.md"),
            
            fileInput("csvs", label="Cloupe file(s)", multiple = TRUE),
            
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

        ## field should not be empty
        shiny::req(input$csvs)
        shiny::req(input$leader)

        ## TODO add dataframe checking here
        ## add package validate
        ## https://cran.r-project.org/web/packages/validate/vignettes/cookbook.html
        # if (2 < 1) {
        #     output$error <- renderText({"Not correct format for output. Please provide correct file"})
        #     return()
        # }

        main_dir <- getwd()
        
        ## read all probe_id which are presented in 76K platform
        all_76k <- read_tsv("./data/probeid_76k.csv", col_names = T, show_col_types = FALSE)

        ## load leader filename
        leader_fname <- input$leader$datapath
        
        ## read leader file
        leader_file <- read_csv(leader_fname, col_names = T) %>% 
            select(current = 1, cond2 = 2, cond1 = 3, feature = 4) %>% 
            # mutate(number = row_number()) %>%
            mutate(number = str_extract(current,"(^[:digit:]+)(?=_)")) %>% 
            rowwise() %>% 
            mutate(new_name = str_interp("${number}_scLoupe_${cond2}_vs_${cond1}_DiffExp_${feature}.tsv"),
                   old_name = str_interp("${current}.csv")) %>% 
            select(-number) %>% 
            select(old_name, new_name)
        
        ## get list of filenames
        input_files <- input$csvs %>% 
            select(old_name = name, datapath)
        
        leader_file <- inner_join(leader_file, input_files, by = "old_name") %>% 
            select(-old_name) %>% 
            select(old_name = datapath, new_name)
   
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
                         
                         walk2(leader_file$old_name, leader_file$new_name, function (x,y) {
                             incProgress(1/n)    
                             convert_cloupe_to_segex(x,y,all_76k)
                             }
                         )
                         
                         setProgress(detail = "Zipping files...")
                         ## zip files
                         zip("data.zip", leader_file$new_name)
                         #combined_name <- paste0(zipdir,"/test1.zip")
                         
                     }
        )
        
        resulted_zip(paste0(zipdir,"/data.zip"))
        #print(combined_name)
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