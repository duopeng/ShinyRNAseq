#Shiny app for EcR RNA seq dataset (Catteruccia lab)
#This tool can perform Differential Expressed Gene (DEG) Analysis and GO Enrichemnt Analysis
#version 3.3 (2023-7-17), Dropped Entrez in the BiomaRt request to avoid an erroor 
#by Duo Peng

library(BiocManager) #needed for deployment at Shinapp.io
options(repos = BiocManager::repositories()) #needed for deployment at Shinapp.io
library(shiny)
library(DT)
library(htmltools)

##################################################
# load data files                                #
##################################################

# read in metadata
metadata <- read.table( "data/metadata.tab", header=T, sep = "\t")
metadata$Replicate <- as.factor( metadata$Replicate )
metadata$Timepoint <- as.factor( metadata$Timepoint )

#read in albi-to-gamb annotation
anno <- read.delim( "data/Gamb2Albi.txt_Albi2Gamb.txt.combined", header= FALSE, sep = "\t" )
anno <- as.matrix(anno)

#read in GO to Gene conversion table
load("data/GO2Gene_anno.RData") #variable loaded is a list named "GO"

##################################################
# Begin Shiny server code                        #
##################################################

#############################
# Shiny server: UI          #
#############################

parameter_tabs <- tagList(
    tags$style("#params { display:none; }"),
    tabsetPanel(id = "params",
                tabPanel("dsRNA Treatment",
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h5("Your contrast will be EcR vs GFP"),
                                  h5("log2FoldChange = EcR/GFP"),

                         ),
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to select tissue(s)"),
                                  checkboxInput(inputId="tr_Midgut", label = "Midgut", value = TRUE ),
                                  checkboxInput(inputId="tr_FatBody", label = "FatBody", value = FALSE ),
                                  checkboxInput(inputId="tr_Ovary", label = "Ovary", value = FALSE),
                         ),
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to select timpoint(s)"),
                                  checkboxInput(inputId="tr_h0", label = "0hr", value = FALSE ),
                                  checkboxInput(inputId="tr_h24", label = "24hr", value = TRUE ),
                                  checkboxInput(inputId="tr_h36", label = "36hr", value = FALSE),
                                  checkboxInput(inputId="tr_h48", label = "48hr", value = FALSE),
                         )
                ),
                tabPanel("Tissue", 
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to define your contrast:"),
                                  splitLayout(cellWidths = c("33%", "33%","33%"), 
                                              checkboxInput(inputId="tissuec1_Midgut", label = "Midgut", value = TRUE ),
                                              checkboxInput(inputId="tissuec1_FatBody", label = "FatBody", value = FALSE ),
                                              checkboxInput(inputId="tissuec1_Ovary", label = "Ovary", value = FALSE)
                                             
                                  ),
                                  splitLayout(cellWidths = c(rep("14%",7)),
                                              h4(""),h4(""),h4(""),
                                              h4("vs"),
                                              h4(""),h4(""),h4("")
                                              
                                  ),
                                  splitLayout(cellWidths = c("33%", "33%","33%"), 
                                              checkboxInput(inputId="tissuec2_Midgut", label = "Midgut", value = FALSE ),
                                              checkboxInput(inputId="tissuec2_FatBody", label = "FatBody", value = TRUE ),
                                              checkboxInput(inputId="tissuec2_Ovary", label = "Ovary", value = FALSE)
                                              
                                  ),                                 
                         ),
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to select timpoint(s)"),
                                  checkboxInput(inputId="tis_h0", label = "0hr", value = FALSE ),
                                  checkboxInput(inputId="tis_h24", label = "24hr", value = FALSE ),
                                  checkboxInput(inputId="tis_h36", label = "36hr", value = TRUE),
                                  checkboxInput(inputId="tis_h48", label = "48hr", value = FALSE),
                         ),
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to select treatment(s)"),
                                  checkboxInput(inputId="tis_GFP", label = "dsGFP", value = TRUE ),
                                  checkboxInput(inputId="tis_EcR", label = "dsEcR", value = TRUE ),
                                  checkboxInput(inputId="tis_GFP_EcR", label = "report *differential* DEGs between dsEcR and dsGFP", value = TRUE )
                         )
                ),
                tabPanel("Timepoint",
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to define your contrast:"),
                                  splitLayout(cellWidths = c("25%", "25%","25%","25%"), 
                                              checkboxInput(inputId="timec1_0", label = "0hr", value = FALSE ),
                                              checkboxInput(inputId="timec1_24", label = "24hr", value = TRUE ),
                                              checkboxInput(inputId="timec1_36", label = "36hr", value = FALSE),
                                              checkboxInput(inputId="timec1_48", label = "48hr", value = FALSE)
                                              
                                  ),
                                  splitLayout(cellWidths = c(rep("14%",7)),
                                              h4(""),h4(""),h4(""),
                                              h4("vs"),
                                              h4(""),h4(""),h4("")
                                              
                                  ),
                                  splitLayout(cellWidths = c("25%", "25%","25%","25%"), 
                                              checkboxInput(inputId="timec2_0", label = "0hr", value = TRUE ),
                                              checkboxInput(inputId="timec2_24", label = "24hr", value = FALSE ),
                                              checkboxInput(inputId="timec2_36", label = "36hr", value = FALSE),
                                              checkboxInput(inputId="timec2_48", label = "48hr", value = FALSE)
                                              
                                  ),                                 
                         ),
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to select treatment(s)"),
                                  checkboxInput(inputId="time_GFP", label = "dsGFP", value = TRUE ),
                                  checkboxInput(inputId="time_EcR", label = "dsEcR", value = TRUE ),
                                  checkboxInput(inputId="time_GFP_EcR", label = "report *differential* DEGs between dsEcR and dsGFP", value = TRUE )
                         ),
                         fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                                  h4("check to select tissue(s)"),
                                  checkboxInput(inputId="time_Midgut", label = "Midgut", value = FALSE ),
                                  checkboxInput(inputId="time_FatBody", label = "FatBody", value = TRUE ),
                                  checkboxInput(inputId="time_Ovary", label = "Ovary", value = FALSE),
                         )
                )
    )
)

ui <- fluidPage(
  tags$head(tags$style('
     #my_tooltip {
      position: absolute;
      width: 600px;
      z-index: 100;
      padding: 0;
     }
  ')),
  
  tags$script('
    $(document).ready(function() {
      // id of the plot
      $("#DEGplot").mousemove(function(e) { 

        // ID of uiOutput
        $("#my_tooltip").show();         
        $("#my_tooltip").css({             
          top: (e.pageY + 5) + "px",             
          left: (e.pageX + 5 -700) + "px"         
        });     
      });     
    });
  '),  
  
    tags$head(
        tags$style(HTML("
                        .shiny-output-error-myclass {
                          color: red;
                          font-size: 30px;}
                        ")),
        tags$style(type="text/css", "select { max-width: 270px; }"),
        tags$style(type="text/css", ".span4 { max-width: 290px; }"),
        tags$style(type="text/css", ".well { max-width: 280px; }")
    ),
    
    sidebarLayout(
        sidebarPanel(
            h3("Please select various options and click on the update button to generate plots"),
            selectInput("constrast", "Please select a contrast", 
                        choices = c("dsRNA Treatment", "Tissue", "Timepoint")
            ),
            parameter_tabs,
            fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                     h4("Differential gene expression analysis parameters"),
                     textInput( "DEG_p_thres", value="0.05","p-val cutoff", width = "100px" ),
                     checkboxInput(inputId="DEG_multi_correction", label = "multiple test correction", value = TRUE)

            ),
            fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                     h4("GO enrichment parameters"),
                     checkboxInput(inputId="GO_flag", label = "Perform GO enrichment (requires internect connection)", value = TRUE),
                     textInput( "GO_p_thres", value="0.05", "p-val cutoff", width = "100px" ),
                     checkboxInput(inputId="GO_multi_correction", label = "multiple test correction", value = TRUE)

            ),
            fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:2px',
                     h4("Plotting parameters"),
                     splitLayout(cellWidths = c("33%", "33%","33%"), 
                                 textInput("dotsize", "point size", value="3", width = "100px" ),
                                 textInput("plotwidth", "plot width", value="800", width = "100px" ),
                                 textInput("plotheight", "plot height", value="800", width = "100px" )
                     ),
            ),
            fluidRow(style='display:block;padding:7.5px;margin:0 0 6px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;',
                     actionButton("updDatSel", "Update")
            ),
            
        ),
        mainPanel(
            
            fluidRow(
                tabsetPanel(type = "tabs",
                            tabPanel("Selected samples", 
                                    #fluidRow(column(1,offset=0),column(12,offset=0, tableOutput(outputId = 'metatable' ))),
                                    fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput("metatable"))),
                                    fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'summary')))
                            ),
                            tabPanel("Volcano plot", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput('DEGplot',hover = hoverOpts(id="plot_hover",delay=250, delayType= "throttle")))),
                                     uiOutput("my_tooltip")
                            ),
                            tabPanel("Go enrichment - Up genes", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'Goupgenesmf' , height="auto"))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'placeholer2',height="40px" ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'Goupgenesbp' , height="auto")))
                            ),
                            tabPanel("Go enrichment - Down genes", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'Godowngenesmf', height="auto" ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'placeholer3' ,height="40px"))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'Godowngenesbp' , height="auto")))
                            ),
                            tabPanel("Go enrichment - Up&Down genes", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'Goallgenesmf', height="auto" ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'placeholer4' ,height="40px"))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'Goallgenesbp', height="auto" )))
                            ),
                            tabPanel("Go details - Up genes", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, span(textOutput("DTGOUpgenes_caption"), style="color:black; font-size:30px ;font-family:Arial" ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput('DTGOUpgenes' ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'placeholer6' ))),
                                     
                            ),
                            tabPanel("Go details - Down genes", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, span(textOutput("DTGODowngenes_caption"), style="color:black; font-size:30px ;font-family:Arial" ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput('DTGODowngenes' ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'placeholer7' ))),
                                     
                            ),
                            tabPanel("Go details - Up&Down genes", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, span(textOutput("DTGOUpDowngenes_caption"), style="color:black; font-size:30px ;font-family:Arial" ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput('DTGOUpDowngenes' ))),
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'placeholer8' ))),
                                     
                            ),
                             tabPanel("DEG details - Up genes", 
                                     #fluidRow(column(1,offset=0),column(12,offset=0, tableOutput(outputId = 'Upgenes' )))
                                     fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput('Upgenes' )))
                            ),
                            tabPanel("DEG details - Down  genes", 
                                     #fluidRow(column(1,offset=0),column(12,offset=0, tableOutput(outputId = 'Downgenes' )))
                                     fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput('Downgenes' )))
                            ),
                            tabPanel("DEG details - Up&Down genes", 
                                     #fluidRow(column(1,offset=0),column(12,offset=0, tableOutput(outputId = 'Siggenes' )))
                                     fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput('Siggenes' )))
                            ),
                            tabPanel("DEG details - All genes", 
                                     #fluidRow(column(1,offset=0),column(12,offset=0, tableOutput(outputId = 'Allgenes' )))
                                     downloadButton("downloadData", "Download"),
                                     fluidRow(column(1,offset=0),column(12,offset=0, DT::dataTableOutput('Allgenes' )))
                                     
                            )

                            # tabPanel("Heatmap",plotOutput(outputId = 'heatmap', width= paste0(lg_width,'px'), height = paste0(lg_height,'px'))
                            
                            # ),
                            # tabPanel("Clustermap", plotOutput(outputId = 'clusterplot',width= paste0(lg_width,'px'), height = paste0(lg_height,'px'))
                            # )
                )            
            )
            
        )
    )
)

#############################
# Shiny server: server code #
#############################

server <- function(input, output,session) {


  observeEvent(input$constrast, {
      updateTabsetPanel(session, "params", selected = input$constrast)
  }) 
  
  volcanodata4tooltip = reactiveValues()
  volcanodata4tooltip$df = NULL
  
  observe({ #UPDATE TISSUE SELECTION CHECK BOX
    tis_GFP_EcR_flag = input$tis_GFP_EcR == "TRUE"
    if (tis_GFP_EcR_flag == "TRUE")
    {
      updateCheckboxInput(session, "tis_GFP", value = tis_GFP_EcR_flag)
      updateCheckboxInput(session, "tis_EcR", value = tis_GFP_EcR_flag)
    }
    else{
      updateCheckboxInput(session, "tis_GFP", value = !tis_GFP_EcR_flag)
      updateCheckboxInput(session, "tis_EcR", value = tis_GFP_EcR_flag)
    }
    
    time_GFP_EcR_flag = input$time_GFP_EcR == "TRUE"
    if (time_GFP_EcR_flag == "TRUE")
    {
      updateCheckboxInput(session, "time_GFP", value = time_GFP_EcR_flag)
      updateCheckboxInput(session, "time_EcR", value = time_GFP_EcR_flag)
    }
    else{
      updateCheckboxInput(session, "time_GFP", value = !time_GFP_EcR_flag)
      updateCheckboxInput(session, "time_EcR", value = time_GFP_EcR_flag)
    }
    
  })
  
  
  observeEvent(input$updDatSel, {
    #print(input$constrast)
    #print(input$color)
    #print(FPKM_GOI_metadata)
    #print("begin subset")

    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Progress: ", value = 0)
    progress$inc(1/8, detail = paste("Loading required R packages"))
    
    library( DESeq2 )
    library( ggplot2 )
    library( biomaRt )
    library(EnhancedVolcano)
    library(data.table)
    library(tidyverse)
    library(httr)
    library(jsonlite)
    library(curl)
    library(RCurl)

    
    progress$inc(1/8, detail = paste("Subsetting data according to user selection"))
    #####################
    #subset (meta) data #
    #####################
    if (input$constrast=="dsRNA Treatment"){
      
      df = subset_time(metadata,"Treatment")
      df = subset_tissue(df,"Treatment")      
      
    }
    if (input$constrast=="Tissue"){
      
      #validate contrast selection
      print(validate_contrast(c(input$tissuec1_Midgut, input$tissuec1_FatBody, input$tissuec1_Ovary),
                              c(input$tissuec2_Midgut, input$tissuec2_FatBody, input$tissuec2_Ovary)))
            
      shiny::validate(
        shiny::need(validate_contrast(c(input$tissuec1_Midgut, input$tissuec1_FatBody, input$tissuec1_Ovary),
                                      c(input$tissuec2_Midgut, input$tissuec2_FatBody, input$tissuec2_Ovary)),
                    "Please check your contrast definition: only one sample can be selected on each side of 'vs', and cannot be the same"),
        errorClass = "myclass"
      )
      
      df = subset_tissue_contrast(metadata)
      df = subset_time(df,"Tissue")
      df = subset_treatment(df,"Tissue")

    }
    if (input$constrast=="Timepoint"){
      
      #validate contrast selection
      print(validate_contrast(c(input$timec2_0, input$timec2_24, input$timec2_36, input$timec2_48),
                              c(input$timec1_0, input$timec1_24, input$timec1_36, input$timec1_48)))
      
      shiny::validate(
        shiny::need(validate_contrast(c(input$timec2_0, input$timec2_24, input$timec2_36, input$timec2_48),
                                      c(input$timec1_0, input$timec1_24, input$timec1_36, input$timec1_48)),
                    "Please check your contrast definition: only one sample can be selected on each side of 'vs', and cannot be the same"),
        errorClass = "myclass"
      )
      
      df = subset_time_contrast(metadata)
      df = subset_tissue(df,"Timepoint")
      df = subset_treatment(df,"Timepoint")
    }
    
    #print(df)
    progress$inc(1/8, detail = paste("Retrieving lastest Anopheles gambiae annotation"))
    
    # create database for annotation before running loops
    species.db <- "agambiae_eg_gene"
    mart <- useMart( host="https://metazoa.ensembl.org","metazoa_mart", dataset = species.db )
    filters <- "ensembl_gene_id"
    attributes <- c( "ensembl_gene_id", "description","external_gene_name")
    
    progress$inc(1/8, detail = paste("Retrieving lastest Anopheles gambiae annotation"))
   
    #####################
    #DEG analysis       #
    #####################
    dds <- NULL
    dds_GFP <- NULL
    dds_EcR <- NULL
    contrasts <- NULL
    Tissue_differential_DEG_Flag = input$constrast=="Tissue" && input$tis_GFP_EcR==TRUE
    Time_differential_DEG_Flag = input$constrast=="Timepoint" && input$time_GFP_EcR==TRUE
    
    ###################################
    
    if (input$constrast=="dsRNA Treatment"){
      contrasts <- list( c("Treatment", "EcR", "GFP"))
    }
    if (input$constrast=="Tissue"){
      selection = c(input$tissuec1_Midgut, input$tissuec1_FatBody, input$tissuec1_Ovary,input$tissuec2_Midgut, input$tissuec2_FatBody, input$tissuec2_Ovary)
      tissue = c("Midgut", "FatBody", "Ovary", "Midgut", "FatBody", "Ovary")
      tissue_subset = tissue[selection]
      contrasts <- list( c("Tissue", tissue_subset))
    }
    if (input$constrast=="Timepoint"){
      selection = c(input$timec1_0, input$timec1_24, input$timec1_36, input$timec1_48, input$timec2_0, input$timec2_24, input$timec2_36, input$timec2_48)
      timpoints = c(0,24,36,48,0,24,36,48)
      timpoints_subset = timpoints[selection]
      contrasts <- list( c("Timepoint", timpoints_subset))     
    }
    
    print(contrasts)
    
    #print out metadata after subsetting
    output$metatable <- DT::renderDataTable({
        DT::datatable(
                      df[,-2],
                      caption = htmltools::tags$caption("Data selected and used in current analyses", style="color:black; font-size:30px ;font-family:Arial"),
                      rownames=FALSE
                      ) 
          })

    progress$inc(1/8, detail = paste("Performing differentialy gene expression analysis, this can take a few minutes"))
    
    # generate object for DESeq2 , cds stands for CountDataSet
    cds <- NULL
    dds <- NULL

    if ( Tissue_differential_DEG_Flag || Time_differential_DEG_Flag ) # report differential DEG
    {
        ####################################
        print("reporting differential DEGs") 
        if (input$constrast=="Tissue"){
          selection = c(input$tissuec1_Midgut, input$tissuec1_FatBody, input$tissuec1_Ovary,input$tissuec2_Midgut, input$tissuec2_FatBody, input$tissuec2_Ovary)
          tissue = c("Midgut", "FatBody", "Ovary", "Midgut", "FatBody", "Ovary")
          tissue_subset = tissue[selection]
          contrasts <- list( c("Tissue", tissue_subset))
        }
        if (input$constrast=="Timepoint"){
          selection = c(input$timec1_0, input$timec1_24, input$timec1_36, input$timec1_48, input$timec2_0, input$timec2_24, input$timec2_36, input$timec2_48)
          timpoints = c(0,24,36,48,0,24,36,48)
          timpoints_subset = timpoints[selection]
          contrasts <- list( c("Timepoint", timpoints_subset))     
        }
        
        print(contrasts)
        
        #print out metadata after subsetting
        output$metatable <- DT::renderDataTable({
          DT::datatable(
            df[,-2],
            caption = htmltools::tags$caption("Data selected and used in current analyses", style="color:black; font-size:30px ;font-family:Arial"),
            rownames=FALSE
          ) 
        })
        
        #progress$inc(1/8, detail = paste("Performing differentialy gene expression analysis"))
        
        # subset metadata
        df_GFP = df[which(df$Treatment=="GFP"),]
        df_EcR = df[which(df$Treatment=="EcR"),]
        
        print("df_GFP")
        print(df_GFP)
        print("df_EcR")
        print(df_EcR)
        # generate object for DESeq2 , cds stands for CountDataSet
        cds_GFP <- NULL
        cds_EcR <- NULL
        
        if (contrasts[[1]][1] == "Timepoint")
        {
          cds_GFP <- DESeqDataSetFromHTSeqCount( sampleTable = df_GFP,
                                             directory = ".",
                                             design = ~Timepoint)
          cds_EcR <- DESeqDataSetFromHTSeqCount( sampleTable = df_EcR,
                                                 directory = ".",
                                                 design = ~Timepoint)
        }
        if (contrasts[[1]][1] == "Tissue")
        {
          cds_GFP <- DESeqDataSetFromHTSeqCount( sampleTable = df_GFP,
                                             directory = ".",
                                             design = ~Tissue)
          cds_EcR <- DESeqDataSetFromHTSeqCount( sampleTable = df_EcR,
                                             directory = ".",
                                             design = ~Tissue)
        }
        
        #progress$inc(1/8, detail = paste("Performing differentialy gene expression analysis: estimating dispersions and fitting models"))
        # NBM fit, dispersion estimate, dds stands for Deseq2DataSeq
        dds_GFP <- DESeq(cds_GFP)
        dds_EcR <- DESeq(cds_EcR)
        #plotDispEsts( dds, main = "RNASeq" )
        
    }
    else{
      
      if (contrasts[[1]][1] == "Treatment")
      {
        cds <- DESeqDataSetFromHTSeqCount( sampleTable = df,
                                           directory = ".",
                                           design = ~Treatment)
        
      }
      if (contrasts[[1]][1] == "Timepoint")
      {
        cds <- DESeqDataSetFromHTSeqCount( sampleTable = df,
                                           directory = ".",
                                           design = ~Timepoint)
      }
      if (contrasts[[1]][1] == "Tissue")
      {
        cds <- DESeqDataSetFromHTSeqCount( sampleTable = df,
                                           directory = ".",
                                           design = ~Tissue)
      }
      
      #progress$inc(1/8, detail = paste("Performing differentialy gene expression analysis: estimating dispersions and fitting models"))
      # NBM fit, dispersion estimate, dds stands for Deseq2DataSeq
      dds <- DESeq(cds)
      #plotDispEsts( dds, main = "RNASeq" )
      
    }
    
    
    ###########
    ##DEG data#
    ###########
    p_val_sel = NULL
    p_val_sel_thres = as.numeric(input$DEG_p_thres)
    #print(p_val_sel_thres)
    if (input$DEG_multi_correction == FALSE){p_val_sel = "pvalue"}
    else {p_val_sel = "padj"}
    
    
    ######## calculate *differential* DEG ########
    #variables to store DEG results
    table1N2_dedup_sig = NULL
    table1N2_dedup_up = NULL
    table1N2_dedu_down = NULL
    #variables of additions to title of data table
    EcR_GFP_updown_title = NULL
    EcR_GFP_up_title = NULL
    EcR_GFP_down_title = NULL
    EcR_GFP_all_title = NULL
    
    if ( Tissue_differential_DEG_Flag || Time_differential_DEG_Flag ) 
    {
        currcon <- contrasts[[1]]
        #print(paste0("currcon:",currcon))
        res_GFP <- results( dds_GFP, contrast = currcon )
        res_EcR <- results( dds_EcR, contrast = currcon )
        
        #######the data table 
        EcR_GFP_updown_title = "This table consists of EcR-only DEGs and GFP-only DEGs, separable by the 'Sig-source' column"
        EcR_GFP_up_title = "Up regulated genes consists of (1) up regulated EcR-only DEGs, (2)and down regulated GFP-only DEGs, separable by the 'Sig-source' column"
        EcR_GFP_down_title = "Down regulated genes consists of (1) down regulated EcR-only DEGs, (2)and up regulated GFP-only DEGs, separable by the 'Sig-source' column"
        EcR_GFP_all_title = "This table contains 2 sets of results for each gene, GFG and EcR, separable by the 'source' column"
        
        #######all genes GFP ######
        all <- subset( res_GFP )
        all.genes <- rownames( all )
        all.anno <- getBM( attributes = attributes,filters = filters,  values = all.genes,
                           mart = mart, uniqueRows = T)
        
        table1=as.data.frame(all) 
        table2=as.data.frame(all.anno)
        
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table based on by.x and by.y
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedup_ALL=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedup_ALL= table1N2_dedup_ALL[order(table1N2_dedup_ALL$padj,table1N2_dedup_ALL$pvalue,decreasing=FALSE),]
        table1N2_dedup_ALL = table1N2_dedup_ALL[, -which(names(table1N2_dedup_ALL) %in% c("V2", "stat"))]
        colnames(table1N2_dedup_ALL)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol",	"A. albi gene ID")
        table1N2_dedup_ALL_GFP = table1N2_dedup_ALL
        table1N2_dedup_ALL_GFP$source = rep("GFP",nrow(table1N2_dedup_ALL_GFP))
        
        #######all genes EcR ######
        all <- subset( res_EcR )
        all.genes <- rownames( all )
        all.anno <- getBM( attributes = attributes,filters = filters,  values = all.genes,
                           mart = mart, uniqueRows = T)
        
        table1=as.data.frame(all) 
        table2=as.data.frame(all.anno)
        
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table based on by.x and by.y
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedup_ALL=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedup_ALL= table1N2_dedup_ALL[order(table1N2_dedup_ALL$padj,table1N2_dedup_ALL$pvalue,decreasing=FALSE),]
        table1N2_dedup_ALL = table1N2_dedup_ALL[, -which(names(table1N2_dedup_ALL) %in% c("V2", "stat"))]
        colnames(table1N2_dedup_ALL)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol",	"A. albi gene ID")
        table1N2_dedup_ALL_EcR = table1N2_dedup_ALL       
        table1N2_dedup_ALL_EcR$source = rep("EcR",nrow(table1N2_dedup_ALL_EcR))
        
        ##merge GFP and EcR
        table1N2_dedup_ALL_EcR_GFP  = rbind.data.frame(table1N2_dedup_ALL_GFP,table1N2_dedup_ALL_EcR)
        
        output$Allgenes <- DT::renderDataTable({
          DT::datatable(table1N2_dedup_ALL_EcR_GFP, 
                        extensions = 'Buttons',
                        options = list(lengthMenu = c(100, 250, 500,1000), 
                                       pageLength = 100,
                                       paging = TRUE,
                                       searching = TRUE,
                                       fixedColumns = TRUE,
                                       autoWidth = TRUE,
                                       ordering = TRUE,
                                       dom = 'Bfrtip',
                                       buttons = c('copy', 'csv', 'excel')
                                       ),
                        class = "display",
                        caption = htmltools::tags$caption(paste0(length(table1N2_dedup_ALL_EcR_GFP$GeneID), " genes (all genes in the genome)"), htmltools::tags$br(), paste0(EcR_GFP_all_title), style="color:black; font-size:30px ;font-family:Arial"),
                        rownames=FALSE
          )%>% 
            formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
        })
        datasetInput1 <- reactive({
          table1N2_dedup_ALL_EcR_GFP
        })
        
        myData <- reactive({table1N2_dedup_ALL_EcR_GFP})
        output$downloadData <- downloadHandler(
          filename = function() {
            paste("data-", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            write.csv(myData(), file)
          }
        )
        
        
        #####Sig genes GFP##############
        idx = which(res_GFP[,p_val_sel]<p_val_sel_thres)
        Sig <- res_GFP[idx,]
        Sig.genes <- rownames( Sig )
        Sig.anno <- getBM( attributes = attributes, filters = filters, values = Sig.genes,
                           mart = mart)
        table1=as.data.frame(Sig) 
        table2=as.data.frame(Sig.anno)
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedup_sig=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedup_sig = table1N2_dedup_sig[order(table1N2_dedup_sig$padj,table1N2_dedup_sig$pvalue,decreasing=FALSE),]
        table1N2_dedup_sig = table1N2_dedup_sig[, -which(names(table1N2_dedup_sig) %in% c("V2", "stat"))]
        colnames(table1N2_dedup_sig)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol", "A. albi gene ID")
        table1N2_dedup_sig_GFP = table1N2_dedup_sig
        table1N2_dedup_sig_GFP$Sig_source = rep("GFP",nrow(table1N2_dedup_sig_GFP))

        #####Sig genes EcR##############
        idx = which(res_EcR[,p_val_sel]<p_val_sel_thres)
        Sig <- res_EcR[idx,]
        Sig.genes <- rownames( Sig )
        Sig.anno <- getBM( attributes = attributes, filters = filters, values = Sig.genes,
                           mart = mart)
        table1=as.data.frame(Sig) 
        table2=as.data.frame(Sig.anno)
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedup_sig=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedup_sig = table1N2_dedup_sig[order(table1N2_dedup_sig$padj,table1N2_dedup_sig$pvalue,decreasing=FALSE),]
        table1N2_dedup_sig = table1N2_dedup_sig[, -which(names(table1N2_dedup_sig) %in% c("V2", "stat"))]
        colnames(table1N2_dedup_sig)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol", "A. albi gene ID")
        table1N2_dedup_sig_EcR = table1N2_dedup_sig
        table1N2_dedup_sig_EcR$Sig_source = rep("EcR",nrow(table1N2_dedup_sig_EcR))
        
        ##merge GFP and EcR
        table1N2_dedup_sig = rbind.data.frame(table1N2_dedup_sig_GFP,table1N2_dedup_sig_EcR)
  
        output$Siggenes <- DT::renderDataTable({
          DT::datatable(table1N2_dedup_sig, 
                        options = list(lengthMenu = c(100, 250, 500,1000), pageLength = 100),
                        caption = htmltools::tags$caption(paste0(length(table1N2_dedup_sig$GeneID), " up&down regulated genes"), htmltools::tags$br(), paste0(EcR_GFP_updown_title), style="color:black; font-size:30px ;font-family:Arial"),
                        rownames=FALSE
          )%>% 
            formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
        })
        
        #####get GFP_only_DEGs and EcR_only_DEGs##############
        GFP_only_GeneID = setdiff(table1N2_dedup_sig_GFP$GeneID,table1N2_dedup_sig_EcR$GeneID)
        EcR_only_GeneID = setdiff(table1N2_dedup_sig_EcR$GeneID,table1N2_dedup_sig_GFP$GeneID)
        
        GFP_only_DEG = subset(table1N2_dedup_sig_GFP, GeneID %in% GFP_only_GeneID )
        EcR_only_DEG = subset(table1N2_dedup_sig_EcR, GeneID %in% EcR_only_GeneID )
        #####Up genes EcR_only_DEG:up or GFP_only_DEG:down##############
        idx = which(EcR_only_DEG[,"Log2FoldChange"] > 0)
        up_EcR <- EcR_only_DEG[idx,]
        
        idx = which(GFP_only_DEG[,"Log2FoldChange"] < 0)
        down_GFP <- GFP_only_DEG[idx,]
        
        up = rbind.data.frame(up_EcR,down_GFP)
        table1N2_dedup_up = up
        
        output$Upgenes <- DT::renderDataTable({
          DT::datatable(table1N2_dedup_up, options = list(lengthMenu = c(100, 250, 500,1000), pageLength = 100),
                        caption = htmltools::tags$caption(paste0(length(table1N2_dedup_up$GeneID), " up regulated genes"), htmltools::tags$br(), paste0(EcR_GFP_up_title), style="color:black; font-size:25px ;font-family:Arial"),
                        rownames=FALSE
          )%>% 
            formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
        })
        
        
        #####Down genes EcR_only_DEG:down or GFP_only_DEG:up##############
        idx = which(EcR_only_DEG[,"Log2FoldChange"] < 0)
        down_EcR <- EcR_only_DEG[idx,]
        
        idx = which(GFP_only_DEG[,"Log2FoldChange"] > 0)
        up_GFP <- GFP_only_DEG[idx,]
        
        down = rbind.data.frame(down_EcR,up_GFP)
        table1N2_dedu_down = down
        
        output$Downgenes <- DT::renderDataTable({ DT::datatable(table1N2_dedu_down, options = list(lengthMenu = c(100, 250, 500,1000), pageLength = 100),
                                                                caption = htmltools::tags$caption(paste0(length(table1N2_dedu_down$GeneID), " down regulated genes"),htmltools::tags$br(),paste0(EcR_GFP_down_title), style="color:black; font-size:25px ;font-family:Arial"),
                                                                rownames=FALSE
        ) %>% 
            formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
        }) 
    }
    else{

        ######## calculate *normal* DEGs ######## 
        #######all genes######
        currcon <- contrasts[[1]]
        #print(paste0("currcon:",currcon))
        res <- results( dds, contrast = currcon )
        ## Order by adjusted p-value
        #res <- res[order(res$padj), ]
        all <- subset( res )
        all.genes <- rownames( all )
        all.anno <- getBM( attributes = attributes,filters = filters,  values = all.genes,
                           mart = mart, uniqueRows = T)
        
        table1=as.data.frame(all) 
        table2=as.data.frame(all.anno)
        
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table based on by.x and by.y
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedup_ALL=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedup_ALL= table1N2_dedup_ALL[order(table1N2_dedup_ALL$padj,table1N2_dedup_ALL$pvalue,decreasing=FALSE),]
        table1N2_dedup_ALL = table1N2_dedup_ALL[, -which(names(table1N2_dedup_ALL) %in% c("V2", "stat"))]
        colnames(table1N2_dedup_ALL)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol", "A. albi gene ID")
        #output$Allgenes <- renderTable({table1N2_dedup_ALL},      caption = paste0(length(table1N2_dedup_ALL$GeneID), " genes (all genes in the genome)"),
                                       #caption.placement = getOption("xtable.caption.placement", "top"), 
                                       #caption.width = getOption("xtable.caption.width", NULL))
        
        output$Allgenes <- DT::renderDataTable({
          DT::datatable(table1N2_dedup_ALL, 
                        options = list(lengthMenu = c(100, 250, 500,1000), 
                                       pageLength = 100,
                                       paging = TRUE,
                                       searching = TRUE,
                                       fixedColumns = TRUE,
                                       autoWidth = TRUE,
                                       ordering = TRUE,
                                       dom = 'Bfrtip',
                                       buttons = c('copy', 'csv', 'excel')
                        ),
                        class = "display",
                        caption = htmltools::tags$caption(paste0(length(table1N2_dedup_ALL$GeneID), " genes (all genes in the genome)"), style="color:black; font-size:25px ;font-family:Arial"),
                        rownames=FALSE
                        )%>% 
            formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
        })
        
        myData = reactive({table1N2_dedup_ALL})
        output$downloadData <- downloadHandler(
          filename = function() {
            paste("data-", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            write.csv(myData(), file)
          }
        )
        
        #####Sig genes##############
        idx = which(res[,p_val_sel]<p_val_sel_thres)
        Sig <- res[idx,]
        Sig.genes <- rownames( Sig )
        Sig.anno <- getBM( attributes = attributes, filters = filters, values = Sig.genes,
                           mart = mart)
        table1=as.data.frame(Sig) 
        table2=as.data.frame(Sig.anno)
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedup_sig=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedup_sig = table1N2_dedup_sig[order(table1N2_dedup_sig$padj,table1N2_dedup_sig$pvalue,decreasing=FALSE),]
        table1N2_dedup_sig = table1N2_dedup_sig[, -which(names(table1N2_dedup_sig) %in% c("V2", "stat"))]
        colnames(table1N2_dedup_sig)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol", "A. albi gene ID")
        #output$Siggenes <- renderTable({table1N2_dedup_sig},      caption = paste0(length(table1N2_dedup_sig$GeneID), " up&down regulated genes"),
                                                                #caption.placement = getOption("xtable.caption.placement", "top"), 
                                                                #caption.width = getOption("xtable.caption.width", NULL))
        output$Siggenes <- DT::renderDataTable({
          DT::datatable(table1N2_dedup_sig, 
                        options = list(lengthMenu = c(100, 250, 500,1000), pageLength = 100),
                        caption = htmltools::tags$caption(paste0(length(table1N2_dedup_sig$GeneID), " up&down regulated genes"), style="color:black; font-size:25px ;font-family:Arial"),
                        rownames=FALSE
                        )%>% 
            formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
        })
        
        #####Up genes############
        idx = which(Sig[,"log2FoldChange"] > 0)
        up <- Sig[idx,]
        #up <- subset(res, p_val_sel < p_val_sel_thres & log2FoldChange > 0 )
        up.genes <- rownames( up )
        up.anno <- getBM( attributes = attributes, filters = filters, values = up.genes,
                          mart = mart)
        table1=as.data.frame(up) 
        table2=as.data.frame(up.anno)
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedup_up=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedup_up = table1N2_dedup_up[order(table1N2_dedup_up$padj,table1N2_dedup_up$pvalue,decreasing=FALSE),]
        table1N2_dedup_up = table1N2_dedup_up[, -which(names(table1N2_dedup_up) %in% c("V2", "stat"))]
        colnames(table1N2_dedup_up)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol", "A. albi gene ID")
        #output$Upgenes <- renderTable({table1N2_dedup_up},      caption = paste0(length(table1N2_dedup_up$GeneID), " up regulated genes"),
                                                                #caption.placement = getOption("xtable.caption.placement", "top"), 
                                                                #caption.width = getOption("xtable.caption.width", NULL))
        output$Upgenes <- DT::renderDataTable({
          DT::datatable(table1N2_dedup_up, options = list(lengthMenu = c(100, 250, 500,1000), pageLength = 100),
                        caption = htmltools::tags$caption(paste0(length(table1N2_dedup_up$GeneID), " up regulated genes"), style="color:black; font-size:30px ;font-family:Arial"),
                        rownames=FALSE
                        )%>% 
            formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
        })
        
        #####Down genes############
        idx = which(Sig[,"log2FoldChange"] < 0)
        down <- Sig[idx,]
        #down <- subset(res, p_val_sel < p_val_sel_thres & log2FoldChange < 0 )
        down.genes <- rownames( down )
        down.anno <- getBM( attributes = attributes, filters = filters, values = down.genes,
                            mart = mart)
        table1=as.data.frame(down) 
        table2=as.data.frame(down.anno)
        table1=setDT(table1, keep.rownames = TRUE)[] # convert rownames into the first row
        x=colnames(table1) # fix table1 's column names
        x[1]="geneID" # add "geneID" to column names
        colnames(table1)=x # fix table 1 's column names
        table1N2=merge(table1,table2,by.x="geneID",by.y="ensembl_gene_id",all.x=TRUE) # merge annotation and results table
        table1N2_dedup=table1N2[!duplicated(table1N2[,"geneID"]),] # remove duplicates entrys (a result of multiple hits per query in the bioMart DB)
        idx=match(table1N2_dedup$geneID,anno[,1])##annotate genes using custom lookup list
        
        table1N2_dedu_down=as.data.frame(cbind(table1N2_dedup,anno[idx,3],anno[idx,2]))
        table1N2_dedu_down = table1N2_dedu_down[order(table1N2_dedu_down$padj,table1N2_dedu_down$pvalue,decreasing=FALSE),]
        table1N2_dedu_down = table1N2_dedu_down[, -which(names(table1N2_dedu_down) %in% c("V2", "stat"))]
        colnames(table1N2_dedu_down)= c("GeneID",	"Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val",	"description",	"gene_symbol", "A. albi gene ID")
        
        output$Downgenes <- DT::renderDataTable({ DT::datatable(table1N2_dedu_down, options = list(lengthMenu = c(100, 250, 500,1000), pageLength = 100),
                                                                caption = htmltools::tags$caption(paste0(length(table1N2_dedu_down$GeneID), " down regulated genes"), style="color:black; font-size:30px ;font-family:Arial"),
                                                                rownames=FALSE
                                                                ) %>% 
                                                  formatRound(c("Avg. expr",	"Log2FoldChange",	"SE(log2FC)",	"p-val",	"adjusted p-val"),3)
                                        }) 
                                    
    
    }
        
    
    ################
    ##plot volcano #
    ################
    
    if ( Tissue_differential_DEG_Flag || Time_differential_DEG_Flag ) 
    {
      output$DEGplot <- renderPlot({isolate({
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Volcano plot is disabled for 'differential DEG' analysis, since there are two sets of DEGs (one from GFP and one from EcR).\n",
                                     "\n",
                                     ""), 
             cex = 1.6, col = "black")
      })})
    }
    else
    {
        p_val_sel = NULL
        if (input$DEG_multi_correction == FALSE){p_val_sel = "pvalue"}
        else {p_val_sel = "padj"}
        fac = 1
        currcon <- contrasts[[fac]]
        res_for_eVolcano <- results( dds, contrast = currcon )
        
        #print(p_val_sel)
        res_for_eVolcano$geneID = rownames(res_for_eVolcano)
        #print(res_for_eVolcano[which(res_for_eVolcano$geneID=="AGAP004519"),])
        #print(res_for_eVolcano[which(res_for_eVolcano$geneID=="AGAP002593"),])
        #print(res_for_eVolcano)
        
        output$DEGplot <- renderPlot({
          isolate({
            EnhancedVolcano(res_for_eVolcano,
                                 #selectLab = c("AGAP029539", "AGAP001826", "AGAP004203") ,
                                 lab = rownames(res_for_eVolcano),
                                 x = 'log2FoldChange',
                                 y = p_val_sel,
                                 #xlim = c(-5, 8)
                                 title=paste0( currcon[2], "_vs_", currcon[3]),
                                 subtitle = "",
                                 pCutoff = as.numeric(input$DEG_p_thres),
                                 FCcutoff=0.5849625, # fold change = 1.5
                                 #transcriptPointSize= as.numeric(input$dotsize),
                                 #transcriptLabSize = 4.0,
                                 col = c("grey30", "forestgreen", "royalblue", "red2"),
                                 axisLabSize = 24,
                                 titleLabSize = 24,
                                 subtitleLabSize = 24,
                                 captionLabSize = 24,
                                 legendLabSize = 24,
            )

          })
        },
          height= as.numeric(input$plotheight), width= as.numeric(input$plotwidth)
        )
        
        x = res_for_eVolcano$log2FoldChange
        y = -log(res_for_eVolcano[,p_val_sel], base = 10)
        volcano_data  = data.frame("Log2foldchange" = x, "neg_Log10P"= y, "geneID" = rownames(res_for_eVolcano))
       # print(table1N2_dedup_ALL)
        annotation_df = table1N2_dedup_ALL[,c("GeneID","description")]
        colnames(annotation_df) = c("geneID","description")
        volcano_data = merge(volcano_data, annotation_df, by="geneID", all.x = TRUE)
        #print(volcano_data)
        volcanodata4tooltip$df=volcano_data
    }
    
    
    output$DEGtab <- renderTable({
      isolate({})
    })
    
    progress$inc(1/8, detail = paste("Conducting GO enrichment analysis"))
    
    #print(table1N2_dedup_ALL[which(table1N2_dedup_ALL$GeneID=="AGAP004880"),])
    #print(table1N2_dedup_sig[which(table1N2_dedup_sig$GeneID=="AGAP004880"),])
    #print(table1N2_dedup_up[which(table1N2_dedup_up$GeneID=="AGAP004880"),])
    #print(table1N2_dedu_down[which(table1N2_dedu_down$GeneID=="AGAP004880"),])
    #####################
    #GO enrichment      #
    #####################
    
    if (input$GO_flag==TRUE)
    {
         #get significant gene IDs
        sigGenes = unique(table1N2_dedup_sig$GeneID)
        upGenes = unique(table1N2_dedup_up$GeneID)
        downGenes = unique(table1N2_dedu_down$GeneID)
        if (length(sigGenes)>5000) # if  significant genes >10000, do a hard cutoff, since the limit for GO enrichment API is 10000 gene IDs
        {
          sigGenes <- sigGenes[1:5000]
        } 
        if (length(upGenes)>5000)
        {
          upGenes = upGenes[1:5000]
        }
        if (length(downGenes)>5000)
        {
          downGenes = downGenes[1:5000]
        }
      
        sig_GeneIDs = paste(sigGenes, collapse = ",")
        up_GeneIDs = paste(upGenes, collapse = ",")
        down_GeneIDs = paste(downGenes, collapse = ",")
        ############## molecular fucntion #############
        path <- "http://pantherdb.org/services/oai/pantherdb/enrich/overrep"
        AnnoDataSet = "GO:0003674" #molecular function
        
        print(paste0("enrichment analysis with ", length(sigGenes), " genes"))
        results_df_sig_mf =  PantherGOenrich(sig_GeneIDs,path, AnnoDataSet )
        Sys.sleep(1)
        print(paste0("enrichment analysis with ", length(upGenes), " genes"))
        results_df_up_mf  =  PantherGOenrich(up_GeneIDs,path, AnnoDataSet )
        Sys.sleep(1)
        print(paste0("enrichment analysis with ", length(downGenes), " genes"))
        results_df_down_mf =  PantherGOenrich(down_GeneIDs,path, AnnoDataSet )
    
        ############## END molecular fucntion #############
        
        
        ############## biological processes #############
        AnnoDataSet = "GO:0008150" #biological processes
        
        print(paste0("enrichment analysis with ", length(sigGenes), " genes"))
        results_df_sig_bp =  PantherGOenrich(sig_GeneIDs,path, AnnoDataSet )
        Sys.sleep(1)
        print(paste0("enrichment analysis with ", length(upGenes), " genes"))
        results_df_up_bp  =  PantherGOenrich(up_GeneIDs,path, AnnoDataSet )
        Sys.sleep(1)
        print(paste0("enrichment analysis with ", length(downGenes), " genes"))
        results_df_down_bp =  PantherGOenrich(down_GeneIDs,path, AnnoDataSet )
        
        ############## END biological processes #############
        
        ############################
        #render GO enrichment plots#
        ############################
        
        Goterm_num_sig_mf=nrow(results_df_sig_mf %>% filter(plus_minus=="+"))
        plot_height=cal_plot_height(Goterm_num_sig_mf)
        
        if(Goterm_num_sig_mf>=1){
          output$Goallgenesmf <- renderPlot({
            isolate({
              title_size=30
              axis_text_size=25
              axix_title_size=25
              legend_title_size=axix_title_size
              legend_text_size=axis_text_size
      
              plot(ggplot(results_df_sig_mf %>% filter(plus_minus=="+"), aes(x=hit_perc, y=term.label.wrapped,  colour=fdr, size=number_in_list)) +
                     geom_point() +
                     expand_limits(x=0) +
                     labs(x="Hits (%)", y="GO term", colour="p value (adjusted)", size="Count")+
                     ggtitle(paste0("Molecular function enrichment", " (",Goterm_num_sig_mf,")")) +  
                     theme(
                       axis.title.x=element_text(size=axix_title_size),
                       axis.title.y=element_text(size=axix_title_size), 
                       axis.text.x=element_text(size=axix_title_size),
                       axis.text.y=element_text(size=axix_title_size),     
                       plot.title = element_text(hjust=0.5, size=title_size),
                       legend.text=element_text(size=legend_text_size),
                       legend.title=element_text(size=legend_title_size),
                     )+
                     scale_color_gradient(low = "red", high = "blue")#+
                     #scale_size(range=c(pngheight/400,pngheight/100))
              )
            })
          }, height = plot_height)
        }else{
          output$Goallgenesmf <- renderPlot({
            isolate({
              par(mar = c(0,0,0,0))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, paste("There are no significantly enriched GO terms of molecular functions\n",
                                           "\n",
                                           ""), 
                   cex = 1.6, col = "black")
            })
          })
        }
        
        Goterm_num_sig_b=nrow(results_df_sig_bp %>% filter(plus_minus=="+"))
        plot_height=cal_plot_height(Goterm_num_sig_b)
        
        if(Goterm_num_sig_b>=1){
          output$Goallgenesbp <- renderPlot({
            isolate({
              title_size=30
              axis_text_size=25
              axix_title_size=25
              legend_title_size=axix_title_size
              legend_text_size=axis_text_size
              plot(ggplot(results_df_sig_bp %>% filter(plus_minus=="+"), aes(x=hit_perc, y=term.label.wrapped,  colour=fdr, size=number_in_list)) +
                     geom_point() +
                     expand_limits(x=0) +
                     labs(x="Hits (%)", y="GO term", colour="p value (adjusted)", size="Count")+
                     ggtitle(paste0("Biological process enrichment", " (",Goterm_num_sig_b,")")) +    
                     theme(
                       axis.title.x=element_text(size=axix_title_size),
                       axis.title.y=element_text(size=axix_title_size), 
                       axis.text.x=element_text(size=axix_title_size),
                       axis.text.y=element_text(size=axix_title_size),     
                       plot.title = element_text(hjust=0.5, size=title_size),
                       legend.text=element_text(size=legend_text_size),
                       legend.title=element_text(size=legend_title_size),
                     )+
                     scale_color_gradient(low = "red", high = "blue")#+
                   #scale_size(range=c(pngheight/400,pngheight/100))
              )
            })
          }, height = plot_height)
        }else{
          output$Goallgenesbp <- renderPlot({
            isolate({
              par(mar = c(0,0,0,0))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, paste("There are no significantly enriched GO terms of biological processes\n",
                                           "\n",
                                           ""), 
                   cex = 1.6, col = "black")
            })
          })
        }
        
        ########GO plot up genes########
        Goterm_num_up_mf=nrow(results_df_up_mf %>% filter(plus_minus=="+"))
        plot_height=cal_plot_height(Goterm_num_up_mf)
        
        if(Goterm_num_up_mf>=1){
          output$Goupgenesmf <- renderPlot({
            isolate({
              title_size=30
              axis_text_size=25
              axix_title_size=25
              legend_title_size=axix_title_size
              legend_text_size=axis_text_size
              
              plot(ggplot(results_df_up_mf %>% filter(plus_minus=="+"), aes(x=hit_perc, y=term.label.wrapped,  colour=fdr, size=number_in_list)) +
                     geom_point() +
                     expand_limits(x=0) +
                     labs(x="Hits (%)", y="GO term", colour="p value (adjusted)", size="Count")+
                     ggtitle(paste0("Molecular function enrichment", " (",Goterm_num_up_mf,")")) + 
                     theme(
                       axis.title.x=element_text(size=axix_title_size),
                       axis.title.y=element_text(size=axix_title_size), 
                       axis.text.x=element_text(size=axix_title_size),
                       axis.text.y=element_text(size=axix_title_size),     
                       plot.title = element_text(hjust=0.5, size=title_size),
                       legend.text=element_text(size=legend_text_size),
                       legend.title=element_text(size=legend_title_size),
                     )+
                     scale_color_gradient(low = "red", high = "blue")#+
                   #scale_size(range=c(pngheight/400,pngheight/100))
              )
            })
          }, height = plot_height)
        }else{
          output$Goupgenesmf <- renderPlot({
          isolate({
            par(mar = c(0,0,0,0))
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, paste("There are no significantly enriched GO terms of molecular functions\n",
                                         "\n",
                                         ""), 
                 cex = 1.6, col = "black")
          })
          })
        }
        
        
        
        Goterm_num_up_bp=nrow(results_df_up_bp %>% filter(plus_minus=="+"))
        plot_height=cal_plot_height(Goterm_num_up_bp)
        
        if(Goterm_num_up_bp>=1){
          output$Goupgenesbp <- renderPlot({
            isolate({
              title_size=30
              axis_text_size=25
              axix_title_size=25
              legend_title_size=axix_title_size
              legend_text_size=axis_text_size
              plot(ggplot(results_df_up_bp %>% filter(plus_minus=="+"), aes(x=hit_perc, y=term.label.wrapped,  colour=fdr, size=number_in_list)) +
                     geom_point() +
                     expand_limits(x=0) +
                     labs(x="Hits (%)", y="GO term", colour="p value (adjusted)", size="Count")+
                     ggtitle(paste0("Biological process enrichment", " (",Goterm_num_up_bp,")")) +  
                     theme(
                       axis.title.x=element_text(size=axix_title_size),
                       axis.title.y=element_text(size=axix_title_size), 
                       axis.text.x=element_text(size=axix_title_size),
                       axis.text.y=element_text(size=axix_title_size),     
                       plot.title = element_text(hjust=0.5, size=title_size),
                       legend.text=element_text(size=legend_text_size),
                       legend.title=element_text(size=legend_title_size),
                     )+
                     scale_color_gradient(low = "red", high = "blue")#+
                   #scale_size(range=c(pngheight/400,pngheight/100))
              )
            })
          }, height = plot_height)
        }else{
          output$Goupgenesbp <- renderPlot({
            isolate({
              par(mar = c(0,0,0,0))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, paste("There are no significantly enriched GO terms of biological processes\n",
                                           "\n",
                                           ""), 
                   cex = 1.6, col = "black")
            })
          })
        }
    
        ########GO plot down genes########
        Goterm_num_down_mf=nrow(results_df_down_mf %>% filter(plus_minus=="+"))
        plot_height=cal_plot_height(Goterm_num_down_mf)
        
        if(Goterm_num_down_mf>=1){
          output$Godowngenesmf <- renderPlot({
            isolate({
              title_size=30
              axis_text_size=25
              axix_title_size=25
              legend_title_size=axix_title_size
              legend_text_size=axis_text_size
              
              plot(ggplot(results_df_down_mf %>% filter(plus_minus=="+"), aes(x=hit_perc, y=term.label.wrapped,  colour=fdr, size=number_in_list)) +
                     geom_point() +
                     expand_limits(x=0) +
                     labs(x="Hits (%)", y="GO term", colour="p value (adjusted)", size="Count")+
                     ggtitle(paste0("Molecular function enrichment", " (",Goterm_num_down_mf,")")) + 
                     theme(
                       axis.title.x=element_text(size=axix_title_size),
                       axis.title.y=element_text(size=axix_title_size), 
                       axis.text.x=element_text(size=axix_title_size),
                       axis.text.y=element_text(size=axix_title_size),     
                       plot.title = element_text(hjust=0.5, size=title_size),
                       legend.text=element_text(size=legend_text_size),
                       legend.title=element_text(size=legend_title_size),
                     )+
                     scale_color_gradient(low = "red", high = "blue")#+
                   #scale_size(range=c(pngheight/400,pngheight/100))
              )
            })
          }, height = plot_height)
        }else{
          output$Godowngenesmf <- renderPlot({
            isolate({
              par(mar = c(0,0,0,0))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, paste("There are no significantly enriched GO terms of molecular functions\n",
                                           "\n",
                                           ""), 
                   cex = 1.6, col = "black")
            })
          })
        }
        
        Goterm_num_down_bp=nrow(results_df_down_bp %>% filter(plus_minus=="+"))
        plot_height=cal_plot_height(Goterm_num_down_bp)
        
        if(Goterm_num_down_bp>=1){
          output$Godowngenesbp <- renderPlot({
            isolate({
              title_size=30
              axis_text_size=25
              axix_title_size=25
              legend_title_size=axix_title_size
              legend_text_size=axis_text_size
              plot(ggplot(results_df_down_bp %>% filter(plus_minus=="+"), aes(x=hit_perc, y=term.label.wrapped,  colour=fdr, size=number_in_list)) +
                     geom_point() +
                     expand_limits(x=0) +
                     labs(x="Hits (%)", y="GO term", colour="p value (adjusted)", size="Count")+
                     ggtitle(paste0("Biological process enrichment", " (",Goterm_num_down_bp,")")) +   
                     theme(
                       axis.title.x=element_text(size=axix_title_size),
                       axis.title.y=element_text(size=axix_title_size), 
                       axis.text.x=element_text(size=axix_title_size),
                       axis.text.y=element_text(size=axix_title_size),     
                       plot.title = element_text(hjust=0.5, size=title_size),
                       legend.text=element_text(size=legend_text_size),
                       legend.title=element_text(size=legend_title_size),
                     )+
                     scale_color_gradient(low = "red", high = "blue")#+
                   #scale_size(range=c(pngheight/400,pngheight/100))
              )
            })
          }, height = plot_height)
        }else{
          output$Godowngenesbp <- renderPlot({
            isolate({
              par(mar = c(0,0,0,0))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, paste("There are no significantly enriched GO terms of biological processes\n",
                                           "\n",
                                           ""), 
                   cex = 1.6, col = "black")
            })
          })
        }
      
        #############################
        #render GO & genes tables ###
        #############################

        ##render DT tables
        ########GO data up genes########
        Goterm_num_up_mf=nrow(results_df_up_mf %>% filter(plus_minus=="+"))
        Goterm_num_up_bp=nrow(results_df_up_bp %>% filter(plus_minus=="+"))
        Goterm_num_up = Goterm_num_up_mf + Goterm_num_up_bp
        
        if(Goterm_num_up>=1){
          ##get data for DT datatable
          DT_Dat_Goup = format_GOdat_for_DTtable(results_df_up_mf %>% filter(plus_minus=="+"),
                                                 results_df_up_bp %>% filter(plus_minus=="+"),
                                                 table1N2_dedup_up,table1N2_dedup_ALL)
          
          output$DTGOUpgenes_caption <- renderText(paste0(Goterm_num_up," enriched GO terms, click on the plus icon for GeneIDs"))
          
          output$DTGOUpgenes <- DT::renderDataTable({
          
                              ## DT datatable
                              DT::datatable(DT_Dat_Goup, callback = callback, escape = -2,
                                            rownames=TRUE,
                                            options = list(lengthMenu = c(20,50,100,500), pageLength = 100,
                                                           columnDefs = list(
                                                             list(visible = FALSE, targets = ncol(DT_Dat_Goup)),
                                                             list(orderable = FALSE, className = 'details-control', targets = 1),
                                                             list(className = "dt-center", targets = "_all")
                                                           )
                                            )
                              ) %>%
                                formatRound(c("p-val",	"Adjusted.p-val",	"Fold.enrichment","Expected","Hit_perc"),3) 
          })
        }else{
          output$DTGOUpgenes <- DT::renderDataTable({
            
                  ## DT datatable
                  DT::datatable(data.frame(Warning="There are no significantly enriched GO terms"))
            })
        }
        ########GO data down genes########
        Goterm_num_down_mf=nrow(results_df_down_mf %>% filter(plus_minus=="+"))
        Goterm_num_down_bp=nrow(results_df_down_bp %>% filter(plus_minus=="+"))
        Goterm_num_down = Goterm_num_down_mf + Goterm_num_down_bp
        
        if(Goterm_num_down>=1){
          ##get data for DT datatable
          DT_Dat_Godown = format_GOdat_for_DTtable(results_df_down_mf %>% filter(plus_minus=="+"),
                                                   results_df_down_bp %>% filter(plus_minus=="+"),
                                                   table1N2_dedu_down,table1N2_dedup_ALL)
          
          output$DTGODowngenes_caption <- renderText(paste0(Goterm_num_down," enriched GO terms, click on the plus icon for GeneIDs"))
          
          output$DTGODowngenes <- DT::renderDataTable({
            
            ## DT datatable
            DT::datatable(DT_Dat_Godown, callback = callback, escape = -2,
                          rownames=TRUE,
                          options = list(lengthMenu = c(20,50,100,500), pageLength = 100,
                                         columnDefs = list(
                                           list(visible = FALSE, targets = ncol(DT_Dat_Godown)),
                                           list(orderable = FALSE, className = 'details-control', targets = 1),
                                           list(className = "dt-center", targets = "_all")
                                         )
                          )
            ) %>%
              formatRound(c("p-val",	"Adjusted.p-val",	"Fold.enrichment","Expected","Hit_perc"),3) 
          })
        }else{
          output$DTGODowngenes <- DT::renderDataTable({
            
            ## DT datatable
            DT::datatable(data.frame(Warning="There are no significantly enriched GO terms"))
          })
        }
        
        ########GO data sig genes########
        Goterm_num_sig_mf=nrow(results_df_sig_mf %>% filter(plus_minus=="+"))
        Goterm_num_sig_bp=nrow(results_df_sig_bp %>% filter(plus_minus=="+"))
        Goterm_num_sig = Goterm_num_sig_mf + Goterm_num_sig_bp
        
        if(Goterm_num_sig>=1){
          ##get data for DT datatable
          DT_Dat_Goall = format_GOdat_for_DTtable(results_df_sig_mf %>% filter(plus_minus=="+"),
                                                  results_df_sig_bp %>% filter(plus_minus=="+"),
                                                  table1N2_dedup_sig,table1N2_dedup_ALL)
          
          output$DTGOUpDowngenes_caption <- renderText(paste0(Goterm_num_sig," enriched GO terms, click on the plus icon for GeneIDs"))
          
          output$DTGOUpDowngenes <- DT::renderDataTable({
            
            ## DT datatable
            DT::datatable(DT_Dat_Goall, callback = callback, escape = -2,
                          rownames=TRUE,
                          options = list(lengthMenu = c(20,50,100,500), pageLength = 100,
                                         columnDefs = list(
                                           list(visible = FALSE, targets = ncol(DT_Dat_Goall)),
                                           list(orderable = FALSE, className = 'details-control', targets = 1),
                                           list(className = "dt-center", targets = "_all")
                                         )
                          )
            ) %>%
              formatRound(c("p-val",	"Adjusted.p-val",	"Fold.enrichment","Expected","Hit_perc"),3) 
          })
        }else{
          output$DTGOUpDowngenes <- DT::renderDataTable({
            
            ## DT datatable
            DT::datatable(data.frame(Warning="There are no significantly enriched GO terms"))
          })
        }
    }
    else{ #GO enrichment turned off
      
      output$Goallgenesmf <- renderPlot({
        isolate({
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          text(x = 0.5, y = 0.5, paste("GO enrichment analysis turned off\n",
                                       "\n",
                                       ""), 
               cex = 1.6, col = "black")
        })
      })
      
      output$Goupgenesmf <- renderPlot({
        isolate({
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          text(x = 0.5, y = 0.5, paste("GO enrichment analysis turned off\n",
                                       "\n",
                                       ""), 
               cex = 1.6, col = "black")
        })
      })
      
      output$Godowngenesmf <- renderPlot({
        isolate({
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          text(x = 0.5, y = 0.5, paste("GO enrichment analysis turned off\n",
                                       "\n",
                                       ""), 
               cex = 1.6, col = "black")
        })
      })
      output$Goallgenesbp <- renderPlot({
        isolate({
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          text(x = 0.5, y = 0.5, paste("GO enrichment analysis turned off\n",
                                       "\n",
                                       ""), 
               cex = 1.6, col = "black")
        })
      })
      
      output$Goupgenesbp <- renderPlot({
        isolate({
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          text(x = 0.5, y = 0.5, paste("GO enrichment analysis turned off\n",
                                       "\n",
                                       ""), 
               cex = 1.6, col = "black")
        })
      })
      
      output$Godowngenesbp <- renderPlot({
        isolate({
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          text(x = 0.5, y = 0.5, paste("GO enrichment analysis turned off\n",
                                       "\n",
                                       ""), 
               cex = 1.6, col = "black")
        })
      })
    }
  
    #####################
    #GSEA              #
    #####################
#    progress$inc(1/8, detail = paste("Conducting GSEA analysis"))
#    
#    if (input$GSEA_flag==TRUE)
#    {}
#    else{}   
    
    
    
    

    
    
  }) #end of observeEvent(input$updDatSel,{})

  
  
  
  #######################################
  ##define plotting functions          ##
  ####################################### 
  
  output$summary <- renderPlot({ 
    
    ####validate contrast (print warning msg)
    shiny::validate(
      shiny::need(validate_contrast(c(input$tissuec1_Midgut, input$tissuec1_FatBody, input$tissuec1_Ovary),
                                    c(input$tissuec2_Midgut, input$tissuec2_FatBody, input$tissuec2_Ovary)),
                  "Please check your contrast definition: only one sample can be selected on each side of 'vs', and cannot be the same"),
      errorClass = "myclass"
    )
    shiny::validate(
      shiny::need(validate_contrast(c(input$timec2_0, input$timec2_24, input$timec2_36, input$timec2_48),
                                    c(input$timec1_0, input$timec1_24, input$timec1_36, input$timec1_48)),
                  "Please check your contrast definition: only one sample can be selected on each side of 'vs', and cannot be the same"),
      errorClass = "myclass"
    )

  })
  
  output$my_tooltip <- renderUI({
    
    hover <- input$plot_hover
    #print(str(hover))
    req("x" %in% names(hover))
    hover$mapping$x = "Log2foldchange"
    hover$mapping$y = "neg_Log10P"
    #print(hover$x)
    #print(hover$y)
    #y <- nearPoints(volcano_data, hover, threshold = 90, maxpoints = 20)
    #print(y)
    #req(nrow(y) != 0)
    #print(isolate(volcanodata4tooltip$df[which(volcanodata4tooltip$df$geneID=="AGAP004519"),]))
    #print(isolate(volcanodata4tooltip$df[which(volcanodata4tooltip$df$geneID=="AGAP002593"),]))
    verbatimTextOutput("vals")
    
  })
  
  output$vals <- renderPrint({
    hover <- input$plot_hover
    req("x" %in% names(hover))
    hover$mapping$x = "Log2foldchange"
    hover$mapping$y = "neg_Log10P"
    y <- nearPoints(isolate(volcanodata4tooltip$df), hover, threshold = 10, maxpoints = 10)[c("geneID","description")]
    req(nrow(y)>=1)
    #req(nrow(y) != 0)
    print(y, row.names=FALSE)
    #paste0(hover$x," ",hover$y)

  })  
  
  #######################################
  ##define functions                   ##
  #######################################
  
  #validate contrast selection
  validate_contrast <- function(list1,list2)
  {

    if (any(list1)==FALSE || any(list2)==FALSE ) #check if one side of vs is empty
    {
      return(FALSE)
    }
    else if (sum(table(list1)["TRUE"], table(list2)["TRUE"])!=2)
    {
      return(FALSE)
    }
    else{
      for (i in 1:length(list1))
      {
        if(all(sapply(list(list1[i],list2[i]),function(x) x==TRUE))) #check if overlapping conditions are selected
        {
          return(FALSE)
        }
      }
      return(TRUE)
    }

  }
  
  #GO enrichment
  PantherGOenrich <- function(genelist,path,AnnoDataSet)
  {
    #get p-val thres
    GO_p_thres = as.numeric(input$GO_p_thres)
    #send http request
    request <- curl_fetch_memory(paste0(path,"?geneInputList=",genelist,"&annotDataSet=",AnnoDataSet,"&organism=7165"))
    print(paste0("api call response code: ", request$status_code))
    response <- rawToChar(request$content)
    response <- fromJSON(response, flatten = TRUE)
    results <- response$results
    results_df <- results$result  %>%  data.frame()
    
    #filter results
    if (input$GO_multi_correction == TRUE)
    {
      results_df <- results_df %>% filter(fdr<GO_p_thres)
      #print("GO fdr")
    }else
    {
      results_df <- results_df %>% filter(pvalue<GO_p_thres)
      #print("GO pval")
    }
    results_df <- results_df %>% filter(term.label!="UNCLASSIFIED")
    results_df$hit_perc <- results_df$number_in_list/results_df$number_in_reference*100
    
    #string wrap (for plotting)
    results_df$term.label.wrapped <- str_wrap(results_df$term.label,50)
    #print(results_df)
    return(results_df)
  }
  
  cal_plot_height <-function(Goterm_num)
  {
    if ((Goterm_num^(1.05)*38)<350)
    {
      return(350)
    }else{
      return((Goterm_num^(1.05))*40)
    }
  }
  
  #######################################
  ##define function to subset meta data##
  #######################################
  subset_time_contrast <- function(df)
  {
    selection = c(input$timec1_0, input$timec1_24, input$timec1_36, input$timec1_48, input$timec2_0, input$timec2_24, input$timec2_36, input$timec2_48)
    timpoints = c(0,24,36,48,0,24,36,48)
    timpoints_subset = timpoints[selection]
    #print(timpoints_subset)
    df = df[df$Timepoint %in% timpoints_subset,]
    #print(df)
    return(df)
  }
  
  subset_tissue_contrast <- function(df)
  {
    selection = c(input$tissuec1_Midgut, input$tissuec1_FatBody, input$tissuec1_Ovary,input$tissuec2_Midgut, input$tissuec2_FatBody, input$tissuec2_Ovary)
    tissue = c("Midgut", "FatBody", "Ovary", "Midgut", "FatBody", "Ovary")
    tissue_subset = tissue[selection]
    #print(tissue_subset)
    df = df[df$Tissue %in% tissue_subset,]
    #print(df)
    return(df)
  } 
  
  subset_time <- function(df,flag)
  {
    if (flag=="Treatment")
    {   
      selection = c(input$tr_h0,input$tr_h24,input$tr_h36,input$tr_h48)
      timpoints = c(0,24,36,48)
      timpoints_subset = timpoints[selection]
      #print(timpoints_subset)
      df = df[df$Timepoint %in% timpoints_subset,]
      #print(df)
      return(df)
    }
    if (flag=="Tissue")
    {   
      selection = c(input$tis_h0,input$tis_h24,input$tis_h36,input$tis_h48)
      timpoints = c(0,24,36,48)
      timpoints_subset = timpoints[selection]
      #print(timpoints_subset)
      df = df[df$Timepoint %in% timpoints_subset,]
      #print(df)
      return(df)
    }
  }
  
  subset_tissue <- function(df,flag)
  {
    if (flag=="Treatment")
    {
      selection = c(input$tr_Midgut,input$tr_FatBody,input$tr_Ovary)
      tissues = c("Midgut","FatBody","Ovary")
      tissues_subset = tissues[selection]
      #print(tissues_subset)
      df = df[df$Tissue %in% tissues_subset,]
      #print(df)
      return(df)
    }
    if (flag=="Timepoint")
    {
      selection = c(input$time_Midgut,input$time_FatBody,input$time_Ovary)
      tissues = c("Midgut","FatBody","Ovary")
      tissues_subset = tissues[selection]
      #print(tissues_subset)
      df = df[df$Tissue %in% tissues_subset,]
      #print(df)
      return(df)
    }
  }   
  
  subset_treatment <- function(df,flag)
  {
    if (flag=="Tissue")
    {   
      selection = c(input$tis_EcR,input$tis_GFP)
      treatment = c("EcR","GFP")
      treatment_subset = treatment[selection]
      #print(treatment_subset)
      df = df[df$Treatment %in% treatment_subset,]
      #print(df)
      return(df)

    }
    if (flag=="Timepoint")
    {  
      selection = c(input$time_EcR,input$time_GFP)
      treatment = c("EcR","GFP")
      treatment_subset = treatment[selection]
      #print(treatment_subset)
      df = df[df$Treatment %in% treatment_subset,]
      #print(df)
      return(df)    
    }
  }
  #########################################
  ##define a functions used in GO DT table#
  #########################################
  
  #function for fining genes in enriched GO lists
  find_GO_genes_intersect <- function(GOlist, genelist)
  {
    GO_listgenes=list()
    for (GOterm in GOlist)
    {
      Genes = intersect(GO[[GOterm]], genelist)  
      #print(paste0(GOterm," ",length(Genes)))  
      GO_listgenes[[GOterm]]=Genes
    } 
    return(GO_listgenes)
  }
  
  #fucntion for creating DT table for GO data
  format_GOdat_for_DTtable <-function(results_mf, results_bp, table1N2_dedup_x,table1N2_dedup_ALL) # only vary the first 3 parameters
  {
    results_mf$GO.type = rep("Molecular function",length(results_mf$term.id))
    results_bp$GO.type = rep("Biological process",length(results_bp$term.id))  
    GO_results_for_table_mf = results_mf[c("term.id","GO.type","term.label" , "pValue","fdr","fold_enrichment","number_in_list","number_in_reference","expected","hit_perc")]
    GO_results_for_table_bp = results_bp[c("term.id","GO.type","term.label" , "pValue","fdr","fold_enrichment","number_in_list","number_in_reference","expected","hit_perc")]
    
    GO_results_for_table = rbind(GO_results_for_table_mf,GO_results_for_table_bp)
    colnames(GO_results_for_table) = c("GO.ID","Go.type","GO.label" , "p-val","Adjusted.p-val","Fold.enrichment","Number.in.list","Number.in.reference","Expected","Hit_perc")
    #create expandable subtables
    subtable_list = list()
    
    #prep annotations
    annotation_df = table1N2_dedup_ALL[,c("GeneID","description")]
    colnames(annotation_df) = c("GeneID","Description")
    #mf
    GO_listgenes = find_GO_genes_intersect(results_mf$term.id,table1N2_dedup_x$GeneID) #usage: find_GO_genes_intersect(GOlist, genelist)
    for (cur_GOterm in names(GO_listgenes))
    {
      subtable = data.frame(GeneID = GO_listgenes[[cur_GOterm]], stringsAsFactors = FALSE)
      subtable = merge(subtable,annotation_df,by="GeneID", all.x = TRUE) #add annotation
      subtable_list[[length(subtable_list)+1]] = subtable # append current subtable to the list
    }
    #bp
    GO_listgenes = find_GO_genes_intersect(results_bp$term.id,table1N2_dedup_x$GeneID) #usage: find_GO_genes_intersect(GOlist, genelist)
    for (cur_GOterm in names(GO_listgenes))
    {
      subtable = data.frame(GeneID = GO_listgenes[[cur_GOterm]], stringsAsFactors = FALSE)
      subtable = merge(subtable,annotation_df,by="GeneID", all.x = TRUE) #add annotation
      subtable_list[[length(subtable_list)+1]] = subtable # append current subtable to the list
    }
    
    subtables <- lapply(subtable_list, purrr::transpose)
    
    #merge subtables with DT table
    DT_Dat <- cbind(" " = "&oplus;", GO_results_for_table, "_details" = I(subtables))
    
    return(DT_Dat)
  }
  
  #############define a function-like variable used in GO DT table
  ## the callback
  callback = JS(
    "table.column(1).nodes().to$().css({cursor: 'pointer'});",
    "",
    "// make the table header of the nested table",
    "var format = function(d, childId){",
    "  if(d != null){",
    "    var html = ", 
    "      '<table class=\"display compact hover\" id=\"' + childId + '\"><thead><tr>';",
    "    for (var key in d[d.length-1][0]) {",
    "      html += '<th>' + key + '</th>';",
    "    }",
    "    html += '</tr></thead></table>'",
    "    return html;",
    "  } else {",
    "    return '';",
    "  }",
    "};",
    "",
    "// row callback to style the rows of the child tables",
    "var rowCallback = function(row, dat, displayNum, index){",
    "  if($(row).hasClass('odd')){",
    "    $(row).css('background-color', 'papayawhip');",
    "    $(row).hover(function(){",
    "      $(this).css('background-color', '#E6FF99');",
    "    }, function() {",
    "      $(this).css('background-color', 'papayawhip');",
    "    });",
    "  } else {",
    "    $(row).css('background-color', 'lemonchiffon');",
    "    $(row).hover(function(){",
    "      $(this).css('background-color', '#DDFF75');",
    "    }, function() {",
    "      $(this).css('background-color', 'lemonchiffon');",
    "    });",
    "  }",
    "};",
    "",
    "// header callback to style the header of the child tables",
    "var headerCallback = function(thead, data, start, end, display){",
    "  $('th', thead).css({",
    "    'border-top': '3px solid indigo',", 
    "    'color': 'indigo',",
    "    'background-color': '#fadadd'",
    "  });",
    "};",
    "",
    "// make the datatable",
    "var format_datatable = function(d, childId){",
    "  var dataset = [];",
    "  var n = d.length - 1;",
    "  for(var i = 0; i < d[n].length; i++){",
    "    var datarow = $.map(d[n][i], function (value, index) {",
    "      return [value];",
    "    });",
    "    dataset.push(datarow);",
    "  }",
    "  var id = 'table#' + childId;",
    "  if (Object.keys(d[n][0]).indexOf('_details') === -1) {",
    "    var subtable = $(id).DataTable({",
    "                 'data': dataset,",
    "                 'autoWidth': true,",
    "                 'deferRender': true,",
    "                 'info': false,",
    "                 'lengthChange': false,",
    "                 'ordering': d[n].length > 1,",
    "                 'order': [],",
    "                 'paging': false,",
    "                 'scrollX': false,",
    "                 'scrollY': false,",
    "                 'searching': false,",
    "                 'sortClasses': false,",
    "                 'rowCallback': rowCallback,",
    "                 'headerCallback': headerCallback,",
    "                 'columnDefs': [{targets: '_all', className: 'dt-center'}]",
    "               });",
    "  } else {",
    "    var subtable = $(id).DataTable({",
    "            'data': dataset,",
    "            'autoWidth': true,",
    "            'deferRender': true,",
    "            'info': false,",
    "            'lengthChange': false,",
    "            'ordering': d[n].length > 1,",
    "            'order': [],",
    "            'paging': false,",
    "            'scrollX': false,",
    "            'scrollY': false,",
    "            'searching': false,",
    "            'sortClasses': false,",
    "            'rowCallback': rowCallback,",
    "            'headerCallback': headerCallback,",
    "            'columnDefs': [", 
    "              {targets: -1, visible: false},", 
    "              {targets: 0, orderable: false, className: 'details-control'},", 
    "              {targets: '_all', className: 'dt-center'}",
    "             ]",
    "          }).column(0).nodes().to$().css({cursor: 'pointer'});",
    "  }",
    "};",
    "",
    "// display the child table on click",
    "table.on('click', 'td.details-control', function(){",
    "  var tbl = $(this).closest('table'),",
    "      tblId = tbl.attr('id'),",
    "      td = $(this),",
    "      row = $(tbl).DataTable().row(td.closest('tr')),",
    "      rowIdx = row.index();",
    "  if(row.child.isShown()){",
    "    row.child.hide();",
    "    td.html('&oplus;');",
    "  } else {",
    "    var childId = tblId + '-child-' + rowIdx;",
    "    row.child(format(row.data(), childId)).show();",
    "    td.html('&CircleMinus;');",
    "    format_datatable(row.data(), childId);",
    "  }",
    "});")
  
  #########finished defined function
  
}

# Run the application 
shinyApp(ui = ui, server = server)
