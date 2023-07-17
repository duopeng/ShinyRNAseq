#Shiny app for the EcR RNA seq dataset (Catteruccia lab)
#This tool can retrieve FPKM and TPM data and make plots
#version 3 (2020-5-25)
#by Duo Peng


library(shiny)
library(ggplot2)
library(ggbeeswarm)

df <- data.frame(x=1:8, y=1, col=letters[1:8])
g <- ggplot(df, aes(x=x, y=y, color=col)) + geom_point(size=5) +
  scale_color_brewer(palette="Set1")
colors <- ggplot_build(g)$data[[1]]$colour

##################
##read count data#
##################

#metadata
#########
metadata<-read.delim("data/metadata.tab",header=T,stringsAsFactors = FALSE)
#fix variable name of randomized sample#
metadata$Randomized_Num = paste0("X",metadata$Randomized_Num)

#FPKM
#####
FPKMtabfile="data/stringtie_abundance_Agambiae_parsed_FPKM.tab"
dat_FPKM<-read.delim(FPKMtabfile,header=T,stringsAsFactors = FALSE)


#TPM
####
TPMtabfile="data/stringtie_abundance_Agambiae_parsed_TPM.tab"
dat_TPM<-read.delim(TPMtabfile,header=T,stringsAsFactors = FALSE)

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
                         fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                                  h4("check to select tissue(s)"),
                                  checkboxInput(inputId="tr_Midgut", label = "Midgut", value = TRUE ),
                                  checkboxInput(inputId="tr_FatBody", label = "FatBody", value = TRUE ),
                                  checkboxInput(inputId="tr_Ovary", label = "Ovary", value = TRUE),
                         ),
                         fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                                  h4("check to select timpoint(s)"),
                                  checkboxInput(inputId="tr_h0", label = "0hr", value = TRUE ),
                                  checkboxInput(inputId="tr_h24", label = "24hr", value = TRUE ),
                                  checkboxInput(inputId="tr_h36", label = "36hr", value = TRUE),
                                  checkboxInput(inputId="tr_h48", label = "48hr", value = TRUE),
                         )
                ),
                tabPanel("Tissue", 
                         fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                                  h4("check to select timpoint(s)"),
                                  checkboxInput(inputId="tis_h0", label = "0hr", value = FALSE ),
                                  checkboxInput(inputId="tis_h24", label = "24hr", value = FALSE ),
                                  checkboxInput(inputId="tis_h36", label = "36hr", value = TRUE),
                                  checkboxInput(inputId="tis_h48", label = "48hr", value = TRUE),
                         ),
                         fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                                                         h4("check to select treatment(s)"),
                                                         checkboxInput(inputId="tis_EcR", label = "dsEcR", value = TRUE ),
                                                         checkboxInput(inputId="tis_GFP", label = "dsGFP", value = TRUE ),
                         )
                ),
                tabPanel("Time",
                         fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                                                         h4("check to select treatment(s)"),
                                                         checkboxInput(inputId="time_EcR", label = "dsEcR", value = TRUE ),
                                                         checkboxInput(inputId="time_GFP", label = "dsGFP", value = FALSE ),
                         ),
                         fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                                  h4("check to select tissue(s)"),
                                  checkboxInput(inputId="time_Midgut", label = "Midgut", value = FALSE ),
                                  checkboxInput(inputId="time_FatBody", label = "FatBody", value = TRUE ),
                                  checkboxInput(inputId="time_Ovary", label = "Ovary", value = FALSE),
                         )
                )
    )
)

ui <- fluidPage(
  
    tags$head(
      tags$style(HTML("
        .shiny-output-error-myclass {
          color: red;
          font-size: 30px;
          
        }
      "))
    ),
  
    sidebarLayout(
         sidebarPanel(
            h3("Please select various options and click on the update button to generate plots"),
            selectInput("constrast", "Please select a contrast", 
                        choices = c("dsRNA Treatment", "Tissue", "Time")
            ),
            parameter_tabs,
            fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                     h4("gene to plot (AGAP geneID)"),
                     textInput("geneID", "geneID", value="AGAP029539", width = "200px" )
            ),
            fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;color:#333333;background-color:#f5f5f5;border:1px solid #cccccc;border-radius:4px',
                     h4("Plotting parameters"),
                     splitLayout(cellWidths = c("33%", "33%","33%"), 
                         textInput("dotsize", "point size", value="4", width = "100px" ),
                         textInput("plotwidth", "plot width", value="1000", width = "100px" ),
                         textInput("plotheight", "plot height", value="600", width = "100px" )
                     ),

                     selectInput("color", "Points colored by:", 
                                 choices = c( "Tissue","Timepoint","Rep_Num", "Treatment")),   
                     selectInput("shape", "Points shaped by:", 
                                 choices = c("Timepoint", "Rep_Num", "Tissue","Treatment" )),


                                             ),
            fluidRow(style='display:block;padding:9.5px;margin:0 0 10px;font-size:13px;line-height:1.42857143;word-break:normal;word-wrap:break-word;',
                     actionButton("updDatSel", "Update")
            )
        ),
        mainPanel(

            fluidRow(
                tabsetPanel(type = "tabs",
                            tabPanel("FPKM", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'FPKMplot')
                                     )), 
                            ),
                            tabPanel("FPKM data", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, tableOutput('FPKMtable')
                                     )),
                            ),
                            tabPanel("TPM", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, plotOutput(outputId = 'TPMplot' )
                                     ))
                            ),
                            tabPanel("TPM data", 
                                     fluidRow(column(1,offset=0),column(12,offset=0, tableOutput(outputId = 'TPMtable' )
                                     ))
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
# Shiny server: Server      #
#############################

server <- function(input, output,session) {

  axis_text_size=20
  axix_title_size=20
  title_size=20
  legend_title_size=axix_title_size
  legend_text_size=axis_text_size
  
  
    observeEvent(input$constrast, {
        updateTabsetPanel(session, "params", selected = input$constrast)
    }) 
    
    observeEvent(input$updDatSel, {
        #print(input$constrast)
        #print(input$color)
        #print(FPKM_GOI_metadata)
        #print("begin subset")
      

        
        #check if gene ID is in the table (the first validate prints the message, the second validate stop the code below from executing)
        output$FPKMplot <- renderPlot({ 
          shiny::validate(
            shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
            errorClass = "myclass"
          )
        })
        
        shiny::validate(
          shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
          errorClass = "myclass"
        )
      
        if (input$constrast=="dsRNA Treatment")
        {
            #subsetdata
            ###########
            df_list = subset_gene()
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df_FPKM)
            df_list = subset_time(df_FPKM,df_TPM,"Treatment")
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df_FPKM)
            df_list = subset_tissue(df_FPKM,df_TPM,"Treatment")
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df_FPKM)
            
            #generate output plots and tables
            #################################
            
            #FPKM
            output$FPKMplot <- renderPlot({ 
              
              #check if gene ID is in the table
              shiny::validate(
                shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
                errorClass = "myclass"
              )
              
                isolate({
                  df_FPKM$FPKM = as.double(as.character(df_FPKM$FPKM))
                  df_FPKM$Rep_Num = as.factor(df_FPKM$Rep_Num)
                    ggplot(df_FPKM, aes_string(x="Treatment", y="FPKM")) +
                      geom_violin(trim=FALSE, adjust = 1)+
                      geom_quasirandom(aes_string(color = as.name(input$color), shape = as.name(input$shape), size =as.integer(input$dotsize)))+
                        theme(
                            axis.title.x=element_text(size=axix_title_size, margin = margin(t = 15, r = 0, b = 0, l = 0)),
                            axis.title.y=element_text(size=axix_title_size), 
                            axis.text.x=element_text(size=axix_title_size),
                            axis.text.y=element_text(size=axix_title_size),
                            plot.title = element_text(hjust=0.5, size=title_size),
                            legend.text=element_text(size=legend_text_size),
                            legend.title=element_text(size=legend_title_size),
                        )+
                      scale_color_manual(values=c( "red","black", colors[-1]))+
                      #scale_colour_discrete(size = as.integer(input$dotsize))
                      guides(shape = guide_legend(override.aes = list(size = as.integer(input$dotsize))),
                             color = guide_legend(override.aes = list(size = as.integer(input$dotsize)))
                             )+
                      scale_size(guide= FALSE)+
                      ggtitle(input$geneID)
                })
            },
            height = as.numeric(input$plotheight), width = as.numeric(input$plotheight))
            
            df2_FPKM=df_FPKM
            df2_FPKM$Randomized_Num = sub('.', '', df2_FPKM$Randomized_Num)
            output$FPKMtable <- renderTable(df2_FPKM)
            
            #TPM
            output$TPMplot <- renderPlot({ 
              
              #check if gene ID is in the table
              shiny::validate(
                shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
                errorClass = "myclass"
              )
              
              isolate({
                df_TPM$TPM = as.double(as.character(df_TPM$TPM))
                df_TPM$Rep_Num = as.factor(df_TPM$Rep_Num)
                ggplot(df_TPM, aes_string(x="Treatment", y="TPM")
                ) +
                  geom_violin(trim=FALSE, adjust = 1)+
                  geom_quasirandom(aes_string(color = as.name(input$color), shape = as.name(input$shape), size =as.integer(input$dotsize)))+
                  theme(
                    axis.title.x=element_text(size=axix_title_size, margin = margin(t = 15, r = 0, b = 0, l = 0)),
                    axis.title.y=element_text(size=axix_title_size), 
                    axis.text.x=element_text(size=axix_title_size),
                    axis.text.y=element_text(size=axix_title_size),
                    plot.title = element_text(hjust=0.5, size=title_size),
                    legend.text=element_text(size=legend_text_size),
                    legend.title=element_text(size=legend_title_size),
                  )+
                  #scale_colour_discrete(size = as.integer(input$dotsize))
                  guides(shape = guide_legend(override.aes = list(size = as.integer(input$dotsize))),
                         color = guide_legend(override.aes = list(size = as.integer(input$dotsize)))
                  )+
                  scale_size(guide= FALSE)+
                  ggtitle(input$geneID)
              })
            },
            height = as.numeric(input$plotheight), width = as.numeric(input$plotheight))
            
            df2_TPM=df_TPM
            df2_TPM$Randomized_Num = sub('.', '', df2_TPM$Randomized_Num)
            output$TPMtable <- renderTable(df2_TPM)
            
        }
        if (input$constrast=="Tissue")
        {
            #subsetdata
            ###########
            df_list = subset_gene()
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df)
            df_list = subset_treatment(df_FPKM,df_TPM,"Tissue")
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df)
            df_list = subset_time(df_FPKM,df_TPM,"Tissue")
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df)
            
            #generate output plots and tables
            #################################
            
            #FPKM
            output$FPKMplot <- renderPlot({ 
              
              #check if gene ID is in the table
              shiny::validate(
                shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
                errorClass = "myclass"
              )
              
              isolate({
                df_FPKM$FPKM = as.double(as.character(df_FPKM$FPKM))
                df_FPKM$Rep_Num = as.factor(df_FPKM$Rep_Num)
                ggplot(df_FPKM, aes_string(x="Tissue", y="FPKM")
                ) +
                  geom_violin(trim=FALSE, adjust = 1)+
                  geom_quasirandom(aes_string(color = as.name(input$color), shape = as.name(input$shape), size =as.integer(input$dotsize)))+
                  theme(
                    axis.title.x=element_text(size=axix_title_size, margin = margin(t = 15, r = 0, b = 0, l = 0)),
                    axis.title.y=element_text(size=axix_title_size), 
                    axis.text.x=element_text(size=axix_title_size),
                    axis.text.y=element_text(size=axix_title_size),
                    plot.title = element_text(hjust=0.5, size=title_size),
                    legend.text=element_text(size=legend_text_size),
                    legend.title=element_text(size=legend_title_size),
                  )+
                  #scale_colour_discrete(size = as.integer(input$dotsize))
                  guides(shape = guide_legend(override.aes = list(size = as.integer(input$dotsize))),
                         color = guide_legend(override.aes = list(size = as.integer(input$dotsize)))
                  )+
                  scale_size(guide= FALSE)+
                  ggtitle(input$geneID)
              })
            },
            height = as.numeric(input$plotheight), width = as.numeric(input$plotheight))
            
            df_FPKM2=df_FPKM
            df_FPKM2$Randomized_Num = sub('.', '', df_FPKM2$Randomized_Num)
            output$FPKMtable <- renderTable(df_FPKM2)            
            
            #TPM
            output$TPMplot <- renderPlot({ 
              
              #check if gene ID is in the table
              shiny::validate(
                shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
                errorClass = "myclass"
              )
              
              isolate({
                df_TPM$TPM = as.double(as.character(df_TPM$TPM))
                df_TPM$Rep_Num = as.factor(df_TPM$Rep_Num)
                ggplot(df_TPM, aes_string(x="Tissue", y="TPM", )
                ) +
                  geom_violin(trim=FALSE, adjust = 1)+
                  geom_quasirandom(aes_string(color = as.name(input$color), shape = as.name(input$shape), size =as.integer(input$dotsize)))+
                  theme(
                    axis.title.x=element_text(size=axix_title_size, margin = margin(t = 15, r = 0, b = 0, l = 0)),
                    axis.title.y=element_text(size=axix_title_size), 
                    axis.text.x=element_text(size=axix_title_size),
                    axis.text.y=element_text(size=axix_title_size),
                    plot.title = element_text(hjust=0.5, size=title_size),
                    legend.text=element_text(size=legend_text_size),
                    legend.title=element_text(size=legend_title_size),
                  )+
                  #scale_colour_discrete(size = as.integer(input$dotsize))
                  guides(shape = guide_legend(override.aes = list(size = as.integer(input$dotsize))),
                         color = guide_legend(override.aes = list(size = as.integer(input$dotsize)))
                  )+
                  scale_size(guide= FALSE)+
                  ggtitle(input$geneID)
              })
            },
            height = as.numeric(input$plotheight), width = as.numeric(input$plotheight))
            
            df_TPM2=df_TPM
            df_TPM2$Randomized_Num = sub('.', '', df_TPM2$Randomized_Num)
            output$TPMtable <- renderTable(df_TPM2)       
            
        }
        if (input$constrast=="Time")
        {
            #subsetdata
            ###########
            df_list = subset_gene()
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df_FPKM)
            df_list = subset_treatment(df_FPKM,df_TPM,"Time")
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df)
            df_list = subset_tissue(df_FPKM,df_TPM,"Time")
            df_FPKM = df_list[[1]]
            df_TPM = df_list[[2]]
            #print(df)
 
            #generate output plots and tables
            #################################

            ##FPKM
            output$FPKMplot <- renderPlot({ 
              
              #check if gene ID is in the table
              shiny::validate(
                shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
                errorClass = "myclass"
              )
              
              isolate({
                df_FPKM$FPKM = as.double(as.character(df_FPKM$FPKM))
                df_FPKM$Rep_Num = as.factor(df_FPKM$Rep_Num)
                ggplot(df_FPKM, aes_string(x="Timepoint", y="FPKM")
                ) +
                  geom_violin(trim=FALSE, adjust = 1)+
                  geom_quasirandom(aes_string(color = as.name(input$color), shape = as.name(input$shape), size =as.integer(input$dotsize)))+
                  theme(
                    axis.title.x=element_text(size=axix_title_size, margin = margin(t = 15, r = 0, b = 0, l = 0)),
                    axis.title.y=element_text(size=axix_title_size), 
                    axis.text.x=element_text(size=axix_title_size),
                    axis.text.y=element_text(size=axix_title_size),
                    plot.title = element_text(hjust=0.5, size=title_size),
                    legend.text=element_text(size=legend_text_size),
                    legend.title=element_text(size=legend_title_size),
                  )+
                  
                  #scale_colour_discrete(size = as.integer(input$dotsize))
                  guides(shape = guide_legend(override.aes = list(size = as.integer(input$dotsize))),
                         color = guide_legend(override.aes = list(size = as.integer(input$dotsize)))
                  )+
                  scale_size(guide= FALSE)+
                  ggtitle(input$geneID)
              })
            },
            height = as.numeric(input$plotheight), width = as.numeric(input$plotheight))
            
            df_FPKM2=df_FPKM
            df_FPKM2$Randomized_Num = sub('.', '', df_FPKM2$Randomized_Num)
            output$FPKMtable <- renderTable(df_FPKM2)  
            
            ##TPM
            output$TPMplot <- renderPlot({ 
              isolate({
                
                #check if gene ID is in the table
                shiny::validate(
                  shiny::need(input$geneID %in% dat_TPM$gene_ID, "Please check your gene ID, it is invalid or obsolete"),
                  errorClass = "myclass"
                )
                
                df_TPM$TPM = as.double(as.character(df_TPM$TPM))
                df_TPM$Rep_Num = as.factor(df_TPM$Rep_Num)
                ggplot(df_TPM, aes_string(x="Timepoint", y="TPM")
                ) +
                  geom_violin(trim=FALSE, adjust = 1)+
                  geom_quasirandom(aes_string(color = as.name(input$color), shape = as.name(input$shape), size =as.integer(input$dotsize)))+
                  theme(
                    axis.title.x=element_text(size=axix_title_size, margin = margin(t = 15, r = 0, b = 0, l = 0)),
                    axis.title.y=element_text(size=axix_title_size), 
                    axis.text.x=element_text(size=axix_title_size),
                    axis.text.y=element_text(size=axix_title_size),
                    plot.title = element_text(hjust=0.5, size=title_size),
                    legend.text=element_text(size=legend_text_size),
                    legend.title=element_text(size=legend_title_size),
                  )+
                  #scale_colour_discrete(size = as.integer(input$dotsize))
                  guides(shape = guide_legend(override.aes = list(size = as.integer(input$dotsize))),
                         color = guide_legend(override.aes = list(size = as.integer(input$dotsize)))
                  )+
                  scale_size(guide= FALSE)+
                  ggtitle(input$geneID)
              })
            },
            height = as.numeric(input$plotheight), width = as.numeric(input$plotheight))
            
            df_TPM2=df_TPM
            df_TPM2$Randomized_Num = sub('.', '', df_TPM2$Randomized_Num)
            output$TPMtable <- renderTable(df_TPM2)  
        }
        

    })
    
    #######################################
    ##define function to subset FPKM data##
    #######################################
    
    #subset gene
    subset_gene <- function()
    {
        #print(dat_FPKM)
        FPKM_GOI = dat_FPKM[dat_FPKM$gene_ID==input$geneID,]
        #print(FPKM_GOI)
        t_FPKM_GOI = data.frame(FPKM = t(FPKM_GOI)[-1,],
                            Randomized_Num = names(t(FPKM_GOI)[-1,]) ) # remove gene name row , then transpose
        df_FPKM = merge(metadata,t_FPKM_GOI, by = "Randomized_Num")
        df_FPKM$Timepoint = as.character(df_FPKM$Timepoint)
        #print(df_FPKM)
        TPM_GOI = dat_TPM[dat_TPM$gene_ID==input$geneID,]
        t_TPM_GOI = data.frame(TPM = t(TPM_GOI)[-1,],
                                Randomized_Num = names(t(TPM_GOI)[-1,]) ) # remove gene name row , then transpose
        df_TPM = merge(metadata,t_TPM_GOI, by = "Randomized_Num")
        df_TPM$Timepoint = as.character(df_TPM$Timepoint)
        
        return(list(df_FPKM,df_TPM))
    }
      
    subset_time <- function(df_FPKM,df_TPM,flag)
    {
        if (flag=="Treatment")
        {   
            selection = c(input$tr_h0,input$tr_h24,input$tr_h36,input$tr_h48)
            timpoints = c(0,24,36,48)
            timpoints_subset = timpoints[selection]
            #print(timpoints_subset)
            index = df_FPKM$Timepoint %in% timpoints_subset
            df_FPKM = df_FPKM[index,]
            
            selection = c(input$tr_h0,input$tr_h24,input$tr_h36,input$tr_h48)
            timpoints = c(0,24,36,48)
            timpoints_subset = timpoints[selection]
            #print(timpoints_subset)
            index = df_TPM$Timepoint %in% timpoints_subset
            df_TPM = df_TPM[index,]
        }
        if (flag=="Tissue")
        {   
            selection = c(input$tis_h0,input$tis_h24,input$tis_h36,input$tis_h48)
            timpoints = c(0,24,36,48)
            timpoints_subset = timpoints[selection]
            #print(timpoints_subset)
            index = df_FPKM$Timepoint %in% timpoints_subset
            df_FPKM = df_FPKM[index,]
            
            selection = c(input$tis_h0,input$tis_h24,input$tis_h36,input$tis_h48)
            timpoints = c(0,24,36,48)
            timpoints_subset = timpoints[selection]
            #print(timpoints_subset)
            index = df_TPM$Timepoint %in% timpoints_subset
            df_TPM = df_TPM[index,]
        }        
      return(list(df_FPKM, df_TPM))
    }
    
    subset_tissue <- function(df_FPKM,df_TPM,flag)
    {
        if (flag=="Treatment")
        {
            selection = c(input$tr_Midgut,input$tr_FatBody,input$tr_Ovary)
            tissues = c("Midgut","FatBody","Ovary")
            tissues_subset = tissues[selection]
            #print(tissues_subset)
            index = df_FPKM$Tissue %in% tissues_subset
            df_FPKM = df_FPKM[index,]
            
            selection = c(input$tr_Midgut,input$tr_FatBody,input$tr_Ovary)
            tissues = c("Midgut","FatBody","Ovary")
            tissues_subset = tissues[selection]
            #print(tissues_subset)
            index = df_TPM$Tissue %in% tissues_subset
            df_TPM = df_TPM[index,]
        }
        if (flag=="Time")
        {
            selection = c(input$time_Midgut,input$time_FatBody,input$time_Ovary)
            tissues = c("Midgut","FatBody","Ovary")
            tissues_subset = tissues[selection]
            #print(tissues_subset)
            index = df_FPKM$Tissue %in% tissues_subset
            df_FPKM = df_FPKM[index,]
            
            selection = c(input$time_Midgut,input$time_FatBody,input$time_Ovary)
            tissues = c("Midgut","FatBody","Ovary")
            tissues_subset = tissues[selection]
            #print(tissues_subset)
            index = df_TPM$Tissue %in% tissues_subset
            df_TPM = df_TPM[index,]
        }
      return(list(df_FPKM, df_TPM))
    }
    
    subset_treatment <- function(df_FPKM,df_TPM,flag)
    {
        if (flag=="Tissue")
        {   
            selection = c(input$tis_EcR,input$tis_GFP)
            treatment = c("EcR","GFP")
            treatment_subset = treatment[selection]
            #print(treatment_subset)
            index = df_FPKM$Treatment %in% treatment_subset
            df_FPKM = df_FPKM[index,]
            
            selection = c(input$tis_EcR,input$tis_GFP)
            treatment = c("EcR","GFP")
            treatment_subset = treatment[selection]
            #print(treatment_subset)
            index = df_TPM$Treatment %in% treatment_subset
            df_TPM = df_TPM[index,]
        }
        if (flag=="Time")
        {  
          selection = c(input$time_EcR,input$time_GFP)
          treatment = c("EcR","GFP")
          treatment_subset = treatment[selection]
          #print(treatment_subset)
          index = df_FPKM$Treatment %in% treatment_subset
          df_FPKM = df_FPKM[index,]
   
          selection = c(input$time_EcR,input$time_GFP)
          treatment = c("EcR","GFP")
          treatment_subset = treatment[selection]
          #print(treatment_subset)
          index = df_TPM$Treatment %in% treatment_subset
          df_TPM = df_TPM[index,]       
        }
      return(list(df_FPKM, df_TPM))
    }
    
}

# Run the application 
shinyApp(ui = ui, server = server)






