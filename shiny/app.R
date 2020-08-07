library(shiny)
library(plotly)
library(shinyjs)
library(jsonlite)

ui <- fluidPage(
    useShinyjs(),
    tabsetPanel(type="tabs",
                tabPanel("Variant",
                         fluidRow(style = "border-style: solid;border-width:  1px; border-color: lightblue;",
                             column(3,
                                    fileInput("file1", "Choose json.gz file",accept = c("application/gzip",".json.gz"))
                                    ),
                             column(3, offset = 1,
                                    selectInput("mapping_selector","Mappings", choices=c())
                                    )
                         ),                         
                         tags$div(id="main_div",width="100%",
                                  actionButton("variant_btn", "Variant info",width="100%"),
                                  hidden(tags$div(id="variant_div",align="center",tableOutput("variant"))),

                                  actionButton("gnomad_btn", "GnomAD",width="100%"),
                                  hidden(tags$div(id="gnomad_div",align="center",plotlyOutput("gnomad"))),
                                  
                                  actionButton("regulation_btn","ENSEMBL Regulation",width="100%"),
                                  hidden(tags$div(id="regulation_div",align="center",tableOutput("regulation"))),
                                  
                                  actionButton("gwas_btn","GWAS signals nearby",width="100%"),
                                  hidden(tags$div(id="gwas_div",align="center",tableOutput("gwas"))),
                                  
                                  actionButton("vep_btn","VEP annotation",width="100%"),
                                  hidden(tags$div(id="vep_div",align="center",tableOutput("vep"))),
                                  
                                  actionButton("population_btn","1KG",width="100%"),
                                  hidden(tags$div(id="population_div",align="center",plotlyOutput("population"))),
                                  
                                  actionButton("pubmed_btn","PubMed",width="100%"),
                                  hidden(tags$div(id="pubmed_div",align="center",tableOutput("pubmed"))),

                                  actionButton("phenotype_btn","Nearby variants with phenotypes",width="100%"),
                                  hidden(tags$div(id="phenotype_div",align="center",tableOutput("phenotype"))),

                                  actionButton("gene_btn","Nearby genes",width="100%"),
                                  hidden(tags$div(id="gene_div",align="center",tableOutput("gene"))),

                                  actionButton("gtex_genes_btn","GTEx",width="100%"),
                                  hidden(tags$div(id="gtex_genes_div",align="center",tableOutput("gtex_genes")))
                                  )
                         ),
                tabPanel("Genes",
                         sidebarPanel(
                             selectInput("gene_selector","Genes", choices=c()),
                             width=3
                         ),
                         mainPanel(
                             tags$div(id="main_genes_div",width="100%",
                                          actionButton("gene2_btn", "Gene info",width="100%"),
                                          hidden(tags$div(id="gene2_div",align="center",tableOutput("gene2"))),

                                          actionButton("uniprot_btn", "UniProt",width="100%"),
                                          hidden(tags$div(id="uniprot_div",align="center",tableOutput("uniprot"))),
                                          
                                          actionButton("gwas2_btn","GWAS signals",width="100%"),
                                          hidden(tags$div(id="gwas2_div",align="center",tableOutput("gwas2"))),
                                          
                                          actionButton("gxa_btn","Gene expression",width="100%"),
                                          hidden(tags$div(id="gxa_div",align="center",plotlyOutput("gxa"))),
                                          
                                          actionButton("gtex_variants_btn","GTEx eQTLs",width="100%"),
                                          hidden(tags$div(id="gtex_variants_div",align="center",tableOutput("gtex_variants"))),
                                          
                                          actionButton("mouse_btn","Mouse homologs",width="100%"),
                                          hidden(tags$div(id="mouse_div",align="center",tableOutput("mouse"))),
                                          
                                          actionButton("go_btn","GO terms",width="100%"),
                                          hidden(tags$div(id="go_div",align="center",tableOutput("go"))),
                                          )
                             )
                         )
                )
)

server <- function(input, output,session) {

    mapping_lookup <- reactive({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        L2<-unlist(lapply(lapply(json_data()[,grepl("variant_table",colnames(json_data()))==T],function (x) fromJSON(x,flatten=T)),function(x) x[x$Key %in% c("Location"),]["Value"]),use.names=F)
        L1<-unlist(lapply(colnames(json_data())[grepl("variant_table",colnames(json_data()))],function(x) substring(x,14)),use.names=F)
        mapping_lookup <- data.frame(as.character(L2),as.character(L1))
        colnames(mapping_lookup) <- c("mapping","suffix")

        mapping_lookup
    })

    json_data<-reactive({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- stream_in(file(fname$datapath))

        L <- lapply(colnames(data[,grepl("go_table_ENSG",colnames(data))==T]),function (x) substring(x,10))
        updateSelectInput(session,"gene_selector",label="Genes",choices=L,selected=L[1])
        
        data
    })
    
    ##----------------------------------------------------------------------------------------------------
    observeEvent(input$gene2_btn,{
        toggle("gene2_div")
    })
    
    observeEvent(input$uniprot_btn,{
        toggle("uniprot_div")
    })
    
    observeEvent(input$gwas2_btn,{
        toggle("gwas2_div")
    })
    
    observeEvent(input$gtex_variants_btn,{
        toggle("gtex_variants_div")
    })
    
    observeEvent(input$mouse_btn,{
        toggle("mouse_div")
    })
    
    observeEvent(input$go_btn,{
        toggle("go_div")
    })
    
    observeEvent(input$gxa_btn,{
        toggle("gxa_div")
    })
    
    ##----------------------------------------------------------------------------------------------------

    observeEvent(input$variant_btn,{
        toggle("variant_div")
    })
    
    observeEvent(input$gnomad_btn,{
        toggle("gnomad_div")
    })
    
    observeEvent(input$regulation_btn,{
        toggle("regulation_div")
    })
    
    observeEvent(input$gwas_btn,{
        toggle("gwas_div")
    })
    
    observeEvent(input$population_btn,{
        toggle("population_div")
    })
    
    observeEvent(input$vep_btn,{
        toggle("vep_div")
    })
    
    observeEvent(input$pubmed_btn,{
        toggle("pubmed_div")
    })
    
    observeEvent(input$phenotype_btn,{
        toggle("phenotype_div")
    })
    
    observeEvent(input$gene_btn,{
        toggle("gene_div")
    })
    
    observeEvent(input$gtex_genes_btn,{
        toggle("gtex_genes_div")
    })

    ##----------------------------------------------------------------------------------------------------

    output$gene2<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        gene <- input$gene_selector
        tname <- paste("gene2_table_",gene,sep="")
        as.data.frame(lapply(data[,colnames(data) %in% c(tname)],function(x) fromJSON(x, flatten=T)))
    })

    output$uniprot<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        gene <- input$gene_selector
        tname <- paste("uniprot_table_",gene,sep="")
        as.data.frame(lapply(data[,colnames(data) %in% c(tname)],function(x) fromJSON(x, flatten=T)))
    })

    output$gwas2<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        gene <- input$gene_selector
        tname <- paste("gwas2_table_",gene,sep="")
        as.data.frame(lapply(data[,colnames(data) %in% c(tname)],function(x) fromJSON(x, flatten=T)))
    }, sanitize.text.function = function(x) x)

    output$gtex_variants<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        gene <- input$gene_selector
        tname <- paste("gtex_variants_table_",gene,sep="")
        as.data.frame(lapply(data[,colnames(data) %in% c(tname)],function(x) fromJSON(x, flatten=T)))
    })

    output$mouse<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        gene <- input$gene_selector
        tname <- paste("mouse_table_",gene,sep="")
        as.data.frame(lapply(data[,colnames(data) %in% c(tname)],function(x) fromJSON(x, flatten=T)))
    })

    output$go<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        gene <- input$gene_selector
        tname <- paste("go_table_",gene,sep="")
        as.data.frame(lapply(data[,colnames(data) %in% c(tname)],function(x) fromJSON(x, flatten=T)))
    })

    output$gxa <- renderPlotly({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        gene <- input$gene_selector
        tname <- paste("gxa_table_",gene,sep="")
        df <- as.data.frame(lapply(data[,colnames(data) %in% c(tname)],function(x) fromJSON(x, flatten=T)))
        rownames(df)<-df$Experiment
        df2<-df[,!colnames(df) %in% c("Experiment")]
        plot_ly(x=colnames(as.matrix(df2)),y=rownames(as.matrix(df2)),z=as.matrix(df2),type="heatmap")
    })

    ##----------------------------------------------------------------------------------------------------

    output$variant<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        sfx <- unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F)
        tname <- paste0("variant_table",sfx)
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,16))
        df
    })

    output$gnomad <- renderPlotly({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        sfx <- unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F)
        tname <- paste0("gnomad_table",sfx)
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,15))
        df[,!names(df) %in% c("Population")]<-sapply(df[,!names(df) %in% c("Population")],as.numeric)
        plot_ly(x=df$Population,y=df[,2],type ='bar',name=names(df)[2],marker=list(color='red')) %>% add_trace(y = df[,3],name=names(df)[3],marker=list(color='blue')) %>% layout(barmode='stack',yaxis = list(title='AF'),xaxis=list(title="Population",showticklabels=T))
    })
    
    output$population <- renderPlotly({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        sfx <- unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F)
        tname <- paste0("population_table",sfx)        
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,19))
        df[,!names(df) %in% c("Population")]<-sapply(df[,!names(df) %in% c("Population")],as.numeric)
        plot_ly(x=df$Population,y=df[,2],type ='bar',name=names(df)[2],marker=list(color='red')) %>% add_trace(y = df[,3],name=names(df)[3],marker=list(color='blue')) %>% layout(barmode='stack',yaxis = list(title='AF'),xaxis=list(title="Population",showticklabels=T))
    })
    
    output$regulation<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        tname <- paste0("regulation_table",unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F))        
        as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
    })

    output$gwas<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        tname <- paste0("gwas_table",unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F))
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,13))
        df
    })

    output$pubmed<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        tname <- paste0("pubmed_table",unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F))
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,15))
        df
    })

    output$vep<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        tname <- paste0("vep_table",unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F))
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,12))
        df
    })

    output$phenotype<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        tname <- paste0("phenotype_table",unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F))
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,18))
        df
    })

    output$gene<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        tname <- paste0("gene_table",unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F))
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,13))
        df
    })

    output$gtex_genes<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        tname <- paste0("gtex_genes_table",unlist(mapping_lookup()[mapping_lookup()["mapping"]==input$mapping_selector,]["suffix"],use.names=F))
        df <- as.data.frame(lapply(data[tname],function(x) fromJSON(x, flatten=T)))
        colnames(df)<-lapply(colnames(df),function(x) substring(x,19))
        df
    })

# ====================================================== MAPPING SELECTOR =========================================================
    
    observe({        
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        updateSelectInput(session,"mapping_selector",label="Mappings",choices=as.character(mapping_lookup()$mapping))
    })

    observeEvent(input$mapping_selector,{
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        mapping <- input$mapping_selector
        sfx <- unlist(mapping_lookup()[mapping_lookup()["mapping"]==mapping,]["suffix"],use.names=F)
        print(paste(mapping,sfx))

        for (tname in c("variant","regulation","gnomad","gwas","vep","population","pubmed","phenotype","gene","gtex_genes")){
            tname2 <- paste0(tname,"_table",sfx)
            btn_name <- paste0(tname,"_btn")
            print(tname2)
            if (nrow(as.data.frame(lapply(json_data()[tname2],function(x) fromJSON(x, flatten=T))))==0){
                disable(btn_name)
            }
            else{
                enable(btn_name)
            }
        }
    })
# ======================================================== GENE SELECTOR ==========================================================
    
    observeEvent(input$gene_selector,{
        if (nrow(as.data.frame(lapply(json_data()[,colnames(json_data()) %in% c(paste("uniprot_table_",input$gene_selector,sep=""))],function(x) fromJSON(x, flatten=T))))==0)
            disable("uniprot_btn")
        else
            enable("uniprot_btn")

        if (nrow(as.data.frame(lapply(json_data()[,colnames(json_data()) %in% c(paste("gwas2_table_",input$gene_selector,sep=""))],function(x) fromJSON(x, flatten=T))))==0)
            disable("gwas2_btn")
        else
            enable("gwas2_btn")

        if (nrow(as.data.frame(lapply(json_data()[,colnames(json_data()) %in% c(paste("gtex_variants_table_",input$gene_selector,sep=""))],function(x) fromJSON(x, flatten=T))))==0)
            disable("gtex_variants_btn")
        else
            enable("gtex_variants_btn")

        if (nrow(as.data.frame(lapply(json_data()[,colnames(json_data()) %in% c(paste("mouse_table_",input$gene_selector,sep=""))],function(x) fromJSON(x, flatten=T))))==0)
            disable("mouse_btn")
        else
            enable("mouse_btn")

        if (nrow(as.data.frame(lapply(json_data()[,colnames(json_data()) %in% c(paste("go_table_",input$gene_selector,sep=""))],function(x) fromJSON(x, flatten=T))))==0)
            disable("go_btn")
        else
            enable("go_btn")

        if (nrow(as.data.frame(lapply(json_data()[,colnames(json_data()) %in% c(paste("gxa_table_",input$gene_selector,sep=""))],function(x) fromJSON(x, flatten=T))))==0)
            disable("gxa_btn")
        else
            enable("gxa_btn")        
    })
# =============================================================================================================================    

}

shinyApp(ui, server)

