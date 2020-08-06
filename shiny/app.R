library(shiny)
library(plotly)
library(shinyjs)
library(jsonlite)

ui <- fluidPage(
    useShinyjs(),
    tabsetPanel(type="tabs",
                tabPanel("Variant",
                         sidebarLayout(
                             sidebarPanel(
                                 fileInput("file1", "Choose json.gz file",accept = c("application/gzip",".json.gz")),
                                 selectInput("mapping_selector","Mappings", choices=c()),
                                 width=3
                             ),
                             mainPanel(
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

                                          actionButton("phenotype_btn","Variants with phenotypes",width="100%"),
                                          hidden(tags$div(id="phenotype_div",align="center",tableOutput("phenotype"))),

                                          actionButton("gene_btn","Nearby genes",width="100%"),
                                          hidden(tags$div(id="gene_div",align="center",tableOutput("gene"))),

                                          actionButton("gtex_genes_btn","GTEx",width="100%"),
                                          hidden(tags$div(id="gtex_genes_div",align="center",tableOutput("gtex_genes")))
                                          )
                             )
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
        L1<-unlist(lapply(colnames(json_data()[,grepl("variant_table",colnames(json_data()))==T]),function (x) substring(x,14)),use.names=F)
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

    observeEvent(input$variant0_btn,{
        toggle("variant0_div")
    })
    
    observeEvent(input$gnomad0_btn,{
        toggle("gnomad0_div")
    })
    
    observeEvent(input$regulation0_btn,{
        toggle("regulation0_div")
    })
    
    observeEvent(input$gwas0_btn,{
        toggle("gwas0_div")
    })
    
    observeEvent(input$population0_btn,{
        toggle("population0_div")
    })
    
    observeEvent(input$vep0_btn,{
        toggle("vep0_div")
    })
    
    observeEvent(input$pubmed0_btn,{
        toggle("pubmed0_div")
    })
    
    observeEvent(input$phenotype0_btn,{
        toggle("phenotype0_div")
    })
    
    observeEvent(input$gene0_btn,{
        toggle("gene0_div")
    })
    
    observeEvent(input$gtex_genes0_btn,{
        toggle("gtex_genes0_div")
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
    })

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

    output$variant0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        as.data.frame(lapply(data$variant_table0,function(x) fromJSON(x, flatten=T)))
    })

    output$gnomad0 <- renderPlotly({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        df <- as.data.frame(lapply(data$gnomad_table0,function(x) fromJSON(x, flatten=T)))
        df[,!names(df) %in% c("Population")]<-sapply(df[,!names(df) %in% c("Population")],as.numeric)
        plot_ly(x=df$Population,y=df[,2],type ='bar',name=names(df)[2],marker=list(color='red')) %>% add_trace(y = df[,3],name=names(df)[3],marker=list(color='blue')) %>% layout(barmode='stack',yaxis = list(title='AF'),xaxis=list(title="Population",showticklabels=T))
    })
    
    output$population0 <- renderPlotly({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data<-json_data()
        df <- as.data.frame(lapply(data$population_table0,function(x) fromJSON(x, flatten=T)))
        df[,!names(df) %in% c("Population")]<-sapply(df[,!names(df) %in% c("Population")],as.numeric)
        plot_ly(x=df$Population,y=df[,2],type ='bar',name=names(df)[2],marker=list(color='red')) %>% add_trace(y = df[,3],name=names(df)[3],marker=list(color='blue')) %>% layout(barmode='stack',yaxis = list(title='AF'),xaxis=list(title="Population",showticklabels=T))
    })
    
    output$regulation0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        as.data.frame(lapply(data$regulation_table0,function(x) fromJSON(x, flatten=T)))
    })

    output$gwas0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        as.data.frame(lapply(data$gwas_table0,function(x) fromJSON(x, flatten=T)))
    })

    output$vep0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        as.data.frame(lapply(data$vep_table0,function(x) fromJSON(x, flatten=T)))
    })

    output$uniprot0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        as.data.frame(lapply(data$uniprot_table0,function(x) fromJSON(x, flatten=T)))
    })

    output$phenotype0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        as.data.frame(lapply(data$phenotype_table0,function(x) fromJSON(x, flatten=T)))
    })

    output$gene0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        as.data.frame(lapply(data$gene_table0,function(x) fromJSON(x, flatten=T)))
    })

    output$gtex_genes0<-renderTable({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- json_data()
        as.data.frame(lapply(data$gtex_genes_table0,function(x) fromJSON(x, flatten=T)))
    })

# ====================================================== MAPPING SELECTOR =========================================================
    
    observe({        
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        updateSelectInput(session,"mapping_selector",label="Mappings",choices=as.character(mapping_lookup()$mapping))
    })

    observeEvent(input$mapping_selector,{
        mapping <- input$mapping_selector
        sfx <- unlist(mapping_lookup()[mapping_lookup()["mapping"]==mapping,]["suffix"],use.names=F)
        print(paste(mapping,sfx,sep=" "))

        for (tname in c("regulation","gnomad","gwas","vep","population","pubmed","phenotype","gene","gtex")){
            tname2 <- paste(tname,"_table",sfx)
            btn_name <- paste(tname,"_btn")
            
            if (nrow(as.data.frame(lapply(data[tname2],function(x) fromJSON(x, flatten=T))))==0){
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

