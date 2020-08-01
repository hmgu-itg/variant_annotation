library(shiny)
library(shinyjs)
library(jsonlite)

ui <- fluidPage(
    useShinyjs(),
    sidebarLayout(
        sidebarPanel(
            fileInput("file1", "Choose json.gz file",accept = c("application/gzip",".json.gz")),
            selectInput("mapping_selector","Mappings", choices=c()),
            width=3
            ## tags$hr(),
            ## checkboxInput("header", "Header", TRUE)
        ),
        mainPanel(
            tags$div(id="main_div",width="100%",
            actionButton("variant0_btn", "Variant info",width="100%"),
            hidden(tags$div(id="variant0_div",align="center",tableOutput("variant0"))),

            actionButton("gnomad0_btn", "GnomAD",width="100%"),
            #hidden(tags$div(id="gnomad0_div",align="center",tableOutput("gnomad0"))),
            hidden(tags$div(id="gnomad0_div",align="center",plotlyOutput("gnomad0"))),
            
            actionButton("regulation0_btn","ENSEMBL Regulation",width="100%"),
            hidden(tags$div(id="regulation0_div",align="center",tableOutput("regulation0"))),
            
            actionButton("gwas0_btn","GWAS signals nearby",width="100%"),
            hidden(tags$div(id="gwas0_div",align="center",tableOutput("gwas0"))),
            
            actionButton("vep0_btn","VEP annotation",width="100%"),
            hidden(tags$div(id="vep0_div",align="center",tableOutput("vep0"))),
            
            actionButton("population0_btn","1KG",width="100%"),
            #hidden(tags$div(id="population0_div",align="center",tableOutput("population0")))
            hidden(tags$div(id="population0_div",align="center",plotlyOutput("population0")))
            )
        )
    )
)

server <- function(input, output,session) {
    json_data<-reactive({
        fname<-input$file1
        if (is.null(fname))
            return(NULL)

        data <- stream_in(file(fname$datapath))
        
        if (nrow(as.data.frame(lapply(data$regulation_table0,function(x) fromJSON(x, flatten=T))))==0){
            disable("regulation0_btn")
            }
        
        if (nrow(as.data.frame(lapply(data$gnomad_table0,function(x) fromJSON(x, flatten=T))))==0){
            disable("gnomad0_btn")
            }
        
        if (nrow(as.data.frame(lapply(data$gwas_table0,function(x) fromJSON(x, flatten=T))))==0){
            disable("gwas0_btn")
            }
        
        if (nrow(as.data.frame(lapply(data$vep_table0,function(x) fromJSON(x, flatten=T))))==0){
            disable("vep0_btn")
            }
        
        if (nrow(as.data.frame(lapply(data$population_table0,function(x) fromJSON(x, flatten=T))))==0){
            disable("population0_btn")
            }

        data
    })
    
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
    
    ## output$gnomad0<-renderTable({
    ##     fname<-input$file1
    ##     if (is.null(fname))
    ##         return(NULL)

    ##     data<-json_data()
    ##     as.data.frame(lapply(data$gnomad_table0,function(x) fromJSON(x, flatten=T)))
    ## })

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

    ## output$population0<-renderTable({
    ##     fname<-input$file1
    ##     if (is.null(fname))
    ##         return(NULL)

    ##     data <- json_data()
    ##     as.data.frame(lapply(data$population_table0,function(x) fromJSON(x, flatten=T)))
    ## })

    observe({updateSelectInput(session,"mapping_selector",label="Mappings",choices=unlist(lapply(lapply(json_data()[,grepl("variant_table",colnames(json_data()))==T],function(x) fromJSON(x, flatten=T)),function(x) {z<-as.data.frame(x);z[z["Key"]=="Location","Value"]}),use.names=F))})
    
}

shinyApp(ui, server)

