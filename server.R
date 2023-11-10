# These options are listed in https://shiny.rstudio.com/reference/shiny/0.14/shiny-options.html
options(repos = BiocManager::repositories())
options(shiny.reactlog=TRUE) 
options(shiny.maxRequestSize=2000*1024^3) #to set the maximum file upload size (2GB)
options(shiny.usecairo=FALSE)
#this should limit the file input to 500 Mbye per file - by default Shiny limits file uploads to 5MB per file
#options(shiny.fullstacktrace = TRUE)
#options(shiny.trace = TRUE)
options(shiny.error = browser)
#options(repos = BiocManager::repositories())
#options(repos = c("BioCsoft" = "https://bioconductor.org/packages/3.10/bioc")) #to be run b4 2 deploy (it was 3.8)
storeWarn<- getOption("warn")
options(warn = -1) #this will suppress warning (they were annoying especially from the plotly functions)
platfrom <- .Platform$OS.type #this is to configure the heavier algorithms which could be multi-threaded

source("helpers.R")

################################################### SERVER FUNCTION -----

server <- function(input, output, session) {

  ################################################### session setup -----    
  onStop(function() cat("Session stopped\n"))
  ###  for the onStop and on onStart see https://shiny.rstudio.com/reference/shiny/1.0.5/onStop.html  
  onStart = function() {
    cat("Doing application setup\n")
    onStop(function() {
      cat("Doing application cleanup\n")})}
  
  app_dir <<- getwd()
  data_dir <<- paste0(app_dir, "/tmpdata")
  temp_dir <<- tempdir()
  del_tmpdata(type = "all") #to delete all data from previous sessions
  
  setwd(temp_dir) 
  old_file <- list.files(path = temp_dir)
  unlink(x = old_file, recursive = TRUE) #to delete all data in the temp dir as well
  setwd(app_dir)
  
  withConsoleRedirect <- function(containerId, expr) {
    # Change type="output" to type="message" to catch stderr
    # (messages, warnings, and errors) instead of stdout.
    txt <- capture.output(results <- expr, type = "message")
    if (length(txt) > 0) {
      insertUI(paste0("#", containerId), where = "beforeEnd",
               ui = paste0(txt, "\n", collapse = "")
      )
    }
    results
  }
  
  hideTab(inputId = "preClusteringWorkflow", target = "Clean")
  hideTab(inputId = "preClusteringWorkflow", target = "Scale")
  hideTab(inputId = "preClusteringWorkflow", target = "Align")
  hideTab(inputId = "preClusteringWorkflow", target = "Downsample")
  hideTab(inputId = "preClusteringWorkflow", target = "Concatenate")
  
  ################################################### process checkbox -----    
  
  output$txt <- renderText({
    icons <- paste(input$process, collapse = ", ")
    paste("You chose", icons)})
  output$checktxt <- renderText({
    validate(need(expr = (!is.null(input$process)), "Please select a module!"))
  })
  
  
  ################################################### hide/show tabs preClustering -----    
  
  observeEvent(input$process, {
    if((!("clean" %in% input$process))){hideTab(inputId = "preClusteringWorkflow", target = "Clean")}
    if((!("scale" %in% input$process))){hideTab(inputId = "preClusteringWorkflow", target = "Scale")}
    if((!("align" %in% input$process))){hideTab(inputId = "preClusteringWorkflow", target = "Align")}
    if((!("downsample" %in% input$process))){hideTab(inputId = "preClusteringWorkflow", target = "Downsample")}
    if((!("concatenate" %in% input$process))){hideTab(inputId = "preClusteringWorkflow", target = "Concatenate")}})
  
  observeEvent(input$process, {
    validate(need(expr = !is.null(input$process), message = "select at least one module"))
    if((("clean" %in% input$process))){showTab(inputId = "preClusteringWorkflow", target = "Clean")}
    if((("scale" %in% input$process))){showTab(inputId = "preClusteringWorkflow", target = "Scale")}
    if((("align" %in% input$process))){showTab(inputId = "preClusteringWorkflow", target = "Align")}
    if((("downsample" %in% input$process))){showTab(inputId = "preClusteringWorkflow", target = "Downsample")}
    if((("concatenate" %in% input$process))){showTab(inputId = "preClusteringWorkflow", target = "Concatenate")}})
  
  ################################################### download_manual -----    
  
  output$download_manual <- downloadHandler(
    filename = function() {
      paste("cytoChain_manual", "pdf", sep=".")
    },
    
    content <- function(file) {
      file.copy("./www/cytoChain_manual 2.1.pdf", file)
    },
    contentType = "application/pdf"
  )
  
  ################################################### __________Loading -----    

  dI <- reactiveVal(NULL) #This will be the entry flowSet. 
  #It will change at the end of this procedure with the adding of the cell Id after cell_Id_fS()
  dI.H <- reactiveVal(NULL) #This is the Handled flowSet at the end of each optimization workflow procedure
  dI_val <- reactiveValues(inFile = NULL, pieColor = NULL, fS_table = NULL, fS_plot = NULL, fS_summary = NULL, fS_dim = NULL, fS_p = NULL)
  
  observeEvent(eventExpr = input$runLoad, {
    
    print("start parsing process")
    output$console_output_pre <- renderPrint("start parsing process  - waiting...")
    options(warn = 0)
    png()
    clean_graph()
    
    dI(NULL); dI.c(NULL); dI.t(NULL); dI.a(NULL); dI.d(NULL); dI.co(NULL); dI.H(NULL)
    dI_val$inFile <- NULL; dI_val$pieColor <- NULL; dI_val$fS_table <- NULL; dI_val$fS_plot <- NULL
    dI_val$fS_dim <- NULL; dI_val$fS_p <- NULL
    del_tmpdata(type = "all")

    if(!is.null(input$fSinput)||!is.null(input$zipinput)) {
      
      if (!is.null(input$fSinput)){
        alphafSinput <-  input$fSinput[order(input$fSinput$name),] 
        lista <- loading_fS(inputfile = alphafSinput$datapath, tipo = "flowSet")
        fS <- lista[[1]]
        nomi_files <- alphafSinput$name
        dI(fS)
        if (inherits(x = fS, "flowSet"))
          {fS <- file_fcs(fS, stringhe = nomi_files)
          dI(fS)}}
      
      if (!is.null(input$zipinput)){
        lista <- loading_fS(inputfile = input$zipinput$datapath, tipo = "zip")
        fS <- lista[[1]]
        nomi_files <- lista[[2]]
        dI(fS)
        if(inherits(x = fS, "flowSet")){
          {fS <- file_fcs(fS, stringhe = nomi_files)
          dI(fS)}}}}
      
      if(inherits(x = dI(), "flowSet")){
        set.seed(input$set_seed)
        
        dI_val$fS_dim <- dim_fS(flow_Set  = dI())
        if(!is.null(input$fSinput)){
          dI_val$inFile <- input$fSinput[order(input$fSinput$name),] #was dI_val$inFile <- input$fSinput
          lista <- table_fS(flow_Set  = dI(), file_name = dI_val$inFile$name)}
        if(!is.null(input$zipinput)){
          dI_val$inFile <- input$zipinput
          lista <- table_fS(flow_Set  = dI(), file_name = nomi_files)}
        dI_val$fS_table <- lista[[1]]
        dI_val$fS_plot <- lista[[2]]
        dI_val$fS_summary <- lista[[3]]
        tabDim <- dI_val$fS_dim
        
        name <- pData(parameters(dI()[[1]]))$name
        del_tmpdata(type = "all")
        if (!('cell_Id' %in% name))
        {fS.Id <- cell_Id_fS(flow_Set = dI(), stringhe = nomi_files)
        dI(fS.Id)}
        
        colorset <- set_colors(length(dI()))
        dI_val$fScolor <- colorset
        color.sample <<- colorset
        
        if(identical(x = color.sample, y = "You exceed the maximum number of samples")){
          dI(color.sample)}
        print("Parsing process ends")
        output$console_output_pre <- renderPrint("parsing process ends")
        output$load_panelStatus <- reactive({input$load_panelStatus=="show"})#??? this could be also used to select the sidebars
        outputOptions(output, "load_panelStatus", suspendWhenHidden = FALSE)
        message("System info:  ", Sys.info())
        if (Sys.info()["sysname"] == "Windows") 
        {play(x = bell)}}# else {js$ding()}} #to try to play the sound on server side. Invain
      else {
        flowSet_output <- dI()
        dI(flowSet_output)}
      
      output$checktxt_parsing_fS <- renderText({
        if (!is.null(dI())){
          if (inherits(x = dI(), "character")) 
          {messaggio <- dI()
          if (grepl("wrong dimension or type", messaggio)) {
            msg <- paste0(messaggio, 
                          " - Maybe there are some inconsistencies in the sample's names. Try to fix it with the 'Panel Editor'")}
          else
          {msg <- messaggio}
          validate(need(expr =  (inherits(x = dI(),"flowSet")), message = msg))}}})
  })
  
  ################################################### runLoadShow ----  
  observeEvent(eventExpr = input$runLoadShow,{
    validate(need(expr = !is.null(dI()), message = "please select some samples"))
    if(inherits(x = dI(),"flowSet")) {
      print("start data visualization of the loaded flowSet - waiting...")
      output$console_output_pre <- renderPrint("start data visualization of the loaded flowSet - waiting...")
      validate(need(expr = (inherits(x = dI(),"flowSet")), message = "please select some samples"))
      set.seed(input$set_seed)
      
      dI_val$pieColor <- pie_color(color.sample)
      dI_val$fS_p <- visualize_rangefS(flow_Set  = dI())
      output$contents <- renderDT(DT::datatable(dI_val$inFile, options = list(pageLength = 20)))
      output$fSsummary <- renderDT(DT::datatable(dI_val$fS_summary, options = list(pageLength = 20))%>% 
                                     DT::formatStyle(columns = c(1,2), `font-size` = '16pt'))
      output$samples <- renderDT(DT::datatable(dI_val$fS_table, options = list(pageLength = 20)))
      output$sample_plot <- renderPlot({dI_val$fS_plot})
      output$markers <- renderDT(DT::datatable(dI_val$fS_dim, options = list(pageLength = 20)))
      output$marker_range <- renderPlotly({dI_val$fS_p})
      output$samples_color <- renderPlot({dI_val$pieColor})
      
      output$entryEvents <- renderValueBox({
        valueBox(value = dI_val$fS_summary[2],subtitle =  "Total entry events", icon = icon("list"), color = "green")})
      print("Visualization process ends")
      output$console_output_pre <- renderPrint("visualization process ends")}
    else{
      result <- "please load your FCS files in order to be parsed from cytoChain"
      dI(result)}
    
    if(!is.null(input$fSinput)){
      table_fS(flow_Set  = dI(), file_name = dI_val$inFile$name)}
    if(!is.null(input$zipinput)){
      nomi_files <- dI_val$fS_table$`sample name`
      table_fS(flow_Set  = dI(), file_name = nomi_files)}
    
    output$downloadfS_table_csv <- downloadHandler(
      filename = function() {paste("fS_table","csv", sep = ".")},
      content <- function(file) {file.copy("./tmpdata/fS_table.csv", file)},
      contentType = "text/csv")
    
    dim_fS(flow_Set  = dI())
    output$downloadfS_dim_csv <- downloadHandler(
      filename = function() {paste("fS_dim","csv", sep = ".")},
      content <- function(file) {file.copy("./tmpdata/fS_dim.csv", file)},
      contentType = "text/csv")
    
    output$checktxt_loading_fS <- renderText({
      if (!is.null(dI())){
        if (inherits(x = dI(),"character")) 
        {messaggio <- dI()
        validate(need(expr =  (inherits(x = dI(),"flowSet")), message = messaggio))}}})
  })
  
  
  ################################################### ____________Clean -----    
  ################################################### runClean  -----    
  dI.c <- reactiveVal(NULL)
  dI.c_val <- reactiveValues(fS_table = NULL, dim_table = NULL, plot_dim_table = NULL, fS_summary = NULL, 
                             plot_comp = NULL, plot_visualize = NULL, QC_table = NULL)
  
  observeEvent(eventExpr = input$runClean, {
  
    del_tmpdata(type = "all")
    dI.c(NULL)
    dI.c_val$fS_table <- NULL; dI.c_val$dim_table <- NULL; dI.c_val$plot_dim_table <- NULL; 
    dI.c_val$fS_summary[2] <- "refresh...";
    dI.c_val$fS_summary[3] <- "refresh...";
    dI.c_val$plot_comp <- NULL; dI.c_val$plot_visualize <- NULL;  dI.c_val$QC_table = NULL
  
    if(!("clean" %in% input$process)){
      result <- "please go back to the 'Loading & parsing samples' tab and select the 'clean' module"
      dI.c(result)}
    else
      if (!(inherits(x = dI(),"flowSet")))
      {
        result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
        dI.c(result)}
    else {
      dI.c(dI())
      
    print("start cleaning process")
      output$console_output_pre <- renderPrint("start cleaning process - waiting...")
      
      fFvectorName <- vector(mode = "character", length = length(dI.c()))
      for (i in seq_along(1:length(dI.c()))){
        fF <- dI.c()[[i]]
        no_time <- NULL
        if(!("time" %in% colnames(exprs(fF)))&&!("Time" %in% colnames(exprs(fF)))){ #!!! try to undestand why it does not show
          no_time <- paste0("the 'Time' information of sample nr.",i,"is not present inside the expression matrix. ",
                            "Probably other samples as well are without the 'time' channel: ",
                            "The flow rate quality check, the signal acquisition check and their related event cleaning actions", 
                            "are not possible on those samples")}
        fFvectorName[i] <- basename(fF@description$FILENAME)
        p1 = paste0("C_p1_",as.character(input$clean_param1))
        p2 = paste0("p2_",as.character(input$clean_param2))
        fFvectorName[i] <- paste0("./tmpdata/",p1,"_",p2,"_",fFvectorName[i])
        fF@description$FILENAME <- basename(fFvectorName[i])
       
        nome_fF.c <- cleaning_fS(flow_Frame = fF, fFname = fFvectorName[i], 
                                 par1 = input$clean_param1, par2 = input$clean_param2, seed = input$set_seed)
        print(paste0("sample nr.",as.character(i) ," cleaned"))
        
        check <- "OK"
        if (nome_fF.c != fFvectorName[i]){
          check <- "fail"
          if (grepl("This does not seem to be a valid FCS2.0, FCS3.0 or FCS3.1 file", nome_fF.c)) {
            result <- paste0("the sample nr.", as.character(i),"cannot be parsed as regular FCS files")
            result <- paste0(result, " - flow_auto_qc routine reports: ", nome_fF.c )} 
          if (grepl("time channel", nome_fF.c)) {
            result <- paste0("Is the time channel present in the expression matrix of the ",as.character(i)," flowFrame?")
            result <- paste0(result, " - flow_auto_qc routine reports: ", nome_fF.c )}
          else{
            result <- paste0("Something went wrong with the flowFrames nr.",as.character(i))
            result <- paste0(result, " - flow_auto_qc routine reports: ", nome_fF.c )}
          dI.c(result)
          break
          }}
      
      if (check != "fail"){
        
        nome = substr(x = fFvectorName[i], start = 11, stop = 15)
        nome = paste0("^", nome, "*")
        fFfiles <- list.files(path = "./tmpdata", pattern = nome, full.names = TRUE) 
        #fFfiles <- str_sort(x = fFfiles, numeric = TRUE) #this is necessary to sort in the correct way the flowFrames in the flowSet
        fS.c <- read.flowSet(file = fFfiles, truncate_max_range = FALSE)
      
        for (i in seq_along(1:length(fS.c)))
          {fS.c[[i]]@description$FILENAME <- basename(fS.c[[i]]@description$FILENAME)}
      
        dI.c(fS.c)
        dI.H(fS.c)}
      
      if (inherits(x = dI.c(),"flowSet")){
        
        output$clean_panelStatus <- reactive({input$clean_panelStatus=="show"})
        outputOptions(output, "clean_panelStatus", suspendWhenHidden = FALSE)
        dI.c_val$QC_table <- table_QC(dir_QC = "./resultsQC")
        if (Sys.info()["sysname"] == "Windows") play(x = bell)
      }
      
      print("cleaning process ends")
      output$console_output_pre <- renderPrint("cleaning process ends")}
    
    output$checktxt_cleaning_fS <- renderText({
      if (!(inherits(x = dI.c(),"flowSet"))){
        if (inherits(x = dI.c(),"character")) 
        {messaggio <- dI.c()
        validate(need(expr =  (inherits(x = dI.c(),"flowSet")), message = messaggio))}}
        if (inherits(x = no_time,"character")) 
        {validate(need(expr =  (!inherits(x = no_time, "character")), message = no_time))}
      })
    
  })
  ################################################### runCleanShow  -----    
  observeEvent(eventExpr = input$runCleanShow, {
  
    if(inherits(x = dI.c(),"flowSet")) {
      
      print("start data visualization of the cleaned flowSet")
      output$console_output_pre <- renderPrint("start data visualization of the cleaned flowSet - waiting...")
      validate(need(expr = !is.null(dI.c()), message = "press run to clean FCS files"))
      dI.c_val$fS_table <- table_fS(flow_Set  = dI.c())[[1]]
      lista <- dim_comparison(dI(), dI.c())
      dI.c_val$dim_table <- lista[[1]]
      dI.c_val$plot_dim_table <- lista[[2]]
      dI.c_val$fS_summary <- lista[[3]]
      dI.c_val$plot_comp <- plot_comparison(dI(), dI.c())
      if (inherits(x = dI.c(),"flowSet")){ 
        cdata <- session$clientData
        dI.c_val$plot_visualize <- visualize_expr_fS(flow_Set = dI.c(), colonne = 3,
                                            width = cdata$output_marker_density_clean_width, 
                                            height = cdata$output_marker_density_clean_height)}
      else
      {dI.c_val$plot_visualize <- visualize_expr_fF(flow_Frame = dI.c())}
      
      output$cleaned <- renderDT(DT::datatable(dI.c_val$fS_table, options = list(pageLength = 20)))
      output$cleanedSummary <- renderDT(DT::datatable(dI.c_val$fS_summary, options = list(pageLength = 20)) %>% 
                                          DT::formatStyle(columns = c(1,2), `font-size` = '16pt'))
      output$clean_comparison <- renderDT(DT::datatable(dI.c_val$dim_table, options = list(pageLength = 20)))
      output$clean_comparison_plot <- renderPlot({dI.c_val$plot_dim_table})
      output$plot_clean_comparison <- renderPlotly({dI.c_val$plot_comp})
      output$marker_density_clean <- renderPlotly({dI.c_val$plot_visualize}) #rendered without plotly for problems
      output$QC_clean <- renderDT(DT::datatable(dI.c_val$QC_table, options = list(pageLength = 20)))
      
      output$downloadCleanedfS_table_csv <- downloadHandler(
        filename = function() {paste("fS_table","csv", sep = ".")},
        content <- function(file) {file.copy("./tmpdata/fS_table.csv", file)},
        contentType = "text/csv")
      
      output$downloadComparison_table_csv <- downloadHandler(
        filename = function() {paste("fS_comparison_table","csv", sep = ".")},
        content <- function(file) {file.copy("./tmpdata/fS_comparison_table.csv", file)},
        contentType = "text/csv")
      
      output$downloadQC_table_csv <- downloadHandler(
        filename = function() {paste("_QCmini","txt", sep = ".")},
        content <- function(file) {
          
          file.copy("./resultsQC/_QCmini.txt", file)},
        contentType = "text/csv")
      
      output$c.entryEvents <- renderValueBox({
        valueBox(value = dI_val$fS_summary[2],subtitle =  "Total entry events", icon = icon("list"), color = "green")})
      output$c.handledEvents <- renderValueBox({
        valueBox(value = dI.c_val$fS_summary[3],subtitle =  "Total output events", icon = icon("list"), color = "red")})
      output$c.percEvents <- renderValueBox({
        measure <- 1 - (as.numeric(dI.c_val$fS_summary[2]) - as.numeric(dI.c_val$fS_summary[3]))/as.numeric(dI.c_val$fS_summary[2])
        valueBox(value = paste0(round(100 * measure, digits = 1),"%"),
                 subtitle =  "Percentage of remaining events", icon = icon("list"), color = "red")})
      
      print("Visualization of the cleaning process ends")
      output$console_output_pre <- renderPrint("visualization of the cleaning process ends")
     
      output$downloadCleanfS <- downloadHandler(
        filename = function() {
          paste("cfS", "zip", sep = ".")
        },
        content = function(file){
          del_tmpdata()
          file_path <- file
          tmpdir <- tempdir()
          setwd(tempdir())
          fSlength <- length(dI.c())
          nr_sample <- 1:fSlength
          for (i in seq_along(nr_sample)) {
            
            fFvectorName <- dI.c()[[i]]@description$FILENAME
            write.FCS(dI.c()[[i]], filename = fFvectorName)
          }
          
          Zip_Files <- list.files(path = getwd(), pattern = "^C_p")
          zip(zipfile = file_path, files = Zip_Files)
          setwd(app_dir)
        },
        contentType = "application/zip")
    } 
    else {
      result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
      dI.c(result)}
    
    output$checktxt_cleaning_fS <- renderText({
      if (!is.null(dI.c())){
        if (inherits(x = dI.c(),"character")) 
        {messaggio <- dI.c()
        validate(need(expr =  (inherits(x = dI.c(),"flowSet")), message = messaggio))}}})
  })
  
  
  ################################################### ____________Scale -----    
  
  dI.t <- reactiveVal(NULL)
  dI.t_val <- reactiveValues(prange = NULL, pdens = NULL)
  
  observeEvent(eventExpr = input$runScale, {
    
    del_tmpdata(type = "all")
    dI.t(NULL)
    dI.t_val$prange <- NULL; dI.t_val$pdens <- NULL
    
    if(!("scale" %in% input$process)){
      result <- "please go back to the 'Loading & parsing samples' tab and select the 'scale' module"
      dI.t(result)}
    else
      if (!(inherits(x = dI(),"flowSet")))
      {
        result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
        dI.t(result)
      }
    else{ 
      dI.t(dI())
      if (inherits(x = dI.c(),"flowSet")) {
        dI.t(dI.c())}
      
      dI.H(dI.t())
      
      print("start transforming process")
      output$console_output_pre <- renderPrint("start transforming process - waiting...")
      output$scale_panelStatus <- reactive({input$scale_panelStatus=="show"})
      outputOptions(output, "scale_panelStatus", suspendWhenHidden = FALSE)
      
      fFvectorName <- vector(mode = "character", length = length(dI.t()))
      fS.t <- dI.t()
      
      for (i in seq_along(1:length(fS.t))){
        fF <- fS.t[[i]]
        fFvectorName[i] <- basename(fF@description$FILENAME)
        p = paste0("T_",as.character(input$trans),"_")
        fFvectorName[i] <- paste0(p,fFvectorName[i])}
      
      fS.t <- file_fcs(flow_Set = fS.t, stringhe = fFvectorName)
      fS.t <- trans_fS(flow_Set = fS.t, type_of = input$trans, 
                       param = input$trans_param_arcsinh, 
                       min_quantile = input$trans_param_min, max_quantile = input$trans_param_max)
      
      if (!(inherits(x = fS.t,"flowSet"))){
        fS.t <- "Cytochain reported issues in calculating the scaled flowSet"}
      
      if (inherits(x = fS.t,"flowSet")){
        for (i in seq_along(1:length(fS.t)))
          {fS.t[[i]]@description$FILENAME <- basename(fS.t[[i]]@description$FILENAME)}}
      
      dI.t(fS.t)
      dI.H(fS.t)
      
      print("transforming process ends")
      output$console_output_pre <- renderPrint("transforming process ends")
      if (Sys.info()["sysname"] == "Windows") play(x = bell)}
    
      output$checktxt_transforming_fS <- renderText({
        if (!is.null(dI.t())){
          if (inherits(x = dI.t(),"character")) 
          {messaggio <- dI.t()
          validate(need(expr =  (inherits(x = dI.t(),"flowSet")), message = messaggio))}}})
      
  })
  
  ################################################### runScaleShow -----    
  observeEvent(eventExpr = input$runScaleShow, {
    
    if(inherits(x = dI.t(),"flowSet")) { 

      validate(need(expr = !is.null(dI.t()), message = "press run to transform FCS files"))
      print("start data visualization of the scaled flowSet")
      output$console_output_pre <- renderPrint("start data visualization of the scaled flowSet - waiting...")
      
      dI.t_val$prange <- visualize_rangefS(flow_Set = dI.t())
      # to relay the height/width of the plot's container, we'll query this session's client data 
      #http://shiny.rstudio.com/articles/client-data.html
      cdata <- session$clientData
      dI.t_val$pdens <- visualize_expr_fS(flow_Set = dI.t(), colonne = 2,
                                          width = cdata$output_marker_density_trans_width, 
                                          height = cdata$output_marker_density_trans_height)
      print("Visualization process ends")
      output$console_output_pre <- renderPrint("visualization of the scaled flowSet ends") 
      output$marker_range_trans <- renderPlotly({dI.t_val$prange})
      
      res.expr <- exprs_sub(flow.Set = dI.t())
      desc <- res.expr[[3]]
      numero <-length(desc)
      output$marker_density_trans <- renderPlotly({dI.t_val$pdens})}
    else {
        result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
        dI.t(result)}
    
  output$checktxt_transforming_fS <- renderText({
      if (!is.null(dI.t())){
        if (inherits(x = dI.t(),"character")) 
        {messaggio <- dI.t()
        validate(need(expr =  (inherits(x = dI.t(),"flowSet")), message = messaggio))}}})
  })
  
  #################################################################### saveScale -----
  
  observeEvent(eventExpr = input$saveScale , {
   
    if(inherits(x = dI.t(),"flowSet")) {
      del_tmpdata("all")
      print("save the whole transformed flowSet with all the original channels")
      output$console_output_pre <- renderPrint("save the transformed flowSet with all the original channels - waiting...")  
      fFvectorName <- vector(mode = "character", length(dI.t()))
      for (i in seq_along(1:length(dI.t()))){
        fFvectorName[[i]] <- dI.t()[[i]]@description$FILENAME}
      
      save_fS(in_flow_Set = dI(), handled_flow_Set = dI.t(), stringhe=fFvectorName)
      fFvectorName_temp <- paste0("./tmpdata/complete_",fFvectorName)
      fS.complete <- read.flowSet(files = fFvectorName_temp, truncate_max_range = FALSE)
      fFvectorName_comp <- basename(fFvectorName_temp)
      
      output$saveScale_panelStatus <- reactive({input$saveScale_panelStatus=="show"})
      outputOptions(output, "saveScale_panelStatus", suspendWhenHidden = FALSE)
      print("samples produced")
      output$console_output_pre <- renderPrint("samples produced")
      if (Sys.info()["sysname"] == "Windows") play(x = bell)
    }
    
    output$downloadTransfS <- downloadHandler(
      filename = function() {
        paste("tfS", "zip", sep = ".")
      },
      content = function(file){
        del_tmpdata()
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        fSlength <- length(dI.t())
        nr_sample <- 1:fSlength
        for (i in seq_along(nr_sample)){
          fFvectorName <- dI.t()[[i]]@description$FILENAME
          write.FCS(dI.t()[[i]], filename = fFvectorName)}
        
        Zip_Files <- list.files(path = getwd(), pattern = "^T_")
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")
    
    output$downloadTransfS_complete <- downloadHandler(
      filename = function() {
        paste("tfS_complete", "zip", sep = ".")
      },
      content = function(file){
        del_tmpdata()
        
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        fSlength <- length(dI.t())
        nr_sample <- 1:fSlength
        for (i in seq_along(nr_sample)){
          write.FCS(fS.complete[[i]], filename = fFvectorName_comp[[i]])}
        Zip_Files <- list.files(path = getwd(), pattern = "^complete_")
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")
  })
  
  ################################################### _____________Align ------    
  dI.a <- reactiveVal(NULL)
  dI.a_val <- reactiveValues(align_list = NULL, densb4 = NULL, densafter = NULL, dens = NULL)
  
  output$marker2align <- renderDT(DT::datatable(dI_val$fS_dim, options = list(pageLength = 30))) 
  #dI_val$fS_dim viene popolata all'interno della procedura di loading come prima operazione attraverso dim_fS. 
  #Valutare se spostare all'interno di dim_fS quella parta della funzione cell_Id_fS che scrive la descrizione del 
  #marker là dove non è presente. Indagare anche per l'input dei campioni
  #dal metadata workflow. Vedi anche https://github.com/RGLab/flowStats/issues/24 
  output$mark2align.res = renderPrint(input$marker2align_rows_selected)
  
  
  observeEvent(eventExpr = input$runAlign, {
   
    del_tmpdata(type = "all")
    dI.a(NULL)
    dI.a_val$align_list <- NULL; dI.a_val$densb4 <- NULL; dI.a_val$densafter <- NULL; dI.a_val$dens <- NULL;
    validate(need(expr = !is.null(input$runAlign), message = "please select some samples to be aligned"))
    transformed <- FALSE
    
    if(!("align" %in% input$process)){
      result <- "please go back to the 'Loading & parsing samples' tab and select the 'align' module"
      dI.a(result)}
    else
      if (!(inherits(x = dI(),"flowSet")))
      {
        result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
        dI.a(result)}
    if (is.null(input$marker2align_rows_selected)){
      result <- "please select some markers to align"
      dI.a(result)}
    else {
      if (length(dI()) > 1){
        
        dI.a(dI())
        if (inherits(x = dI.c(),"flowSet")) {
          dI.a(dI.c())} #dI()->dI.c->dI.t->dI.d
        if (inherits(x = dI.t(),"flowSet")) {
          transformed <- TRUE
          dI.a(dI.t()) #dI()->dI.c->dI.t->dI.d
          dI.c(NULL)}
      
        dI.H(dI.a())
      
        print("start aligning process")
        output$console_output_pre <- renderPrint("start aligning process - waiting...")
        output$align_panelStatus <- reactive({input$align_panelStatus=="show"})
        outputOptions(output, "align_panelStatus", suspendWhenHidden = FALSE)
        fS.a <- dI.a()
        fFvectorName <- vector(mode = "character", length = length(fS.a))
      
        for (i in seq_along(1:length(fS.a))){
          fF <- fS.a[[i]]
          fFvectorName[i] <- basename(fF@description$FILENAME)
          fFvectorName[i] <- paste0("A_",fFvectorName[i])}
      
        fS.a <- file_fcs(flow_Set = fS.a, stringhe = fFvectorName)
      
        dI.a_val$align_list <- align_flowStats(flow_Set = fS.a, marker_list = input$marker2align_rows_selected,
                                             trans = transformed)
        #align_flowStats strip-off the "cell_Id" from the samples
        fS.a <- dI.a_val$align_list[[1]]
      
        if (inherits(x = fS.a,"flowSet")){
          for (i in seq_along(1:length(fS.a)))
          {fS.a[[i]]@description$FILENAME <- basename(fS.a[[i]]@description$FILENAME)}
          dI.H(fS.a)}
      
        dI.a(fS.a)
      
        print("aligning process ends")
        output$console_output_pre <- renderPrint("aligning process ends")
        
        if (Sys.info()["sysname"] == "Windows") play(x = bell)}
      
      else{
        dI.a("It needs at least two flowFrames in your flowSet to run the alignment process")
        dI.a_val$align_list[[1]] <- "It needs at least two flowFrames in your flowSet to run the alignment process"}}
      
    output$checktxt_aligning_fS <- renderText({
      if (!is.null(dI.a())){
        if (inherits(x = dI.a(),"character")) 
        {messaggio <- dI.a()
        validate(need(expr =  (inherits(x = dI.a(),"flowSet")), message = messaggio))}}})
  })
  
  
  ################################################### runAlignShow ---- 
  observeEvent(eventExpr = input$runAlignShow, {
    
    if(inherits(x = dI.a(),"flowSet")) {
      validate(need(expr = !is.null(dI.a()), message = "press run to align FCS files"))
      print("start data visualization of the aligned flowSet")
      output$console_output_pre <- renderPrint("start data visualization of the aligned flowSet - waiting...")
      validate(need(expr = !is.null(dI.a()), message = "It needs valid FCS files"))
      dI.a_val$densb4 <- dI.a_val$align_list[[2]]
      dI.a_val$densafter <- dI.a_val$align_list[[3]]

      cdata <- session$clientData
      dI.a_val$dens <- visualize_expr_fS(flow_Set = dI.a(), colonne = 2,
                                          width = cdata$output_plot_align_after2_width, 
                                          height = cdata$output_plot_align_after2_height)
      
      if (inherits(x = dI.a_val$densb4,"trellis")){
      output$plot_align_b4 <- renderImage({
        outfile <- "./tmpdata/PlotAlignB4.png"
        list(src = outfile, contentType = 'image/png', alt = "marker density plot b4 the alignment")}, deleteFile = F)}
      
      if (inherits(x = dI.a_val$densafter,"trellis")){
      output$plot_align_after <- renderImage({
        outfile <- "./tmpdata/PlotAlignAfter.png"
        list(src = outfile, contentType = 'image/png', alt = "marker density plot b4 the alignment")}, deleteFile = F)}
      
      #if (inherits(dI.a_val$densb4,"ggplot")){
      #output$plot_align_b4 <- renderPlot({dI.a_val$densb4})}
      #if (inherits(dI.a_val$densafter,"ggplot")){
      #output$plot_align_after <- renderPlot({dI.a_val$densafter})}
      if (inherits(x = dI.a_val$dens,"plotly")){
        output$plot_align_after2 <- renderPlotly({dI.a_val$dens})}
      print("Visualization process ends")
      output$console_output_pre <- renderPrint("visualization of the aligned flowSet ends")
    }
    else{
      result <- dI.a_val$align_list[[1]]
      dI.a(result)}
    
    output$checktxt_aligning_fS <- renderText({
    if (!is.null(dI.a())){
      validate(need(expr =  (inherits(x = dI.a(),"flowSet")), message = dI.a()))}})
  })
  
  #################################################################### saveAlign -----
  
  observeEvent(eventExpr = input$saveAlign , {
    
    if(inherits(x = dI.a(),"flowSet")) {
      del_tmpdata("all")
      print("save the whole aligned flowSet with all the original channels")
      output$console_output_pre <- renderPrint("save the aligned flowSet with all the original channels - waiting...")  
      fFvectorName <- vector(mode = "character", length(dI.a()))
      for (i in seq_along(1:length(dI.a()))){
        fFvectorName[[i]] <- dI.a()[[i]]@description$FILENAME}
      save_fS(in_flow_Set = dI(), handled_flow_Set = dI.a(), stringhe = fFvectorName)
      
      fFvectorName_temp <- paste0("./tmpdata/complete_",fFvectorName)
      fS.complete <- read.flowSet(files = fFvectorName_temp, truncate_max_range = FALSE)
      fFvectorName_comp <- basename(fFvectorName_temp)
    
      output$saveAlign_panelStatus <- reactive({input$saveAlign_panelStatus=="show"})
      outputOptions(output, "saveAlign_panelStatus", suspendWhenHidden = FALSE)
      if (Sys.info()["sysname"] == "Windows") play(x = bell)
      print("samples produced")
      output$console_output_pre <- renderPrint("samples produced")}
    
    output$downloadAlignedfS <- downloadHandler(
      filename = function() {
        paste("afS", "zip", sep = ".")
      },
      content = function(file){
        
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        fSlength <- length(dI.a())
        nr_sample <- 1:fSlength
        for (i in seq_along(nr_sample)) {
          fFvectorName <- dI.a()[[i]]@description$FILENAME
          write.FCS(dI.a()[[i]], filename = fFvectorName)}
    
        Zip_Files <- list.files(path = getwd(), pattern = "^A_")  
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")
    
    output$downloadAlignedfS_complete <- downloadHandler(
      filename = function() {
        paste("afS_complete", "zip", sep = ".")
      },
      content = function(file){
        del_tmpdata()
        
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        fSlength <- length(dI.a())
        nr_sample <- 1:fSlength
        for (i in seq_along(nr_sample)){
          write.FCS(fS.complete[[i]], filename = fFvectorName_comp[[i]])}
        Zip_Files <- list.files(path = getwd(), pattern = "^complete_")
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")
  })
  
  ################################################### ______Downsample -----    
  
  dI.d <- reactiveVal(NULL)
  dI.d.b4 <- reactiveVal(NULL)
  dI.d_val <- reactiveValues(fS_table = NULL, fS_comp = NULL, fS_comp_plot = NULL, fS_summary = NULL, pd = NULL, dd = NULL, 
                             b4_score_plot = NULL, b4_dens_plot = NULL, after_score_plot = NULL, after_dens_plot = NULL,
                             comb_flag = NULL, no_down_flag = NULL)
  
  observeEvent(eventExpr = input$runDown, {
    
    validate(need(expr = !is.null(input$runDown), message = "please select some samples to be down-sampled"))
    dI.d(NULL); dI.d.b4(NULL);
    
    dI.d_val$fS_table <- NULL; dI.d_val$fS_comp <- NULL; dI.d_val$fS_comp_plot <- NULL; dI.d_val$pd <- NULL; 
    dI.d_val$fS_summary[2] <- "refresh...";
    dI.d_val$fS_summary[3] <- "refresh...";
    dI.d_val$dd <- NULL; dI.d_val$b4_score_plot <- NULL; dI.d_val$b4_dens_plot <- NULL; dI.d_val$after_score_plot <- NULL; 
    dI.d_val$after_dens_plot <- NULL; dI.d_val$comb_flag <- NULL; dI.d_val$no_down_flag <- NULL
    
    if(!("downsample" %in% input$process)){
      result <- "please go back to the 'Loading & parsing samples' tab and select the 'downsample' module"
      dI.d(result)}
        
    if (!(inherits(x = dI(),"flowSet"))){
          result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
          dI.d(result)}
        else{
          dI.d(dI())
          dI.d.b4(dI())
          if (inherits(x = dI.c(),"flowSet")) {
            dI.d(dI.c()) #dI()->dI.c->dI.t->dI.a 
            dI.d.b4(dI.c())}
          if (inherits(x = dI.t(),"flowSet")) {
            dI.d(dI.t()) #dI()->dI.c->dI.t->dI.a 
            dI.d.b4(dI.t())
            dI.c(NULL)} #I can't reset the the involved flowSet -> dI.d.b4() (I just reset the last in the chain not involved
          if (inherits(x = dI.a(),"flowSet")) {
            dI.d(dI.a()) #dI()->dI.c->dI.t->dI.a 
            dI.d.b4(dI.a())
            dI.t(NULL)}
          
            dI.H(dI.d())
            dI.d(NULL)
            print("start downsampling process")
            output$console_output_pre <- renderPrint("start downsampling process - waiting...")
            del_tmpdata(type = "all")
            fFvectorName <- vector(mode = "character", length = length(dI.H()))
            fS.d <- dI.H()
            for (i in seq_along(1:length(fS.d))){
              fFvectorName[i] <- fS.d[[i]]@description$FILENAME}
              
            dI.d_val$comb_flag <- FALSE
            dI.d_val$no_down_flag <- FALSE
            for (i in seq_along(1:length(fS.d))){
              fFvectorName[i] <- basename(fS.d[[i]]@description$FILENAME)}
            fS <- file_fcs(flow_Set = dI.H(), stringhe = fFvectorName)
            
            for (i in seq_along(1:length(dI.H()))){
              if (dim(exprs(dI.H()[[i]]))[1] < MIN_SAMPLE_LENGTH){ #MIN_SAMPLE_LENGTH
                dI.d_val$no_down_flag <- TRUE
                fS.d <- paste0("Impossible to perform downsampling process if the number of events of sample nr.",as.character(i),
                               " is too small (", as.character(dim(exprs(dI.H()[[i]]))[1]), " < MIN_SAMPLE_LENGTH = ", 
                               as.character(MIN_SAMPLE_LENGTH),")")
                break}}
              
            if ((input$down_algo == "spade")&&(dI.d_val$no_down_flag==FALSE)){
              dI.d_val$comb_flag <- TRUE
              fS.d <- spade_fS(flow_Set = fS, stringhe = fFvectorName, type = input$down_type, p1 = input$down_percentile, 
                               p2 = 100, seed = input$set_seed)
#in case type = input$down_type = down_val, the flowSet returned, have the density dim: no downsampling is executed 
              if (inherits(x = fS.d,"flowSet")){
                
                fFfiles <- list.files(path = "./tmpdata", pattern="^downdens_", full.names = TRUE) 
                #fFfiles <- str_sort(x = fFfiles, numeric = TRUE) 
                #this is necessary to sort in the correct way the flowFrames in the flowSet
                fS <- read.flowSet(files = fFfiles, truncate_max_range = FALSE)
                score <- plot_score(flow_Set = fS, type = input$down_algo, sample_color = color.sample)
                if (inherits(x = score,"list")){
                  dI.d_val$b4_score_plot <- score[[2]]
                  dI.d_val$b4_dens_plot <- score[[3]]}
                else
                  {fS.d <- score}
                  #trying to spade downsampling per value
                if (input$down_type == "down_val"){
                  prefix <- downsampling_fS(flow_Set = fS, stringhe = fFvectorName, algo = input$down_algo, 
                                            type = input$down_type, perc = input$down_percentile, 
                                            value = input$down_value, seed = input$set_seed)
                  if ((inherits(x = prefix,"character"))&&unique((substr(prefix, start = 1, stop = 6)!="System"))) 
                    {fS.d <- read.flowSet(files = prefix, truncate_max_range = FALSE)} else {fS.d <- prefix}}
                
                if (input$down_type == "equalize"){
                  prefix <- downsampling_fS(flow_Set = fS, stringhe = fFvectorName, algo = input$down_algo, 
                                            type = input$down_type, perc = input$down_percentile, 
                                           value = input$down_value, seed = input$set_seed, cut_off = input$down_number)
                if ((inherits(x = prefix, "character"))&&unique((substr(prefix, start = 1, stop = 6)!="System"))) 
                  {fS.d <- read.flowSet(files = prefix, truncate_max_range = FALSE)} else {fS.d <- prefix}}
                
                if (inherits(x = fS.d,"flowSet")){
                  score <- plot_score(flow_Set = fS.d, type = input$down_algo, sample_color = color.sample)
                  if (inherits(x = score,"list")){
                    dI.d_val$after_score_plot <- score[[2]]
                    dI.d_val$after_dens_plot <- score[[3]]}}}} #end of spade part

            if ((input$down_algo == "Random")&&(dI.d_val$no_down_flag==FALSE)){
              
              #in case of random selection, there is no need to produce a score dimension: -> downsampling_fS directly
              dI.d_val$comb_flag <- TRUE
              dI.d_val$b4_score_plot <- NULL
              dI.d_val$b4_dense_plot <- NULL
              dI.d_val$after_score_plot <- NULL
              dI.d_val$after_dense_plot <- NULL
              if (input$down_type=="down_perc"){
                prefix <- downsampling_fS(flow_Set = fS, stringhe = fFvectorName, algo = input$down_algo, 
                                          type = input$down_type, perc = input$down_percentile, 
                                        value = input$down_value, seed = input$set_seed)
                if ((inherits(x = prefix,"character"))&&unique((substr(prefix, start = 1, stop = 6)!="System"))) 
                {fS.d <- read.flowSet(files = prefix, truncate_max_range = FALSE)} else {fS.d <- prefix}}
              if (input$down_type=="down_val"){
                fS.d <- "You cannot perform a random downsampling with the type 'by value'"}
              if (input$down_type=="equalize"){
                prefix <- downsampling_fS(flow_Set = fS, stringhe = fFvectorName, algo = input$down_algo, 
                                          type = input$down_type, perc = input$down_percentile, 
                                          value = input$down_value, seed = input$set_seed, cut_off = input$down_number)
                if ((inherits(x = prefix,"character"))&&(unique(substr(prefix, start = 1, stop = 6)!="System")))
                  {fS.d <- read.flowSet(files = prefix, truncate_max_range = FALSE)} else {fS.d <- prefix}}} 

            if (input$down_algo %in% c("RKOF", "LOOP", "LOF", "LDF", "KNN_SUM", "KNN_AGG", "LDE")&&(dI.d_val$no_down_flag==FALSE)){
              dI.d_val$comb_flag <- TRUE
              fS <- outlier_score(flow_Set = dI.H(), stringhe = fFvectorName, algo = input$down_algo, kappa = input$down_kappa, 
                                    kappa_max = input$down_kappa_max, par = input$down_par, seed = input$set_seed)
#the outlier_score is the routine producing the score with the selected algorithm (no spade)                
              if (inherits(x = fS,"flowSet")){
                score <- plot_score(flow_Set = fS, type = input$down_algo, sample_color = color.sample)
                if (inherits(x = score,"list")){
                  dI.d_val$b4_score_plot <- score[[2]]
                  dI.d_val$b4_dens_plot <- score[[3]]}
                else
                  {fS.d <- score}
                if (input$down_type=="equalize"){
                prefix <- downsampling_fS(flow_Set = fS, stringhe = fFvectorName, algo = input$down_algo, 
                                          type = input$down_type, perc = input$down_percentile, 
                                         value = input$down_value, seed = input$set_seed, cut_off = input$down_number)}
                else {prefix <- downsampling_fS(flow_Set = fS, stringhe = fFvectorName, algo = input$down_algo, 
                                                type = input$down_type, perc = input$down_percentile, 
                                               value = input$down_value, seed = input$set_seed)}
                if ((inherits(x = prefix, "character"))&&unique((substr(prefix, start = 1, stop = 6)!="System")))
                  {fS.d <- read.flowSet(files = prefix, truncate_max_range = FALSE)} else {fS.d <- prefix}}
                if (inherits(x = fS.d,"flowSet")){
                  score <- plot_score(flow_Set = fS.d, type = input$down_algo, sample_color = color.sample)
                  if (inherits(x = score,"list")){
                    dI.d_val$after_score_plot <- score[[2]]
                    dI.d_val$after_dens_plot <- score[[3]]}}}
      
      if (inherits(x = fS.d,"flowSet")){
        for (i in seq_along(1:length(fS.d)))
        {fS.d[[i]]@description$FILENAME <- basename(fS.d[[i]]@description$FILENAME)}}
            
      print("downsampling process ends")
      output$console_output_pre <- renderPrint("downsampling process ends")
      if (Sys.info()["sysname"] == "Windows") play(x = bell)
      if (!(inherits(x = fS.d,"flowSet"))){
        dI.d_val$comb_flag == FALSE
        dI.H(fS.d)}
      else{
      output$down_panelStatus <- reactive({input$down_panelStatus=="show"})
      outputOptions(output, "down_panelStatus", suspendWhenHidden = FALSE)
      #at the end of the procedure dI.d()=fS.d; dI.d.b4() is the entry flowSet
      dI.H(fS.d)}}
    
    output$checktxt_downsampling_fS <- renderText({
      if (!is.null(dI.H())){
        if (inherits(x = dI.H(),"character")) 
          {messaggio <- dI.H()
          validate(need(expr =  (inherits(x = dI.H(),"flowSet")), message = messaggio))}}})
    
  })
  
  ################################################### runDownShow ---- 
  observeEvent(eventExpr = input$runDownShow, {
    
    if((inherits(x = dI.H(),"flowSet"))&&(dI.d_val$comb_flag == TRUE)) {
  
      validate(need(expr = !is.null(dI.H()), message = "press run to downsample the FCS files"))
      print("start data visualization process of the downsampled flowSet")
      output$console_output_pre <- renderPrint("start data visualization process of the downsampled flowSet - waiting...")
      
      dI.d_val$fS_table <- table_fS(flow_Set  = dI.H())[[1]]
      
      lista <- dim_comparison(dI(), dI.H())
      dI.d_val$fS_summary <- lista[[3]]
      lista <- dim_comparison(dI.d.b4(), dI.H())
      dI.d_val$fS_comp_plot <- lista[[2]]
      dI.d_val$fS_comp <- lista[[1]]
      dI.d_val$pd <- plot_comparison(dI.d.b4(), dI.H())
      
      cdata <- session$clientData
      dI.d_val$dd <- visualize_expr_fS(flow_Set = dI.H(), colonne = 3,
                                          width = cdata$output_marker_down_width, 
                                          height = cdata$output_marker_down_height)
      
      output$downsampled <- renderDT(DT::datatable(dI.d_val$fS_table, options = list(pageLength = 20)))
      output$downsampledSummary <- renderDT(DT::datatable(dI.d_val$fS_summary, options = list(pageLength = 20)) %>% 
                                              DT::formatStyle(columns = c(1,2), `font-size` = '16pt'))
      output$down_comparison_plot <- renderPlot({dI.d_val$fS_comp_plot})
      output$down_comparison <- renderDT(DT::datatable(dI.d_val$fS_comp, options = list(pageLength = 20)))
      
      output$plot_down_comparison <- renderPlotly({dI.d_val$pd})
      output$marker_density_down <- renderPlotly({dI.d_val$dd})
      
      if (input$down_algo!="Random"){
        output$down_algoStatus <- reactive({input$down_algoStatus=="show"})
        outputOptions(output, "down_algoStatus", suspendWhenHidden = FALSE)
        output$downsampled_b4Score <- renderPlotly(expr = dI.d_val$b4_score_plot)
        output$downsampled_afterScore <- renderPlotly(expr = dI.d_val$after_score_plot)
        output$downsampled_b4Dens <- renderPlotly(expr = dI.d_val$b4_dens_plot)
        output$downsampled_afterDens <- renderPlotly(expr = dI.d_val$after_dens_plot)}
      
      output$downloadDownsampled <- downloadHandler(
        filename = function() {paste("fS_table","csv", sep = ".")},
        content <- function(file) {file.copy("./tmpdata/fS_table.csv", file)},
        contentType = "text/csv")
      
      output$downloadComparisonTable <- downloadHandler(
        filename = function() {paste("fS_comparison_table","csv", sep = ".")},
        content <- function(file) {file.copy("./tmpdata/fS_comparison_table.csv", file)},
        contentType = "text/csv")
      
      output$d.entryEvents <- renderValueBox({
        valueBox(value = dI_val$fS_summary[2],subtitle =  "Total entry events", icon = icon("list"), color = "green")})
      output$d.handledEvents <- renderValueBox({
        valueBox(value = dI.d_val$fS_summary[3],subtitle =  "Total output events", icon = icon("list"), color = "red")})
      output$d.percEvents <- renderValueBox({
        measure <- 1 - (as.numeric(dI.d_val$fS_summary[2]) - as.numeric(dI.d_val$fS_summary[3]))/as.numeric(dI.d_val$fS_summary[2])
        valueBox(value = paste0(round(100*measure, digits = 1),"%"),
                 subtitle =  "Percentage of remaining events", icon = icon("list"), color = "red")})
      print("Visualization process ends")
      output$console_output_pre <- renderPrint("visualization process ends")}
    
    if(!(inherits(x = dI.H(),"flowSet"))){
      if (dI.d_val$comb_flag == FALSE)
        {dI.H("impossible to perform downsampling process with such a combination of parameters entries")}
      if ((dI.d_val$comb_flag == TRUE)&&(dI.d_val$no_down_flag == TRUE))
      {dI.H(paste0("Impossible to perform downsampling process if the number of events of one of the sample is too small 
                   ( < MIN_SAMPLE_LENGTH = ",
                   as.character(MIN_SAMPLE_LENGTH)))}}
    
    output$checktxt_downsampling_fS <- renderText({
      if (!is.null(dI.H())){
        if (inherits(x = dI.H(),"character")) 
        {messaggio <- dI.H()
        validate(need(expr =  (inherits(x = dI.H(),"flowSet")), message = messaggio))}}})
    
  })
  
  
  #################################################################### saveDown -----
  
  observeEvent(eventExpr = input$saveDown , {
    
    if(inherits(x = dI.H(),"flowSet")) {
      del_tmpdata("all")
      print("save the whole downsampled flowSet with all the original channels")
      output$console_output_pre <- renderPrint("save the dowsampled flowSet with all the original channels - waiting...")  
      fFvectorName <- vector(mode = "character", length(dI.H()))
      for (i in seq_along(1:length(dI.H()))){
        fFvectorName[[i]] <- dI.H()[[i]]@description$FILENAME}
      
      save_fS(in_flow_Set = dI(), handled_flow_Set = dI.H(), stringhe = fFvectorName)
      fFvectorName_temp <- paste0("./tmpdata/complete_",fFvectorName)
      fS.complete <- read.flowSet(files = fFvectorName_temp, truncate_max_range = FALSE)
      fFvectorName_comp <- basename(fFvectorName_temp)
      
      output$saveDown_panelStatus <- reactive({input$saveDown_panelStatus=="show"})
      outputOptions(output, "saveDown_panelStatus", suspendWhenHidden = FALSE)
      if (Sys.info()["sysname"] == "Windows") play(x = bell)
      print("samples produced")
      output$console_output_pre <- renderPrint("samples produced")}
    
    output$downloadDownfS <- downloadHandler(
      filename = function() {
        paste("dfS", "zip", sep = ".")
      },
      content = function(file){
        del_tmpdata()
        
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        fSlength <- length(dI.H())
        nr_sample <- 1:fSlength
        for (i in seq_along(nr_sample)){
          fFvectorName <- dI.H()[[i]]@description$FILENAME
          write.FCS(dI.H()[[i]], filename = fFvectorName)}
        
        Zip_Files <- list.files(path = getwd(), pattern = "^D_*") 
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")
    
    output$downloadDownfS_complete <- downloadHandler(
      filename = function() {
        paste("dfS_complete", "zip", sep = ".")
      },
      content = function(file){
        del_tmpdata()
        
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        fSlength <- length(dI.H())
        nr_sample <- 1:fSlength
        for (i in seq_along(nr_sample)){
          write.FCS(fS.complete[[i]], filename = fFvectorName_comp[[i]])}
        Zip_Files <- list.files(path = getwd(), pattern = "^complete_")
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")
    })
  
  ################################################### ______Concatenate -----    
  
  dI.co <- reactiveVal(NULL)
  dI.co.b4 <- reactiveVal(NULL)
  dI.co_val <- reactiveValues(fF_table = NULL, plot_visualize = NULL)
  
  observeEvent(eventExpr = input$runConc, {
    
    dI.co(NULL)
    dI.co.b4(NULL)
    dI.co_val$fF_table <- NULL; dI.co_val$plot_visualize <- NULL;
    dI.d.b4(NULL)
    
    validate(need(expr = !is.null(input$runConc), message = "please select some samples to be concatenated"))
    if(!("concatenate" %in% input$process)){
      result <- "please go back to the 'Loading & parsing samples' tab and select the 'concatenate' module"
      dI.co(result)}
    else{
      if (!(inherits(x = dI(),"flowSet")))
      {
        result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
        dI.co(result)}
      else{
        dI.co(dI())
        if (inherits(x = dI.H(),"flowSet")) {dI.co(dI.H())}

        print("start concatenating process")
        output$console_output_pre <- renderPrint("start concatenating process - waiting...")
        output$conc_panelStatus <- reactive({input$conc_panelStatus=="show"})
        outputOptions(output, "conc_panelStatus", suspendWhenHidden = FALSE)
        fS <- concatenating_fS(flow_Set = dI(), stringa = "conc_sample")
        fS.co <- concatenating_fS(flow_Set = dI.co(), stringa = "conc_sample")
        dI.co.b4(fS)
        dI.co(fS.co)
    
        print("concatenating process ends")
        output$console_output_pre <- renderPrint("concatenating process ends")
        if (Sys.info()["sysname"] == "Windows") play(x = bell)
      }}
    
    output$checktxt_concatenating_fS <- renderText({
      if (!is.null(dI.co())){
        if (inherits(x = dI.co(),"character")) 
        {messaggio <- dI.co()
        validate(need(expr =  (inherits(x = dI.co(),"flowFrame")), message = messaggio))}}})
  })
  
  ################################################### runConcShow ---- 
  observeEvent(eventExpr = input$runConcShow , {
    
    if(inherits(x = dI.co(),"flowFrame")) {
      print("start data visualization process of the concatenated flowSet")
      output$console_output_pre <- renderPrint("start data visualization process of the concatenated flowSet - waiting...")
      
      del_tmpdata()
      validate(need(expr = !is.null(dI.co()), message = "press run to concatenate FCS files"))
      dI.co_val$fF_table <- table_fF(flow_Frame = dI.co())
      dI.co_val$plot_visualize <- visualize_expr_fF(flow_Frame  = dI.co())
      
      output$concatenated <- renderTable({dI.co_val$fF_table})
      output$marker_density_conc <- renderPlotly({dI.co_val$plot_visualize})
      
      print("Visualization process of the concatenated ends")
      output$console_output_pre <- renderPrint("data visualization of the concatenated flowSet ends")
      
      output$downloadConcfS <- downloadHandler(
        filename = function() {
          paste("cofF", "zip", sep = ".")
        },
        content = function(file){
          del_tmpdata()
          file_path <- file
          tmpdir <- tempdir()
          setwd(tempdir())
          write.FCS(dI.co(), filename = "conc_fF.fcs")
          Zip_Files <- list.files(path = getwd(), pattern = "^conc_fF*")
          zip(zipfile = file_path, files = Zip_Files)
          setwd(app_dir)
        },
        contentType = "application/zip")
    }else{
      result <- "please go back to the 'Loading & parsing samples' tab and select some valid FCS files"
      dI.co(result)}
    
    output$checktxt_concatenating_fS <- renderText({
      if (!is.null(dI.co())){
        if (inherits(x = dI.co(),"character")) 
        {messaggio <- dI.co()
        validate(need(expr =  (inherits(x = dI.co(),"flowSet")), message = messaggio))}}})
    
  })
  

  ################################################### _____Edit metadata-----    
  ################################################### runMeta -----    
  meta.csv <- reactiveVal(NULL)
  meta.table <- reactiveVal(NULL)
  fS.in <- reactiveVal(NULL)
  fSmeta <- reactiveVal(NULL)
  metaResume <- reactiveVal(NULL)
  
  observeEvent(eventExpr = input$runMeta, handlerExpr =  {
    
    fSmeta(NULL); dI.c(NULL); dI.t(NULL); dI.a(NULL); dI.d(NULL); dI.d.b4(NULL); dI.co(NULL); dI.co.b4(NULL) 
    #the only active flowSet are dI() and dI.H()
    clean_graph()
    del_tmpdata(type = "all")
  
    if(!is.null(input$fSmeta_in)||!is.null(input$zipmeta_in)) {
      if (!is.null(input$fSmeta_in)){
        alphafSinput <-  input$fSmeta_in[order(input$fSmeta_in$name),] 
        lista <- loading_fS(inputfile = alphafSinput$datapath, tipo = "flowSet")
        fSmeta(lista[[1]])
        if (inherits(x = lista[[1]],"flowSet")){
          fS <- file_fcs(lista[[1]], stringhe = alphafSinput$name)
          fSmeta(fS)}}
      if (!is.null(input$zipmeta_in)){
        lista <- loading_fS(inputfile = input$zipmeta_in$datapath, tipo = "zip")
        fSmeta(lista[[1]])
        array_nomi <- lista[[2]]
        if (inherits(x = lista[[1]],"flowSet")){
          fS <- file_fcs(lista[[1]], stringhe = array_nomi)
          fSmeta(fS)}}}
      else{
        if (inherits(x = dI(),"flowSet")) {fSmeta(dI())}
        if (inherits(x = dI.H(),"flowSet")) {fSmeta(dI.H())}}

    validate(need(expr = !is.null(fSmeta()), message = "It needs valid FCS files"))
   
    if (inherits(x = fSmeta(),"flowSet")){ 
      mt_ok <- FALSE
      meta.csv(input$meta_file)
      toggle <- ifelse(test = !is.null(meta.csv()), yes = T, no = F)
      if (toggle){
        options(warn = 2) # Turn warnings into errors so they can be trapped
        result <- try(expr = read.csv(file=input$meta_file$datapath, header=TRUE, sep=",", check.names = TRUE, 
                                      #with check names all spaces in names are mutated in "."
                                    colClasses = c("character", "character", "character", "character", "character", "character"))) 
                                                #??? to set the last elemetn to "Date"
        if (inherits(x = result,"data.frame")){
          options(warn = 0)
          mt_in <- result
          good <- TRUE}
        else {  # Process any error messages
          # Ignore warnings while processing errors
          options(warn = -1)
          msg <- geterrmessage()
          mt_in <- paste("The table you loaded seems not to be in the correct format. System reports: ", msg, sep = " ")
        # Restore default warning reporting
          good <- FALSE
          options(warn=0)}
        
        if (good == FALSE){ # let's try with csv separated with ";"
          options(warn = 2) # Turn warnings into errors so they can be trapped
          result <- try(expr = read.csv(file=input$meta_file$datapath, header=TRUE, sep=";", check.names = TRUE, 
                                        #with check names all spaces in names are mutated in "."
                                        colClasses = c("character", "character", "character", "character", "character", "character"))) 
                                                #??? to set the last elemtn to "Date"
          if (inherits(x = result,"data.frame")){
            options(warn = 0)
            mt_in <- result
            good <- TRUE} 
          else {
            options(warn = -1)
            msg <- geterrmessage()
            mt_in <- paste("The table you loaded seems not to be in the correct format. System reports: ", msg, sep = " ")
            # Restore default warning reporting
            good <- FALSE
            options(warn=0)}}
        
        if ((inherits(x = mt_in, "data.frame"))&&(ncol(mt_in)==7)&&(nrow(mt_in)==length(fSmeta()))){
          names(mt_in) <- trimws(names(mt_in))
          names(mt_in) <- gsub("\\.", "_", names(mt_in)) #to change "." to "_"
          names(mt_in) <- make.names(names = names(mt_in), unique = TRUE) #to make syntactically valid names out of character vectors
          mt_in[,3] <-  trimws(mt_in[,3]) #to remove spaces 
          mt_in[,4] <-  trimws(mt_in[,4]) #to remove spaces 
          mt_in[,5] <-  trimws(mt_in[,5]) #to remove spaces 
          mt_in[,6] <-  trimws(mt_in[,6]) #to remove spaces 
          mt_in[,3] <-  gsub("\\.", "_", mt_in[,3]) #to change "." to "_"
          mt_in[,4] <-  gsub("\\.", "_", mt_in[,4]) #to change "." to "_"
          mt_in[,5] <-  gsub("\\.", "_", mt_in[,5]) #to change "." to "_"
          mt_in[,6] <-  gsub("\\.", "_", mt_in[,6]) #to change "." to "_"
          mt_in[,3] <- make.names(names = mt_in[,3], unique = FALSE)
          mt_in[,4] <- make.names(names = mt_in[,4], unique = FALSE)
          mt_in[,5] <- make.names(names = mt_in[,5], unique = FALSE)
          mt_in[,6] <- make.names(names = mt_in[,6], unique = FALSE)
          mt_in[,7] <- as.Date(x = mt_in[,7],format="%m/%d/%Y")
          set.seed(input$set_seed)
          color.sample <<- set_colors(nrow(mt_in))
          
          mt <<- mt_in
          if ((length(unique(mt_in[,6]))==(length(unique(mt_in[,7]))))&&(input$time_step==TRUE))
            { output$meta_panelStatus <- reactive({input$meta_panelStatus=="show"})
              outputOptions(output, "meta_panelStatus", suspendWhenHidden = FALSE)
              mt_ok <- TRUE}
          if (input$time_step==FALSE)
          { output$meta_panelStatus <- reactive({input$meta_panelStatus=="show"})
          outputOptions(output, "meta_panelStatus", suspendWhenHidden = FALSE)
          mt_ok <- TRUE}
          }} 
      else{ #toggle
        if (!is.null(input$fSmeta_in)){
          alphafSinput <- input$fSmeta_in[order(input$fSmeta_in$name),]
          nomi <- forbidden_char(alphafSinput$name) #was nomi <- forbidden_char(input$fSmeta_in$name)
          mt <<- edit_table(flow_Set = fSmeta(), file_name = nomi)
          color.sample <<- set_colors(nrow(mt))
          output$meta_panelStatus <- reactive({input$meta_panelStatus=="show"})
          outputOptions(output, "meta_panelStatus", suspendWhenHidden = FALSE)
          mt_ok <- TRUE}
      if ((!is.null(input$zipmeta_in))&&(length(array_nomi)>0)){
        nomi <- forbidden_char(array_nomi)
        mt <<- edit_table(flow_Set = fSmeta(), file_name = nomi)
        color.sample <<- set_colors(nrow(mt))
        output$meta_panelStatus <- reactive({input$meta_panelStatus=="show"})
        outputOptions(output, "meta_panelStatus", suspendWhenHidden = FALSE)
        mt_ok <- TRUE}
      if(is.null(input$fSmeta_in)&&is.null(input$zipmeta_in)){
        mt <<- edit_table(flow_Set = fSmeta())
        color.sample <<- set_colors(nrow(mt))
        output$meta_panelStatus <- reactive({input$meta_panelStatus=="show"})
        outputOptions(output, "meta_panelStatus", suspendWhenHidden = FALSE)
        mt_ok <- TRUE}}
      if(mt_ok)
        {meta.table(mt)}
      else
      {if (!(inherits(x = mt_in,"data.frame"))){
        output$checktxt_loading_meta_fS <- renderText({
        {messaggio <- paste0("something is wrong with the loaded meta-table. System reports: ", mt_in, 
                             ". You may start downloading meta-table automatically ",         
                             "generated when you do not provide any table in input. 
                             Try to manipulate that one according to your use and try with this new version")
        validate(need(expr =  (mt_ok), message = messaggio))}})}
        else{
          output$checktxt_loading_meta_fS <- renderText({
          {messaggio <- paste0("something is wrong with the loaded meta-table.", "You may start downloading meta-table automatically ",         
                               "generated when you do not provide any table in input. 
                               Try to manipulate that one according to your use and try with this new version")
        validate(need(expr =  (mt_ok), message = messaggio))}})}}}
    else
      {if (!(inherits(x = fSmeta(),"flowSet"))){
        output$checktxt_loading_meta_fS <- renderText({
        {messaggio <- fSmeta()
          validate(need(expr =  (inherits(x = fSmeta(),"flowSet")), message = messaggio))}})
      }}
  })
  
  ################################################### runMetaLoad -----    
  observeEvent(eventExpr = input$runMetaload, handlerExpr = {
    validate(need(expr = !is.null(fSmeta()), message = "It needs valid FCS files"))
    validate(need(expr = (inherits(x = mt,"data.frame")), message = "It needs valid dataframe"))
    validate(need(expr = !is.null(meta.table()), message = "It needs a valid meta.table"))
    hotmeta <- input$hot_meta
    if (!is.null(hotmeta)) 
    {#
      #mt <<- hot_to_r(hotmeta) #commented because hot_to_r transforms the dates in 'NA'
    date_s <- sort(mt$date)
    unique.date <- length(unique(date_s))
    uni.tag1 <<- length(unique(mt[,3]))
    uni.tag2 <<- length(unique(mt[,4]))
    uni.tag3 <<- length(unique(mt[,5]))
    uni.tag4 <<- length(unique(mt[,6]))
    
    if (input$time_step == TRUE){
      validate(need(expr = (length(date_s)==length(fSmeta())), message = "It needs some congruent meta-data"))
      validate(need(expr = ((date_s[1])==mt$date[1] &&  #check the date value
                              date_s[length(date_s)]==mt$date[length(mt$date)] &&  #check the date length
                              uni.tag4 == unique.date), message = "It needs some congruent meta-data"))}
    
    #totale <- uni.tag1 + uni.tag2 + uni.tag3 + input$MetaClustNum
    #since the color.metaclust is much more important to set than the rest of the colors, let's put it as a first color set in order to assign
    #always the same colors when you need to perform more than one analysis with different tags having different uni.tag
    totale <- uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4 + 50
    #the input$MetaClustNum parameter is in another tab (the Hign dim analysis tag). 
    #Let's assume to always assigne the maximum value (50) for input$MetaClustNum
    
    colori <- color_set[1:totale]
    #in order to have the same set of colors for the meta_clusters each time you perform the meta-clustering
    color.metaclust <<- colori[1:50] 
    # it was color.metaclust <<- colori[(uni.tag1 + uni.tag2 + uni.tag3 + uni.marker + uni.pheno + 1):length(colori)]
    color.tag1 <<- colori[(50 + 1):(50 + uni.tag1)]
    color.tag2 <<- colori[(50 + uni.tag1 + 1):(50 + uni.tag1 + uni.tag2)]
    color.tag3 <<- colori[(50 + uni.tag1 + uni.tag2 + 1):(50 + uni.tag1 + uni.tag2 + uni.tag3)]
    color.tag4 <<- colori[(50 + uni.tag1 + uni.tag2 + uni.tag3 + 1):(50 + uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4)]}
    
    v.tag1 <- rep.int(x = 1, times = length(fSmeta()))
    v.tag2 <- rep.int(x = 1, times = length(fSmeta()))
    v.tag3 <- rep.int(x = 1, times = length(fSmeta()))
    v.tag4 <- rep.int(x = 1, times = length(fSmeta()))
    names(v.tag1) <- mt[,3]
    names(v.tag2) <- mt[,4]
    names(v.tag3) <- mt[,5]
    names(v.tag4) <- mt[,6]
    col.tag1 <- vector(mode = "character",length = length(fSmeta()))
    col.tag2 <- vector(mode = "character",length = length(fSmeta()))
    col.tag3 <- vector(mode = "character",length = length(fSmeta()))
    col.tag4 <- vector(mode = "character",length = length(fSmeta()))
    
    for(i in 1:length(fSmeta())){
      char.1 <- mt[,3][i]
      char.2 <- mt[,4][i]
      char.3 <- mt[,5][i]
      char.4 <- mt[,6][i]
      for (j in 1:uni.tag1){
        fattore <- levels(factor(mt[,3], as.character(unique(mt[,3]))))[j]
        if (char.1 == fattore){
          col.tag1[i] <- color.tag1[j]}}
      for (j in 1:uni.tag2){
        fattore <- levels(factor(mt[,4], as.character(unique(mt[,4]))))[j]
        if (char.2 == fattore){
          col.tag2[i] <- color.tag2[j]}}
      for (j in 1:uni.tag3){
        fattore <- levels(factor(mt[,5], as.character(unique(mt[,5]))))[j]
        if (char.3 == fattore){
          col.tag3[i] <- color.tag3[j]}}
      for (j in 1:uni.tag4){
        fattore <- levels(factor(mt[,6], as.character(unique(mt[,6]))))[j]
        if (char.4 == fattore){
          col.tag4[i] <- color.tag4[j]}}}
    
    tit_col.1 = paste0("'", colnames(mt)[3],"' ", "tag colors along your ", length(fSmeta()), " samples")
    tit_col.2 = paste0("'", colnames(mt)[4],"' ", "tag colors along your ", length(fSmeta()), " samples")
    tit_col.3 = paste0("'", colnames(mt)[5],"' ", "tag colors along your ", length(fSmeta()), " samples")
    tit_col.4 = paste0("'", colnames(mt)[6],"' ", "tag colors along your ", length(fSmeta()), " samples")
    
    tit_q.1 = paste0("'", colnames(mt)[3],"' ", "tag event distribution")
    tit_q.2 = paste0("'", colnames(mt)[4],"' ", "tag event distribution")
    tit_q.3 = paste0("'", colnames(mt)[5],"' ", "tag event distribution")
    tit_q.4 = paste0("'", colnames(mt)[6],"' ", "tag event distribution")
    
    sample_event <- vector(mode = "numeric", length = length(fSmeta()))
    tag1_event <- vector(mode = "numeric", length = uni.tag1)
    tag2_event <- vector(mode = "numeric", length = uni.tag2)
    tag3_event <- vector(mode = "numeric", length = uni.tag3)
    tag4_event <- vector(mode = "numeric", length = uni.tag4)
    
    for (i in 1:length(fSmeta())){sample_event[i] <- dim(exprs(fSmeta()[[i]]))[1]}
    med <- median(sample_event)
    mad <- mad(x = sample_event, center = median(sample_event)) #median absolute deviation
    # In univariate statistics, the Median Absolute Deviation is the most robust dispersion/scale measure in presence of 
    #outliers -> median +|- 2.5*MAD 
    # is a good method for outlier detection
    sum <- 0
    for (i in 1:uni.tag1){
      subs_tag1<- subset(mt, mt[,3] == unique(mt[,3])[i], select=c(file_name,sample_id))
      for(j in 1:length(subs_tag1$sample_id)){sum <- sum + dim(exprs(fSmeta()[[as.numeric(rownames(subs_tag1))[j]]]))[1]}
      tag1_event[i] <- sum
      names(tag1_event) <- unique(mt[,3])
      sum <- 0}

    for (i in 1:uni.tag2){
      subs_tag2<- subset(mt, mt[,4] == unique(mt[,4])[i], select=c(file_name,sample_id))
      for(j in 1:length(subs_tag2$sample_id)){sum <- sum + dim(exprs(fSmeta()[[as.numeric(rownames(subs_tag2))[j]]]))[1]}
      tag2_event[i] <- sum
      names(tag2_event) <- unique(mt[,4])
      sum <- 0}
    
    for (i in 1:uni.tag3){
      subs_tag3<- subset(mt, mt[,5] == unique(mt[,5])[i], select=c(file_name,sample_id))
      for(j in 1:length(subs_tag3$sample_id)){sum <- sum + dim(exprs(fSmeta()[[as.numeric(rownames(subs_tag3))[j]]]))[1]}
      tag3_event[i] <- sum
      names(tag3_event) <- unique(mt[,5])
      sum <- 0}
    
    for (i in 1:uni.tag4){
      subs_tag4<- subset(mt, mt[,6] == unique(mt[,6])[i], select=c(file_name,sample_id))
      for(j in 1:length(subs_tag4$sample_id)){sum <- sum + dim(exprs(fSmeta()[[as.numeric(rownames(subs_tag4))[j]]]))[1]}
      tag4_event[i] <- sum
      names(tag4_event) <- unique(mt[,6])
      sum <- 0}
    if ((uni.tag1>1)||(uni.tag2>1)||(uni.tag3>1)||(uni.tag3>1)){
      output$condition_tag <- reactive({input$condition_tag1=="show"})} #this is the condition where no tag is set 
                                                                        #(typically when only one sample is loaded) )
  
    if(uni.tag1>1){
      output$condition_tag1 <- reactive({input$condition_tag1=="show"})
      outputOptions(output, "condition_tag1", suspendWhenHidden = FALSE)
      output$tag1.color <- renderPlot({barplot(height = v.tag1, col=col.tag1, axes = FALSE, space = 0.1, main = tit_col.1)})
      output$tag1.bar <- renderPlot({ xx <- barplot(height = tag1_event, xlab = "grouped sample", legend = TRUE, 
                                                    col=color.tag1, main = tit_q.1) ; 
        axis(3, at = xx, labels=tag1_event, tick = TRUE, pos = c(1,1), cex.axis=1.2)})}
    if(uni.tag2>1){
      output$condition_tag2 <- reactive({input$condition_tag2=="show"})
      outputOptions(output, "condition_tag2", suspendWhenHidden = FALSE)
      output$tag2.color <- renderPlot({barplot(height = v.tag2, col=col.tag2, axes = FALSE, space = 0.1, main = tit_col.2)})
      output$tag2.bar <- renderPlot({ xx <- barplot(height = tag2_event, xlab = "grouped sample", legend = TRUE, 
                                                    col=color.tag2, main = tit_q.2) ; 
        axis(3, at = xx, labels=tag2_event, tick = TRUE, pos = c(1,1), cex.axis=1.2)})}
    if(uni.tag3>1){
      output$condition_tag3 <- reactive({input$condition_tag3=="show"})
      outputOptions(output, "condition_tag3", suspendWhenHidden = FALSE)
      output$tag3.color <- renderPlot({barplot(height = v.tag3, col=col.tag3, axes = FALSE, space = 0.1, main = tit_col.3)})
      output$tag3.bar <- renderPlot({ xx <- barplot(height = tag3_event, xlab = "grouped sample", legend = TRUE, 
                                                    col=color.tag3, main = tit_q.3) ; 
        axis(3, at = xx, labels=tag3_event, tick = TRUE, pos = c(1,1), cex.axis=1.2)})}
    if(uni.tag4>1){
      output$condition_tag4 <- reactive({input$condition_tag4=="show"})
      outputOptions(output, "condition_tag4", suspendWhenHidden = FALSE)
      output$tag4.color <- renderPlot({barplot(height = v.tag4, col=col.tag4, axes = FALSE, space = 0.1, main = tit_col.4)})
      output$tag4.bar <- renderPlot({ xx <- barplot(height = tag4_event, xlab = "grouped sample", legend = TRUE, 
                                                    col=color.tag4, main = tit_q.4) ; 
        axis(3, at = xx, labels=tag4_event, tick = TRUE, pos = c(1,1), cex.axis=1.2)})}
  
    metaResume(capture.output(list(sample_event, tag1_event, tag2_event, tag3_event, tag4_event)))
    meta.table(mt)
    output$console_output_meta <- renderPrint("sample metadata loaded")
    if (Sys.info()["sysname"] == "Windows") play(x = dingdong)
  })
  
  output$hot_meta <- renderRHandsontable(
    expr = {
      if (!is.null(meta.table())){
        rhandsontable(data = mt)  %>% 
          hot_context_menu(
            allowColEdit = FALSE,
            useTypes = FALSE,
            customOpts = list(
              csv = list(name = "Download to CSV",
                         callback = htmlwidgets::JS(
                           "function (key, options) {
                           var csv = csvString(this);
                           var link = document.createElement('a');
                           link.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(csv));
                           link.setAttribute('download', 'data.csv');
                           document.body.appendChild(link);
                           link.click();
                           document.body.removeChild(link);
                           }")
              ))) %>%
          hot_col(c(1,2), readOnly = TRUE)
      }
      else 
        rhandsontable(data = mt)
    })
  
  output$downloadmetaResume <- downloadHandler(
    filename = function() {paste("meta_resume","csv", sep = ".")},
    content <- function(file) {
      capture.output(metaResume(), file = "./tmpdata/meta_resume.csv")
      file.copy("./tmpdata/meta_resume.csv",file)},
    contentType = "text/csv")
  
  
  ################################################### runMetaMarkers -----    
  
  meta.markers <- reactiveVal(NULL)
  markers.csv <- reactiveVal(NULL)
  
  observeEvent(eventExpr = input$runMetaMarkers, handlerExpr =  {
    
    validate(need(expr = !is.null(fSmeta()), message = "It needs valid FCS files"))
    markers.csv(input$marker_file)
    toggle <- ifelse(test = !is.null(markers.csv()), yes = T, no = F)
    if (toggle) {
      options(warn = 2) # Turn warnings into errors so they can be trapped
      result <- try(expr = read.csv(file=input$marker_file$datapath, header=TRUE, sep=",", 
                                    colClasses = c("character", "character", "character", "logical")))
      if (inherits(x = result,"data.frame")){
        options(warn = 0)
        mk <<- result
        good <- TRUE
        mk_ok <- TRUE}
      else {  # Process any error messages
        # Ignore warnings while processing errors
        options(warn = -1)
        msg <- geterrmessage()
        mk <<- paste("The table you loaded seems not to be in the correct format. System reports: ", msg, sep = " ")
        # Restore default warning reporting
        options(warn=0)
        good <- FALSE
        mk_ok <- FALSE}
      if (good == FALSE){
        options(warn = 2) # Let's try with semicolon (;)
        result <- try(expr = read.csv(file=input$marker_file$datapath, header=TRUE, sep=";", 
                                      colClasses = c("character", "character", "character", "logical")))
        if (inherits(x = result,"data.frame")){
          options(warn = 0)
          mk <<- result
          good <- TRUE
          mk_ok <- TRUE}
        else {  # Process any error messages
          # Ignore warnings while processing errors
          options(warn = -1)
          msg <- geterrmessage()
          mk <<- paste("The table you loaded seems not to be in the correct format. System reports: ", msg, sep = " ")
          # Restore default warning reporting
          options(warn=0)
          good <- FALSE
          mk_ok <- FALSE }}}
    else {
      mk <<- edit_markers(flow_Set = fSmeta())
      mk_ok <- TRUE}
    
    if (mk_ok){
      meta.markers(mk)}
    else
    {if (!(inherits(x = mk,"data.frame"))){
      output$checktxt_loading_meta_marker <- renderText({
        {messaggio <- paste0("something is wrong with the loaded marker-table. System reports: ", mk, 
                             ". You may start downloading meta-table automatically ",         
                             "generated when you do not provide any table in input. 
                             Try to manipulate that one according to your use and try with this new version")
        validate(need(expr =  (mk_ok), message = messaggio))}})}
      else{
        output$checktxt_loading_meta_marker <- renderText({
          {messaggio <- paste0("something is wrong with the loaded marker-table.", "You may start downloading meta-table automatically ",         
                               "generated when you do not provide any table in input. 
                               Try to manipulate that one according to your use and try with this new version")
          validate(need(expr =  (mk_ok), message = messaggio))}})}}
  })
  
  ################################################### runMetaMarkersLoad -----    
  
  observeEvent(eventExpr = input$runMetaMarkersLoad, handlerExpr = {
    
    validate(need(expr = !is.null(meta.markers()), message = "It needs valid marker entries"))
    hotmarker <- input$hot_markers
    if (!is.null(hotmarker))
   
    {mk <<- hot_to_r(hotmarker)
    
     fF <- fSmeta()[[1]] 
     nomi <- fF@parameters@data$name
     
     nomi_proibiti <- grep("FSC-|SSC-|cell_Id|score|Time|SampleID|clustering", nomi, value = FALSE, ignore.case = TRUE) 
     # there is no "density", in case this is present in the list of channels
     if (length(nomi_proibiti)>0) {nomi <- nomi[-nomi_proibiti]}
     
     nomi_marcatori <- mk$dimension_name
     m <- match(nomi_marcatori, nomi)
     
     validate(need(expr = (all(!is.na(m))), message = "please check channels names"))
     
    uni.marker <<- length(mk$selected[mk$selected==TRUE])
    if (uni.marker == 0){
      mk$selected <<- as.logical(toupper(mk$selected))
      uni.marker <<- length(mk$selected[mk$selected==TRUE])}
    
    validate(need(expr = (uni.marker > 1), message = "please select at least two markers"))
    markers_desc <- exprs_sub(flow.Set = fSmeta())[[3]] #marker description
    
    marcatori_in_tabella <- mk$marker_description
    n_p <- grep("density",marcatori_in_tabella, value = FALSE, ignore.case = TRUE) 
    if (length(n_p)>0) {marcatori_in_tabella <- marcatori_in_tabella[-n_p]}
    #this is because I did not check for density before. The density does not provide any description when generated by the spade routine
    
    check_in <- identical(x = marcatori_in_tabella, y = markers_desc) 
    validate(need(expr = (check_in), message = "the entered table does not match the marker's description of the entry flowSet")) 
    #???to be further studied
    
    totale <- uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4 + 50 + uni.marker
    #the input$MetaClustNum parameter is in another tab (the Hign dim analysis tag). 
    #Let's assume to always assigne the maximum value (50) for input$MetaClustNum
    
    colori <- color_set[1:totale]
    #in order to have the same set of colors for the meta_clusters each time you perform the meta-clustering
    color.metaclust <<- colori[1:50] 
    # it was color.metaclust <<- colori[(uni.tag1 + uni.tag2 + uni.tag3 + uni.marker + uni.pheno + 1):length(colori)]
    
    color.marker <<- colori[(50 + uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4 + 1):(50 + uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4 + uni.marker)]}
    
    meta.markers(mk)
    output$console_output_meta <- renderPrint("marker's metadata loaded")
    if (Sys.info()["sysname"] == "Windows") play(x = dingdong)
  })
  
  output$hot_markers <- renderRHandsontable(
    expr = {
      if (!is.null(meta.markers())){
        rhandsontable(data = mk)  %>% 
          hot_context_menu(
            allowColEdit = FALSE,
            useTypes = FALSE,
            customOpts = list(
              csv = list(name = "Download to CSV",
                         callback = htmlwidgets::JS(
                           "function (key, options) {
                         var csv = csvString(this);
                         var link = document.createElement('a');
                         link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                           encodeURIComponent(csv));
                         link.setAttribute('download', 'data.csv');
                       document.body.appendChild(link);
                         link.click();
                         document.body.removeChild(link);
                       }"
                         ))))}
      else 
        rhandsontable(data = mk)
    })
  
  
  ################################################### runPhenoTable -----    
  
  pheno.table <- reactiveVal(NULL)
  pheno.csv <- reactiveVal(NULL)
  pheno_state <- reactiveValues(state = NULL) #$state fileInput(inputId = "pheno_file", label = h4("pheno-data table")...)
  
  observeEvent(eventExpr = input$runPhenoTable, handlerExpr =  {
   
    validate(need(expr = !is.null(fSmeta()), message = "It needs valid FCS files"))
    validate(need(expr = !is.null(meta.markers()), message = "It needs valid marker entries"))
    
    pheno.csv(input$pheno_file)
    toggle <- ifelse(test = !is.null(pheno.csv()), yes = T, no = F)
    if (toggle) {
      options(warn = 2) # Turn warnings into errors so they can be trapped
      
      result <- try(expr = read.csv(file=input$pheno_file$datapath, header=TRUE, sep=",", check.names = FALSE))
      if (inherits(x = result,"data.frame")){
        options(warn = 0)
        ph <<- result
        good <- TRUE
        ph_ok <- TRUE
        pheno_state$state <- "uploaded"} #state
      else {  # Process any error messages
        # Ignore warnings while processing errors
        options(warn = -1)
        msg <- geterrmessage()
        ph <<- paste("The table you loaded seems not to be in the correct format. System reports: ", msg, sep = " ")
        # Restore default warning reporting
        options(warn=0)
        ph_ok <- FALSE
        good <- FALSE}
      if (good == FALSE){
        options(warn = 2) # Let's try with semicolon (;)
        result <- try(expr = read.csv(file=input$pheno_file$datapath, header=TRUE, sep=";", check.names = FALSE))
        if (inherits(x = result,"data.frame")){
          options(warn = 0)
          ph <<- result
          good <- TRUE
          ph_ok <- TRUE
          pheno_state$state <- "uploaded"} #state}
        else {  # Process any error messages
          # Ignore warnings while processing errors
          options(warn = -1)
          msg <- geterrmessage()
          ph <<- paste("The table you loaded seems not to be in the correct format. System reports: ", msg, sep = " ")
          # Restore default warning reporting
          options(warn=0)
          ph_ok <- FALSE
          good <- FALSE}}}
    else {
      ph <<- edit_pheno(flow_Set = fSmeta(), selected = meta.markers()$selected)
      pheno.table(ph)
      ph_ok <- TRUE}
    
    if (ph_ok){
      pheno.table(ph)}
    else
    {if (!(inherits(x = ph,"data.frame"))){
      output$checktxt_loading_meta_pheno <- renderText({
        {messaggio <- paste0("something is wrong with the loaded pheno-table. System reports: ", ph, 
                             ". You may start downloading meta-table automatically ",         
                             "generated when you do not provide any table in input. 
                             Try to manipulate that one according to your use and try with this new version")
        validate(need(expr =  (ph_ok), message = messaggio))}})}
      else{
        output$checktxt_loading_meta_pheno <- renderText({
          {messaggio <- paste0("something is wrong with the loaded pheno-table.", "You may start downloading meta-table automatically ",         
                               "generated when you do not provide any table in input. 
                               Try to manipulate that one according to your use and try with this new version")
          validate(need(expr =  (ph_ok), message = messaggio))}})}}
    
  })

  ################################################### runPhenoTableLoad -----    
  observeEvent(eventExpr = input$runPhenoTableload, handlerExpr = {
    
    validate(need(expr = !is.null(meta.markers()), message = "It needs valid marker entries"))
    hotpheno <- input$hot_pheno
    ph <<- hot_to_r(hotpheno)
    if (!is.null(hotpheno))
    {
      validate(need(expr = !is.null(pheno.table()), message = "It needs a valid pheno.table"))
      marker_selected <- mk[mk$selected,]$marker_description
      marker_phenotable <- colnames(ph)[-1]
      validate(need(expr = (setequal(marker_selected,marker_phenotable)), message = "It needs a valid pheno.table")) #to check the markers
      good_char <- as.vector(as.matrix(ph[,-1]))
      elementi <- length(good_char)
      plus <- length(which(good_char == '+'))
      menus <- length(which(good_char == '-'))
      mult <- length(which(good_char == '*'))
      validate(need(expr = ((plus + menus + mult) == elementi), message = "It needs a valid pheno.table")) 
      #to check the correctness of the entries
      
      uni.pheno <<- nrow(ph) + 1 + 1 #you could have also the unassigned color and the overlap color
      validate(need(expr = (uni.pheno > 1), message = "please enter some phenotypes"))
      
      totale <- uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4 + 50 + uni.marker + uni.pheno
      #the input$MetaClustNum parameter is in another tab (the Hign dim analysis tag). 
      #Let's assume to always assigne the maximum value (50) for input$MetaClustNum
      
      colori <- color_set[1:totale]
      #in order to have the same set of colors for the meta_clusters each time you perform the meta-clustering
      color.metaclust <<- colori[1:50] 
      # it was color.metaclust <<- colori[(uni.tag1 + uni.tag2 + uni.tag3 + uni.marker + uni.pheno + 1):length(colori)]
      
      color.pheno <<- colori[(50 + uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4 + uni.marker + 1):
                               (50 + uni.tag1 + uni.tag2 + uni.tag3 + uni.tag4 + uni.marker + uni.pheno)]
      
      v.tag <- rep.int(x = 1, times = uni.pheno)
      name.tag <- sort(append(as.character(ph[,1]), "UNASSIGNED"))
      #name.tag <- append(as.character(ph[,1]), "unassigned")
      names(v.tag) <- name.tag
      
      #color.pheno[length(color.pheno) - 1] <<- '#D3D3D3' #lightgrey
      pos <-  which(name.tag %in% "UNASSIGNED")
      color.pheno[pos] <<- '#D3D3D3' #lightgrey
      names(color.pheno) <<- name.tag
      
      name.tag <- append(name.tag, "OVERLAP")
      names(v.tag) <- name.tag
      
      color.pheno[length(color.pheno)] <<- '#000000' #black
      names(color.pheno) <<- name.tag
      
      if(uni.pheno>1){
        output$pheno.color <- renderPlot({barplot(height = v.tag, col=color.pheno, axes = FALSE, space = 0.1, main = "phenotype colors")})}}
      
      output$meta_phenoStatus <- reactive({input$meta_phenoStatus=="show"})
      outputOptions(output, "meta_phenoStatus", suspendWhenHidden = FALSE)
      
    pheno.table(ph)
    output$console_output_meta <- renderPrint("phenotype metadata loaded")
    if (Sys.info()["sysname"] == "Windows") play(x = dingdong)
  })
  
  output$hot_pheno <- renderRHandsontable(
    expr = {
      colwidth <- rep(x=70, ncol(ph))
      colwidth[1] <- 200 
      if (!is.null(pheno.table())){
        rhandsontable(data = ph)  %>% 
          hot_context_menu(
            allowColEdit = FALSE, 
            useTypes = FALSE,
            customOpts = list(csv = list(name = "Download to CSV",
                                                        callback = htmlwidgets::JS(
                                                          "function (key, options) {
                                                        var csv = csvString(this);
                                                        var link = document.createElement('a');
                                                        link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                                                        encodeURIComponent(csv));
                                                        link.setAttribute('download', 'data.csv');
                                                        document.body.appendChild(link);
                                                        link.click();
                                                        document.body.removeChild(link);}")))) %>% 
          hot_cols(colWidths = colwidth, manualColumnResize = TRUE) %>% 
          hot_col(col = names(ph)[-1], allowInvalid = FALSE, type = "dropdown")
        #hot_row(c(1,2,3), readOnly = TRUE)
      }else 
        {rhandsontable(data = ph)}
    })
  
  ################################################### resetMetaTable -----    
  
  observeEvent(eventExpr = input$resetMetaTable, handlerExpr = {
    
    pheno_state$state <- "reset" #$state
    ph <<- edit_pheno()
    pheno.table(ph)
    output$hot_pheno <- renderRHandsontable(
      expr = {
        colwidth <- rep(x=70, ncol(ph))
        colwidth[1] <- 200 
        if (!is.null(pheno.table())){
          rhandsontable(data = ph)  %>% 
            hot_context_menu(allowRowEdit = TRUE, allowColEdit = TRUE) %>% 
            hot_context_menu(customOpts = list(csv = list(name = "Download to CSV",
                                                          callback = htmlwidgets::JS(
                                                            "function (key, options) {
                                                        var csv = csvString(this);
                                                        var link = document.createElement('a');
                                                        link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                                                        encodeURIComponent(csv));
                                                        link.setAttribute('download', 'data.csv');
                                                        document.body.appendChild(link);
                                                        link.click();
                                                        document.body.removeChild(link);}")))) %>% 
            hot_cols(colWidths = colwidth, manualColumnResize = TRUE) %>% 
            hot_col(col = names(ph)[-1], allowInvalid = FALSE, type = "dropdown")
          #hot_row(c(1,2,3), readOnly = TRUE)
        }
        else {rhandsontable(data = ph)}})
  })
  
  pheno_file_input <- reactive({
    if (is.null(pheno_state$state)) {
      return(NULL)
    } else if (pheno_state$state  == 'uploaded') {
      return(input$pheno_file)
    } else if (pheno_state$state  == 'reset') {
      return(NULL)
    }
  })
  
  output$pheno_file_status <- renderText({
    return(paste("Uploaded file:", pheno_file_input()$name))})
 
  ################################################### ____flowSet Analysis -----    
  ################################################### runflowSetAna -----    
  
  dI.ana <- reactiveValues(MDS.tag1 = NULL, MDS.tag2 = NULL, MDS.tag3 =  NULL, MDS.tag4 =  NULL,
                           dia.tag1 = NULL, dia.tag2 = NULL, dia.tag3 = NULL, dia.tag4 = NULL)
  
  dI.ana_val <- reactiveVal(NULL)
  observeEvent(eventExpr = input$runflowSetAna, handlerExpr = {
  
    if((inherits(x = fSmeta(),"flowSet"))&&(!is.null(meta.table()))&&(any(mk$selected, na.rm = TRUE))) {
      print("start showing MDS plots")
      output$console_output_meta <- renderPrint("start showing MDS plots - waiting...")
      
      if(length(fSmeta())>2){
        output$flowSetAna_panelStatus <- reactive({input$flowSetAna_panelStatus=="show"})
        outputOptions(output, "flowSetAna_panelStatus", suspendWhenHidden = FALSE)
        dI.ana$MDS.sample_id <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), type = "sample_id", 
                                            selected = mk$selected, metodo = input$mdsMethod) 
        if (uni.tag1>1){dI.ana$MDS.tag1 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                       type = "tag1", selected = mk$selected, metodo = input$mdsMethod)}
        if (uni.tag2>1){dI.ana$MDS.tag2 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                       type = "tag2", selected = mk$selected, metodo = input$mdsMethod)}
        if (uni.tag3>1){dI.ana$MDS.tag3 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                       type = "tag3", selected = mk$selected, metodo = input$mdsMethod)}
        if (uni.tag4>1){dI.ana$MDS.tag4 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                       type = "tag4", selected = mk$selected, metodo = input$mdsMethod)}}
      else
        {dI.ana_val("To get an MDS plot you need at least 3 samples. Load three '.fcs' files or skip this type of analysis")}
    
      if (uni.tag1>1) {dI.ana$dia.tag1 <- diagnostic_plot(flow_Set = fSmeta(), meta = meta.table(), type = "tag1", selected = mk$selected)}
      if (uni.tag2>1) {dI.ana$dia.tag2 <- diagnostic_plot(flow_Set = fSmeta(), meta = meta.table(), type = "tag2", selected = mk$selected)}
      if (uni.tag3>1) {dI.ana$dia.tag3 <- diagnostic_plot(flow_Set = fSmeta(), meta = meta.table(), type = "tag3", selected = mk$selected)}
      if (uni.tag4>1) {dI.ana$dia.tag4 <- diagnostic_plot(flow_Set = fSmeta(), meta = meta.table(), type = "tag4", selected = mk$selected)}
     
      if(inherits(x = (dI.ana$MDS.sample_id)[[1]],"plotly")){
      output$MDSsample_id <- renderPlotly({dI.ana$MDS.sample_id[[1]]})
      output$MDS_HD <- renderPlot({dI.ana$MDS.sample_id[[2]]})
      if (uni.tag1>1) {output$MDStag1 <- renderPlotly({dI.ana$MDS.tag1[[1]]})}
      if (uni.tag2>1) {output$MDStag2 <- renderPlotly({dI.ana$MDS.tag2[[1]]})}
      if (uni.tag3>1) {output$MDStag3 <- renderPlotly({dI.ana$MDS.tag3[[1]]})}
      if (uni.tag4>1) {output$MDStag4 <- renderPlotly({dI.ana$MDS.tag4[[1]]})}
      
      if (uni.tag1>1) {output$dia_tag1 <- renderPlot({dI.ana$dia.tag1})}
      if (uni.tag2>1) {output$dia_tag2 <- renderPlot({dI.ana$dia.tag2})}
      if (uni.tag3>1) {output$dia_tag3 <- renderPlot({dI.ana$dia.tag3})}
      if (uni.tag4>1) {output$dia_tag4 <- renderPlot({dI.ana$dia.tag4})}
      
      print("MDS process ends")
      output$console_output_meta <- renderPrint("showing MDS plots ends")}
      else{
          result <- dI.ana$MDS.sample_id
          output$checktxt_plotting_mds <- renderText({
            if (!is.null(dI.ana$MDS.sample_id)) 
              {validate(need(expr =  (inherits(x = (dI.ana$MDS.sample_id)[[1]],"plotly")), message = dI.ana$MDS.sample_id))}})
          output$checktxt_plotting_mds <- renderText({if (inherits(x = dI.ana_val(),"character"))
            validate(need(expr =  (!(inherits(x = dI.ana_val(),"character"))), message = dI.ana_val()))})}}
    else {
      dI.ana_val("please load your FCS files, the sample's and the marker's metatables in order to be parsed from cytoChain")
      
      output$checktxt_plotting_mds <- renderText({
        if (!is.null(dI.ana_val())){
          if (inherits(x = dI.ana_val(),"character")) 
          {messaggio <- dI.ana_val()
          validate(need(expr =  (!(inherits(dI.ana_val(),"character"))), message = messaggio))}}})}
  })
  
  
  ################################################### ____flowSet evaluation  -----    
  #################################################################### runPCAEva ----  
  
  dI.PCA.eva <- reactiveValues(pca.res = NULL, pca.df = NULL, pca2.df = NULL, cos2_plot = NULL, cos2_var = NULL)
  
  observeEvent(eventExpr = input$runPCAEva, handlerExpr = {
   
    if((inherits(x = fSmeta(),"flowSet"))&&(!is.null(meta.table()))&&(any(mk$selected, na.rm = TRUE))) {
    
    print("start showing PCA evaluation plots")
    output$console_output_meta <- renderPrint("start showing PCA evaluation plots - waiting...")
    output$flowSetEva_panelStatus <- reactive({input$flowSetEva_panelStatus=="show"})
    outputOptions(output, "flowSetEva_panelStatus", suspendWhenHidden = FALSE)
    
    pca_list <- evaluation_plotPCA(flow_Set = fSmeta(), seed = input$set_seed, selected = mk$selected)
    dI.PCA.eva$pca.res <- pca_list[[1]]
    dI.PCA.eva$pca.df <- pca_list[[2]]
    
    dI.PCA.eva$pca2.df <- pca_list[[3]]
    dI.PCA.eva$cos2_plot <- pca_list[[4]]
    dI.PCA.eva$cos2_var <- pca_list[[5]]
    
    print("evaluation PCA plot process ends")
    
    output$console_output_meta <- renderPrint("PCA evaluation plots ends")
    output$pca_res <- renderPlot({dI.PCA.eva$pca.res})
    output$pca_df <- renderTable({dI.PCA.eva$pca.df})
    
    output$pca2_df <- renderTable({dI.PCA.eva$pca2.df})
    output$cos2plot <- renderPlot({dI.PCA.eva$cos2_plot})
    output$cos2var <- renderPlot({dI.PCA.eva$cos2_var})}

    else{
      dI.PCA.eva$pca.res <- "please load your FCS files, the sample's and the marker's metatables in order to be parsed from cytoChain"
      output$checktxt_runPCAEva <- renderText({
        validate(need(expr =  (!(inherits(x = dI.PCA.eva$pca.res,"character"))), message = dI.PCA.eva$pca.res))})}
  })
  
  
  #################################################################### runKmeansEva ----
  dI.kmeans.eva <- reactiveValues(kmeans.out = NULL, kmeans.df = NULL)
  
  observeEvent(eventExpr = input$runKmeansEva, handlerExpr = {
    
    if(inherits(x = fSmeta(),"flowSet")) {
    validate(need(expr = !is.null(fSmeta()), message = "It needs valid FCS files"))
    
    print("start showing kmeans evaluation plots")
    output$console_output_meta <- renderPrint("start showing kmeans evaluation plots - waiting...")
    output$runKmeansEva_panelStatus <- reactive({input$runKmeansEva_panelStatus=="show"})
    outputOptions(output, "runKmeansEva_panelStatus", suspendWhenHidden = FALSE)
    
    kmeans_list <- evaluation_plotKmeans(data_Set = fSmeta(),  k_max = input$kmeans_cluster, selected = mk$selected)
    dI.kmeans.eva$kmeans.out <- kmeans_list[[1]]
    dI.kmeans.eva$kmeans.df <- kmeans_list[[2]]
    
    output$kmeans_out <- renderPlot({dI.kmeans.eva$kmeans.out})
    output$kmeansTable_out <- renderTable({dI.kmeans.eva$kmeans.df})
    
    print("evaluation plot kmeans process ends")
    output$console_output_meta <- renderPrint("evaluation kmeans plot process ends")
    if (Sys.info()["sysname"] == "Windows") play(x = gong)}
  })
  
  
  
  #################################################################### runMclustEva ----
  
  dI.mclust.eva <- reactiveValues(mclust.out = NULL, mclust.df = NULL)
  
  observeEvent(eventExpr = input$runMclustEva, handlerExpr = {
    
    validate(need(expr = !is.null(fSmeta()), message = "It needs valid FCS files"))
    
    print("start showing mclust evaluation plots")
    output$console_output_meta <- renderPrint("start showing mclust evaluation plots - waiting...")
    output$runMclustEva_panelStatus <- reactive({input$runMclustEva_panelStatus=="show"})
    outputOptions(output, "runMclustEva_panelStatus", suspendWhenHidden = FALSE)
    
    mclust_list <- evaluation_Mclust(flow_Set = fSmeta(),  g_max = input$max_mclust, selected = mk$selected)
    dI.mclust.eva$mclust.out <- mclust_list[[1]]
    dI.mclust.eva$mclust.df <- mclust_list[[2]]
    
    output$mclust_out <- renderPlot({dI.mclust.eva$mclust.out})
    output$mclustTable_out <- renderTable({dI.mclust.eva$mclust.df})
    
    print("evaluation plot with mclust process ends")
    output$console_output_meta <- renderPrint("evaluation mclust plots ends")
    if (Sys.info()["sysname"] == "Windows") play(x = gong)
  })
  
  
  
  #################################################################### _____Map evaluation ----
  #################################################################### runfStSNEinit ----
  dI.tSNE.eva.init <- reactiveValues(mapres_init =  NULL, mapplot_init = NULL)
  observeEvent(eventExpr = input$runfStSNEinit, handlerExpr = {
 
    if (!(inherits(x = dI(),"flowSet"))&&(any(mk$selected, na.rm = TRUE))&&(inherits(x = fSmeta(),"flowSet"))) 
      #case of samples loaded from metadata workflow
    {fS.init <- fSmeta()}
    else
    {fS.init <- dI()}
    
    if(inherits(x = fS.init,"flowSet")){
      print("start generating tSNE map for the selected flowSet")
      output$console_output_meta <- renderPrint("start showing tSNE plots for the selected flowSet - waiting...")
      output$tSNEinit_panelStatus <- reactive({input$tSNEunique_panelStatus=="show"})
      outputOptions(output, "tSNEinit_panelStatus", suspendWhenHidden = FALSE)
        
      dI.tSNE.eva.init$mapres_init <- map_gen(flow_Set = fS.init, type = "tSNE", selected_markers = mk$selected, 
                                              seed_nr = input$set_seed, limit = FALSE, maxCell = input$mapMax,  
                                              Rtsne_pca = input$tSNE_PCA, Rtsne_perp = input$perplexity, Rtsne_theta = input$theta, 
                                              Rtsne_iter = input$tSNEiter, Rtsne_eta = input$tSNEeta)
    
      if(inherits(x = dI.tSNE.eva.init$mapres_init,"character")){
          dI.tSNE.eva.init$mapres_init <- paste0("map generation function says: ", dI.tSNE.eva.init$mapres_init)}
      
      if(inherits(x = dI.tSNE.eva.init$mapres_init,"list")){  
      dI.tSNE.eva.init$mapplot_init <- map_plot_comp(map_df = dI.tSNE.eva.init$mapres_init[[1]], 
                                                       map_inds = dI.tSNE.eva.init$mapres_init[[2]], type = "tortazza")
      #output$mapunique_out <- renderPlot({dI.tSNE.eva.init$mapplot_init})
      #print("tSNE map generation for selected flowset process ends")
      #output$console_output_meta <- renderPrint("tSNE evaluation plot for selected flowset process ends")
    
      output$mapinit_out <- renderPlot({dI.tSNE.eva.init$mapplot_init})
      print("tSNE map generation for entry flowset process ends")
      output$console_output_meta <- renderPrint("tSNE evaluation plot for entry flowset process ends")
      if (Sys.info()["sysname"] == "Windows") play(x = gong)}}
    
    if ((!(inherits(x = fS.init,"flowSet")))||!(any(mk$selected, na.rm = TRUE)))
      {dI.tSNE.eva.init$mapres_init <- "Please select a valid flowSet and the related metadata"}
    
    output$checktxt_eva_init <- renderText({
      if (!(inherits(x = fS.init,"flowSet"))&&(!(inherits(x = dI(),"flowSet"))))
        {validate(need(expr =  (inherits(x = fS.init,"flowSet")), message = dI.tSNE.eva.init$mapres_init))}
      if (inherits(x = dI.tSNE.eva.init$mapres_init,"character"))
        {validate(need(expr =  !(inherits(x = dI.tSNE.eva.init$mapres_init,"character")), message = dI.tSNE.eva.init$mapres_init))}})
  })
  
  #################################################################### runfStSNEfinal ----
  dI.tSNE.eva.final <- reactiveValues(mapres_final = NULL, mapplot_final = NULL)
  
  observeEvent(eventExpr = input$runfStSNEfinal, handlerExpr = {
      
    if (!(inherits(x = dI(),"flowSet")))
      fS.init <- "please go back to the 'Loading & parsing samples' menu and select a flowSet to handle"
    else{ 
      fS.final <- NULL
      if (inherits(x = dI(),"flowSet")) {fS.final <- dI()}    
      if (inherits(x = dI.H(),"flowSet")) {fS.final <- dI.H()}    
    
      if(inherits(x = fS.final,"flowSet")){
        print("start generating tSNE map for the handled flowSet")
        output$console_output_meta <- renderPrint("start showing tSNE evaluation plots for the final flowSet - waiting...")
      
        output$tSNEfinal_panelStatus <- reactive({input$tSNEfinal_panelStatus=="show"})
        outputOptions(output, "tSNEfinal_panelStatus", suspendWhenHidden = FALSE)
      
        dI.tSNE.eva.final$mapres_final <- map_gen(flow_Set = fS.final, type = "tSNE", two_three = FALSE, selected_markers = mk$selected, 
                                                  seed_nr = input$set_seed, limit = FALSE, maxCell = input$mapMax,  
                                                  Rtsne_pca = input$tSNE_PCA,
                                                  Rtsne_perp = input$perplexity, Rtsne_theta = input$theta, 
                                                  Rtsne_iter = input$tSNEiter, Rtsne_eta = input$tSNEeta)
      
        dI.tSNE.eva.final$mapplot_final <- map_plot_comp(map_df = dI.tSNE.eva.final$mapres_final[[1]], 
                                                     map_inds = dI.tSNE.eva.final$mapres_final[[2]], type = "tortazza")
        output$mapfinal_out <- renderPlot({dI.tSNE.eva.final$mapplot_final})
        print("tSNE evaluation plot for the final flowSet process ends")
        output$console_output_meta <- renderPrint("tSNE evaluation plot for the final flowSet process ends")
        if (Sys.info()["sysname"] == "Windows") play(x = gong)}}
    
  output$checktxt_eva_final <- renderText({
    if (is.null(dI())) 
      {validate(need(expr = (inherits(x = dI(),"flowSet")), 
                    message = "please go back to the 'Loading & parsing samples' tab and select a flowSet to handle"))}})
  })
  
  
  #################################################################### runKtSne -----
  
  dI.tSNE.K <- reactiveValues(tSNEx_init =  NULL, tSNEy_init = NULL, fF_init = NULL, clusters_init = NULL)
  
  observeEvent(eventExpr = input$runKtSne, handlerExpr = {
    
    if (!(inherits(x = dI(),"flowSet"))&&(any(mk$selected, na.rm = TRUE))&&(inherits(x = fSmeta(),"flowSet"))) 
      #case of samples loaded from metadata workflow
    {fS.init <- fSmeta()}
    else
    {fS.init <- dI()}
    
    if (inherits(x = fS.init,"flowSet"))
    {
      print("start flowSet concatenation process for the tSNE map of the handled samples")
      output$console_output_meta <- renderPrint("start flowSet concatenation process for handled samples - waiting...")
      fF_init <- concatenating_fS(flow_Set = fS.init, stringa = "conc_sample")
      print("concatenation process ends")
      output$console_output_meta <- renderPrint("concatenation process for the handled samples ends")
      output$runKtSne_panelStatus <- reactive({input$runKtSne_panelStatus=="show"})
      outputOptions(output, "runKtSne_panelStatus", suspendWhenHidden = FALSE)
      
      dI.tSNE.K$tSNEx_init <- dI.tSNE.eva.init$mapres_init[[1]]$map_x
      dI.tSNE.K$tSNEy_init <- dI.tSNE.eva.init$mapres_init[[1]]$map_y
      
      fF_init <- add_dim(flow_Frame = fF_init, dim_name = "tSNEx", dim_vect = dI.tSNE.K$tSNEx_init)
      dI.tSNE.K$fF_init <- add_dim(flow_Frame = fF_init, dim_name = "tSNEy", dim_vect = dI.tSNE.K$tSNEy_init)
      
      print("tSNEx and tSNEy dimensions added")
      output$console_output_meta <- renderPrint("tSNEx and tSNEy dimensions added")
      
      exprsub <- exprs(dI.tSNE.K$fF_init)
      col_select <- c("tSNEx","tSNEy") 
      exprsub <- subset(exprsub, select=col_select)
      
      if (input$tsne_algo == "k_means"){
        resKmeans_init <- kmeans(exprsub, input$NumClust, nstart=10,iter.max = 200, algorithm = "MacQueen")
        #resKmeans_init <- meanShift(queryData = exprsub, algorithm="KDTREE")
      ### try also kmeans(x = exprsub, centers = input$NumClust, nstart=input$NumClust, iter.max = 200, algorithm = "Hartigan-Wong")
        dI.tSNE.K$clusters_init <- resKmeans_init[[1]]} #with kmeans
        #dI.tSNE.K$clusters_init <- resKmeans_init$assignment} #with meanShift
      
      
      if (input$tsne_algo == "SOM"){
        dim1 <- floor(sqrt(input$NumClust))
        dim2 <- ceiling(sqrt(input$NumClust))
        somres <- som(scale(exprsub), grid = somgrid(xdim = dim1, ydim = dim2, "rectangular"))
        dI.tSNE.K$clusters_init <- somres$unit.classif}
      
      colorCluster <- color_set[1:input$NumClust]
      tSNEdf_init <- data.frame(map_x = dI.tSNE.K$tSNEx_init, map_y = dI.tSNE.K$tSNEy_init)
      
      dI.tSNE.K$tSNE_init <- map_plot_comp(map_df = tSNEdf_init, clust = dI.tSNE.K$clusters_init, 
                                            color_clusters = colorCluster, type = "tortazza")
      output$tSNEkmeans_clust_init <- renderPlot({dI.tSNE.K$tSNE_init})
      output$tSNEkmeans_clust_init_info <- renderText({paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)})
      print("tSNE cluster evaluation plot process ends")
      output$console_output_meta <- renderPrint("tSNE cluster evaluation plot process ends")
      if (Sys.info()["sysname"] == "Windows") play(x = dingdong)

      fF_out <- add_dim(flow_Frame = dI.tSNE.K$fF_init, dim_name = "kmeans", dim_vect = dI.tSNE.K$clusters_init)
      if (inherits(x = dI.clust$mapping,"integer")){
        fF_out <- add_dim(flow_Frame = fF_out, dim_name = "clustering", dim_vect = dI.clust$mapping)
        if (inherits(x = dI.metaClust(),"list"))
        {
          code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass
          cell_clustering <- code_clustering[dI.clust$mapping]
          fF_out <- add_dim(flow_Frame = fF_out, dim_name = "meta_clustering", dim_vect = cell_clustering)}}
      
      dI.tSNE.K$fF_final <- fF_out}
      
      output$downloadKmeansTsne <- downloadHandler(
        filename = function() {
          paste("kTsne", "zip", sep = ".")
        },
        content = function(file){
          
          file_path <- file
          tmpdir <- tempdir()
          setwd(tempdir())
          write.FCS(dI.tSNE.K$fF_final, filename = "kTsne.fcs")
          Zip_Files <- list.files(path = getwd(), pattern = "^kTsne*")
          #Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
          #this is necessary to sort in the correct way the flowFrames in the flowSet
          zip(zipfile = file_path, files = Zip_Files)
          setwd(app_dir)
        },
        contentType = "application/zip")  
  })
  
  
  #################################################################### runKtSne_unique -----
  
  dI.tSNE_unique.K <- reactiveValues(tSNEx_final =  NULL, tSNEy_final = NULL, fF_final = NULL, clusters_final = NULL)
  
  observeEvent(eventExpr = input$runKtSne_unique, handlerExpr = {
    
    if (!(inherits(x = dI(),"flowSet")))
      {fS.init <- "please go back to the 'Loading & parsing samples' menu and select a flowSet to handle"}
    else
      {
        fS.final <- NULL
        if (inherits(x = dI(),"flowSet")) {fS.final <- dI()}    
        if (inherits(x = dI.H(),"flowSet")) {fS.final <- dI.H()}    

        if (inherits(x = fS.final,"flowSet")){
          print("start flowSet concatenation process for the tSNE map of the selected samples")
          output$console_output_meta <- renderPrint("start flowSet concatenation process for selected samples - waiting...")
          output$runKtSneunique_panelStatus <- reactive({input$runKtSneunique_panelStatus=="show"})
          outputOptions(output, "runKtSneunique_panelStatus", suspendWhenHidden = FALSE)
          
          dI.tSNE_unique.K$tSNEx_final <- dI.tSNE.eva.final$mapres_final[[1]]$map_x
          dI.tSNE_unique.K$tSNEy_final <- dI.tSNE.eva.final$mapres_final[[1]]$map_y
          fF_final <- concatenating_fS(flow_Set = fS.final, stringa = "conc_sample")
        
          fF_final <- add_dim(flow_Frame = fF_final, dim_name = "tSNEx", dim_vect = dI.tSNE_unique.K$tSNEx_final)
          dI.tSNE_unique.K$fF_final <- add_dim(flow_Frame = fF_final, dim_name = "tSNEy", dim_vect = dI.tSNE_unique.K$tSNEy_final)
        
          print("tSNEx and tSNEy dimensions added")
          output$console_output_meta <- renderPrint("tSNEx and tSNEy dimensions added")
        
          exprsub <- exprs(dI.tSNE_unique.K$fF_final)
          col_select <- c("tSNEx","tSNEy") 
          exprsub <- subset(exprsub, select=col_select)
          
          if (input$tsne_algo == "k_means"){
    ### try also kmeans(x = exprsub, centers = input$NumClust, nstart=input$NumClust, iter.max = 200, algorithm = "Hartigan-Wong")
          resKmeans_final <- kmeans(exprsub, input$NumClust, nstart=10,iter.max = 200, algorithm = "MacQueen")
          #resKmeans_final <- meanShift(queryData = exprsub, algorithm="KDTREE")
             
          dI.tSNE_unique.K$clusters_final <- resKmeans_final[[1]]} #with kmeans
          #dI.tSNE_unique.K$clusters_final <- resKmeans_final$assignment} #with meanShift
           
          if (input$tsne_algo == "SOM"){
            dim1 <- floor(sqrt(input$NumClust))
            dim2 <- ceiling(sqrt(input$NumClust))
            somres <- som(scale(exprsub), grid = somgrid(xdim = dim1, ydim = dim2, "rectangular"))
            dI.tSNE_unique.K$clusters_final <- somres$unit.classif}
          
    
          colorCluster <- color_set[1:input$NumClust]
          tSNEdf_final <- data.frame(map_x = dI.tSNE_unique.K$tSNEx_final, map_y = dI.tSNE_unique.K$tSNEy_final)
    
          dI.tSNE_unique.K$tSNE_final <- map_plot_comp(map_df = tSNEdf_final, clust = dI.tSNE_unique.K$clusters_final, 
                                                 color_clusters = colorCluster, type = "tortazza")
    
          print("tSNE cluster evaluation plot process ends")
          output$console_output_meta <- renderPrint("tSNE cluster evaluation plot process ends")
          if (Sys.info()["sysname"] == "Windows") play(x = gong)
    
          fF_out <- add_dim(flow_Frame = dI.tSNE_unique.K$fF_final, dim_name = "kmeans", dim_vect = dI.tSNE_unique.K$clusters_final)
        
          if (inherits(x = dI.clust$mapping,"integer")){
            {fF_out <- add_dim(flow_Frame = fF_out, dim_name = "clustering", dim_vect = dI.clust$mapping)}
          
          if (inherits(x = dI.metaClust(),"list")){
            code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass
            cell_clustering <- code_clustering[dI.clust$mapping]
            fF_out <- add_dim(flow_Frame = fF_out, dim_name = "meta_clustering", dim_vect = cell_clustering)}}}}
    
    dI.tSNE_unique.K$fF_final <- fF_out
    output$tSNEkmeans_clust_unique <- renderPlot({
      dI.tSNE_unique.K$tSNE_final})
    output$tSNEkmeans_clust_unique_info <- renderText({paste0("x=", input$uniqe_plot_click$x, "\ny=", input$unique_plot_click$y)})
    
    output$downloadKmeansunique <- downloadHandler(
      filename = function() {
        paste("kTsne", "zip", sep = ".")
      },
      content = function(file){
        
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        write.FCS(dI.tSNE_unique.K$fF_final, filename = "kTsne.fcs")
        Zip_Files <- list.files(path = getwd(), pattern = "^kTsne*")
        #Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) #this is necessary to sort in the correct way the flowFrames in the flowSet
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")  
  })
  
  ################################################### ________Clustering -----    

  dI.clust<- reactiveValues(mapping = NULL, codes = NULL, mst = NULL, NumCluster = NULL, sil_plot = NULL, sil_summary = NULL)
  dI.flow <- reactiveValues(labelI = NULL, meta_clust = NULL, mapping = NULL, labelII = NULL, 
                           pheno_mapsI = NULL, pheno_mapsII = NULL, quantity = NULL, time_step = NULL)
  
  observeEvent(eventExpr = input$runClust, handlerExpr = {
    
    clean_graph()
    del_tmpdata(type = "all")
    
    if (inherits(x = fSmeta(),"flowSet")){
      color.clust <<- color_set[1:input$k] #Warning max ClustNum == 538! max(input$k) is currently 200
      
      if (any(mk$selected, na.rm = TRUE)){
        
        dI.metaClust(NULL) #this is to reset the metaclustering results
        dI.labelClust$label.res <- NULL #this is to reset the labelling I results
        dI.labelClust$hm_cluster <- NULL #this is to reset the labelling I results
        dI.label(NULL) #this is to reset the labelling II results
        dI.labelClust$label.res <- NULL; #dI.map_val$showmapcomp <- NULL; #this to try to avoid the re-generation of the tsne map
        dI.label_val$hm_label <- NULL; dI.label_val$hm_label_plus <- NULL
        dI.label_val$hm_pheno <- NULL; dI.label_val$hm_pheno_plus <- NULL
        
        dI.label_p_val$gg_map <- NULL; dI.label_p_val$gg_map_sample <- NULL;  
        dI.label_p_val$gg_map_tag1 <- NULL;  dI.label_p_val$gg_map_tag2 <- NULL;  dI.label_p_val$gg_map_tag3 <- NULL
        dI.label_p_val$gg_map_tag4 <- NULL
        
        dI.quant_val$bar_plot_sample <- NULL; 
        dI.quant_val$bar_plot_tag1 = NULL; dI.quant_val$bar_plot_tag2 <- NULL; dI.quant_val$bar_plot_tag3 <- NULL
        dI.quant_val$bar_plot_tag4 <- NULL
        dI.quant_val$med_expr_tag1 = NULL; dI.quant_val$med_expr_tag2 <- NULL; dI.quant_val$med_expr_tag3 <- NULL
        dI.quant_val$med_expr_tag4 <- NULL 
        
        dI.flow$labelI <- TRUE; dI.flow$meta_clust <- TRUE; dI.flow$labelII <- TRUE; 
        dI.flow$pheno_mapsI <- TRUE; dI.flow$pheno_mapsII <- TRUE; dI.flow$quantity <- TRUE; dI.flow$time_step <- TRUE
        
        print("start clustering process")
        output$console_output_clust <- renderPrint("start clustering process - waiting...")
       
        if (input$ClustAlgo == "flowSOM"){
          flowSOM.res <- flowSOM_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, seed = input$set_seed)
          if (inherits(x = flowSOM.res[[1]],"FlowSOM")){
            dI.clust$mapping <- as.integer(flowSOM.res[[2]]$map$mapping[,1]) # vector [flowSet nr. of cell = cluster nr.]
            dI.clust$codes <- flowSOM.res[[2]]$map$codes # matrix [metaclust nr. x selected_markers]
            dI.clust$mst = flowSOM.res[[3]]
            dI.clust$NumCluster <- 100 #flowSOM uses always 100 clusters}
            color.clust <<- color_set[1:dI.clust$NumCluster]} else{
              dI.clust$mapping <- flowSOM.res}}
      
      if (input$ClustAlgo == "KMeans"){
        kmeans.res <- kmeans_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, Algo = input$KmeansAlgo,
                                   K = input$k, seed = input$set_seed)
        if (inherits(x = kmeans.res,"list")){
          dI.clust$mapping <- kmeans.res[[1]]
          dI.clust$codes <- kmeans.res[[2]]
          dI.clust$NumCluster <- input$k
          color.clust <<- color_set[1:dI.clust$NumCluster]}else{
            dI.clust$mapping <- kmeans.res}}
        
      if ((input$ClustAlgo == "Rphenograph")||(input$ClustAlgo == "FastPG")){
        if (input$ClustAlgo == "Rphenograph"){
          phenograph.res <- phenograph_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, 
                                             K = input$k, seed = input$set_seed, implementation = "Rphenograph")}
        if (input$ClustAlgo == "FastPG"){
         phenograph.res <- phenograph_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, 
                                              K = input$k, seed = input$set_seed, implementation = "FastPG")}
        if (inherits(x = phenograph.res,"list")){
          dI.clust$mapping <- phenograph.res[[1]]
          dI.clust$codes <- phenograph.res[[2]]
          dI.clust$modularity <- phenograph.res[[3]] #this is introduced with fastPG
          dI.clust$NumCluster <- nlevels(factor(dI.clust$mapping))
          color.clust <<- color_set[1:dI.clust$NumCluster]}else{
          dI.clust$mapping <- phenograph.res}}
        
        if (input$ClustAlgo == "DepecheR"){
          depeche.res <- DepecheR_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, 
                                               K = input$k, seed = input$set_seed)
          if (inherits(x = depeche.res,"list")){
            dI.clust$mapping <- depeche.res[[1]]
            dI.clust$codes <- depeche.res[[2]]
            dI.clust$NumCluster <- nlevels(factor(dI.clust$mapping))
            color.clust <<- color_set[1:dI.clust$NumCluster]}else{
              dI.clust$mapping <- depeche.res}}
        
        if(inherits(x = dI.clust$mapping,"integer")) {
          output$clust_panelStatus <- reactive({input$clust_panelStatus=="show"})
          outputOptions(output, "clust_panelStatus", suspendWhenHidden = FALSE)
        
        output$console_output_clust <- renderPrint("clustering process ends")
        print("clustering process ends")}
        
        if (input$clusteringPerf){
          print("start clusters evaluation process")
          cell_clustering <- dI.clust$mapping #[1:Number of events]
          
          flow_Frame <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_sample")
          
          if ((input$ClusteringPdistance)=="minkowski"){
            minkowski <- sum(meta.markers()$selected)}
          else{minkowski = NULL}
          
          clust_color <- color.clust
          plot_res <- cluster_eva(fF = flow_Frame, clustering = cell_clustering, color = clust_color, distance = input$ClusteringPdistance, 
                                  power = minkowski, setseed = input$set_seed, type = "clust")
          dI.clust$sil_plot <- plot_res[[1]] 
          dI.clust$sil_summary <- plot_res[[2]]
          
          if(inherits(x = dI.clust$sil_summary,"data.frame")) {
          output$clust_eva_panelStatus <- reactive({input$clust_eva_panelStatus=="show"})
          outputOptions(output, "clust_eva_panelStatus", suspendWhenHidden = FALSE)
          
          output$console_output_clust <- renderPrint("clusters evaluation process end")
          print("clusters evaluation process end")}}
      
      output$console_output_clust <- renderPrint("clustering process ends")
      if (Sys.info()["sysname"] == "Windows") play(x = bell)} 
      else {dI.clust$mapping <- "Please enter the meta-marker's related table in the 'Metadata & Assays' workflow"}
    } #if (inherits(x=fSmeta(),"flowSet"))
    else
    {dI.clust$mapping <- "You have to deal with a flowSet expreriment: please load your '.fcs' files in the pre-Clustering workflow or through the
    Metadata&Assays workflow"}
    
    output$checktxt_clust <- renderText({if (!(inherits(x = dI.clust$mapping,"integer"))) 
      {validate(need(expr =  (inherits(x = dI.clust$mapping,"integer")), message = dI.clust$mapping))}})
  })
  
  #################################################################### plotClust -----
  
  dI.clust_val <- reactiveValues(markerPlot = NULL, compPlot = NULL)
  observeEvent(eventExpr = input$plotClust , {
    
    clean_graph()
    validate(need(expr = (inherits(x = fSmeta(),"flowSet")), message = "It needs valid FCS files"))
 
    if(inherits(x = dI.clust$mapping,"integer")) {
      print("start visualization process for the clustering data")
      output$console_output_clust <- renderPrint("start visualization process for the clustering data - waiting...")
      
      output$plotClust_panelStatus <- reactive({input$plotClust_panelStatus=="show"})
      outputOptions(output, "plotClust_panelStatus", suspendWhenHidden = FALSE)

      if (input$ClustAlgo == "flowSOM"){
        
        flowSOM_plot(mst = dI.clust$mst , markers_df = meta.markers())
        path_to_file <- "./tmpdata"
        Zip_Files <- list.files(path = path_to_file, pattern = "^Plot", full.names = TRUE)
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) #this is necessary to sort in the correct way the flowFrames in the flowSet
        zip(zipfile = "./tmpdata/plotClust.zip", files = Zip_Files)
  
        output$markerPlot <- renderImage({
          outfile <- "./tmpdata/PlotMarker.png"
          list(src = outfile, contentType = 'image/png', alt = "markers plot in a MST graph")}, deleteFile = F)
        
        output$compPlot <- renderImage({
          outfile <- "./tmpdata/PlotComp.png"
          list(src = outfile, contentType = 'image/png', alt = "components of the MST graph")}, deleteFile = F)
        
        output$downloadClust <- downloadHandler(
          filename = function() {
            paste("flowSOMoutput", "zip", sep=".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/plotClust.zip", file)
          },
          contentType = "application/zip"
        )
      }
      
      if ((input$ClustAlgo == "Rphenograph")||(input$ClustAlgo == "FastPG")){
        eventi <- as.character(length(dI.clust$mapping))
        numCluster <- as.character(nlevels(factor(dI.clust$mapping)))
        para <- as.character(input$k)
        modularity <- as.character(dI.clust$modularity)
        if (input$ClustAlgo == "Rphenograph"){
          text_report <- paste0("<h2>RPhenograph (the R implementation of Phenograph clustering algorithm), run with K = ",
                                para , ", grouping the dataset of ",eventi , " events in ",
                                numCluster , " clusters, with ", modularity, " as modularity<h2>")
          output$phenograph_text <- renderText(text_report)
          text_report <- paste0("RPhenograph (the R implementation of Phenograph clustering algorithm), run with K = ",
                                para , ", grouping the dataset of ",eventi , " events in ",
                                numCluster , " clusters, with ", modularity, " as modularity")
          cat(text_report,file="./tmpdata/Rphenograph_report.txt",append=TRUE)}
        if (input$ClustAlgo == "FastPG"){
          text_report <- paste0("<h2>FastPG (an enhanced R implementation Phenograph clustering algorithm), run with K = ",
                                para , ", grouping the dataset of ",eventi , " events in ",
                                numCluster , " clusters, with ", modularity, " as modularity<h2>")
          output$phenograph_text <- renderText(text_report)
          text_report <- paste0("FastPG (an enhanced R implementation Phenograph clustering algorithm), run with K = ",
                                para , ", grouping the dataset of ",eventi , " events in ",
                                numCluster , " clusters, with ", modularity, " as modularity")
          cat(text_report,file="./tmpdata/FastPG_report.txt",append=TRUE)}}
      
      if (input$ClustAlgo == "KMeans"){
        eventi <- as.character(length(dI.clust$mapping))
        numCluster <- as.character(nlevels(factor(dI.clust$mapping)))
        para <- as.character(input$k)
        output$kmeans_text <- renderText(paste0("<h2>Kmeans clustering algorithm, run with K = ",para, 
                                                ", grouped the dataset of ",eventi, 
                                                "events in ", numCluster, " clusters<h2>"))}
      
      if (input$ClustAlgo == "DepecheR"){
        eventi <- as.character(length(dI.clust$mapping))
        numCluster <- as.character(nlevels(factor(dI.clust$mapping)))
        para <- as.character(input$k)
        output$depeche_text <- renderText(paste0("<h2>DEPECHE clustering algorithm, run with K = ",para, 
                                                ", grouped the dataset of ",eventi, 
                                                " events in ", numCluster, " clusters<h2>"))}
      
      if (input$clusteringPerf == TRUE){
        output$sil_plot <- renderPlot({dI.clust$sil_plot})
        output$sil_summary <- renderTable({dI.clust$sil_summary})
        
        output$downloadPerf_csv <- downloadHandler(filename <- function() {paste("silhouette_clust","csv", sep = ".")},
                                                   content  <- function(file) {file.copy("./tmpdata/silhouette_clust.csv", file)},
                                                   contentType = "text/csv" )}
        
      print("visualization process for the clustering data ends")
      output$console_output_clust <- renderPrint("visualization process for the clustering data ends")}
  })
  
  
  #################################################################### saveClust -----
  
  dI.clust_save <- reactiveValues(fFclust = NULL)
  
  observeEvent(eventExpr = input$saveClust , {

    validate(need(expr = (inherits(x = fSmeta(),"flowSet")), message = "It needs valid FCS files"))
    if((inherits(x = dI.clust$mapping,"integer"))&&(inherits(x = fSmeta(),"flowSet"))) {
      print("producing concatenated sample with the clustering data")
      output$console_output_clust <- renderPrint("producing concatenated sample with the clustering data - waiting...")  
      output$clust_save_panelStatus <- reactive({input$clust_save_panelStatus=="show"})
      outputOptions(output, "clust_save_panelStatus", suspendWhenHidden = FALSE)
      
      if (length(dI()) > 1){
        fFvectorName <- vector(mode = "character", length(dI.H()))
        for (i in seq_along(1:length(fSmeta()))){
          fFvectorName[[i]] <- basename(fSmeta()[[i]]@description$FILENAME)}
        save_fS(in_flow_Set = dI(), handled_flow_Set = fSmeta(), stringhe =fFvectorName)
        
        fFvectorName_comp <- paste0("./tmpdata/complete_",fFvectorName)
        fS.complete <- read.flowSet(files = fFvectorName_comp, truncate_max_range = FALSE)
        fF <- concatenating_fS(flow_Set = fS.complete, stringa = "conc_clust_sample")}
      else
        {fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_clust_sample")}
        
      dI.clust_save$fFclust <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)
      print("concatenated sample with clustering data produced")
      output$console_output_clust <- renderPrint("concatenated sample with clustering data produced")
      if (Sys.info()["sysname"] == "Windows") play(x = bell)}
    
    output$downloadFCSClust <- downloadHandler(
      
      filename = function() {
        paste("cluster", "zip", sep = ".")
      },
      content = function(file){
        file_path <- file
        tmpdir <- tempdir()
        setwd(tempdir())
        write.FCS(dI.clust_save$fFclust, filename = "cluster.fcs")
        Zip_Files <- list.files(path = getwd(), pattern = "^cluster*")
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE)
        zip(zipfile = file_path, files = Zip_Files)
        setwd(app_dir)
      },
      contentType = "application/zip")
  })
  
  
  ################################################# ___________Label I -------
  ################################################# labelClust -------
  
  dI.labelClust <- reactiveValues(label.res = NULL, cs.res = NULL, hm_cluster = NULL)

  observeEvent(eventExpr = input$labelClust , {
   
    clean_graph()
    if((inherits(x = dI.clust$codes,"matrix"))&&(inherits(x = fSmeta(),"flowSet"))&&(dI.flow$labelI)){
        if (input$signature_finding_method=="Densities"){
          options(warn = 1)
          dI.labelClust$cs.res <- try(expr = cluster_signature(fS = fSmeta(), selected.marker = mk$selected, clustId = dI.clust$mapping,
                                               color_marker = color.marker, tipo="cluster", central = input$HeatMapCentral, 
                                               min_quantile = input$minQuantile, max_quantile = input$maxQuantile,
                                               pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold), silent = FALSE)
          
          if (inherits(x = dI.labelClust$cs.res,"list")){options(warn = 0)} else errore <- dI.labelClust$cs.res}
      
        if((inherits(x = dI.clust$mapping,"integer"))&&(input$minQuantile<input$maxQuantile)&&(input$pos_threshold>=input$neg_threshold)) {
          print("start visualization process for the clustering heatmap")
          output$console_output_clust <- renderPrint("start visualization process for the clustering data - waiting...")
          output$labelClust_panelStatus <- reactive({input$labelClust_panelStatus=="show"})
          outputOptions(output, "labelClust_panelStatus", suspendWhenHidden = FALSE)
          
          print("start labelling I process")
          output$console_output_clust <- renderPrint("start labelling I process with the clustered metadata - waiting...")
          mapping <- dI.clust$mapping 
          if(nrow(ph)>0){
            
            if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
            label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = mapping, 
                                  phenoquery = ph, central = input$HeatMapCentral, 
                                  min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                  pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                  cluster.signature = dI.labelClust$cs.res[[2]])} else 
            {label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = mapping, 
                                   phenoquery = ph, central = input$HeatMapCentral, 
                                   min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                   pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold)}
            dI.labelClust$label.res <- label.res}
        
            output$console_output_clust <- renderPrint("start producing heatmap")
            print("start producing heatmap")
            output$plot_labelClust_panelStatus <- reactive({input$plot_labelClust_panelStatus=="show"})
            outputOptions(output, "plot_labelClust_panelStatus", suspendWhenHidden = FALSE)
            mapping <- dI.clust$mapping 
            #color_clusters <- color_set[(uni.tag1 + uni.tag2 + uni.tag3 + uni.marker + uni.pheno + 1):length(color_set)]
            color_clusters <- color.clust
          if (is.data.frame(dI.labelClust$label.res)){
            if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
              dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label", 
                                                   clustering = mapping, selected_markers = meta.markers()$selected, 
                                                   cluster_labelling = dI.labelClust$label.res, 
                                                   pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                   central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                   pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                   dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                   color_clusters = color_clusters, color_label = color.pheno,
                                                   cluster.signature = dI.labelClust$cs.res[[1]])}else
              {dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), 
                                                type_HM = "full_label", clustering = mapping, selected_markers = meta.markers()$selected, 
                                                cluster_labelling = dI.labelClust$label.res, 
                                                pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                color_clusters = color_clusters, color_label = color.pheno)}}
            else{
              if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
                dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label", 
                                                            clustering = mapping, selected_markers = meta.markers()$selected, 
                                                            pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                            central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                            pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                            dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                            color_clusters = color_clusters, color_label = color.pheno, 
                                                            cluster.signature = dI.labelClust$cs.res[[1]])}else
                {dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label", 
                                                            clustering = mapping, selected_markers = meta.markers()$selected, 
                                                            pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                            central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                            pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                            dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                            color_clusters = color_clusters, color_label = color.pheno)}}
          dI.flow$labelI <- FALSE
          output$console_output_clust <- renderPrint("labelling process I ends")
          print("labelling process I ends")}
        if (!(inherits(x = dI.clust$mapping,"integer"))){
          dI.labelClust$label.res <- "Please perform the clustering algorithm: without clustering you cannot perfom this procedure"}
        if (!(input$minQuantile < input$maxQuantile)){
          dI.labelClust$label.res <- "Please fix your quantile related setting: it must be minQuantile < input$maxQuantile"}
        if (input$pos_threshold < input$neg_threshold){
          dI.labelClust$label.res <- "Please fix your labelling threshold setting: it must be pos_threshold > neg_threshold"} 
        if (!(inherits(x = dI.clust$mapping,"integer"))){
          dI.labelClust$label.res <- "Please perform the clustering algorithm: without clustering you cannot perfom this procedure"}
        
      #else{dI.labelClust$label.res <- paste0("Please enter a phenotype table in the 'Metadata & Assays' workflow or skip this and go directly to the ", 
      #                                    "'Meta-clustering' step")}
      dI.flow$labelI <- FALSE
      if (Sys.info()["sysname"] == "Windows") play(x = bell)}
    
      output$checktxt_labelClust <- renderText({
        if (!(inherits(x = dI.clust$codes,"matrix")))
          {validate(need(expr =  (inherits(x = dI.clust$codes,"matrix")), message = "Please, perform clustering first"))}
        if (!(inherits(x = dI.clust$codes,"matrix")))
          {validate(need(expr =  (inherits(x = fSmeta(),"matrix")), message = "You have not load any valid '.fcs' file"))}
        if (!(is.data.frame(dI.labelClust$label.res)))
          {validate(need(expr =  (is.data.frame(dI.labelClust$label.res)), message = dI.labelClust$label.res))}
        if (!is.list(dI.labelClust$hm_cluster))
          {validate(need(expr =  (is.list(dI.labelClust$hm_cluster)), message = dI.labelClust$hm_cluster))}})
  })
  
  ################################################# plot_labelClust -------
  
  observeEvent(eventExpr = input$plot_labelClust , {
    
    #if ((inherits(x=dI.clust$codes,"matrix"))&&(is.data.frame(dI.labelClust$label.res))){
    if (inherits(x = dI.clust$codes,"matrix")){
     
      clean_graph()
      output$console_output_clust <- renderPrint("start producing heatmap")
      print("start producing heatmap")
      output$plot_labelClust_panelStatus <- reactive({input$plot_labelClust_panelStatus=="show"})
      outputOptions(output, "plot_labelClust_panelStatus", suspendWhenHidden = FALSE)
      
      mapping <- dI.clust$mapping 
      color_clusters <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
 
      if(inherits(x = dI.labelClust$hm_cluster,'list')){
        path_to_file <- "./tmpdata"
        Zip_Files <- list.files(path = path_to_file, pattern = "^HM_", full.names = TRUE)
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) #this is necessary to sort in the correct way the flowFrames in the flowSet
        zip(zipfile = "./tmpdata/plotClust.zip", files = Zip_Files)
        
        cdata <- session$clientData
        
        output$HM_full_label <- renderImage({
          outfile <- "./tmpdata/HM_full_label.png"
          list(src = outfile, contentType = 'image/png', alt = "cluster labeld heatmap",
               width = cdata$output_HM_full_label_width, height = cdata$output_HM_full_label_height)}, deleteFile = F)
        
        output$console_output_clust <- renderPrint("heatmap process end")
        print("visualization process for the clustering heatmap ends")
        
        output$downloadClust_p <- downloadHandler(
          filename = function() {
            paste("output", "zip", sep=".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/plotClust.zip", file)
          },
          contentType = "application/zip"
        )
        
        output$downloadClust_csv <- downloadHandler(
          filename = function() {
            paste("HM_full_label","csv", sep = ".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/HM_full_label.csv", file)
          },
          contentType = "text/csv"
        )
        
        output$downloadClustPheno_csv <- downloadHandler(
          filename = function() {
            paste("pheno_full_label","csv", sep = ".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/pheno_full_label.csv", file)
          },
          contentType = "text/csv"
        )
      } 
    }
    output$checktxt_plot_labelClust <- renderText({
      if (!(inherits(x = dI.clust$codes,"matrix")))
        {validate(need(expr =  (inherits(x = dI.clust$codes,"matrix")), message = "Please, perform clustering first"))}})
  })
  
  
  ################################################### _____Meta-Clustering -----    
  ################################################### runMetaClust -----    
  
  dI.metaClust <- reactiveVal(NULL)
  dI.metaClust_val <- reactiveValues(hm_cluster = NULL, sil_plot = NULL, sil_summary = NULL, cs.res = NULL)
  
  observeEvent(eventExpr = input$runMetaClust, handlerExpr = {

    if((inherits(x = dI.clust$codes,"matrix"))&&(inherits(x = dI.clust$mapping,"integer"))&&(dI.flow$meta_clust)){
      print("start meta-clustering process")
      output$console_output_clust <- renderPrint("start meta-clustering process - waiting...")
      
      dI.metaClust(NULL) #this is to reset the metaclustering results
      
      if ((input$MetaClustAlgo == "ConsensusClusterPlus")&&(dI.clust$NumCluster>input$MetaClustNum)){
         mc.res <- meta_clustering(codes = dI.clust$codes, maxK = input$MetaClustNum, reps = 100, 
                        pItem = 1.0, clusterAlg = input$MetaClustCriteria, distance = input$MetaClustDist, 
                        seed = input$set_seed, savedir = "./tmpdata", save = TRUE)
          dI.metaClust(mc.res) 
        #Notice, depending of the dataset, the meta_clustering can produce a number of meta_clusters which is less than input$MetaClustNum 
        
        if(inherits(x = dI.metaClust(),"list")) {
          output$metaClust_panelStatus <- reactive({input$metaClust_panelStatus=="show"})
          outputOptions(output, "metaClust_panelStatus", suspendWhenHidden = FALSE)
          
          code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
          cell_clustering <- code_clustering[dI.clust$mapping] #[1:Number of events]
          
          if (input$signature_finding_method=="Densities"){
            dI.metaClust_val$cs.res <- try(expr = cluster_signature(fS = fSmeta(), selected.marker = mk$selected, clustId = cell_clustering,
                                                        color_marker = color.marker, tipo ="meta", central = input$HeatMapCentral,
                                                        min_quantile = input$minQuantile, max_quantile = input$maxQuantile,
                                                        pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold), silent = FALSE)
          if (inherits(x = dI.metaClust_val$cs.res,"list")){options(warn = 0)} else errore <- dI.metaClust_val$cs.res}
        
          if (is.data.frame(dI.labelClust$label.res)){
            if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
            dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering, 
                                                         cluster_labelling = NULL, pheno_table = ph, Nclust = input$MetaClustNum, 
                                                         selected_markers = meta.markers()$selected, 
                                                         central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                                         max_quantile = input$maxQuantile, 
                                                         pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                         dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                         color_clusters = color.metaclust, cluster.signature = dI.metaClust_val$cs.res[[1]])}
            else{
                 dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering, 
                                                              cluster_labelling = NULL, pheno_table = ph, Nclust = input$MetaClustNum, 
                                                              selected_markers = meta.markers()$selected, 
                                                              central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                                              max_quantile = input$maxQuantile, 
                                                              pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                              dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                              color_clusters = color.metaclust)}
            } 
          else 
          {if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
            dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering,
                                                                    cluster_labelling = NULL, pheno_table = NULL, Nclust = input$MetaClustNum, 
                                                                    selected_markers = meta.markers()$selected,
                                                                    central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                                                    max_quantile = input$maxQuantile,
                                                                    pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                                    dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                                    color_clusters = color.metaclust, cluster.signature = dI.metaClust_val$cs.res[[1]])}
            else{
              dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering,
                                                                        cluster_labelling = NULL, pheno_table = NULL, Nclust = input$MetaClustNum, 
                                                                        selected_markers = meta.markers()$selected,
                                                                        central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                                                        max_quantile = input$maxQuantile,
                                                                        pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                                        dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                                        color_clusters = color.metaclust)}}
            if (input$clusteringPerf == TRUE){
              print("start meta-clusters evaluation process")
            
            flow_Frame <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_metasample")
            
            if ((input$ClusteringPdistance)=="minkowski"){
              minkowski <- sum(meta.markers()$selected)}
            else{minkowski = NULL}
            
            clust_color <- color.metaclust
            plot_res <- cluster_eva(fF = flow_Frame, clustering = cell_clustering, color = color.metaclust, distance = input$ClusteringPdistance, 
                                    power = minkowski, setseed = input$set_seed, type = "meta")
            dI.metaClust_val$sil_plot <- plot_res[[1]] 
            dI.metaClust_val$sil_summary <- plot_res[[2]]

          if(inherits(x = dI.metaClust_val$sil_summary,"data.frame")) {
            output$metaClust_eva_panelStatus <- reactive({input$metaClust_eva_panelStatus=="show"})
            outputOptions(output, "metaClust_eva_panelStatus", suspendWhenHidden = FALSE)
            
            output$console_output_clust <- renderPrint("meta-clusters evaluation process end")
            print("clusters evaluation process end")} #this is new
            }
          }}
          
      print("meta-clustering process ends") 
      output$console_output_clust <- renderPrint("meta-clustering process ends") 
      if (Sys.info()["sysname"] == "Windows") play(x = bell)
      dI.flow$meta_clust <- FALSE}
          else {
        if (!(inherits(x = fSmeta(),"flowSet")))
          {(mc.res <- "You have not loaded any '.fcs' files")}
        if(!(inherits(x = dI.clust$mapping,"integer")))
          {(mc.res <- "Please perform the clustering before the meta-clustering!")}
        if(!is.null(dI.clust$NumCluster)){
          if(!(dI.clust$NumCluster>input$MetaClustNum))
          {(dI.metaClust_val$hm_cluster <- 
                "The number of meta-clusters must be lower than the numbers of clusters: Did you check the related paramter's setting?")}}
      }
          
    
    output$checktxt_metaClust <- renderText({
      if (!is.list(dI.metaClust_val$hm_cluster))
        {validate(need(expr =  (is.list(dI.metaClust_val$hm_cluster)), message = dI.metaClust_val$hm_cluster))}})
    
    output$checktxt_metaClust <- renderText({
        if (!(inherits(x = fSmeta(),"flowSet"))){
          validate(need(expr =  (inherits(x = fSmeta(),"flowSet")), message = mc.res))}
      if(!is.null(dI.clust$NumCluster)){
        if (!(dI.clust$NumCluster>input$MetaClustNum)){
          validate(need(expr =  (dI.clust$NumCluster>input$MetaClustNum), 
                        message = "The number of meta-clusters must be less than the number of clusters"))}}
        if (!(inherits(x = dI.metaClust(),"list"))){
              validate(need(expr =  (inherits(x = dI.metaClust(),"list")), message = "Something went wrong with the meta-clustering. Try again..."))}
        if (!(inherits(x = dI.clust$mapping,"integer"))){
          validate(need(expr =  (inherits(x = dI.clust$mapping,"integer")), message = mc.res))}
      })
  })
  
  ################################################### plotMeta -----    
  
  observeEvent(eventExpr = input$plotMeta , {
    
    clean_graph()
    if(inherits(x = dI.metaClust(),"list")) {
      print("start visualization process for the meta-clustering data")
      output$console_output_clust <- renderPrint("start visualization process for the meta-clustering data - waiting...")
      
      if(is.list(dI.metaClust_val$hm_cluster)){
        output$plotMeta_panelStatus <- reactive({input$plotMeta_panelStatus=="show"})
        outputOptions(output, "plotMeta_panelStatus", suspendWhenHidden = FALSE)
        
        path_to_file <- "./tmpdata"
        Zip_Files <- list.files(path = path_to_file, pattern = "^HM_", full.names = TRUE)
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE)
        zip(zipfile = "./tmpdata/plotMeta.zip", files = Zip_Files)
        
        Zip_Files <- list.files(path = path_to_file, pattern = "^tmpdata.", full.names = TRUE)
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
        zip(zipfile = "./tmpdata/consensus_csv.zip", files = Zip_Files)
        
        Zip_Files <- list.files(path = path_to_file, pattern = "^consensus.", full.names = TRUE)
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
        zip(zipfile = "./tmpdata/consensus_png.zip", files = Zip_Files)
        
        output$metaHeatmap_cluster <- renderImage({
          outfile <- "./tmpdata/HM_meta.png"
          list(src = outfile, contentType = 'image/png', alt = "metacluster heatmap")}, deleteFile = F)
        #??? show the traking plots of the ConsensusClusterPlus
        
        print("visualization process for the meta-clustering data ends")
        output$console_output_clust <- renderPrint("visualization process for the meta-clustering data ends")
        
        output$downloadMeta <- downloadHandler(
          filename = function() {
            paste("output", "zip", sep=".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/plotMeta.zip", file)
          },
          contentType = "application/zip"
        )
        
        output$downloadMeta_csv <- downloadHandler(
          filename = function() {
            paste("HM_meta","csv", sep = ".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/HM_meta.csv", file)
          },
          contentType = "text/csv"
        )
        
        output$downloadPhenoMeta_csv <- downloadHandler(
          filename = function() {
            paste("pheno_meta","csv", sep = ".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/pheno_meta.csv", file)
          },
          contentType = "text/csv"
        )
        
        
        output$downloadConsensus_csv <- downloadHandler(
          filename = function() {
            paste("consensus_csv__output", "zip", sep=".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/consensus_csv.zip", file)
          },
          contentType = "application/zip"
        )
        
        output$downloadConsensus_png <- downloadHandler(
          filename = function() {
            paste("consensus_png_output", "zip", sep=".")
          },
          
          content <- function(file) {
            file.copy("./tmpdata/consensus_png.zip", file)
          },
          contentType = "application/zip"
        )
        
        if (input$clusteringPerf == TRUE){
          output$meta_sil_plot <- renderPlotly({dI.metaClust_val$sil_plot})
          output$meta_sil_summary <- renderTable({dI.metaClust_val$sil_summary})
          
          output$downloadMetaPerf_csv <- downloadHandler(filename <- function() {paste("meta_silhouette","csv", sep = ".")},
                                                     content  <- function(file) {file.copy("./tmpdata/meta_silhouette.csv", file)},
                                                     contentType = "text/csv" )}}}
    
    output$checktxt_plotMeta <- renderText({
      if (!is.list(dI.metaClust_val$hm_cluster))
        {validate(need(expr =  (is.list(dI.metaClust_val$hm_cluster)), message = dI.metaClust_val$hm_cluster))}})
    
  })
  
  
  #################################################################### saveMetaClust  ----
  
  dI.metaclust_save <- reactiveValues(fFclust = NULL, fFclust_ori = NULL)
  
  observeEvent(eventExpr = input$saveMetaClust , {
    
    validate(need(expr = (inherits(x = fSmeta(),"flowSet")), message = "It needs valid FCS files"))
    if((inherits(x = dI.metaClust(),"list"))&&(inherits(x = fSmeta(),"flowSet"))&&(inherits(x = dI.clust$mapping,"integer"))){
      print("producing concatenated sample with the meta-clustering data")
      output$console_output_clust <- renderPrint("producing concatenated sample with the meta-clustering data - waiting...")
      output$saveMetaClust_panelStatus <- reactive({input$saveMetaClust_panelStatus=="show"})
      outputOptions(output, "saveMetaClust_panelStatus", suspendWhenHidden = FALSE)
      
      code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
      cell_clustering <- code_clustering[dI.clust$mapping] #[1:Number of events]
      if (length(dI()) > 1){#this because user may enter from the second workflow and thus without any dI()
          
        fFvectorName <- vector(mode = "character", length(fSmeta()))
        for (i in seq_along(1:length(fSmeta()))){
          fFvectorName[[i]] <- basename(fSmeta()[[i]]@description$FILENAME)}
            
        save_fS(in_flow_Set = dI(), handled_flow_Set = fSmeta(), stringhe = fFvectorName)
        fFvectorName_comp <- paste0("./tmpdata/complete_",fFvectorName)
        fS.complete <- read.flowSet(files = fFvectorName_comp, truncate_max_range = FALSE)
        fF <- concatenating_fS(flow_Set = fS.complete, stringa = "conc_meta_sample")}
      else
        {fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_meta_sample")}
  
      dI.metaclust_save$fFclust <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)
      dI.metaclust_save$fFclust <- add_dim(flow_Frame = dI.metaclust_save$fFclust, dim_name = "meta_clusterId", dim_vect = cell_clustering)
      print("concatenated sample with meta-clustering data produced")
      output$console_output_clust <- renderPrint("concatenated sample with meta-clustering data produced")
      if (Sys.info()["sysname"] == "Windows") play(x = bell)
      
      output$downloadFCSMeta_Clust <- downloadHandler(
        
        filename = function() {
          paste("meta_cluster", "zip", sep = ".")
        },
        content = function(file){
          
          file_path <- file
          tmpdir <- tempdir()
          setwd(tempdir())
          write.FCS(dI.metaclust_save$fFclust, filename = "meta_cluster.fcs")
          Zip_Files <- list.files(path = getwd(), pattern = "^meta_cluster*")
          Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
          zip(zipfile = file_path, files = Zip_Files)
          setwd(app_dir)
        },
        contentType = "application/zip")  
    }
    
    output$checktxt_saveMetaClust <- renderText({
      if (!(inherits(x = fSmeta(),"flowSet")))
        {validate(need(expr =  (inherits(x = fSmeta(),"flowSet")), message = "you have not load any '.fcs' file"))}
      if (!(inherits(x = dI.metaClust(),"list")))
        {validate(need(expr =  (inherits(x = dI.metaClust(),"list")), message = "please perform metaclust before"))}
      if (!(inherits(x = dI.clust$mapping,"integer")))
        {validate(need(expr =  (inherits(x = dI.clust$mapping,"integer")), message = "you should have performed at least clustering first"))}})
    
  })
  
  ################################################### __________Mapping -----    
  ################################################### runMapsComp -----    
  
  dI.map <- reactiveVal(NULL)
  
  observeEvent(eventExpr = input$runMapsComp, handlerExpr = {
    
    #if(((is.data.frame(dI.labelClust$label.res))||(class(dI.metaClust())=="list"))&&(class(fSmeta())=="flowSet"))
    if(((is.list(dI.labelClust$hm_cluster))||(inherits(x = dI.metaClust(),"list")))&&(inherits(x = fSmeta(),"flowSet")))
    {
      print("start mapping process")
      output$console_output_clust <- renderPrint("start mapping process - waiting...")
      if ((inherits(x = dI.tSNE.eva.final$mapres_final,"list"))&&(!input$two_three)&&(input$mapType=='tSNE'))
        {mapres <- dI.tSNE.eva.final$mapres_final}
      else
        {mapres <- map_gen(flow_Set = fSmeta(), type = input$mapType, two_three = input$two_three, 
                        selected_markers = mk$selected, seed_nr = input$set_seed, limit = input$map_limit, maxCell = input$mapMax,  
                        #tSNE parameters
                        Rtsne_pca = input$tSNE_PCA, Rtsne_perp = input$perplexity, Rtsne_theta = input$theta,
                        Rtsne_iter = input$tSNEiter, Rtsne_eta = input$tSNEeta,
                        #UMAP paramters
                        UMAP.n_neighbors = input$UMAP.n_neighbors, UMAP.metric = input$UMAP.metric, 
                        UMAP.min_dist = input$UMAP.min_dist, UMAP.spread = input$UMAP.spread, UMAP.random_state = input$set_seed)}
      dI.map(mapres)
      
      if(length(dI.map())==3){
        print("mapping process ends")
        output$console_output_clust <- renderPrint("mapping process ends")
        output$runMapsComp_panelStatus <- reactive({input$runMapsComp_panelStatus=="show"})
        outputOptions(output, "runMapsComp_panelStatus", suspendWhenHidden = FALSE)
        if (Sys.info()["sysname"] == "Windows") play(x = gong)}
     
      }
    else{ 
      if (!(inherits(x = fSmeta(),"flowSet")))
        {dI.map("You have not load any '.fcs' file")}
      if (!(is.data.frame(dI.labelClust$label.res))||(inherits(x = dI.metaClust(),"list")))
        {dI.map("You should perform the labelling action and/or the meta-clustering to produce a colored map")}}
    
    output$checktxt_runMapsComp <- renderText({
      if (inherits(x = dI.map(),"character"))
        {validate(need(expr =  !(inherits(x = dI.map(),"character")), message = dI.map()))}})
  })
  
  ################################################### showmapcomp -----    
  dI.map_val <- reactiveValues(showmapcomp = NULL, jt = NULL, jt_1 = NULL, jt_2 = NULL, jt_3 = NULL, jt_4 = NULL, 
                               jt_meta = NULL, jt_meta_1 = NULL, jt_meta_2 = NULL, jt_meta_3 = NULL, jt_meta_4 = NULL)
  
  observeEvent(eventExpr = input$plotmapcomp , {
    
    clean_graph()
    
    if(inherits(x = dI.map(),"list")) {
      print("start map visualization")
      output$console_output_clust <- renderPrint("start map visualization process - waiting...")
    
      if(inherits(x = dI.metaClust(),"list")){ #case 1
        output$condition_meta <- reactive({input$condition_meta=="show"})
        outputOptions(output, "condition_meta", suspendWhenHidden = FALSE)
        
        code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
        meta_clust <- code_clustering
        cell_clustering <- as.integer(code_clustering[dI.clust$mapping])
        metaClustNum <- input$MetaClustNum
        color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
        showmap <- map_plot_comp(map_df = dI.map()[[1]], type = input$mapType, two_three = input$two_three, map_inds = dI.map()[[2]],
                                 map_reduced = dI.map()[[3]], metadata = mt, clust = dI.clust$mapping, metaclust = meta_clust,
                                 NMC = metaClustNum, color_clusters = color_metaclust)
        fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_meta_sample")
        nomi <- colnames(exprs(fF))
        if (!any(grepl("clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)}
        if (!any(grepl("meta_clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "meta_clusterId", dim_vect = cell_clustering)}
        
        jt <- join_tag(frame = fF, meta = mt, type = "meta", tag1col=color.tag1, tag2col=color.tag2, tag3col=color.tag3, tag4col=color.tag4)} 
      else{
        #if (is.data.frame(dI.labelClust$label.res))
        if (is.list(dI.labelClust$hm_cluster))
        {meta_clust <- dI.labelClust$label.res$numeric
        metaClustNum <- length(dI.labelClust$label.res$numeric)
        
        color_metaclust <- dI.labelClust$hm_cluster$color_clusters
        showmap <- map_plot_comp(map_df = dI.map()[[1]], type = input$mapType, two_three = input$two_three, map_inds = dI.map()[[2]],
                                 map_reduced = dI.map()[[3]], metadata = mt, clust = dI.clust$mapping, 
                                 NMC = metaClustNum, color_clusters = color_metaclust)
        fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_meta_sample")
        nomi <- colnames(exprs(fF))
        if (!any(grepl("clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)}}
        jt <- join_tag(frame = fF, meta = mt, type = "clust", tag1col=color.tag1, tag2col=color.tag2, tag3col=color.tag3, tag4col=color.tag4)}
      
      output$runMapShow_panelStatus <- reactive({input$runMapShow_panelStatus=="show"})
      outputOptions(output, "runMapShow_panelStatus", suspendWhenHidden = FALSE)
      
      dI.map_val$showmapcomp <- showmap
      output$showmapcomp <- renderPlotly({dI.map_val$showmapcomp})
      
      nomi_tag <- colnames(mt)
      nomi_tag <- nomi_tag[-(1:2)]
      if (inherits(x = jt,"list")){
        dI.map_val$jt <- jt[[1]][[1]]
        output$showJtSample <- renderPlot({dI.map_val$jt})
        if(inherits(x = dI.metaClust(),"list")){dI.map_val$jt_meta <- jt[[2]][[1]]; output$showJtMetaSample <- renderPlot({dI.map_val$jt_meta})}
        if (uni.tag1 > 1){
          output$conditionJt_tag1 <- reactive({input$conditionJt_tag1=="show"})
          outputOptions(output, "conditionJt_tag1", suspendWhenHidden = FALSE)
          dI.map_val$jt1 <- jt[[1]][[2]]
          output$showJt1_Sample <- renderPlot({dI.map_val$jt1})
          if(inherits(x = dI.metaClust(),"list")){dI.map_val$jt_meta_1 <- jt[[2]][[2]]; output$showJt1_MetaSample <- renderPlot({dI.map_val$jt_meta_1})}}
        if (uni.tag2 > 1){
          output$conditionJt_tag2 <- reactive({input$conditionJt_tag2=="show"})
          outputOptions(output, "conditionJt_tag2", suspendWhenHidden = FALSE)
          dI.map_val$jt2 <- jt[[1]][[3]]
          output$showJt2_Sample <- renderPlot({dI.map_val$jt2})
          if(inherits(x = dI.metaClust(),"list")){dI.map_val$jt_meta_2 <- jt[[2]][[3]]; output$showJt2_MetaSample <- renderPlot({dI.map_val$jt_meta_2})}}
        if (uni.tag3 > 1){
          output$conditionJt_tag3 <- reactive({input$conditionJt_tag3=="show"})
          outputOptions(output, "conditionJt_tag3", suspendWhenHidden = FALSE)
          dI.map_val$jt3 <- jt[[1]][[4]]
          output$showJt3_Sample <- renderPlot({dI.map_val$jt3})
          if(inherits(x = dI.metaClust(),"list")){dI.map_val$jt_meta_3 <- jt[[2]][[4]]; output$showJt3_MetaSample <- renderPlot({dI.map_val$jt_meta_3})}}
        if (uni.tag4 > 1){
          output$conditionJt_tag4 <- reactive({input$conditionJt_tag4=="show"})
          outputOptions(output, "conditionJt_tag4", suspendWhenHidden = FALSE)
          dI.map_val$jt4 <- jt[[1]][[5]]
          output$showJt4_Sample <- renderPlot({dI.map_val$jt4})
          if(inherits(x = dI.metaClust(),"list")){dI.map_val$jt_meta_4 <- jt[[2]][[5]]; output$showJt4_MetaSample <- renderPlot({dI.map_val$jt_meta_4})}}
        
        print("quantity visualization process ends")
        output$console_output_clust <- renderPrint("quantity plot visualization process ends")
        
        output$downloadJtSample_csv <- downloadHandler(filename = function() {paste0("clusterId_", "sample",".csv")},
                                                          content <- function(file) {file.copy("clusterId_sample.csv", file)}, contentType = "text/csv")
        if(inherits(x = dI.metaClust(),"list")){
          output$downloadJtMetaSample_csv <- downloadHandler(filename = function() {paste0("meta_clasterId_", "sample",".csv")},
                                                         content <- function(file) {file.copy("meta_clusterId_sample.csv", file)}, contentType = "text/csv")}
        
        if (uni.tag1 > 1){
          output$downloadSampleTag1_csv <- downloadHandler(filename = function() {paste0("clusterId_", nomi_tag[1],".csv")},
                                                          content <- function(file) {file.copy(paste0("./tmpdata/clusterId_",nomi_tag[1],".csv"), file)}, contentType = "text/csv")
          if(inherits(x = dI.metaClust(),"list")){
            output$downloadSampleTag1_csv <- downloadHandler(filename = function() {paste0("meta_clusterId_", nomi_tag[1],".csv")},
                                                             content <- function(file) {file.copy(paste0("./tmpdata/meta_clusterId_",nomi_tag[1],".csv"), file)}, contentType = "text/csv")}}
        
        if (uni.tag2 > 1){
          output$downloadSampleTag2_csv <- downloadHandler(filename = function() {paste0("clusterId_", nomi_tag[2],".csv")},
                                                          content <- function(file) {file.copy(paste0("./tmpdata/clusterId_",nomi_tag[1],".csv"), file)}, contentType = "text/csv")
          if(inherits(x = dI.metaClust(),"list")){
            output$downloadSampleTag1_csv <- downloadHandler(filename = function() {paste0("meta_clusterId_", nomi_tag[2],".csv")},
                                                             content <- function(file) {file.copy(paste0("./tmpdata/meta_clusterId_",nomi_tag[2],".csv"), file)}, contentType = "text/csv")}}
        if (uni.tag3 > 1){
          output$downloadSampleTag3_csv <- downloadHandler(filename = function() {paste0("clusterId_", nomi_tag[3],".csv")},
                                                          content <- function(file) {file.copy(paste0("./tmpdata/clusterId_",nomi_tag[1],".csv"), file)}, contentType = "text/csv")
          if(inherits(x = dI.metaClust(),"list")){
            output$downloadSampleTag1_csv <- downloadHandler(filename = function() {paste0("meta_clusterId_", nomi_tag[3],".csv")},
                                                             content <- function(file) {file.copy(paste0("./tmpdata/meta_clusterId_",nomi_tag[3],".csv"), file)}, contentType = "text/csv")}}
        if (uni.tag4 > 1){
          output$downloadSampleTag4_csv <- downloadHandler(filename = function() {paste0("clusterId_", nomi_tag[4],".csv")},
                                                          content <- function(file) {file.copy(paste0("./tmpdata/clusterId_",nomi_tag[1],".csv"), file)}, contentType = "text/csv")}
          if(inherits(x = dI.metaClust(),"list")){
            output$downloadSampleTag1_csv <- downloadHandler(filename = function() {paste0("meta_clusterId_", nomi_tag[4],".csv")},
                                                           content <- function(file) {file.copy(paste0("./tmpdata/meta_clusterId_",nomi_tag[4],".csv"), file)}, contentType = "text/csv")}}
      
      output$console_output_clust <- renderPrint("map visualization ends")
      
      print("map visualization process ends")
      output$console_output_clust <- renderPrint("map visualization process ends")
    }
  })
  
  ################################################### show marker map -----    
  dI.map_marker <- reactiveValues(showmarker = NULL)
  observeEvent(eventExpr = input$plot_marker_map , {
  
    clean_graph()
    if(inherits(x = dI.map(),"list")) {
      print("start marker map visualization")
      output$console_output_clust <- renderPrint("start marker map visualization process - waiting...")
      showmap <- map_plot_comp(map_df = dI.map()[[1]], type = "marker")
      dI.map_marker$showmarker <- showmap}
    
    output$runMapMarker_panelStatus <- reactive({input$runMapMarker_panelStatus=="show"})
    outputOptions(output, "runMapMarker_panelStatus", suspendWhenHidden = FALSE)
    output$showmarkermap <- renderPlot({dI.map_marker$showmarker}, 
                                       width = 1350, height = function() {length(mk$selected[mk$selected==TRUE]) * 200})
    print("map visualization process ends")
    output$console_output_clust <- renderPrint("map visualization process ends")
  })
  
  #################################################################### saveMap  ----
  dI.map_save <- reactiveValues(fFmap = NULL)
  
  observeEvent(eventExpr = input$saveMap , {
  
    if((inherits(x = fSmeta(),"flowSet"))&&(inherits(x = dI.clust$mapping,"integer"))&&(inherits(x = dI.map(),"list"))){
      print("producing concatenated sample with the map data")
      output$console_output_clust <- renderPrint("producing concatenated sample with the map data - waiting...")
      output$saveMap_panelStatus <- reactive({input$saveMap_panelStatus=="show"})
      outputOptions(output, "saveMap_panelStatus", suspendWhenHidden = FALSE)
      
      if(inherits(x = fSmeta(),"flowSet")){
        
        if (length(dI()) > 0){
          fFvectorName <- vector(mode = "character", length(fSmeta()))
          for (i in seq_along(1:length(fSmeta()))){
            fFvectorName[[i]] <- basename(fSmeta()[[i]]@description$FILENAME)}
          
          save_fS(in_flow_Set = dI(), handled_flow_Set = fSmeta(), stringhe = fFvectorName)
          
          fFvectorName_comp <- paste0("./tmpdata/complete_",fFvectorName)
          fS.complete <- read.flowSet(files = fFvectorName_comp, truncate_max_range = FALSE)
          fF <- concatenating_fS(flow_Set = fS.complete, stringa = "conc_meta_sample")}
          else
            {fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_meta_sample")}
        
        #dI.map_save$fFmap <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)
        nomi <- colnames(exprs(fF))
        if (!any(grepl("clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)}
        
        if (inherits(x = dI.metaClust(),"list")){
          code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass 
          cell_clustering <- as.integer(code_clustering[dI.clust$mapping])
          if (!any(grepl("meta_clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "meta_clusterId", dim_vect = cell_clustering)}}
        
        meta <- mt
        nomi_tag <- colnames(meta)
        nomi_tag <- nomi_tag[-(1:2)]
        for (i in 1:4){meta[,i+2] <- (as.numeric(factor(meta[,i+2])))}
          if ((inherits(x = dI.metaClust(),"list"))&&(!any(grepl("meta_clusterId", nomi, fixed = TRUE)))){
            df <- as.data.frame(exprs(fF)[,c("SampleID", "clusterId", "meta_clusterId")])} else {df <- as.data.frame(exprs(fF)[,c("SampleID", "clusterId")])}
          
          meta$SampleID <- as.numeric(1:nrow(meta))
          rj <- right_join(meta,df)
          rj <- rj[,-(1:2)]
          rj$date <- NULL
          
          rj$SampleID <- factor(rj$SampleID)
          rj$clusterId <- factor(rj$clusterId)
          if ((inherits(x = dI.metaClust(),"list"))&&(!any(grepl("meta_clusterId", nomi, fixed = TRUE)))){rj$meta_clusterId <- factor(rj$meta_clusterId)}
      
          if (uni.tag1 > 1){
            tag <-  nomi_tag[1]
            tagdata <- as.integer(factor(rj[,1]))
            fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
          if (uni.tag2 > 2){
            tag <-  nomi_tag[2]
            tagdata <- as.integer(factor(rj[,2]))
            fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
          if (uni.tag3 > 3){
            tag <-  nomi_tag[3]
            tagdata <- as.integer(factor(rj[,3]))
            fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
          if (uni.tag4 > 1){
            tag <-  nomi_tag[4]
            tagdata <- as.integer(factor(rj[,4]))
            fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
          
          fF <- add_dim(flow_Frame = fF, dim_name = "map_X", dim_vect = dI.map()[[1]]$map_x)
          fF <- add_dim(flow_Frame = fF, dim_name = "map_Y", dim_vect = dI.map()[[1]]$map_y)
          if(input$two_three==TRUE)
          {fF <- add_dim(flow_Frame = fF, dim_name = "map_Z", dim_vect = dI.map()[[1]]$map_z)}}
      
      dI.map_save$fFmap <- fF
      print("concatenated sample with the map data produced")
      output$console_output_clust <- renderPrint("concatenated sample with map data produced")
      if (Sys.info()["sysname"] == "Windows") play(x = bell)
      
      output$downloadFCSMap <- downloadHandler(
        
        filename = function() {
          paste("map_cluster", "zip", sep = ".")
        },
        content = function(file){
          
          file_path <- file
          tmpdir <- tempdir()
          setwd(tempdir())
          write.FCS(dI.map_save$fFmap, filename = "map_cluster.fcs")
          write.csv(x =  exprs(dI.map_save$fFmap), file = "map_cluster.csv")
          Zip_Files <- list.files(path = getwd(), pattern = "^map_cluster*")
          Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
          zip(zipfile = file_path, files = Zip_Files)
          setwd(app_dir)
        },
        contentType = "application/zip")}
    
    output$checktxt_saveMap <- renderText({
      if (!(inherits(x = fSmeta(),"flowSet")))
        {validate(need(expr =  (inherits(x = fSmeta(),"flowSet")), message = "you have not load any '.fcs' file"))}
      if (!(inherits(x = dI.metaClust(),"list")))
        {validate(need(expr =  (inherits(x = dI.metaClust(),"list")), message = "please perform metaclust before"))}
      if (!(inherits(dI.clust$mapping)=="integer"))
        {validate(need(expr =  (inherits(x = dI.clust$mapping,"integer")), message = "you should have performed at least clustering first"))}
      if (!(inherits(x = dI.map(),"list")))
        {validate(need(expr =  (inherits(x = dI.map(),"list")), message = "you should have performed the mapping first"))}})
    
  })
  
  ################################################### __________Label II -----
  ################################################### runLabel -----    
  
  dI.label <- reactiveVal(NULL)
  dI.label_val <- reactiveValues(hm_label = NULL, hm_label_plus = NULL, hm_pheno = NULL, hm_pheno_plus = NULL)
  
  observeEvent(eventExpr = input$runLabel, handlerExpr = {
    clean_graph()
    
    if((is.data.frame(dI.labelClust$label.res))&&(inherits(x = fSmeta(),"flowSet"))&&(dI.flow$labelII)){ 
      # This hold only when you pass by the Label I procedure
      print("start labelling II process")
      output$console_output_clust <- renderPrint("start labelling II process - waiting...")
      if(inherits(x = dI.metaClust(),"list")){  
        # case 3: Labelling + Metaclustering
        output$label_panelStatus <- reactive({input$label_panelStatus=="show"})
        outputOptions(output, "label_panelStatus", suspendWhenHidden = FALSE)
        
        code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
        cell_clustering <- code_clustering[dI.clust$mapping]
        metaClustNum <- input$MetaClustNum #[1:Number of events]
        color_label <- color.pheno
        
        if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
          label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                phenoquery = ph, central = input$HeatMapCentral, min_quantile = input$minQuantile,  max_quantile = input$maxQuantile, 
                                pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                cluster_signature = dI.metaClust_val$cs.res[[2]])} else
          {label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                 phenoquery = ph, central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                 max_quantile = input$maxQuantile, pos_threshold = input$pos_threshold, 
                                 neg_threshold = input$neg_threshold)}

        dI.label(label.res)
        print("labelling process II ends")
        output$console_output_clust <- renderPrint("labelling process II ends")
        if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
        dI.label_val$hm_label <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "label", 
                                          clustering = cell_clustering, cluster_labelling = dI.label(), 
                                          Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                          central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                          pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                          dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                          color_clusters = color.metaclust, color_label = color_label, 
                                          cluster.signature = dI.metaClust_val$cs.res[[1]])}
        else{
          dI.label_val$hm_label <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "label", 
                                          clustering = cell_clustering, cluster_labelling = dI.label(), 
                                          Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                          central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                          pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                          dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                          color_clusters = color.metaclust, color_label = color_label)} 
        
        mm <- match(cell_clustering, dI.label()$original_cluster)
        cell_clustering2 <- dI.label()$new_cluster[mm]
        color_label <- dI.label_val$hm_label$color_label
        
        if ((nlevels(factor(cell_clustering2))) > 1){
          if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
          dI.label_val$hm_pheno <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                                            clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                                            Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                                            central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                            pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                            dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                            color_clusters = color_label, cluster_labelling = dI.metaClust_val$cs.res[[1]])}
          else{
            dI.label_val$hm_pheno <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                                              clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                                              Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                                              central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                              pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                              dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                              color_clusters = color_label)}}
        
        else{dI.label_val$hm_pheno <- "The second plot cannot be built: there must be at least two different labelled phenotypes"}}
      
      else{  # case 2 #only labelling
        output$label_panelStatus <- reactive({input$label_panelStatus=="show"})
        outputOptions(output, "label_panelStatus", suspendWhenHidden = FALSE)
       
        cell_clustering <- dI.clust$mapping 
        label.res <- dI.labelClust$label.res
        color_label <- color.pheno
        metaClustNum <- length(dI.labelClust$label.res$new_cluster)
        color_clusters <- dI.labelClust$hm_cluster$color_clusters
        if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
        dI.label_val$hm_label_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label_plus", 
                                                      clustering = cell_clustering, selected_markers = meta.markers()$selected, 
                                                      cluster_labelling = label.res, 
                                                      pheno_table = ph, Nclust = metaClustNum, 
                                                      central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                      pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                      dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                      color_clusters = color_clusters, color_label = color_label, 
                                                      cluster.signature = dI.labelClust$cs.res[[1]])}
        else{dI.label_val$hm_label_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label_plus", 
                                                                  clustering = cell_clustering, selected_markers = meta.markers()$selected, 
                                                                  cluster_labelling = label.res, 
                                                                  pheno_table = ph, Nclust = metaClustNum, 
                                                                  central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                                  pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                                  dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                                  color_clusters = color_clusters, color_label = color_label)}
        
        if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities"))
          {label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                  phenoquery = ph, central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                  pos_threshold =input$pos_threshold, neg_threshold = input$neg_threshold, 
                                  cluster_signature = dI.labelClust$cs.res[[1]])}else{
            label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                      phenoquery = ph, central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                      max_quantile = input$maxQuantile, pos_threshold = input$pos_threshold, 
                                      neg_threshold = input$neg_threshold)}
        dI.label(label.res)
        mm <- match(cell_clustering, dI.label()$original_cluster)
        cell_clustering2 <- dI.label()$new_cluster[mm]
        color_label <- dI.label_val$hm_label_plus$color_label
        
        if ((nlevels(factor(cell_clustering2))) > 1){
          if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
          dI.label_val$hm_pheno_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                                          clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                                          Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                                          central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                          pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                          dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                          color_clusters = color_label, cluster_labelling = dI.clust_val$cs.res[[1]])}
          else{
            dI.label_val$hm_pheno_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                                          clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                                          Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                                          central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                          pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                          dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                          color_clusters = color_label)}}
        
        else{dI.label_val$hm_pheno_plus <- "The second plot cannot be built: there must be at least two different labelled phenotypes"}}
    
    print("labelling process ends")
    output$console_output_clust <- renderPrint("labelling process ends")}
    
      if (!(inherits(x = fSmeta(),"flowSet")))
        {dI.label("You have not loaded any flowSet")}
      if (!(nrow(ph)>0))
        {dI.label("You cannot perform this step without setting the phenotype definition: edit and load the table in 'Metadata & Assays' workflow")}
      if(!(inherits(x = dI.clust$mapping,"integer")))
        {dI.label("You should cluster your flowSet first...")}
      if(!is.data.frame(dI.labelClust$label.res))
        {dI.label("Please perform the labelling in the 'Labelling I' tab")}
    
    output$checktxt_runLabel <- renderText({
      if (!is.null(dI.label()))
        {validate(need(expr = (is.null(dI.label())), message = dI.label()))}})
    
    output$checktxt_runLabel <- renderText({
      if (is.null(dI.label_val$hm_pheno_plus))
      {validate(need(expr = (!is.null(dI.label_val$hm_pheno_plus)), message = dI.label_val$hm_pheno_plus))}})
    })
  
  ################################################# runLabelPlot -------
 
  observeEvent(eventExpr = input$runLabelPlot, {
    
    clean_graph()
    if(is.data.frame(dI.label())) {
      print("start labelled map visualization")
      output$console_output_clust <- renderPrint("start labelled map visualization - waiting...")
      if(inherits(x = dI.metaClust(),"list")){
        output$labelPlot_panelStatus <- reactive({input$labelPlot_panelStatus=="show"})
        outputOptions(output, "labelPlot_panelStatus", suspendWhenHidden = FALSE)
        dI.label_val$hm_label_plus <- NULL
        dI.label_val$hm_pheno_plus <- NULL
        }
      else{
        output$labelPlot_panelStatusPlus <- reactive({input$labelPlot_panelStatusPlus=="show"})
        outputOptions(output, "labelPlot_panelStatusPlus", suspendWhenHidden = FALSE)
        dI.label_val$hm_label <- NULL
        dI.label_val$hm_pheno <- NULL
        } 

      if(is.list(dI.label_val$hm_label)){
        cdata <- session$clientData
        output$metaHeatmap_label <- renderImage({
          outfile <- "./tmpdata/HM_label.png"
          list(src = outfile, contentType = 'image/png', alt = "label heatmap",
          width = cdata$output_metaHeatmap_label_width, height = cdata$output_metaHeatmap_label_height)}, deleteFile = F)
      
        output$downloadLabel_csv <- downloadHandler(
        filename = function() {paste("HM_label","csv", sep = ".")},
        content <- function(file) {file.copy("./tmpdata/HM_label.csv", file)},
        contentType = "text/csv")
        
        path_to_file <- "./tmpdata"}
    
      if(is.list(dI.label_val$hm_pheno)){
        cdata <- session$clientData
        output$metaHeatmap_pheno <- renderImage({
          outfile <- "./tmpdata/HM_pheno.png"
          list(src = outfile, contentType = 'image/png', alt = "phenocluster heatmap",
               width = cdata$output_metaHeatmap_pheno_width, height = cdata$output_metaHeatmap_pheno_height)}, deleteFile = F)
        
        output$downloadPheno_csv <- downloadHandler(
          filename = function() {paste("HM_pheno","csv", sep = ".")},
          content <- function(file) {file.copy("./tmpdata/HM_pheno.csv", file)},
          contentType = "text/csv")
        
        path_to_file <- "./tmpdata"
        Zip_Files <- list.files(path = path_to_file, pattern = "^HM_", full.names = TRUE)
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
        zip(zipfile = "./tmpdata/plotLabel.zip", files = Zip_Files)
        
        output$downloadLabel <- downloadHandler(
          filename = function() {paste("output", "zip", sep=".")},
          content <- function(file) {file.copy("./tmpdata/plotLabel.zip", file)},
          contentType = "application/zip")}
      
      if(is.list(dI.label_val$hm_label_plus)){
        cdata <- session$clientData
        output$metaHeatmap_labelPlus <- renderImage({
          outfile <- "./tmpdata/HM_full_label_plus.png"
          list(src = outfile, contentType = 'image/png', alt = "label heatmap",
               width = cdata$output_metaHeatmap_labelPlus_width, height = cdata$output_metaHeatmap_labelPlus_height)}, deleteFile = F)
        
        output$downloadLabel_csvPlus <- downloadHandler(
          filename = function() {paste("HM_full_label_plus","csv", sep = ".")},
          content <- function(file) {file.copy("./tmpdata/HM_full_label_plus.csv", file)},
          contentType = "text/csv")}
      
      if(is.list(dI.label_val$hm_pheno_plus)){
        cdata <- session$clientData
        output$metaHeatmap_phenoPlus <- renderImage({
          outfile <- "./tmpdata/HM_pheno.png"
          list(src = outfile, contentType = 'image/png', alt = "phenocluster heatmap",
               width = cdata$output_metaHeatmap_phenoPlus_width, height = cdata$output_metaHeatmap_phenoPlus_height)}, deleteFile = F)
        
        output$downloadPheno_csvPlus <- downloadHandler(
          filename = function() {paste("HM_pheno","csv", sep = ".")},
          content <- function(file) {file.copy("./tmpdata/HM_pheno.csv", file)},
          contentType = "text/csv")
        
        path_to_file <- "./tmpdata"
        Zip_Files <- list.files(path = path_to_file, pattern = "^HM_", full.names = TRUE)
        Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
        zip(zipfile = "./tmpdata/plotLabel.zip", files = Zip_Files)
        output$downloadLabelPlus <- downloadHandler(
          filename = function() {paste("output", "zip", sep=".")},
          
          content <- function(file) {file.copy("./tmpdata/plotLabel.zip", file)},
          contentType = "application/zip")}}
      
      print("labelled map visualization process ends")
      output$console_output_clust <- renderPrint("labelled map visualization process ends")
      
      output$checktxt_runLabelPlot <- renderText({
        if (!is.null(dI.label_val$hm_label))
        {validate(need(expr =  (is.data.frame(dI.label_val$hm_label)), message = dI.label_val$hm_label))}
        if (!is.null(dI.label_val$hm_label_plus))
        {validate(need(expr =  (is.data.frame(dI.label_val$hm_label_plus)), message = dI.label_val$hm_label_plus))}
        })
      
      output$checktxt_runLabelPlot <- renderText({
        if (!is.null(dI.label_val$hm_pheno))
          {validate(need(expr =  (is.list(dI.label_val$hm_pheno)), message = dI.label_val$hm_pheno))}
        if (!is.null(dI.label_val$hm_pheno_plus))
          {validate(need(expr =  (is.list(dI.label_val$hm_pheno_plus)), message = dI.label_val$hm_pheno_plus))}
        })
    })
  
  
  #################################################_____Pheno maps I  -------
  
  dI.label_p_val <- reactiveValues(gg_map = NULL, gg_map_sample = NULL,  
                                   gg_map_tag1 = NULL,  gg_map_tag2 = NULL,  gg_map_tag3 = NULL, gg_map_tag4 = NULL)
  
  observeEvent(eventExpr = input$runLabelPlot_p, {
    showMaps <- "Recall to perform the analysis starting from clustering first and then the rest of the desidered workflow"
    clean_graph()
    
    if(((is.data.frame(dI.labelClust$label.res))||(inherits(x = dI.metaClust(),"list")))&&(inherits(x = dI.map(),"list"))&&(dI.flow$pheno_mapsI))
    { 
      #three cases:
      #1) Clustering + Metaclust dI.clust()&&dI.metaClust()
      #2) Clustering + Label (dI.clust()&&dI.label()) Not visibile anymore
      #3) Clustering + Label + MetaClust (dI.clust()&&dI.label()&&dI.metaClust())
      if ((inherits(x = dI.map(),"list"))&&(inherits(x = dI.clust$mapping,"integer")))
      {
        print("start labelled_p map visualization")
        output$console_output_clust <- renderPrint("start labelled map visualization - waiting...")
        output$labelPlotp_panelStatus <- reactive({input$labelPlotp_panelStatus=="show"})
        outputOptions(output, "labelPlotp_panelStatus", suspendWhenHidden = FALSE)
      
        if((is.list(dI.labelClust$hm_cluster))&&((inherits(x = dI.metaClust(),"list")))&&(is.list(dI.label_val$hm_pheno))){
            # case 3)
            code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
            cell_clustering <- as.integer(code_clustering[dI.clust$mapping])
            metaClustNum <- input$MetaClustNum
            mm <- match(cell_clustering, dI.label()$original_cluster)
            cell_label <- dI.label()$new_cluster[mm]
            map_cell <- factor(cell_label[dI.map()[[2]]])
            color_pheno <- dI.label_val$hm_pheno$color_clusters
            showMaps <- plot_labelling(flow_Set = fSmeta(), phenoclust = dI.label(), clustering = dI.clust$mapping, 
                                   m_clustering = cell_clustering , map_res = dI.map()[[1]], map_cell = map_cell, 
                                   selected_markers = mk$selected, metadata = mt, color_clusters = color_pheno)}
        
        if((is.data.frame(dI.labelClust$label.res))&&(!(inherits(x = dI.metaClust(),"list")))&&(inherits(x = fSmeta(),"flowSet"))&&(dI.flow$labelII)){ #case 2
            # case 2)
            
            cell_clustering <- dI.clust$mapping 
            mm <- match(cell_clustering, dI.label()$original_cluster)
            cell_clustering2 <- dI.label()$new_cluster[mm]
            color_label <- dI.label_val$hm_label_plus$color_label
            
            cell_label <- dI.label()$new_cluster[mm]
            map_cell <- factor(cell_label[dI.map()[[2]]])
            
            showMaps <- plot_labelling(flow_Set = fSmeta(), phenoclust = dI.label(), clustering = dI.clust$mapping, 
                                    m_clustering = cell_clustering2 , map_res = dI.map()[[1]], map_cell = map_cell, 
                                   selected_markers = mk$selected, metadata = mt, color_clusters = color_label)}
        
        if(!(is.data.frame(dI.labelClust$label.res))&&(inherits(x = dI.metaClust(),"list"))){
          # case 1)
          code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
          cell_clustering <- code_clustering[dI.clust$mapping]
          metaClustNum <- input$MetaClustNum
          original_cluster = seq(1, metaClustNum)
          new_cluster <- original_cluster
          for (i in 1:metaClustNum){new_cluster[i] <- as.character(original_cluster[i])}
        
          # metaClustNum <- input$MetaClustNum
          pheno_meta <- data.frame("original_cluster" = original_cluster, "new_cluster" = new_cluster)
          mm <- match(cell_clustering, pheno_meta$original_cluster)
          cell_label <- pheno_meta$new_cluster[mm]
          map_cell <- factor(cell_label[dI.map()[[2]]], levels = 1:metaClustNum) # add cell_label to map_df
          color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
          showMaps <- plot_labelling(flow_Set = fSmeta(), phenoclust = pheno_meta, clustering = dI.clust$mapping, 
                                   m_clustering = cell_clustering , map_res = dI.map()[[1]], map_cell = map_cell, 
                                   selected_markers = mk$selected, metadata = mt, color_clusters = color_metaclust)}
      
        if (inherits(x = showMaps,"list")){
          dI.label_p_val$gg_map <- showMaps[[1]]
          dI.label_p_val$gg_map_sample <- showMaps[[2]]
          if (uni.tag1 > 1)
            {dI.label_p_val$gg_map_tag1 <- showMaps[[3]]}
          if (uni.tag2 > 1)
            {dI.label_p_val$gg_map_tag2 <- showMaps[[4]]}
          if (uni.tag3 > 1)
            {dI.label_p_val$gg_map_tag3 <- showMaps[[5]]}
          if (uni.tag4 > 1)
            {dI.label_p_val$gg_map_tag4 <- showMaps[[6]]}
      
          output$showLabeldMap <- renderPlotly({dI.label_p_val$gg_map})
          output$showSampleMap <- renderPlotly({dI.label_p_val$gg_map_sample})
          print("labelled map_p visualization process ends")
          output$console_output_clust <- renderPrint("labelled map visualization process ends")}}
      
      dI.flow$pheno_mapsI <- FALSE}
    
      output$checktxt_runPhenomapsI <- renderText({
        if (!((is.data.frame(dI.labelClust$label.res))||(inherits(x = dI.metaClust(),"list")))){
          testo <- paste0("You should perform the labelling action and/or the meta-clustering to produce a colored map. ",
                          "Recall to always start from the clustering step first")
        validate(need(expr =  (is.list(dI.label_val$hm_label)), message = testo))}
        if (!(inherits(x = dI.map(),"list"))){
          testo <- paste0("You should run the mapping algoritm first")
          validate(need(expr =  ((inherits(x = dI.map(),"list"))), message = testo))}
        if (!(inherits(x = showMaps,"list"))){
          validate(need(expr =  (inherits(x = showMaps,"list")), message = showMaps))}})
  })
  
  
  #################################################_____Pheno maps II  -------
  
  observeEvent(eventExpr = input$runLabelPlot_pp, {
    clean_graph()
    
    if(((is.data.frame(dI.labelClust$label.res))||(inherits(x = dI.metaClust(),"list")))&&(inherits(x = dI.map(),"list"))&&(dI.flow$pheno_mapsII)){
      if((inherits(x = dI.map(),"list"))&&(inherits(x = dI.clust$mapping,"integer"))) {
    
        print("start labelled_++ map visualization")
        output$console_output_clust <- renderPrint("start labelled map visualization - waiting...")
      
        output$labelPlotpp_panelStatus <- reactive({input$labelPlotpp_panelStatus=="show"})
        outputOptions(output, "labelPlotpp_panelStatus", suspendWhenHidden = FALSE)
        
        if (uni.tag1 > 1) {output$showTag1Map <- renderPlotly({dI.label_p_val$gg_map_tag1})}
        if (uni.tag2 > 1) {output$showTag2Map <- renderPlotly({dI.label_p_val$gg_map_tag2})}
        if (uni.tag3 > 1) {output$showTag3Map <- renderPlotly({dI.label_p_val$gg_map_tag3})}
        if (uni.tag4 > 1) {output$showTag4Map <- renderPlotly({dI.label_p_val$gg_map_tag4})}
        
        print("labelled map_++ visualization process ends")
        output$console_output_clust <- renderPrint("labelled map++ visualization process ends")}
    }
    
    dI.flow$pheno_mapsII <- FALSE
    
    output$checktxt_runPhenomapsII <- renderText({
      if (!((is.data.frame(dI.labelClust$label.res))||(inherits(x = dI.metaClust(),"list")))){
        testo <- paste0("You should perform the labelling action and/or the meta-clustering to produce a colored map. ",
                        "Recall to always start from the clustering step first")
        validate(need(expr =  (is.list(dI.label_val$hm_label)), message = testo))}
      if (!(inherits(x = dI.map(),"list"))){
        testo <- paste0("You should run the mapping algoritm first")
        validate(need(expr =  (inherits(x = dI.map(),"list")), message = testo))}
      if (is.null(dI.label_p_val$gg_map_sample)){
        testo <- paste0("Please go to the 'Pheno maps I' first")
        validate(need(expr =  (!is.null(dI.label_p_val$gg_map_sample)), message = testo))}})
  })  
  
  #################################################_____Quantity plots  -------
  
  dI.quant_val <- reactiveValues(bar_plot_sample = NULL, bar_plot_tag1 = NULL, bar_plot_tag2 = NULL, bar_plot_tag3 = NULL, bar_plot_tag4 = NULL, 
                                 dodge_plot_sample = NULL, dodge_plot_tag1 = NULL, dodge_plot_tag2 = NULL, dodge_plot_tag3 = NULL, dodge_plot_tag4 = NULL, 
                                 med_expr_tag1 = NULL, med_expr_tag2 = NULL, med_expr_tag3 = NULL, med_expr_tag4 = NULL)
  
  observeEvent(eventExpr = input$runQuantPlot, {
    
    clean_graph()
    if((inherits(x = dI.clust$mapping,"integer"))&&(dI.flow$quantity))  {
      #three cases:
      #1) Clustering + Metaclust dI.clust()&&dI.metaClust()
      #2) Clustering + Label (dI.clust()&&dI.label())
      #3) Clustering + Label + MetaClust (dI.clust()&&dI.label()&&dI.metaClust())
      print("start quantity plot visualization")
      output$console_output_clust <- renderPrint("start quantity plot visualization - waiting...")
      
      showBars <- NULL
      showMedian <- NULL
    
      if(is.list(dI.labelClust$hm_cluster)){
        if((inherits(x = dI.metaClust(),"list"))&&(is.list(dI.label_val$hm_pheno))){
          output$QuantPlot_panelStatus <- reactive({input$QuantPlot_panelStatus=="show"})
          outputOptions(output, "QuantPlot_panelStatus", suspendWhenHidden = FALSE)
          #case 3)
          code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
          cell_clustering <- code_clustering[dI.clust$mapping]
          metaClustNum <- input$MetaClustNum #[1:Number of events]
          
          mm <- match(cell_clustering, dI.label()$original_cluster)
          cell_label <- dI.label()$new_cluster[mm]
          
          color_pheno <- dI.label_val$hm_pheno$color_clusters
          showBars <- plot_labellingbar(flow_Set = fSmeta(), flag = FALSE, cell_label = cell_label, phenoclust = dI.label(),
                                      NMC = metaClustNum, metadata =  mt, color_pheno = color_pheno) 
          
          showMedian <- plot_median_expr(flow_Set = fSmeta(), phenoclust =  dI.label(), cell_clustering = cell_clustering, 
                                         metadata = mt, NMC = metaClustNum,
                                         color_tag1 = color.tag1, color_tag2 = color.tag2, color_tag3 = color.tag3, color_tag4 = color.tag4)}
        else{
          #if ((is.list(dI.label_val$hm_label_plus))){
          if ((is.list(dI.labelClust$hm_cluster))){
            
            #case 2) with label II performed
            if ((is.list(dI.label_val$hm_label_plus))){
              output$QuantPlot_panelStatus <- reactive({input$QuantPlot_panelStatus=="show"})
              outputOptions(output, "QuantPlot_panelStatus", suspendWhenHidden = FALSE)
              cell_clustering <- dI.clust$mapping 
              mm <- match(cell_clustering, dI.label()$original_cluster)
              cell_clustering <- dI.label()$new_cluster[mm]
              cell_label <- dI.label()$new_cluster[mm]
              dI_label <- dI.label()
              
              metaClustNum <- nlevels(factor(dI.label()$new_cluster))
              color_label <- dI.label_val$hm_label_plus$color_label
              cell_clustering <- dI.clust$mapping } else 
                #case 2) without pheno-data 
              {
                output$QuantPlot_panelStatus <- reactive({input$QuantPlot_panelStatus=="show"})
                outputOptions(output, "QuantPlot_panelStatus", suspendWhenHidden = FALSE)
                
                cell_label <- dI.clust$mapping
                dI_label <- substring(dI.labelClust$hm_cluster$labels_row, first=1, last=3)
                dI_label <-  gsub(pattern = "\\(", x = dI_label, replacement = "")
                dI_label <-  as.numeric(gsub(pattern = "\\s", "", x = dI_label))
                
                dI_label <- data.frame(original_cluster = dI_label, new_cluster = dI_label)
                metaClustNum <- nlevels(factor(dI.clust$mapping)) #metaClustNum <- nlevels(factor(dI.label()$new_cluster))
                color_label <- dI.labelClust$hm_cluster$color_clusters #color_label <- dI.label_val$hm_label_plus$color_label}  
                
                showBars <- plot_labellingbar(flow_Set = fSmeta(), flag = FALSE, cell_label = cell_label, phenoclust = dI_label(),
                                              NMC = metaClustNum, metadata =  mt, color_pheno = color_label) 
                
                showMedian <- plot_median_expr(flow_Set = fSmeta(), phenoclust =  dI_label, cell_clustering = cell_label, 
                                               metadata = mt, NMC = metaClustNum, 
                                               color_tag1 = color.tag1, color_tag2 = color.tag2, color_tag3 = color.tag3, color_tag4 = color.tag4)}}}}
      
      if(!(is.data.frame(dI.labelClust$label.res))&&(inherits(x = dI.metaClust(),"list"))){
        #case 1)
          output$QuantPlot_panelStatus <- reactive({input$QuantPlot_panelStatus=="show"})
          outputOptions(output, "QuantPlot_panelStatus", suspendWhenHidden = FALSE)
          code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
          cell_clustering <- code_clustering[dI.clust$mapping]
          #metaClustNum <- input$MetaClustNum 
          metaClustNum <- length(levels(factor(cell_clustering))) 
          original_cluster = seq(1, metaClustNum)
          new_cluster <- original_cluster
          for (i in 1:metaClustNum){new_cluster[i] <- as.character(original_cluster[i])}
          pheno_meta <- data.frame("original_cluster" = original_cluster, "new_cluster" = new_cluster)
          
          mm <- match(cell_clustering, pheno_meta$original_cluster)
          cell_label <- factor(pheno_meta$new_cluster[mm], levels = 1:metaClustNum)
          
          color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
          showBars <- plot_labellingbar(flow_Set = fSmeta(), flag = TRUE, cell_label = cell_label, phenoclust = pheno_meta,
                                      NMC = metaClustNum, metadata =  mt, color_pheno = color_metaclust) 
          
          showMedian <- plot_median_expr(flow_Set = fSmeta(), phenoclust =  pheno_meta, cell_clustering = cell_clustering, 
                                         metadata = mt, NMC = metaClustNum, 
                                         color_tag1 = color.tag1, color_tag2 = color.tag2, color_tag3 = color.tag3, color_tag4 = color.tag4)}
      
      if (inherits(x = showBars,"list")){
        dI.quant_val$bar_plot_sample <- showBars[[1]]
        dI.quant_val$dodge_plot_sample <- showBars[[6]]
        if (uni.tag1 > 1){
          dI.quant_val$bar_plot_tag1 <- showBars[[2]]
          dI.quant_val$dodge_plot_tag1 <- showBars[[7]]}
        if (uni.tag2 > 1){
          dI.quant_val$bar_plot_tag2 <- showBars[[3]]
          dI.quant_val$dodge_plot_tag2 <- showBars[[8]]}
        if (uni.tag3 > 1){
          dI.quant_val$bar_plot_tag3 <- showBars[[4]]
          dI.quant_val$dodge_plot_tag3 <- showBars[[9]]}
        if (uni.tag4 > 1){
          dI.quant_val$bar_plot_tag4 <- showBars[[5]]
          dI.quant_val$dodge_plot_tag4 <- showBars[[10]]}
      
        output$showBarSample <- renderPlotly({dI.quant_val$bar_plot_sample})
        output$showDodgeSample <- renderPlotly({dI.quant_val$dodge_plot_sample})
        output$showBarTag1 <- renderPlotly({dI.quant_val$bar_plot_tag1})
        output$showDodgeTag1 <- renderPlotly({dI.quant_val$dodge_plot_tag1})
        output$showBarTag2 <- renderPlotly({dI.quant_val$bar_plot_tag2})
        output$showDodgeTag2 <- renderPlotly({dI.quant_val$dodge_plot_tag2})
        output$showBarTag3 <- renderPlotly({dI.quant_val$bar_plot_tag3})
        output$showDodgeTag3 <- renderPlotly({dI.quant_val$dodge_plot_tag3})
      
        if (uni.tag1 > 1)
          {dI.quant_val$med_expr_tag1 <- showMedian[[1]]}
        if (uni.tag2 > 1)
          {dI.quant_val$med_expr_tag2 <- showMedian[[2]]}
        if (uni.tag3 > 1)
          {dI.quant_val$med_expr_tag3 <- showMedian[[3]]}
        if (uni.tag4 > 1)
          {dI.quant_val$med_expr_tag4 <- showMedian[[4]]}
        
        if (uni.tag1 > 1) {output$med_expr_tag1 <- renderPlotly({dI.quant_val$med_expr_tag1})}
        if (uni.tag2 > 1) {output$med_expr_tag2 <- renderPlotly({dI.quant_val$med_expr_tag2})}
        if (uni.tag3 > 1) {output$med_expr_tag3 <- renderPlotly({dI.quant_val$med_expr_tag3})}
        if (uni.tag4 > 1) {output$med_expr_tag4 <- renderPlotly({dI.quant_val$med_expr_tag4})}
        
        print("quantity visualization process ends")
        output$console_output_clust <- renderPrint("quantity plot visualization process ends")
        
        output$downloadQuantSample_csv <- downloadHandler(filename = function() {paste("df_sample","csv", sep = ".")},
                                              content <- function(file) {file.copy("./tmpdata/df_sample.csv", file)}, contentType = "text/csv")
        
        if (uni.tag1 > 1){
          output$downloadQuantTag1_csv <- downloadHandler(filename = function() {paste0("df_",names(mt)[3],".csv")},
                                              content <- function(file) {file.copy(paste0("./tmpdata/df_",names(mt)[3],".csv"), file)}, contentType = "text/csv")}
        if (uni.tag2 > 1){
          output$downloadQuantTag2_csv <- downloadHandler(filename = function() {paste0("df_",names(mt)[4],".csv")},
                                              content <- function(file) {file.copy(paste0("./tmpdata/df_",names(mt)[4],".csv"), file)}, contentType = "text/csv")}
        if (uni.tag3 > 1){
          output$downloadQuantTag3_csv <- downloadHandler(filename = function() {paste0("df_",names(mt)[5],".csv")},
                                              content <- function(file) {file.copy(paste0("./tmpdata/df_",names(mt)[5],".csv"), file)}, contentType = "text/csv")}
        if (uni.tag4 > 1){
          output$downloadQuantTag4_csv <- downloadHandler(filename = function() {paste0("df_",names(mt)[6],".csv")},
                                              content <- function(file) {file.copy(paste0("./tmpdata/df_",names(mt)[6],".csv"), file)}, contentType = "text/csv")}}
      else{ showBars <- "Impossible to perform the quantity analysis"}
      
      dI.flow$quantity <- FALSE
      if (Sys.info()["sysname"] == "Windows") play(x = bell)}
    
    output$checktxt_runQuantPlot <- renderText({
      if (!(inherits(x = dI.clust$mapping,"integer"))){
        testo <- "Please perform the clustering operation firts"
        validate(need(expr =  (inherits(x = dI.clust$mapping,"integer")), message = testo))}
      if (!(inherits(x = showBars,"list"))){
        validate(need(expr =  (inherits(x = showBars,"list")), message = showBars))}
      if (!(inherits(x = dI.labelClust$hm_cluster,"list"))){
        validate(need(expr =  (inherits(x = dI.labelClust$hm_cluster,"list")), message = "Please perform the label I process first"))}
      })
}) #input
  
  
  #################################################____Time-step plots  -------
  
  
  dI.quant_val_p <- reactiveValues(showStream = NULL)
  
  observeEvent(eventExpr = input$runQuantPlot_p, {
    
    clean_graph()
    if((inherits(x = dI.clust$mapping,"integer"))&&(uni.tag4 > 1)&&(dI.flow$time_step)) {
    if (!(dI.flow$quantity)){
      
      #three cases:
      #1) Clustering + Metaclust dI.clust()&&dI.metaClust()
      #2) Clustering + Label (dI.clust()&&dI.label())
      #3) Clustering + Label + MetaClust (dI.clust()&&dI.label()&&dI.metaClust())
      
      print("start quantity plot visualization")
      output$console_output_clust <- renderPrint("start quantity plot visualization - waiting...")
      
      if(is.list(dI.labelClust$hm_cluster)){
      
        if((inherits(x = dI.metaClust(),"list"))&&((is.list(dI.label_val$hm_pheno_plus)))){
          #case 3)
          output$QuantPlotp_panelStatus <- reactive({input$QuantPlotp_panelStatus=="show"})
          outputOptions(output, "QuantPlotp_panelStatus", suspendWhenHidden = FALSE)
          
          code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
          cell_clustering <- code_clustering[dI.clust$mapping]
          metaClustNum <- input$MetaClustNum
          mm <- match(cell_clustering, dI.label()$original_cluster)
          cell_label <- dI.label()$new_cluster[mm]
          
          color_pheno <- dI.label_val$hm_pheno$color_clusters
          if (input$time_step == TRUE){
            dI.quant_val_p$showStream <- plot_stream(flow_Set = fSmeta(), flag = FALSE, phenoclust = dI.label(), cell_label = cell_label, 
                                                     metadata = mt, NMC = metaClustNum, color_pheno = color_pheno)}}
        
          #case 2) with label II performed
          if ((is.list(dI.labelClust$hm_cluster))){
            if ((is.list(dI.label_val$hm_label_plus))){
          output$QuantPlotp_panelStatus <- reactive({input$QuantPlotp_panelStatus=="show"})
          outputOptions(output, "QuantPlotp_panelStatus", suspendWhenHidden = FALSE)
          
          code_clustering <- dI.labelClust$label.res$numeric
          cell_clustering <- code_clustering[dI.clust$mapping]
          metaClustNum <- length(dI.labelClust$label.res$numeric)
          mm <- match(cell_clustering, dI.label()$original_cluster)
          cell_label <- dI.label()$new_cluster[mm]
          if (input$time_step == TRUE){
            dI.quant_val_p$showStream <- plot_stream(flow_Set = fSmeta(), flag = FALSE, phenoclust = dI.label(), cell_label = cell_label, 
                                                     metadata = mt, NMC = metaClustNum, color_pheno = color.pheno)}}
          else{
            #case 2) without pheno-data 
            output$QuantPlotp_panelStatus <- reactive({input$QuantPlotp_panelStatus=="show"})
            outputOptions(output, "QuantPlotp_panelStatus", suspendWhenHidden = FALSE)
            
            cell_label <- dI.clust$mapping
            dI_label <- substring(dI.labelClust$hm_cluster$labels_row, first=1, last=3)
            dI_label <-  gsub(pattern = "\\(", x = dI_label, replacement = "")
            dI_label <-  as.numeric(gsub(pattern = "\\s", "", x = dI_label))
            
            dI_label <- data.frame(original_cluster = dI_label, new_cluster = dI_label)
            metaClustNum <- nlevels(factor(dI.clust$mapping)) #metaClustNum <- nlevels(factor(dI.label()$new_cluster))
            color_label <- dI.labelClust$hm_cluster$color_clusters #color_label <- dI.label_val$hm_label_plus$color_label}  
            
            if (input$time_step == TRUE){
              dI.quant_val_p$showStream <- plot_stream(flow_Set = fSmeta(), flag = FALSE, phenoclust = dI_label(), cell_label = cell_label, 
                                                       metadata = mt, NMC = metaClustNum, color_pheno = color_label)}}}} 
      if(!(is.data.frame(dI.labelClust$label.res))&&(inherits(x = dI.metaClust(),"list"))){
          #case 1)
          output$QuantPlotp_panelStatus <- reactive({input$QuantPlotp_panelStatus=="show"})
          outputOptions(output, "QuantPlotp_panelStatus", suspendWhenHidden = FALSE)
          
          code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
          cell_clustering <- code_clustering[dI.clust$mapping]
          #metaClustNum <- input$MetaClustNum
          metaClustNum <- length(levels(factor(cell_clustering))) 
          original_cluster = seq(1, metaClustNum)
          new_cluster <- original_cluster
          for (i in 1:metaClustNum){new_cluster[i] <- as.character(original_cluster[i])}
          #new_cluster <- as.character(original_cluster)
          pheno_meta <- data.frame("original_cluster" = original_cluster, "new_cluster" = new_cluster)
          mm <- match(cell_clustering, pheno_meta$original_cluster)
          cell_label <- factor(pheno_meta$new_cluster[mm])
          color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
          if (input$time_step == TRUE){
            dI.quant_val_p$showStream <- plot_stream(flow_Set = fSmeta(), flag = TRUE, phenoclust = pheno_meta, cell_label = cell_label, 
                                                     metadata = mt, NMC = metaClustNum, color_pheno = color_metaclust)}}
      
      output$showBarTag4 <- renderPlotly({dI.quant_val$bar_plot_tag4})
      output$showDodgeTag4 <- renderPlotly({dI.quant_val$dodge_plot_tag4})
      
      output$med_expr_tag4 <- renderPlotly({dI.quant_val$med_expr_tag4})
      
      if (inherits(x = dI.quant_val_p$showStream,"list")){
        output$stream_plot_prop <- renderStreamgraph({dI.quant_val_p$showStream[[1]]})
        output$stream_plot_count <- renderStreamgraph({dI.quant_val_p$showStream[[2]]})
        
        print("quantity visualization process ends")
        output$console_output_clust <- renderPrint("quantity plot visualization process ends")}
      else
      {dI.quant_val_p$showStream ="maybe you forgot something..."}
      dI.flow$time_step <- FALSE
    } }
    
  output$checktxt_runQuantPlot_p <- renderText({
      if (!(inherits(x = dI.quant_val_p$showStream,"list")))
      {validate(need(expr =  (inherits(x = dI.quant_val_p$showStream,"list")), message =  dI.quant_val_p$showStream))}
      if (!(uni.tag4>1))
      {validate(need(expr =  (uni.tag4 > 1), message = "It needs at least two different steps to produce the stream plots!"))}
      if (!(inherits(x = dI.clust$mapping,"integer")))
      {validate(need(expr =  (inherits(x = dI.clust$mapping,"integer")), message = "You have to perform the clustering operation first!"))}
      if (dI.flow$quantity)
      {validate(need(expr = (dI.flow$quantity == FALSE), message = "Please , perform the Quantity plots process first..."))}
      if (!(is.list(dI.labelClust$hm_cluster)))
      {validate(need(expr = (inherits(x = dI.labelClust$hm_cluster,"list")), message = "Please , perform the label I process first..."))}
    })
  }) # end of observeEvent
  
  
  
    #################################################____Panel editor  -------
    dI.PE <- reactiveVal(NULL)
    dI.PE_val <- reactiveValues(fF_names = NULL)  
    
    observeEvent(eventExpr = input$panelEdit, {
      
      validate(need(expr = !(is.null(input$FCS_input)&&is.null(input$zip_input)), message = "please select some samples"))
      print("start parsing process")
      output$console_output_pre <- renderPrint("start parsing process  - waiting...")
      options(warn = 0)
      png()
      clean_graph()
      del_tmpdata(type = "all")
      panel.table <- NULL
      
      if(!is.null(input$FCS_input)) {
        alphafSinput <-  input$FCS_input[order(input$FCS_input$name),] 
        lista <- loading_fS(inputfile =  alphafSinput$datapath, tipo = "flowSet")
        fS <- lista[[1]]
        if (inherits(x = fS,"flowSet")){
          array_nomi <- alphafSinput$name #was array_nomi <- lista[[2]]
          dI.PE_val$fF_names <- array_nomi
          fFvectorName <- vector(mode = "character", length = length(alphafSinput$name))
          for (i in seq_along(1:length(alphafSinput$name))){
            fFvectorName[i] = paste0("./tmpdata/", dI.PE_val$fF_names[[i]])
            nome <- panel_fcs(flow_Frame = alphafSinput$datapath[i], stringa = fFvectorName[i])}
            entry_check <- "OK"}
        else{
          array_nomi <- "no_flowSet"
          entry_check <- "fail"}}
      
      if(!is.null(input$zip_input)) {
        lista <- loading_fS(inputfile = input$zip_input$datapath, tipo = "zip")
        fS <- lista[[1]]
        if (inherits(x = fS,"flowSet")){
          array_nomi <- lista[[2]]
          dI.PE_val$fF_names <- array_nomi
          fFvectorName <- vector(mode = "character", length = length(array_nomi))
          for (i in seq_along(1:length(array_nomi))){
            fFvectorName[i] = paste0("./tmpdata/", dI.PE_val$fF_names[[i]])
            nome <- panel_fcs(flow_Frame = fFvectorName[i], stringa = fFvectorName[i])}
            entry_check <- "OK"}
        else{
          array_nomi <- "no_flowSet"
          entry_check <- "fail"}}
        
      if (entry_check == "OK"){
        {
          lista_files <- list.files("./tmpdata", pattern = "*.fcs$", ignore.case = T)
          #lista_files <- str_sort(x = lista_files, numeric = TRUE) 
          files.list <- file.path("./tmpdata", lista_files)
          if (length(files.list) > 0){
            panel.table <- get_panel_table(files.list)}
          output$checktxt_parsing_panel <- renderText({
            if (entry_check=="OK"){
              {messaggio <- "All the selected files seem to be have the correct FCS3.0 format. 
                Check weather they could represents a flowSet with a common names & descriptions for every single flowFrame.
                If this is not the case, then change the table below to fix the mismatching names, process files and dowload your new sample set"
              validate(need(expr =  (entry_check=="fail"), message = messaggio))}}})
          }
        
        if (length(panel.table) > 0){
          observe({
            
            if(!is.null(input$paneleditorui_process_files) && input$paneleditorui_process_files != 0) {
              
              isolate({
                
                df <- rhandsontable::hot_to_r(input$paneleditorui_panel_table)
                for(i in 1:ncol(df)){
                  df[, i] <- gsub("absent", NA, df[, i])}
                df$Remove <- as.logical(df$Remove)
                panel.table$"Most common" <- df$"Most common" <- NULL
                
                showModal(modalDialog(title = "Panel editor report", "File processing started, please wait..."))
                
                rename_parameters_in_files(working.dir = "./tmpdata", out.dir = "renamed", df)
                showModal(modalDialog(title = "Panel editor report",
                                      sprintf("File processing completed. Download your new samples")))
                
                
                handled_file <- list.files(path = "./tmpdata/renamed", full.names = TRUE)
                lista <- loading_fS(inputfile = handled_file, tipo = "flowSet")
                if (inherits(x = lista[[1]],"flowSet")){
                  dI.PE(lista[[1]])} else
                  {
                    output$checktxt_parsing_panel <- renderText({
                    messaggio <- "Something went wrong on handling the entered flowFrames"
                    validate(need(expr =  (entry_check=="fail"), message = messaggio))})}
                
                output$PE_panelStatus <- reactive({input$PE_panelStatus=="show"})
                outputOptions(output, "PE_panelStatus", suspendWhenHidden = FALSE)})}
          })}
        
        if (!is.null(panel.table)){
          output$paneleditorui_panel_table <- rhandsontable::renderRHandsontable({
            temp <- panel.table
            
            for(i in 1:ncol(temp))
            {temp[, i][is.na(temp[, i])] <- "absent"}
            df <- data.frame(Remove = FALSE, temp, check.names = F, stringsAsFactors = F)
            hot <- rhandsontable::rhandsontable(df, rowHeaderWidth = 100)
            hot <- rhandsontable::hot_cols(hot, fixedColumnsLeft = 3, renderer = "
            function(instance, td, row, col, prop, value, cellProperties) {
                if(col == 0)
                    Handsontable.renderers.CheckboxRenderer.apply(this, arguments)
                else {
                Handsontable.renderers.TextRenderer.apply(this, arguments)
                    if(instance.params != null) {
                        if(instance.params.data[row][0])
                            td.style.background = 'lightgrey'
                        else {
                            if(value == 'absent')
                                td.style.background = 'orange'
                            else if(value != instance.params.data[row][2] && col > 2)
                                td.style.background = 'lightpink'}
                    }
                }
                return(td)
            }")
            hot})
        }
      }
      else{
        messaggio <- paste0("Something is wrong with your entered file: System reported: ",fS)
        output$checktxt_parsing_panel <- renderText({
          validate(need(expr =  (!inherits(x = messaggio,"character")), message = messaggio))})}
      
      output$download_FCS <- downloadHandler(
        filename = function() {
          paste("newFCS", "zip", sep = ".")
        },
        content = function(file){
          
          file_path <- file
          tmpdir <- tempdir()
          setwd(tempdir())
          fSlength <- length(dI.PE())
          for (i in seq_along(1:fSlength)) {
            fFvectorName[[i]] <- paste0("new_", dI.PE_val$fF_names[[i]])
            write.FCS(dI.PE()[[i]], filename = fFvectorName[[i]])}
          
          Zip_Files <- list.files(path = getwd(), pattern = "^new*")
          #Zip_Files <- str_sort(x = Zip_Files, numeric = TRUE) 
          zip(zipfile = file_path, files = Zip_Files)
          setwd(app_dir)
        },
        contentType = "application/zip")  
    })
    
    #################################################____Auto Complete  -------
    observeEvent(eventExpr = input$runAuto, {
 
        if((inherits(x = fSmeta(), "flowSet"))&&(!is.null(meta.table()))&&(any(mk$selected, na.rm = TRUE))) {
        output$console_output_clust <- renderPrint("analysis started")
        dI.metaClust(NULL) #this is to reset the metaclustering results
        dI.labelClust$label.res <- NULL #this is to reset the labelling I results
        dI.labelClust$hm_cluster <- NULL #this is to reset the labelling I results
        dI.label(NULL) #this is to reset the labelling II results
        dI.labelClust$label.res <- NULL; #dI.map_val$showmapcomp <- NULL; #this to try to avoid the re-generation of the tsne map
        dI.label_val$hm_label <- NULL; dI.label_val$hm_label_plus <- NULL
        dI.label_val$hm_pheno <- NULL; dI.label_val$hm_pheno_plus <- NULL
        
        dI.label_p_val$gg_map <- NULL; dI.label_p_val$gg_map_sample <- NULL;  
        dI.label_p_val$gg_map_tag1 <- NULL;  dI.label_p_val$gg_map_tag2 <- NULL;  dI.label_p_val$gg_map_tag3 <- NULL; dI.label_p_val$gg_map_tag4 <- NULL
        
        dI.quant_val$bar_plot_sample <- NULL; 
        dI.quant_val$dodge_plot_sample <- NULL; 
        dI.quant_val$bar_plot_tag1 = NULL; dI.quant_val$bar_plot_tag2 <- NULL; dI.quant_val$bar_plot_tag3 <- NULL; dI.quant_val$bar_plot_tag4 <- NULL
        dI.quant_val$dodge_plot_tag1 = NULL; dI.quant_val$dodge_plot_tag2 <- NULL; dI.quant_val$dodge_plot_tag3 <- NULL; dI.quant_val$dodge_plot_tag4 <- NULL
        dI.quant_val$med_expr_tag1 = NULL; dI.quant_val$med_expr_tag2 <- NULL; dI.quant_val$med_expr_tag3 <- NULL; dI.quant_val$med_expr_tag4 <- NULL 
        
        dI.flow$labelI <- TRUE; dI.flow$meta_clust <- TRUE; dI.flow$labelII <- TRUE; 
        dI.flow$pheno_mapsI <- TRUE; dI.flow$pheno_mapsII <- TRUE; dI.flow$quantity <- TRUE; dI.flow$time_step <- TRUE
        
        del_tmpdata("all")
        ## MDS plot ----
        
        clean_graph()
        print("start saving MDS plots")
        if(length(fSmeta())>2){
          dI.ana$MDS.sample_id <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                              type = "sample_id", selected = mk$selected, save = TRUE, metodo = input$mdsMethod) 
          if (uni.tag1>1){dI.ana$MDS.tag1 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                         type = "tag1", selected = mk$selected, save = TRUE, metodo = input$mdsMethod)}
          if (uni.tag2>1){dI.ana$MDS.tag2 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                         type = "tag2", selected = mk$selected, save = TRUE, metodo = input$mdsMethod)}
          if (uni.tag3>1){dI.ana$MDS.tag3 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                         type = "tag3", selected = mk$selected, save = TRUE, metodo = input$mdsMethod)}
          if (uni.tag4>1){dI.ana$MDS.tag4 <- MDSplot_new(flow_Set = fSmeta(), meta = meta.table(), 
                                                         type = "tag4", selected = mk$selected, save = TRUE, metodo = input$mdsMethod)}}
        
        #if(class(dI.ana$MDS.sample_id[[1]])[1]=="plotly"){print("MDS plots saved")}
        if(inherits(x = dI.ana$MDS.sample_id[[1]], "plotly")){print("MDS plots saved")}
        ### PCA plot ----
        print("evaluation PCA plot process start")
        pca_list <- evaluation_plotPCA(flow_Set = fSmeta(), seed = input$set_seed, selected = mk$selected)
        
        ggsave(filename = "PCA.png", plot = pca_list[[1]], device = "png", path = "./tmpdata/")
        write.csv(x = pca_list[[2]],file = "./tmpdata/PCA.csv", row.names = FALSE)
        write.csv(x = pca_list[[3]],file = "./tmpdata/PCA_contribution.csv", row.names = FALSE)
        ggsave(filename = "PCAcos2_plot.png", plot = pca_list[[4]], device = "png", path = "./tmpdata/")
        ggsave(filename = "PCAcos2_var.png", plot = pca_list[[5]], device = "png", path = "./tmpdata/")
        capture.output(metaResume(), file = "./tmpdata/meta_resume.csv")
        
        print("evaluation PCA plot process ends")
        
        ############################################################ 1. Workflow #1 & #3-----
        ### Clustering ----
        clean_graph()
        color.clust <<- color_set[1:input$k]
        print("start clustering process")
      
        if (input$ClustAlgo == "flowSOM"){
          flowSOM.res <- flowSOM_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, seed = input$set_seed)
          if (inherits(flowSOM.res[[1]],"FlowSOM")){
            dI.clust$mapping <- as.integer(flowSOM.res[[2]]$map$mapping[,1]) # vector [flowSet nr. of cell = cluster nr.]
            dI.clust$codes <- flowSOM.res[[2]]$map$codes # matrix [metaclust nr. x selected_markers]
            dI.clust$mst = flowSOM.res[[3]]
            dI.clust$NumCluster <- 100 #flowSOM uses always 100 clusters}
            color.clust <<- color_set[1:dI.clust$NumCluster]} else{
              dI.clust$mapping <- flowSOM.res}}
        
        if (input$ClustAlgo == "KMeans"){
          kmeans.res <- kmeans_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, Algo = input$KmeansAlgo,
                                     K = input$k, seed = input$set_seed)
          if (inherits(x = kmeans.res,"list")){
            dI.clust$mapping <- kmeans.res[[1]]
            dI.clust$codes <- kmeans.res[[2]]
            dI.clust$NumCluster <- input$k
            color.clust <<- color_set[1:dI.clust$NumCluster]}else{
              dI.clust$mapping <- kmeans.res}}
        
        if ((input$ClustAlgo == "Rphenograph")||(input$ClustAlgo == "FastPG")){
          if (input$ClustAlgo == "Rphenograph"){
            phenograph.res <- phenograph_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, 
                                               K = input$k, seed = input$set_seed, implementation = "Rphenograph")
            if (inherits(x = phenograph.res, "list")){
              eventi <- as.character(phenograph.res[[1]])
              numCluster <- as.character(nlevels(phenograph.res[[1]]))
              para <- as.character(input$k)
              modularity <- phenograph.res[[3]] #this is introduced with fastPG
            text_report <- paste0("RPhenograph (the R implementation of Phenograph clustering algorithm), run with K = ",
                                  para , ", grouping the dataset of ",eventi , " events in ",
                                  numCluster , " clusters, with ", modularity, " as modularity")
            cat(text_report,file="./tmpdata/Rphenograph_report.txt",append=TRUE)
            print("Rphenograph report produced")}}
          if (input$ClustAlgo == "FastPG"){
            phenograph.res <- phenograph_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, 
                                               K = input$k, seed = input$set_seed, implementation = "FastPG")
            if (inherits(x = phenograph.res,"list")){
              eventi <- as.character(phenograph.res[[1]])
              numCluster <- as.character(nlevels(phenograph.res[[1]]))
              para <- as.character(input$k)
              modularity <- phenograph.res[[3]]
            text_report <- paste0("FastPG (an enhanced R implementation Phenograph clustering algorithm), run with K = ",
                                  para , ", grouping the dataset of ",eventi , " events in ",
                                  numCluster , " clusters, with ", modularity, " as modularity")
            cat(text_report,file="./tmpdata/FastPG_report.txt",append=TRUE)
            print("FastPG report produced")}}
          if (inherits(x = phenograph.res,"list")){
            dI.clust$mapping <- phenograph.res[[1]]
            dI.clust$codes <- phenograph.res[[2]]
            dI.clust$modularity <- phenograph.res[[3]] #this is introduced with fastPG
            dI.clust$NumCluster <- nlevels(factor(dI.clust$mapping))
            color.clust <<- color_set[1:dI.clust$NumCluster]}else{
              dI.clust$mapping <- phenograph.res}}
        
        if (input$ClustAlgo == "flowClust"){
          flowClust.res <- flowClust_clust(flow_Set = fSmeta(), selected_markers = meta.markers()$selected, 
                                           K = input$k, seed = input$set_seed)
          if (inherits(x = flowClust.res, "list")){
            dI.clust$mapping <- flowClust.res[[1]]
            dI.clust$codes <- flowClust.res[[2]]
            dI.clust$NumCluster <- nlevels(factor(dI.clust$mapping))
            color.clust <<- color_set[1:dI.clust$NumCluster]
            print("flowClust report produced")}else{
              dI.clust$mapping <- flowClust.res}}
        
        output$console_output_clust <- renderPrint("clustering process ends")
        print("clustering process ends")
        
        if (input$clusteringPerf){
          print("start clusters evaluation process")
          cell_clustering <- dI.clust$mapping #[1:Number of events]
          
          flow_Frame <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_sample")
          
          if ((input$ClusteringPdistance)=="minkowski"){
            minkowski <- sum(meta.markers()$selected)}
          else{minkowski = NULL}
          
          clust_color <- color.clust
          plot_res <- cluster_eva(fF = flow_Frame, clustering = cell_clustering, color = clust_color, distance = input$ClusteringPdistance, 
                                  power = minkowski, setseed = input$set_seed, type = "clust")
          dI.clust$sil_plot <- plot_res[[1]] 
          dI.clust$sil_summary <- plot_res[[2]]
          output$console_output_clust <- renderPrint("clusters evaluation process end")
          print("clusters evaluation process end")}
        
        if (input$ClustAlgo == "flowSOM"){
          print("producing flowSOM report")
          flowSOM_plot(mst = dI.clust$mst , markers_df = meta.markers())
          print("flowSOM report produced")}
          
        ### Label I ----
        clean_graph()
        #if((class(dI.clust$codes)=="matrix")&&(class(fSmeta())=="flowSet")&&(dI.flow$labelI))
        #Warning: Error in &&: 'length = 2' in coercion to 'logical(1)' 
        # from release 4.2.0 on:  
        # * Calling && or || with either argument of length greater than one now gives a warning (which it is 
        # intended will become an error).
        # * Calling if() or while() with a condition of length greater than one gives an error rather than a warning.  Consequently,
        #  environment variable _R_CHECK_LENGTH_1_CONDITION_ no longer has any effect   
        if((inherits(dI.clust$codes, "matrix"))&&(inherits(x = fSmeta(),"flowSet"))&&(dI.flow$labelI))  {
          #if(nrow(ph)>0){
            {if (input$signature_finding_method=="Densities"){
              dI.labelClust$cs.res <- try(expr = cluster_signature(fS = fSmeta(), selected.marker = mk$selected, clustId = dI.clust$mapping,
                                          color_marker = color.marker, tipo = "cluster", central = input$HeatMapCentral,
                                          min_quantile = input$minQuantile, max_quantile = input$maxQuantile,
                                          pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold), silent = FALSE)
                if (inherits(x = dI.metaClust_val$cs.res,"list")){options(warn = 0)} else errore <- dI.metaClust_val$cs.res}
              
            if((inherits(x = dI.clust$mapping, "integer"))&&(input$minQuantile<input$maxQuantile)&&(input$pos_threshold>=input$neg_threshold)) {
              print("start labelling process")
              mapping <- dI.clust$mapping 
              if(nrow(ph)>0){
                if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
                label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = mapping, 
                                  phenoquery = ph, central = input$HeatMapCentral, 
                                  min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                  pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                  cluster_signature = dI.labelClust$cs.res[[2]])} else
                {label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = mapping, 
                                  phenoquery = ph, central = input$HeatMapCentral, 
                                  min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                  pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold)}
              
              dI.labelClust$label.res <- label.res}
              
              mapping <- dI.clust$mapping 
              color_clusters <- color.clust
              
              if (is.data.frame(dI.labelClust$label.res)){
                print("start producing heatmap")
                if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
                dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label", 
                                                clustering = mapping, selected_markers = meta.markers()$selected, 
                                                cluster_labelling = dI.labelClust$label.res, 
                                                pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                color_clusters = color_clusters, color_label = color.pheno, 
                                                cluster.signature = dI.labelClust$cs.res[[1]])}
                else{
                  dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label", 
                                                clustering = mapping, selected_markers = meta.markers()$selected, 
                                                cluster_labelling = dI.labelClust$label.res, 
                                                pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                color_clusters = color_clusters, color_label = color.pheno)}
                }
              else
              {
                if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
                dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label", 
                                                  clustering = mapping, selected_markers = meta.markers()$selected, 
                                                  pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                  central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                  pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                  dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                  color_clusters = color_clusters, color_label = color.pheno, 
                                                  cluster_labelling = dI.labelClust$cs.res[[1]])}
                else
                  {dI.labelClust$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label", 
                                                  clustering = mapping, selected_markers = meta.markers()$selected, 
                                                  pheno_table = ph, Nclust = length(dI.labelClust$label.res$new_cluster), 
                                                  central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                  pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                  dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                  color_clusters = color_clusters, color_label = color.pheno)}}
              
           #   dI.flow$labelI <- FALSE}}
              dI.flow$labelI <- FALSE}}
          print("labelling process ends")
          dI.flow$labelI <- FALSE}
        
        ### Meta-clustering -------------------------------------------------------------
        
        clean_graph()
        #if((class(dI.clust$codes)=="matrix")&&(class(dI.clust$mapping)=="integer")&&(dI.flow$meta_clust)){
        if (inherits(x = dI.clust$codes, "matrix")){
          print("start meta-clustering process")
          dI.metaClust(NULL) #this is to reset the metaclustering results
          
          #if ((input$MetaClustAlgo == "ConsensusClusterPlus")&&(dI.clust$NumCluster>input$MetaClustNum)){
          if (dI.clust$NumCluster>input$MetaClustNum){
            mc.res <- meta_clustering(codes = dI.clust$codes, maxK = input$MetaClustNum, reps = 100, 
                                      pItem = 1.0, clusterAlg = input$MetaClustCriteria, distance = input$MetaClustDist, 
                                      seed = input$set_seed, savedir = "./tmpdata", save = FALSE)
            dI.metaClust(mc.res) #Notice, depending of the dataset, the meta_clustering can produce a number of meta_clusters which is less than input$MetaClustNum 
            
            if(inherits(x = dI.metaClust(),"list")) {
              code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
              cell_clustering <- code_clustering[dI.clust$mapping] #[1:Number of events]
              
              if (input$signature_finding_method=="Densities"){
                dI.metaClust_val$cs.res <- try(expr = cluster_signature(fS = fSmeta(), selected.marker = mk$selected, clustId = code_clustering,
                                                  color_marker = color.marker, tipo = "meta", central = input$HeatMapCentral,
                                                  min_quantile = input$minQuantile, max_quantile = input$maxQuantile,
                                                  pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold), silent = FALSE)
                if (inherits(x = dI.metaClust_val$cs.res,"list")){options(warn = 0)} else errore <- dI.metaClust_val$cs.res}
              
              if (is.data.frame(dI.labelClust$label.res)){ 
                if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
                dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering, 
                                                        cluster_labelling = NULL, pheno_table = ph, Nclust = input$MetaClustNum, 
                                                        central = input$HeatMapCentral, selected_markers = meta.markers()$selected, 
                                                        min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                        pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                        dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                        color_clusters = color.metaclust, 
                                                        cluster.signature = dI.metaClust_val$cs.res[[1]])}
                else{
                  dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering, 
                                                        cluster_labelling = NULL, pheno_table = ph, Nclust = input$MetaClustNum, 
                                                        central = input$HeatMapCentral, selected_markers = meta.markers()$selected, 
                                                        min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                        pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                        dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                        color_clusters = color.metaclust)}} 
              else
                {if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
                  dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering,
                                                        cluster_labelling = NULL, pheno_table = NULL, Nclust = input$MetaClustNum, 
                                                        central = input$HeatMapCentral, selected_markers = meta.markers()$selected,
                                                        min_quantile = input$minQuantile, max_quantile = input$maxQuantile,
                                                        pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                        dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                        color_clusters = color.metaclust, 
                                                        cluster.signature = dI.metaClust_val$cs.res[[1]])}
                  else{
                    dI.metaClust_val$hm_cluster <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "meta", clustering = cell_clustering,
                                                        cluster_labelling = NULL, pheno_table = NULL, Nclust = input$MetaClustNum, 
                                                        central = input$HeatMapCentral, selected_markers = meta.markers()$selected,
                                                        min_quantile = input$minQuantile, max_quantile = input$maxQuantile,
                                                        pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                        dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                        color_clusters = color.metaclust)}
                }}
            if (input$clusteringPerf){
              print("start clusters evaluation process")
              
              flow_Frame <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_metasample")
              
              if ((input$ClusteringPdistance)=="minkowski"){
                minkowski <- sum(meta.markers()$selected)}
              else{minkowski = NULL}
              
              clust_color <- color.metaclust
              plot_res <- cluster_eva(fF = flow_Frame, clustering = cell_clustering, color = color.metaclust, distance = input$ClusteringPdistance, 
                                      power = minkowski, setseed = input$set_seed, type = "meta")
              dI.metaClust_val$sil_plot <- plot_res[[1]] 
              dI.metaClust_val$sil_summary <- plot_res[[2]]
              output$console_output_clust <- renderPrint("clusters evaluation process end")
              print("clusters evaluation process end")}} #if(class(dI.metaClust())=="list")
            
          print("meta-clustering process ends")
          dI.flow$meta_clust <- FALSE}
        
        
       ### Mapping -----------------------------------------------------------------------------
        mapres <- NULL
        if(((is.data.frame(dI.labelClust$label.res))||(inherits(x = dI.metaClust(),"list")))&&(inherits(x = fSmeta(),"flowSet"))){
          print("start mapping process")
          if ((inherits(x = dI.tSNE.eva.init$mapres_init,"list"))&&(!input$two_three)&&(input$mapType=='tSNE')){
            mapres <- dI.tSNE.eva.init$mapres_init}
          if ((inherits(x = dI.tSNE.eva.final$mapres_final,"list"))&&(!input$two_three)&&(input$mapType=='tSNE')){
            mapres <- dI.tSNE.eva.final$mapres_final}
          #if (class(mapres)!="list")
          if (!inherits(x = mapres,"list")){
            mapres <- map_gen(flow_Set = fSmeta(), type = input$mapType, two_three = input$two_three, 
            selected_markers = mk$selected, seed_nr = input$set_seed, limit = input$map_limit, maxCell = input$mapMax,  
            #tSNE parameters
            Rtsne_pca = input$tSNE_PCA, Rtsne_perp = input$perplexity, Rtsne_theta = input$theta,
            Rtsne_iter = input$tSNEiter, Rtsne_eta = input$tSNEeta,
            #UMAP paramters
            UMAP.n_neighbors = input$UMAP.n_neighbors, UMAP.metric = input$UMAP.metric, 
            UMAP.min_dist = input$UMAP.min_dist, UMAP.spread = input$UMAP.spread, UMAP.random_state = input$set_seed)}
          #mapres <- readRDS(file = "mapres.RDS") #remember to uncomment the command b4
          dI.map(mapres)
          
        clean_graph()
          if(length(dI.map())==2){print("mapping process ends")}}
        if(inherits(x = dI.map(),"list")) {
          print("start map printing")
          
          if(inherits(x = dI.metaClust(),"list")){ #case 1
            code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
            meta_clust <- code_clustering
            metaClustNum <- input$MetaClustNum
            color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
            showmap <- map_plot_comp(map_df = dI.map()[[1]], type = input$mapType, two_three = input$two_three, map_inds = dI.map()[[2]],
                                     map_reduced = dI.map()[[3]], metadata = mt, clust = dI.clust$mapping, metaclust = meta_clust,
                                     NMC = metaClustNum, color_clusters = color_metaclust, save = TRUE)
            fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_meta_sample")
            nomi <- colnames(exprs(fF))
            if (!any(grepl("clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)}
            if (!any(grepl("meta_clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "meta_clusterId", dim_vect = cell_clustering)}
            
            jt <- join_tag(frame = fF, meta = mt, type = "meta", tag1col=color.tag1, tag2col=color.tag2, tag3col=color.tag3, tag4col=color.tag4)}
          else{
            if (is.data.frame(dI.labelClust$label.res)){
              meta_clust <- dI.labelClust$label.res$numeric
              metaClustNum <- length(dI.labelClust$label.res$numeric)
            
              color_metaclust <- dI.labelClust$hm_cluster$color_clusters
              showmap <- map_plot_comp(map_df = dI.map()[[1]], type = input$mapType, two_three = input$two_three, map_inds = dI.map()[[2]],
                                     map_reduced = dI.map()[[3]], metadata = mt, clust = dI.clust$mapping, 
                                     NMC = metaClustNum, color_clusters = color_metaclust, save = TRUE)}
              fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_meta_sample")
              nomi <- colnames(exprs(fF))
              if (!any(grepl("clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)}
            jt <- join_tag(frame = fF, meta = mt, type = "clust", tag1col=color.tag1, tag2col=color.tag2, tag3col=color.tag3, tag4col=color.tag4)}
          
          print("map printing ends")}
        
        clean_graph()
        if(inherits(x = dI.map(),"list")) {
          print("start marker map printing")
          
          showmap <- map_plot_comp(map_df = dI.map()[[1]], type = "marker", save = TRUE)
          print("marker map printing ends")}
        
        #################################################################### saveMap  ----
        
        if((inherits(x = fSmeta(), "flowSet"))&&(inherits(x = dI.clust$mapping, "integer"))&&(inherits(x = dI.map(),"list"))){
          print("producing concatenated sample with the map data")
          
          if(inherits(x = fSmeta(),"flowSet")){
           
            if (length(dI()) > 0){
              fFvectorName <- vector(mode = "character", length(fSmeta()))
              for (i in seq_along(1:length(fSmeta()))){
                fFvectorName[[i]] <- basename(fSmeta()[[i]]@description$FILENAME)}
              
              save_fS(in_flow_Set = dI(), handled_flow_Set = fSmeta(), stringhe = fFvectorName)
              
              fFvectorName_comp <- paste0("./tmpdata/complete_",fFvectorName)
              fS.complete <- read.flowSet(files = fFvectorName_comp, truncate_max_range = FALSE)
              fF <- concatenating_fS(flow_Set = fS.complete, stringa = "conc_meta_sample")}
            else
            {fF <- concatenating_fS(flow_Set = fSmeta(), stringa = "conc_meta_sample")}
            
            lista_file <- list.files(path = data_dir, pattern = "^conc_meta_sample*")
            if (!(length(lista_file)==0)){
              lista_file <- paste0("./tmpdata/", lista_file)
              unlink(lista_file, recursive = FALSE)}
            
            nomi <- colnames(exprs(fF))
            if (!any(grepl("clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "clusterId", dim_vect = dI.clust$mapping)}
            
            if (inherits(x = dI.metaClust(),"list")){
              code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass 
              cell_clustering <- as.integer(code_clustering[dI.clust$mapping])
              if (!any(grepl("meta_clusterId", nomi, fixed = TRUE))){fF <- add_dim(flow_Frame = fF, dim_name = "meta_clusterId", dim_vect = cell_clustering)}}
            
            if (inherits(x = dI.metaClust(),"list")){
            jt <- join_tag(frame = fF, meta = mt, type = "meta",
                           tag1col=color.tag1, tag2col=color.tag2, tag3col=color.tag3, tag4col=color.tag4)} 
            else {jt <- join_tag(frame = fF, meta = mt, type = "clust",
                                 tag1col=color.tag1, tag2col=color.tag2, tag3col=color.tag3, tag4col=color.tag4)}
            
            meta <- mt
            nomi_tag <- colnames(meta)
            nomi_tag <- nomi_tag[-(1:2)]
            for (i in 1:4){mt[,i+2] <- (as.numeric(factor(meta[,i+2])))}
            if (!any(grepl("meta_clusterId", nomi, fixed = TRUE))){
              df <- as.data.frame(exprs(fF)[,c("SampleID", "clusterId", "meta_clusterId")])} else {df <- as.data.frame(exprs(fF)[,c("SampleID", "clusterId")])}
            
            meta$SampleID <- as.numeric(1:nrow(meta))
            rj <- right_join(meta,df)
            rj <- rj[,-(1:2)]
            rj$date <- NULL
            
            rj$SampleID <- factor(rj$SampleID)
            rj$clusterId <- factor(rj$clusterId)
            if (!any(grepl("meta_clusterId", nomi, fixed = TRUE))){rj$meta_clusterId <- factor(rj$meta_clusterId)}
            
            
            if (uni.tag1 > 1){
              tag <-  nomi_tag[1]
              tagdata <- as.integer(factor(rj[,1]))
              fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
            if (uni.tag2 > 1){
              tag <-  nomi_tag[2]
              tagdata <- as.integer(factor(rj[,2]))
              fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
            if (uni.tag3 > 1){
              tag <-  nomi_tag[3]
              tagdata <- as.integer(factor(rj[,3]))
              fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
            if (uni.tag4 > 1){
              tag <-  nomi_tag[4]
              tagdata <- as.integer(factor(rj[,4]))
              fF <- add_dim(flow_Frame = fF, dim_name = tag, dim_vect = tagdata)}
            
            fF <- add_dim(flow_Frame = fF, dim_name = "map_X", dim_vect = dI.map()[[1]]$map_x)
            fF <- add_dim(flow_Frame = fF, dim_name = "map_Y", dim_vect = dI.map()[[1]]$map_y)
            if(input$two_three==TRUE)
            {fF <- add_dim(flow_Frame = fF, dim_name = "map_Z", dim_vect = dI.map()[[1]]$map_z)}
            
              dI.map_save$fFmap <- fF
              print("concatenated sample with the map data produced")}
          
          write.FCS(dI.map_save$fFmap, filename = "./tmpdata/map_cluster.fcs")
          #write.csv(x = exprs(dI.map_save$fFmap), file = "./tmpdata/map_cluster.csv")
          print("concatenated sample with the map data produced")}
        
        ### Label II ----
        clean_graph()
        if((is.data.frame(dI.labelClust$label.res))&&(inherits(x = fSmeta(),"flowSet"))&&(dI.flow$labelII)){ 
          # This hold only when you pass by the Label I procedure
          print("start labelling II process")
          if(inherits(dI.metaClust(),"list")){  # case 3: Labelling + Metaclustering
            
            code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
            cell_clustering <- code_clustering[dI.clust$mapping]
            metaClustNum <- input$MetaClustNum #[1:Number of events]
            color_label <- color.pheno
            if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
              label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                      phenoquery = ph, central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                      pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                      cluster_signature = dI.metaClust_val$cs.res[[1]])} else
              {label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                      phenoquery = ph, central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                      max_quantile = input$maxQuantile, pos_threshold = input$pos_threshold, 
                                      neg_threshold = input$neg_threshold)}
            dI.label(label.res)
            print("labelling II process ends")
            if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
              dI.label_val$hm_label <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "label", 
                                          clustering = cell_clustering, cluster_labelling = dI.label(), 
                                          Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                          central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                          pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                          dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                          color_clusters = color.metaclust, color_label = color_label, work = 3, 
                                          cluster.signature = dI.metaClust_val$cs.res[[1]])}
            else{dI.label_val$hm_label <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "label", 
                                            clustering = cell_clustering, cluster_labelling = dI.label(), 
                                            Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                            central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                            pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                            dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                            color_clusters = color.metaclust, color_label = color_label, work = 3)}
            
            mm <- match(cell_clustering, dI.label()$original_cluster)
            cell_clustering2 <- dI.label()$new_cluster[mm]
            color_label <- dI.label_val$hm_label$color_label
            
            if ((nlevels(factor(cell_clustering2))) > 1){
              if ((inherits(x = dI.metaClust_val$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
                dI.label_val$hm_pheno <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                              clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                              Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                              central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                              pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                              dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                              color_clusters = color_label, work = 3, 
                                              cluster_labelling = dI.metaClust_val$cs.res[[1]])}
              else
              {dI.label_val$hm_pheno <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                              clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                              Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                              central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                              pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                              dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                                        color_clusters = color_label, work = 3)}}
            else{dI.label_val$hm_pheno <- "The second plot cannot be built: there must be at least two different labelled phenotypes"}}}
          
         #the else statement is removed here: the workflow 2) is performed at the end
        #### Pheno maps -------
        clean_graph()
        if(((is.data.frame(dI.labelClust$label.res))||(inherits(x = dI.metaClust(),"list")))&&(inherits(x = dI.map(),"list"))&&(dI.flow$pheno_mapsI))
        { 
          if ((inherits(x = dI.map(),"list"))&&(inherits(x = dI.clust$mapping,"integer")))
          {
            print("start labelled_p map printing")
            
            if((is.list(dI.labelClust$hm_cluster))&&(is.list(dI.label_val$hm_pheno))){
              if(inherits(x = dI.metaClust(), "list")){
                # case 3)
                code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
                cell_clustering <- code_clustering[dI.clust$mapping]
                metaClustNum <- input$MetaClustNum
                mm <- match(cell_clustering, dI.label()$original_cluster)
                cell_label <- dI.label()$new_cluster[mm]
                map_cell <- factor(cell_label[dI.map()[[2]]])
                
                color_pheno <- dI.label_val$hm_pheno$color_clusters
                showMaps <- plot_labelling(flow_Set = fSmeta(), phenoclust = dI.label(), clustering = dI.clust$mapping, 
                                           m_clustering = cell_clustering , map_res = dI.map()[[1]], map_cell = map_cell, 
                                           selected_markers = mk$selected, metadata = mt, color_clusters = color_pheno, save = TRUE, work = 3)}}
            else
            {if(!(is.data.frame(dI.labelClust$label.res))&&(inherits(x = dI.metaClust(),"list"))){
              # case 1)
              code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
              cell_clustering <- code_clustering[dI.clust$mapping]
              metaClustNum <- input$MetaClustNum
              original_cluster = seq(1, metaClustNum)
              new_cluster <- original_cluster
              for (i in 1:metaClustNum){new_cluster[i] <- as.character(original_cluster[i])}
              
              # metaClustNum <- input$MetaClustNum
              pheno_meta <- data.frame("original_cluster" = original_cluster, "new_cluster" = new_cluster)
              mm <- match(cell_clustering, pheno_meta$original_cluster)
              cell_label <- pheno_meta$new_cluster[mm]
              map_cell <- factor(cell_label[dI.map()[[2]]], levels = 1:metaClustNum) # add cell_label to map_df
              color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
              showMaps <- plot_labelling(flow_Set = fSmeta(), phenoclust = pheno_meta, clustering = dI.clust$mapping, 
                                         m_clustering = cell_clustering , map_res = dI.map()[[1]], map_cell = map_cell, 
                                         selected_markers = mk$selected, metadata = mt, color_clusters = color_metaclust, save = TRUE, work = 1)}}
           
              print("labelled map_p printing process ends")}
          
          dI.flow$pheno_mapsI <- FALSE}
        
        #### Quantity plots -------------------------------------------------------------------------
        
        clean_graph()
        if((inherits(x = dI.clust$mapping, "integer"))&&(dI.flow$quantity))  {
          print("start printing quantity plot")
          showBars <- NULL
          
          if(is.list(dI.labelClust$hm_cluster)&&((is.list(dI.label_val$hm_pheno)))){
            if(inherits(x = dI.metaClust(),"list")){
              
              #case 3)
              code_clustering <- as.integer(dI.metaClust()[[input$MetaClustNum]]$consensusClass) #[1:Number of clusters]
              cell_clustering <- code_clustering[dI.clust$mapping]
              metaClustNum <- input$MetaClustNum #[1:Number of events]
              
              mm <- match(cell_clustering, dI.label()$original_cluster)
              cell_label <- dI.label()$new_cluster[mm]
              
              color_pheno <- dI.label_val$hm_pheno$color_clusters #Error in : $ operator is invalid for atomic vectors
              showBars <- plot_labellingbar(flow_Set = fSmeta(), flag = FALSE, cell_label = cell_label, phenoclust = dI.label(),
                                            NMC = metaClustNum, metadata =  mt, color_pheno = color_pheno, save = TRUE, work = 3) 
              
              showMedian <- plot_median_expr(flow_Set = fSmeta(), phenoclust =  dI.label(), cell_clustering = cell_clustering, 
                                             metadata = mt, NMC = metaClustNum,
                                             color_tag1 = color.tag1, color_tag2 = color.tag2, color_tag3 = color.tag3, color_tag4 = color.tag4, 
                                             save = TRUE, work = 3)}}
          # else is computed afterwords: the second worflow at the end
          if(!(is.data.frame(dI.labelClust$label.res))&&(inherits(x = dI.metaClust(), "list"))){
            #case 1)
            code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
            cell_clustering <- code_clustering[dI.clust$mapping]
            #metaClustNum <- input$MetaClustNum 
            metaClustNum <- length(levels(factor(cell_clustering))) 
            original_cluster = seq(1, metaClustNum)
            new_cluster <- original_cluster
            for (i in 1:metaClustNum){new_cluster[i] <- as.character(original_cluster[i])}
            pheno_meta <- data.frame("original_cluster" = original_cluster, "new_cluster" = new_cluster)
            
            mm <- match(cell_clustering, pheno_meta$original_cluster)
            cell_label <- factor(pheno_meta$new_cluster[mm], levels = 1:metaClustNum)
            
            color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
            showBars <- plot_labellingbar(flow_Set = fSmeta(), flag = TRUE, cell_label = cell_label, phenoclust = pheno_meta,
                                          NMC = metaClustNum, metadata =  mt, color_pheno = color_metaclust, save = TRUE, work = 1) 
            
            showMedian <- plot_median_expr(flow_Set = fSmeta(), phenoclust =  pheno_meta, cell_clustering = cell_clustering, 
                                           metadata = mt, NMC = metaClustNum, 
                                           color_tag1 = color.tag1, color_tag2 = color.tag2, color_tag3 = color.tag3, color_tag4 = color.tag4, 
                                           save = TRUE, work = 1)}
          print("quantity plot printing process ends")
          dI.flow$quantity <- FALSE} # if(class(dI.clust$mapping)=="numeric")
        
        #### Time-step plots ----
        
        clean_graph()
        if((inherits(x = dI.clust$mapping, "integer"))&&(uni.tag4 > 1)&&(dI.flow$time_step)){
          print("start printing plot visualization")
          
          if(is.data.frame(dI.label())){
            if(inherits(x = dI.metaClust(), "list")){
              #case 3)
              
              code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
              cell_clustering <- code_clustering[dI.clust$mapping]
              metaClustNum <- input$MetaClustNum
              mm <- match(cell_clustering, dI.label()$original_cluster)
              cell_label <- dI.label()$new_cluster[mm]
              color_pheno <- dI.label_val$hm_pheno$color_clusters
              if (input$time_step == TRUE){
                dI.quant_val_p$showStream <- plot_stream(flow_Set = fSmeta(), flag = FALSE, phenoclust = dI.label(), cell_label = cell_label, 
                                                       metadata = mt, NMC = metaClustNum, color_pheno = color_pheno)}}} 
          else{  
            if(!(is.data.frame(dI.labelClust$label.res))&&(inherits(x = dI.metaClust(),"list"))){
              #case 1)
              
              code_clustering <- dI.metaClust()[[input$MetaClustNum]]$consensusClass #[1:Number of clusters]
              cell_clustering <- code_clustering[dI.clust$mapping]
              #metaClustNum <- input$MetaClustNum
              metaClustNum <- length(levels(factor(cell_clustering))) 
              original_cluster = seq(1, metaClustNum)
              new_cluster <- original_cluster
              for (i in 1:metaClustNum){new_cluster[i] <- as.character(original_cluster[i])}
              #new_cluster <- as.character(original_cluster)
              pheno_meta <- data.frame("original_cluster" = original_cluster, "new_cluster" = new_cluster)
              mm <- match(cell_clustering, pheno_meta$original_cluster)
              cell_label <- factor(pheno_meta$new_cluster[mm])
              color_metaclust <- dI.metaClust_val$hm_cluster$color_clusters
              if (input$time_step == TRUE){
                dI.quant_val_p$showStream <- plot_stream(flow_Set = fSmeta(), flag = TRUE, phenoclust = pheno_meta, cell_label = cell_label, 
                                                       metadata = mt, NMC = metaClustNum, color_pheno = color_metaclust)}}}
          
          dI.flow$time_step <- FALSE}
        
        ################################################################ 2. Workflow #2 -----
       
        if((is.data.frame(dI.labelClust$label.res))&&(inherits(x = fSmeta(),"flowSet"))&&(dI.flow$labelII)){ 
          # This hold only when you pass by the Label I procedure
          print("start labelling process for workflow 2")
          dI.metaClust(NULL)
          
          cell_clustering <- dI.clust$mapping 
          label.res <- dI.labelClust$label.res
          color_label <- color.pheno
          metaClustNum <- length(dI.labelClust$label.res$new_cluster)
          color_clusters <- dI.labelClust$hm_cluster$color_clusters
          if ((inherits(x = dI.labelClust$cs.res, "list"))&&(input$signature_finding_method=="Densities")){
          dI.label_val$hm_label_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label_plus", 
                                              clustering = cell_clustering, selected_markers = meta.markers()$selected, 
                                              cluster_labelling = label.res, 
                                              pheno_table = ph, Nclust = metaClustNum, 
                                              central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                              pos_threshold = dI.labelClust$cs.res[[1]], neg_threshold = dI.labelClust$cs.res[[1]], 
                                              dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                              color_clusters = color_clusters, color_label = color_label, work = 2, 
                                              cluster.signature = dI.labelClust$cs.res[[1]])}
          else{dI.label_val$hm_label_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "full_label_plus", 
                                                clustering = cell_clustering, selected_markers = meta.markers()$selected, 
                                                cluster_labelling = label.res, 
                                                pheno_table = ph, Nclust = metaClustNum, 
                                                central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold, 
                                                dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                color_clusters = color_clusters, color_label = color_label, work = 2)}
          
          if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){  
          label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                      phenoquery = ph, central = input$HeatMapCentral, 
                                      min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                      pos_threshold = dI.labelClust$cs.res[[1]], neg_threshold = dI.labelClust$cs.res[[1]], 
                                      cluster_signature = dI.labelClust$cs.res[[2]])} else                        
          {label.res <- phenocluster(flow_Set = fSmeta(), selected_markers = mk$selected, clustering = cell_clustering, 
                                      phenoquery = ph, central = input$HeatMapCentral, min_quantile = input$minQuantile, 
                                      max_quantile = input$maxQuantile, pos_threshold = input$pos_threshold, 
                                      neg_threshold = input$neg_threshold)}
          dI.label(label.res)
          mm <- match(cell_clustering, dI.label()$original_cluster)
          cell_clustering2 <- dI.label()$new_cluster[mm]
          color_label <- dI.label_val$hm_label_plus$color_label
            
          if ((nlevels(factor(cell_clustering2))) > 1){
            if ((inherits(x = dI.labelClust$cs.res,"list"))&&(input$signature_finding_method=="Densities")){
              dI.label_val$hm_pheno_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                                clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                                Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                                central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                pos_threshold = dI.labelClust$cs.res[[1]], neg_threshold = dI.labelClust$cs.res[[1]],
                                                dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                color_clusters = color_label, work = 2, 
                                                cluster.labelling = dI.labelClust$cs.res[[1]])}
            else
            {dI.label_val$hm_pheno_plus <- plot_clustering_heatmap_wrapper(flow_Set = fSmeta(), type_HM = "pheno", 
                                                clustering = cell_clustering2,  #cluster_labelling = dI.label(), 
                                                Nclust = metaClustNum, selected_markers = meta.markers()$selected, pheno_table = ph,
                                                central = input$HeatMapCentral, min_quantile = input$minQuantile, max_quantile = input$maxQuantile, 
                                                pos_threshold = input$pos_threshold, neg_threshold = input$neg_threshold,
                                                dist_method = input$HeatMapMethod, hclust_method = input$HeatmapHCCriteria, 
                                                color_clusters = color_label, work = 2)}}
          else{dI.label_val$hm_pheno_plus <- "The second plot cannot be built: there must be at least two different labelled phenotypes"}
          
          print("labelling process for workflow 2 ends")}
        
        if((is.list(dI.labelClust$hm_cluster))&&((is.list(dI.label_val$hm_pheno_plus)))){
          print("start labelled_p map printing for workflow 2 ")
          cell_clustering <- dI.clust$mapping 
          mm <- match(cell_clustering, dI.label()$original_cluster)
          cell_clustering2 <- dI.label()$new_cluster[mm]
          color_label <- dI.label_val$hm_label_plus$color_label
            
          cell_label <- dI.label()$new_cluster[mm]
          map_cell <- factor(cell_label[dI.map()[[2]]])
          
          showMaps <- plot_labelling(flow_Set = fSmeta(), phenoclust = dI.label(), clustering = dI.clust$mapping, 
                                     m_clustering = cell_clustering2 , map_res = dI.map()[[1]], map_cell = map_cell, 
                                     selected_markers = mk$selected, metadata = mt, color_clusters = color_label, save = TRUE, work = 2)
          print("labelled_p map printing for workflow 2 ends")}
        
        clean_graph()
        if((inherits(x = dI.clust$mapping, "integer"))&&(dI.flow$quantity)&&(is.list(dI.label_val$hm_label_plus))){
          
          showBars <- NULL
          print("printing quantity plot for workflow2")
          if(is.list(dI.labelClust$hm_cluster)){
            if ((is.list(dI.label_val$hm_label_plus))){
                cell_clustering <- dI.clust$mapping 
                mm <- match(cell_clustering, dI.label()$original_cluster)
                cell_clustering <- dI.label()$new_cluster[mm]
                cell_label <- dI.label()$new_cluster[mm]
                
                metaClustNum <- nlevels(factor(dI.label()$new_cluster))
                color_label <- dI.label_val$hm_label_plus$color_label
                cell_clustering <- dI.clust$mapping 
                
                showBars <- plot_labellingbar(flow_Set = fSmeta(), flag = FALSE, cell_label = cell_label, phenoclust = dI.label(),
                                              NMC = metaClustNum, metadata =  mt, color_pheno = color_label, save = TRUE, work = 2) 
                
                showMedian <- plot_median_expr(flow_Set = fSmeta(), phenoclust =  dI.label(), cell_clustering = cell_clustering, 
                                               metadata = mt, NMC = metaClustNum, 
                                               color_tag1 = color.tag1, color_tag2 = color.tag2, color_tag3 = color.tag3, color_tag4 = color.tag4,
                                               save = TRUE, work = 2)}}
          
          print("quantity print process for workflow 2 ends")
          dI.flow$quantity <- FALSE}
        
        ### ??? is here to add the workflow2 with or without the pheno table?
        
        output$console_output_clust <- renderPrint("analysis finished")
        if (Sys.info()["sysname"] == "Windows") play(x = gong)
        output$runAuto_panelStatus <- reactive({input$runAuto_panelStatus=="show"})
        outputOptions(output, "runAuto_panelStatus", suspendWhenHidden = FALSE)}
    
        output$downloadRunAuto <- downloadHandler(
            filename = function() {
              paste("Analysis_material", "zip", sep=".")
          },
        
          content <- function(file) {
            
            del_tmpdata()
            file_path <- file
            setwd(data_dir)
          
            Zip_Files <- list.files(path = data_dir, pattern = "*.*")
            zip(zipfile = file_path, files = Zip_Files)
            setwd(app_dir)
            },
          contentType = "application/zip")
      })# end observeEent
    
} #end of server

#session$onSessionEnded(stopApp)