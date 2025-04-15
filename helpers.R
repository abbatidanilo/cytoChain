############################################################################# global var -----
#https://www.shinyapps.io/admin/#/login?redirect=%2Fapplications%2Fall
#https://abbatidanilo.shinyapps.io/cytoChain/

RENDER_DATA <- 35000 #with 50000 the plotting of the labelled map is very slow
MIN_SAMPLE_LENGTH <- 30 #with a sample with less than MIN_SAMPLE_LENGTH the downsampling cannot be performed
MIN_SAMPLE_LENGTH_AFTER_DOWNSAMPLING <- 20 
MAX_EVA <- 20000

MIN_SAMPLE_LENGTH_FOR_CLUSTER_SIGNATURE <- 50 
MIN_PEAK_DISTANCE_TO_SHOW_PLOT_REPORT <- 0.02
MIN_EVENTS_TO_TRUST_PLOT_REPORT <- 100

tabDim <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("dimension name", "marker description"))
app_dir <- new(Class = "character")
data_dir <- new(Class = "character")
temp_dir <- new(Class = "character")
mt <- new(Class = "data.frame")
mk <- new(Class = "data.frame")
ph <- new(Class = "data.frame")

uni.tag1 <- 0
uni.tag2 <- 0
uni.tag3 <- 0
uni.tag4 <- 0
uni.marker <- 0
uni.pheno <- 0

color.sample <- new(Class = "character")
color.tag1 <- new(Class = "character")
color.tag2 <- new(Class = "character")
color.tag3 <- new(Class = "character")
color.tag4 <- new(Class = "character")

color_set <- c('#7200dd', '#E8D2CB', '#ffab00', '#aeea00', '#00c853', '#f50057', '#2962ff', '#00b8d4', '#c51162', '#3e2723', '#ff6d00', '#ffd600', 
               '#64ffda', '#00bfa5', '#aa00ff', '#C3CCFF', '#d50000', '#0091ea', '#6D544F', '#FFCA83', '#ffea00', '#76ff03', '#1de9b6', '#00b0ff', 
               '#717ec5', '#F4E5F6', '#fe4065', '#ffc400', '#ff3d00', '#c6ff00', '#00e676', '#00e5ff', '#2979ff', '#6200ea', '#c6ff00', '#570540', 
               '#E7E35F', '#FFEBD2', '#ffff00', '#b2ff59', '#31680E', '#448aff', '#b388ff', '#ff4081', '#37474f', '#ff6e40', '#ffd740', '#eeff41', 
               '#69f0ae', '#85D9FF', '#536dfe', '#e040fb', '#ff5252', '#6d4c41', '#b652ff', '#ffd180', '#ffff8d', '#ccff90', '#a7ffeb', '#82b1ff',
               '#311b92', '#ff80ab', '#455a64', '#ff9e80', '#ffe57f', '#f4ff81', '#b9f6ca', '#80d8ff', '#8c9eff', '#ea80fc', '#ff8a80', '#dea49f',
               '#bf360c', '#ff6f00', '#827717', '#1b5e20', '#18ffff', '#0d47a1', '#4527a0', '#880e4f', '#546e7a', '#e65100', '#f9a825', '#33691e',
               '#004d40', '#01579b', '#1a237e', '#4a148c', '#b71c1c', '#a7b71c', '#ef6c00', '#fbc02d', '#558b2f', '#00695c', '#0277bd', '#283593',
               '#6a1b9a', '#c62828', '#9e9e9e', '#ff8f00', '#9e9d24', '#2e7d32', '#84ffff', '#1565c0', '#512da8', '#ad1457', '#607d8b', '#795548', 
               '#ffa000', '#afb42b', '#388e3c', '#006064', '#1976d2', '#5e35b1', '#c2185b', '#78909c', '#8d6e63', '#fdd835', '#689f38', '#00796b', 
               '#0288d1', '#303f9f', '#7b1fa2', '#d32f2f', '#bdbdbd', '#d84315', '#ffeb3b', '#7cb342', '#00897b', '#039be5', '#3949ab', '#8e24aa', 
               '#e53935', '#e0e0e0', '#e64a19', '#c0ca33', '#43a047', '#00838f', '#1e88e5', '#673ab7', '#d81b60', '#90a4ae', '#a1887f', '#f57c00',
               '#cddc39', '#4caf50', '#0097a7', '#2196f3', '#7e57c2', '#e91e63', '#b0bec5', '#bcaaa4', '#fb8c00', '#8bc34a', '#009688', '#03a9f4',
               '#3f51b5', '#9c27b0', '#f44336', '#eeeeee', '#f4511e', '#ffb300', '#9ccc65', '#26a69a', '#29b6f6', '#5c6bc0', '#ab47bc', '#ef5350',
               '#f5f5f5', '#ff5722', '#ffc107', '#66bb6a', '#00acc1', '#42a5f5', '#9575cd', '#ec407a', '#cfd8dc', '#d7ccc8', '#ff9800', '#ffee58',
               '#81c784', '#00bcd4', '#64b5f6', '#b39ddb', '#f06292', '#eceff1', '#efebe9', '#ffa726', '#fff176', '#4db6ac', '#4fc3f7', '#7986cb',
               '#ba68c8', '#e57373', '#73e575', '#ff7043', '#ffca28', '#d4e157', '#e6ee9c', '#a5d6a7', '#4dd0e1', '#90caf9', '#d1c4e9', '#f48fb1',
               '#ff8a65', '#ffd54f', '#dce775', '#c5e1a5', '#80cbc4', '#81d4fa', '#9fa8da', '#ce93d8', '#ef9a9a', '#ffb74d', '#fff59d', '#aed581',
               '#ffecb3', '#f0f4c3', '#c8e6c9', '#80deea', '#bbdefb', '#ede7f6', '#f8bbd0', '#ffab91', '#ffe082', '#fffde7', '#dcedc8', '#b2dfdb',
               '#b3e5fc', '#c5cae9', '#e1bee7', '#ffcdd2', '#ffcc80', '#fff9c4', '#1965b0', '#dc050c', '#fff3e0', '#f9fbe7', '#e8f5e9', '#b2ebf2',
               '#e3f2fd', '#f3e5f5', '#ffebee', '#fb8072', '#e0f7fa', '#fff8e1', '#f1f8e9', '#e0f2f1', '#e1f5fe', '#e8eaf6', '#fce4ec', '#ffe0b2',
               '#7bafde', '#d5de7b', '#7570b3', '#a6761d', '#55a1b1', '#33a02c', '#e7298a', '#ff7f00', '#882e72', '#999999', '#beaed4', '#e6ab02',
               '#8dd3c7', '#b2df8a', '#e78ac3', '#fdb462', '#b17ba6',
               '#7200dd', '#dd2c00', '#ffab00', '#aeea00', '#00c853', '#6200ea', '#2962ff', '#00b8d4', '#c51162', '#3e2723', '#ff6d00', '#ffd600',
               '#64dd17', '#00bfa5', '#aa00ff', '#304ffe', '#d50000', '#0091ea', '#4e342e', '#ff9100', '#ffea00', '#76ff03', '#1de9b6', '#00b0ff',
               '#717ec5', '#d500f9', '#fe4065', '#ffc400', '#ff3d00', '#c6ff00', '#00e676', '#00e5ff', '#2979ff', '#f50057', '#7c4dff', '#570540',
               '#5d4037', '#ffab40', '#ffff00', '#b2ff59', '#64ffda', '#448aff', '#b388ff', '#ff4081', '#37474f', '#ff6e40', '#ffd740', '#eeff41',
               '#69f0ae', '#40c4ff', '#536dfe', '#e040fb', '#ff5252', '#b652ff', '#6d4c41', '#ffd180', '#ffff8d', '#ccff90', '#a7ffeb', '#82b1ff',
               '#311b92', '#ff80ab', '#455a64', '#ff9e80', '#ffe57f', '#f4ff81', '#b9f6ca', '#80d8ff', '#8c9eff', '#ea80fc', '#ff8a80', '#dea49f',
               '#bf360c', '#ff6f00', '#827717', '#1b5e20', '#18ffff', '#0d47a1', '#4527a0', '#880e4f', '#546e7a', '#e65100', '#f9a825', '#33691e',
               '#004d40', '#01579b', '#1a237e', '#4a148c', '#b71c1c', '#a7b71c', '#ef6c00', '#fbc02d', '#558b2f', '#00695c', '#0277bd', '#283593',
               '#6a1b9a', '#c62828', '#9e9e9e', '#ff8f00', '#9e9d24', '#2e7d32', '#84ffff', '#1565c0', '#512da8', '#ad1457', '#607d8b', '#795548', 
               '#ffa000', '#afb42b', '#388e3c', '#006064', '#1976d2', '#5e35b1', '#c2185b', '#78909c', '#8d6e63', '#fdd835', '#689f38', '#00796b', 
               '#0288d1', '#303f9f', '#7b1fa2', '#d32f2f', '#bdbdbd', '#d84315', '#ffeb3b', '#7cb342', '#00897b', '#039be5', '#3949ab', '#8e24aa', 
               '#e53935', '#e0e0e0', '#e64a19', '#c0ca33', '#43a047', '#00838f', '#1e88e5', '#673ab7', '#d81b60', '#90a4ae', '#a1887f', '#f57c00',
               '#cddc39', '#4caf50', '#0097a7', '#2196f3', '#7e57c2', '#e91e63', '#b0bec5', '#bcaaa4', '#fb8c00', '#8bc34a', '#009688', '#03a9f4',
               '#3f51b5', '#9c27b0', '#f44336', '#eeeeee', '#f4511e', '#ffb300', '#9ccc65', '#26a69a', '#29b6f6', '#5c6bc0', '#ab47bc', '#ef5350',
               '#f5f5f5', '#ff5722', '#ffc107', '#66bb6a', '#00acc1', '#42a5f5', '#9575cd', '#ec407a', '#cfd8dc', '#d7ccc8', '#ff9800', '#ffee58',
               '#81c784', '#00bcd4', '#64b5f6', '#b39ddb', '#f06292', '#eceff1', '#efebe9', '#ffa726', '#fff176', '#4db6ac', '#4fc3f7', '#7986cb',
               '#ba68c8', '#e57373', '#73e575', '#ff7043', '#ffca28', '#d4e157', '#e6ee9c', '#a5d6a7', '#4dd0e1', '#90caf9', '#d1c4e9', '#f48fb1',
               '#ff8a65', '#ffd54f', '#dce775', '#c5e1a5', '#80cbc4', '#81d4fa', '#9fa8da', '#ce93d8', '#ef9a9a', '#ffb74d', '#fff59d', '#aed581',
               '#ffecb3', '#f0f4c3', '#c8e6c9', '#80deea', '#bbdefb', '#ede7f6', '#f8bbd0', '#ffab91', '#ffe082', '#fffde7', '#dcedc8', '#b2dfdb',
               '#b3e5fc', '#c5cae9', '#e1bee7', '#ffcdd2', '#ffcc80', '#fff9c4', '#1965b0', '#dc050c', '#fff3e0', '#f9fbe7', '#e8f5e9', '#b2ebf2',
               '#e3f2fd', '#f3e5f5', '#ffebee', '#fb8072', '#e0f7fa', '#fff8e1', '#f1f8e9', '#e0f2f1', '#e1f5fe', '#e8eaf6', '#fce4ec', '#ffe0b2',
               '#7bafde', '#d5de7b', '#7570b3', '#a6761d', '#55a1b1', '#33a02c', '#e7298a', '#ff7f00', '#882e72', '#999999', '#beaed4', '#e6ab02',
               '#8dd3c7', '#b2df8a', '#e78ac3', '#fdb462', '#b17ba6')
#269 + 269 = 538 HEX Colors code

color.marker <- new(Class = "character")
color.pheno <- new(Class = "character")
color.clust <- new(Class = "character")
color.metaclust <- new(Class = "character")

bell <- load.wave(where = "www/shipsbell.wav")
dingdong <- load.wave(where = "www/dingdong.wav")
gong <- load.wave(where = "www/gong.wav")

options(warn = 0)
align_list <- NULL
png()
dev.off()

############################################################################# service functions  -----

del_tmpdata <- function(type = "temp_dir") {

if (type == "all"){  
  temp_dir = "./tmpdata"
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)} else
      #Delete files if they exist
    {old_file <- list.files(path = temp_dir)
    if (!(length(old_file)==0)){
      old_file <- paste(temp_dir, old_file, sep = "/")
      unlink(x = old_file, recursive = TRUE, force = TRUE)}}
  
  temp_dir = "./tmpdata/renamed"
    old_file <- list.files(path = temp_dir)
    if (!(length(old_file)==0)){
      old_file <- paste(temp_dir, old_file, sep = "/")
      unlink(old_file, recursive = FALSE)}

  temp_dir = "./resultsQC"
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)} else
      #Delete files if they exist
    {old_file <- list.files(path = temp_dir)
    if (!(length(old_file)==0)){
      old_file <- paste(temp_dir, old_file, sep = "/")
      unlink(old_file, recursive = FALSE)}}

  tmpdir <- tempdir()
  old_file <- list.files(path = tmpdir)
  if (!(length(old_file)==0)){
    old_file <- paste(tmpdir, old_file, sep = "/")
    unlink(old_file, recursive = FALSE)}}
  else
  {tmpdir <- tempdir()
  old_file <- list.files(path = tmpdir)
  if (!(length(old_file)==0)){
    old_file <- paste(tmpdir, old_file, sep = "/")
    unlink(old_file, recursive = FALSE)}}
}


clean_graph <- function(){
  if((names(dev.cur()) == "null device")){
    return()
    }else{dev.off()}
}

# This function it to fix the file names. The internal make.names() is too restricted because does not admit even the 
# "+" and the "-" and too many other chars. Here only a subset of dangerous characters in excel and forbidden in the main OSs
# 
forbidden_char <- function(caratteri){
  for (i in 1:length(caratteri)){
    entry_string <- caratteri[i]
    entry_string <- chartr(old = "|", new = "_", x = entry_string)
    entry_string <- chartr(old = "#", new = "_", x = entry_string)
    entry_string <- chartr(old = "?", new = "_", x = entry_string)
    entry_string <- chartr(old = "*", new = "_", x = entry_string)
    entry_string <- chartr(old = "<", new = "_", x = entry_string)
    entry_string <- chartr(old = ">", new = "_", x = entry_string)
    entry_string <- chartr(old = "!", new = "_", x = entry_string)
    entry_string <- chartr(old = "[", new = "_", x = entry_string)
    entry_string <- chartr(old = "]", new = "_", x = entry_string)
    entry_string <- chartr(old = "{", new = "_", x = entry_string)
    entry_string <- chartr(old = "}", new = "_", x = entry_string)
    entry_string <- chartr(old = ";", new = "_", x = entry_string)
    entry_string <- chartr(old = ",", new = "_", x = entry_string)
    entry_string <- chartr(old = "\\", new = "_", x = entry_string)
    entry_string <- chartr(old = "/", new = "_", x = entry_string)
    caratteri[i] <- entry_string}
  new_strings <- caratteri
  return (new_strings)
}


##################################################################### general functions  -----

loading_fS <- function(inputfile, tipo) {
  if (tipo == "flowSet"){

    options(warn = 1) # Turn warnings into errors so they can be trapped
    fSresult <- try(expr = read.flowSet(files = inputfile, truncate_max_range = FALSE), silent = FALSE)
    ## ----------------------------------------------------------------------------------- ##
    # To compare two different samples in case of an error (typically because they have different numbers of channels: #
    # FCS1 <- read.FCS(filename = "..\\..\\sample_name.fcs", truncate_max_range = FALSE)
    # FCS2 <- read.FCS(filename = "..\\..\\sample_name.fcs", truncate_max_range = FALSE)
    # data1 <- FCS1@parameters@data
    # data2 <- FCS2@parameters@data
    # ... and check the names
    ##
    fcs_files_name <- basename(inputfile)
    if (inherits(x = fSresult,"flowSet")){ 
      options(warn = 0)
    
      for (i in 1:length(fSresult)){
        dim <- dim(exprs(fSresult[[i]]))[1]
        if (dim == 0){
          fSresult <- paste0("The sample nr. ", as.character(i), "  does not contain any event")
          break}}
      
      lista <- list(fSresult, fcs_files_name)
      return(lista)}
    else{ # Process any error messages
      # Ignore warnings while processing errors
      options(warn = -1)
      msg <- fSresult #it was msg <- geterrmessage()
      
      fF_dim <- vector(mode = "numeric", length = length(inputfile))
      options(warn = 0)
      for (i in 1:length(inputfile)){
      fFresult <- try(expr = read.FCS(filename = inputfile[i], truncate_max_range = FALSE), silent = FALSE)
      if (inherits(x = fFresult,"flowFrame")){
        fF_dim[i] <- dim(exprs(fFresult))[2]
        no_break <- TRUE} # this is the case in which there is a error at flowSet level: fcs with different dimensions
      else{
        msg <- fFresult
        no_break = FALSE # this is the case in which there are problems at flowFrame level
        break}}
      sort_dim <- sort(fF_dim, decreasing = TRUE)
      if (length(unique(sort_dim)) != 1){
        if (no_break == TRUE){
          livelli <- as.character(sort(as.numeric(levels(factor(sort_dim))), decreasing = TRUE))
          quale1 <- which(fF_dim %in% as.numeric(livelli[1]))
          quale2 <- which(fF_dim %in% as.numeric(livelli[2]))
          msg <- paste0("Some samples have ",livelli[1]," dimensions (aka markers), but sample nr.",as.character(quale2[1]), 
                        " has only ",livelli[2], " dimensions. 
                        Some samples have ",livelli[2]," dimensions, but sample nr.",as.character(quale1[1]), 
                        " has ",livelli[1], " dimensions. Please check if the samples belong to the same experiment")}
      else { msg <- paste0("Samples nr. ",as.character(i)," has some problems. System reports: <<", fFresult, ">>. 
                           Try to remove the sample from the flowSet")}}
      
      if (length(unique(sort_dim))==1){
        msg <- paste0("There are issues on selected flowSet. System reports: <<", msg, ">> ")}
      
      if (grepl("missing value", msg)) {
        result <- paste("USER ERROR: Did you supply any argument?")} 
      else if (grepl("This does not seem to be a valid",msg)&&(no_break == TRUE)) {
        # not valid entry
        result <- paste0("The files selected cannot be parsed as regular FCS files: try with another flowSet. System reports: <<", 
                         fSresult, ">>")}
      else {result <- paste0("Some errors found on parsing the selected files. ", msg)}
      # Restore default warning reporting
      options(warn=0)
      lista <- list(result, fcs_files_name)
      return(lista)}}
  
  if (tipo == "zip"){
    unzip(zipfile = inputfile, list = FALSE, exdir = data_dir)
    options(warn = 1) # Turn warnings into errors so they can be trapped
    fcs_files <- list.files(path = data_dir, full.names = T)
    fcs_files_name <- basename(fcs_files)
    fSresult <- try(expr = read.flowSet(files = fcs_files, truncate_max_range = FALSE), silent = FALSE)
    if (inherits(x = fSresult,"flowSet")){ 
      options(warn = 0)
      
      for (i in 1:length(fSresult)){
        dim <- dim(exprs(fSresult[[i]]))[1]
        if (dim == 0){
          fSresult <- paste0("The sample nr. ", as.character(i), "  does not contain any event")
          break}}
      
      lista <- list(fSresult, fcs_files_name)
      return(lista)}
    else{ # Process any error messages
      # Ignore warnings while processing errors
      options(warn = -1)
      msg <- fSresult #it was msg <- geterrmessage()
      
      fF_dim <- vector(mode = "numeric", length = length(fcs_files))
      options(warn = 0)
      for (i in 1:length(fcs_files)){
        fFresult <- try(expr = read.FCS(filename = fcs_files[[i]], truncate_max_range = FALSE), silent = FALSE)
        if (inherits(x = fFresult,"flowFrame")){
          fF_dim[i] <- dim(exprs(fFresult))[2]
          no_break <- TRUE} # this is the case in which there is a error at flowSet level: fcs with different dimensions
        else{
          msg <- fFresult
          no_break = FALSE # this is the case in which there are problems at flowFrame level
          break}}
      sort_dim <- sort(fF_dim, decreasing = TRUE)
      if (length(unique(sort_dim)) != 1){
        if (no_break == TRUE){
          livelli <- as.character(sort(as.numeric(levels(factor(sort_dim))), decreasing = TRUE))
          quale1 <- which(fF_dim %in% as.numeric(livelli[1]))
          quale2 <- which(fF_dim %in% as.numeric(livelli[2]))
          msg <- paste0("Some samples have ",livelli[1]," dimensions (aka markers), but sample nr.",as.character(quale2[1]), 
                        " has only ",livelli[2], " dimensions. 
                        Some samples have only ",livelli[2]," dimensions, but sample nr.",as.character(quale1[1]), 
                        " has ",livelli[1], " dimensions. Please check if the samples belong to the same experiment")}
        else { msg <- paste0("Samples nr. ",as.character(i)," has some problems. System reports: <<", fFresult, ">>. 
                           Try to remove the sample from the flowSet")}}
      
      if (length(unique(sort_dim))==1){
        msg <- paste0("There are issues on selected flowSet. System reports: <<", msg, ">> ")}
      
      if (grepl("missing value", msg)) {
        result <- paste("USER ERROR: Did you supply any argument?")} 
      else if (grepl("This does not seem to be a valid",msg)&&(no_break == TRUE)) {
        # not valid entry
        result <- paste0("The files selected cannot be parsed as regular FCS files: try with another flowSet. System reports: <<", 
                         fSresult, ">>")}
      else {result <- paste0("Some errors found on parsing the selected files. ", msg)}
      # Restore default warning reporting
      options(warn=0)
      lista <- list(result, fcs_files_name)
      return(lista)}}
  
}


file_fcs <- function(flow_Set, stringhe) {

  #this function is necessary to unlink the original flowSet from the handled one
  write.flowSet(x = flow_Set, outdir = "./tmpdata", filename = stringhe)
  fFfiles <- paste0("./tmpdata/", stringhe)
  fS <- loading_fS(inputfile = fFfiles, tipo = "flowSet")[[1]] 
  return(fS)
}


table_fS <- function(flow_Set, file_name = NULL) {

  fSlength <- length(flow_Set)
  tablefS <- tibble(`sample name` = "name",`sample id` = "sample id", `number of events` = 1, .rows = length(flow_Set))
  
  for (i in seq_along(1:fSlength)){
    if (!is.null(file_name)) tablefS$`sample name`[i] = file_name[[i]]
    else  { 
      if (inherits(x = flow_Set,"flowSet"))
      tablefS$`sample name`[i] = flow_Set@phenoData@data$name[i] 
    else 
      tablefS$`sample name`[i] = flow_Set@description[["$FIL"]]}
    
    if (fSlength < 100){
      if (i < 10) 
        {tablefS$`sample id`[i] = paste0("sample_0", as.character(i))}
      else
        {tablefS$`sample id`[i] = paste0("sample_", as.character(i))}} 
    else # if fSlength >= 100
      {if (i < 10) 
        {tablefS$`sample id`[i] = paste0("sample_00", as.character(i))}
      else { 
        if (i < 100)
          {tablefS$`sample id`[i] = paste0("sample_0", as.character(i))} 
        else 
          {tablefS$`sample id`[i] = paste0("sample_", as.character(i))}}}
    
    if (inherits(x = flow_Set,"flowSet")) 
      {tablefS$`number of events`[i] = dim(exprs(flow_Set[[i]]))[1]}
    else
      {tablefS$`number of events`[i] = dim(exprs(flow_Set))[1]}}
  
  file_name <- paste0("./tmpdata/fS_table.csv")
  write.csv(x = tablefS, file = file_name)
  
  fS_summary <- data.frame("Nr of samples" = length(flow_Set), "Total nr of events" = sum(tablefS$`number of events`),
                           "Minumum number of events through samples" = min(tablefS$`number of events`),
                           "Sample with the minimum number of events" =  tablefS[tablefS$`number of events` == min(tablefS$`number of events`),]$`sample id`[1],
                           "Maximum number of events through samples" = max(tablefS$`number of events`),
                           "Sample with the maximum number of events" =  tablefS[tablefS$`number of events` == max(tablefS$`number of events`),]$`sample id`[1],
                           "mean number of events through samples" = round(mean(tablefS$`number of events`), digits = 1),
                           "standard deviation of events through samples" = round(sd(tablefS$`number of events`),digits = 1), 
                           check.names = FALSE)
                           #row.names = "flowSet  statistic")
  
  fS_summary <- t(fS_summary)
  colnames(fS_summary) <- "Sample summary"
  
  pl <- ggplot(tablefS, aes(x = `sample id`, y=`number of events`))+
    geom_bar(stat="identity", width=0.7, fill="steelblue")+
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45,hjust=1)) 
  
  lista <- list(tablefS, pl, fS_summary)
  return(lista)
}


table_fF <- function(flow_Frame) {

  tablefF <- tibble(`sample name` = "name",`sample id` = "sample id", `number of events` = 1, .rows = length(flow_Frame))
  tablefF$`sample name`[1] = flow_Frame@parameters@data$name[1]
  tablefF$`sample id`[1] = paste("sample_0", as.character(1), sep = "")
  tablefF$`number of events`[1] = dim(exprs(flow_Frame))[1]
  return(tablefF)}


table_QC <- function(dir_QC) {
  lista_file <- list.files(path = dir_QC, pattern = ".txt$", full.names = T) 
  if (length(lista_file) > 0)
    {QC_table <- read_tsv(file = lista_file)}
  return(QC_table)
}


dim_fS <- function(flow_Set) {
for (i in 1:length(flow_Set)){
    pDname <- as.character(pData(parameters(flow_Set[[i]]))$name)
    pDdesc <- as.character(pData(parameters(flow_Set[[i]]))$desc)
    name_mask <- pDdesc[!is.na(pDdesc)]
    mask <- !(pDdesc %in% name_mask)
    chan_undescribed <- which(mask) #gives the position of the un-descripted channel
    #chan_fisici <- grep("FSC-A|FSC-H|SSC-A|SSC|TIME", pDname, value = FALSE, ignore.case = TRUE) #to find the lasers and the time
    chan_fisici <- grep("FSC-|SSC-|TIME", pDname, value = FALSE, ignore.case = TRUE) #to find the lasers and the time
    
    mask <- !(chan_undescribed %in% chan_fisici)
    chan_to_describe <- chan_undescribed[which(mask)] # this gives the position of the undescribed channels only
    
    if (length(chan_to_describe)>0){
      flow_Set[[i]]@parameters@data$desc[chan_to_describe] <- flow_Set[[i]]@parameters@data$name[chan_to_describe]}}
    #This part has been introduced to put a name in the channel without description
    #The description will then match the name of the marker
  res.expr <- exprs_sub(flow.Set = flow_Set)
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
  posizione <- which(nomi == "density")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]
    desc <- desc[-posizione]}
  
  posizione <- which(nomi == "cell_Id")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]
    desc <- desc[-posizione]}
  dimfS <- tibble(`dimension name` = nomi, `marker description` = desc) %>% rowid_to_column("Id") %>% drop_na()
  
  file_name <- paste0("./tmpdata/fS_dim.csv")
  write.csv(x = dimfS, file = file_name)
  
  return(dimfS)
}


set_colors <- function(number_of_colors) {
MAX_NUMBER_COLOR_BREWER = 74 #from brewer.pal
MAX_NUMBER_COLOR_SET = 2*269 #from the color_set
  if (number_of_colors <= MAX_NUMBER_COLOR_BREWER){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  rmd_color <- sample(col_vector, number_of_colors, replace = FALSE)
  return(rmd_color)}
  if (number_of_colors <= MAX_NUMBER_COLOR_SET){
    return(color_set[1:number_of_colors])
  }else
    {return("You exceed the maximum number of samples")}
}


pie_color <- function(color_vector){
  
  color_length <- length(color_vector)    
  fFvectorName <- vector(mode = "character", color_length)
  
  for (i in seq_along(1:color_length)){
    if (color_length < 100){
      if (i < 10) 
        {fFvectorName[i] = paste0("sample_0", as.character(i))}
      else
        {fFvectorName[i] = paste0("sample_", as.character(i))}} 
    else # if color_length >= 100
      {if (i < 10) 
        {fFvectorName[i] = paste0("sample_00", as.character(i))}
      else { 
        if (i < 100)
          {fFvectorName[i] = paste0("sample_0", as.character(i))} 
        else 
          {fFvectorName[i] = paste0("sample_", as.character(i))}}}}
  
  df <- data.frame(group = fFvectorName,value = 100/color_length)
  bp <- ggplot(df, aes(x="", y=value, fill=group)) +
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=color_vector)
  return(pie)
}

exprs_sub <- function(flow.Set, selected.markers = NULL) {

  pDname <- as.character(pData(parameters(flow.Set[[1]]))$name)
  pDdesc <- as.character(pData(parameters(flow.Set[[1]]))$desc)
  name_mask <- pDdesc[!is.na(pDdesc)]
  mask <- pDdesc %in% name_mask
  pDname_subset <- c(unname(pDname[mask]))
  pDdesc_subset <- c(unname(pDdesc[mask]))
  
  #nomi <- grep("FSC-|SSC-|TIME|cell_Id", name_mask, value = FALSE, ignore.case = TRUE)
  nomi <- grep(pattern = "FSC-|SSC-|TIME|DENSITY|SCORE|CELL_ID", x = name_mask, ignore.case = TRUE, value = FALSE)
  
  if (length(nomi)!=0){
    pDname_subset <- pDname_subset[-c(nomi)]
    pDdesc_subset <- pDdesc_subset[-c(nomi)]}
  
  if (!is.null(selected.markers)){
    pDname_subset <- pDname_subset[selected.markers]
    pDdesc_subset <- pDdesc_subset[selected.markers]}
  
  expr <- fsApply(flow.Set, exprs)
  exprsub <- expr[, pDname_subset]
  lista <- list(exprsub, pDname_subset, pDdesc_subset)
  return(lista)
}

# This functions (read_fcs, write_fcs, update_flowFrame_keywords, write_flowFrame, read_parameters, get_panel_table, rename_fcs_parameters_name, rename_fcs_parameters_desc, 
# rename_parameters_in_files, get_common_names, get_problem_idx) 
# are from premessa: https://github.com/ParkerICI/premessa 

read_fcs <- function(f.name) {
  
  fcs <- flowCore::read.FCS(f.name, truncate_max_range = FALSE)
  ret <- list()
  ret$m <- flowCore::exprs(fcs)
  p.names <-  as.character(flowCore::parameters(fcs)$name)
  colnames(ret$m) <- p.names
  ret$desc <- as.character(flowCore::parameters(fcs)$desc)
  names(ret$desc) <- p.names
  ret$keywords <- flowCore::keyword(fcs)
  return(ret)
}


write_fcs <- function(fcs, out.name) {
  
  # Drop parameter keywords
  keys <- fcs$keywords
  keys <- keys[grep("\\$P[0-9]+.", names(keys), invert = T)]
  
  # Exclude keywords that will be written with the new file
  excl.keywords <- c("FCSversion", "$BEGINANALYSIS", "$BEGINSTEXT", "$BYTEORD", "$DATATYPE", "$ENDANALYSIS",
                     "$ENDSTEXT", "$MODE", "$NEXTDATA", "$TOT", "$PAR", "$BEGINDATA", "$ENDDATA")
  keys <- keys[!(names(keys) %in% excl.keywords)]
  flow.frame <- flowCore::flowFrame(fcs$m)
  flow.frame <- update_flowFrame_keywords(flow.frame, fcs$m, fcs$desc, data.range = 262144)
  flowCore::keyword(flow.frame) <- keys
  marker.names <- fcs$desc
  names(marker.names) <- colnames(fcs$m)
  flowCore::markernames(flow.frame) <- marker.names
  write_flowFrame(flow.frame, out.name)
}

#Modified from https://github.com/nolanlab/cytofCore
#Fixes the parameters in the flowFrame, based on information
#from the corresponding exprs matrix
#if desc is different from NULL, also sets the parmeter description
#(i.e. $PnS)
#Set the range from the data or from a fixed value
update_flowFrame_keywords <- function(flowFrame, exprs.m, desc = NULL, data.range = "data") {
  
  params <- flowCore::parameters(flowFrame)
  pdata <- flowCore::pData(params)
  
  if(is.null(desc))
    desc <- colnames(exprs.m)
  
  for (i in 1:ncol(flowFrame)) {
    s <- paste("$P",i,"S",sep="")
    n <- paste("$P",i,"N",sep="")
    r <- paste("$P",i,"R",sep="")
    b <- paste("$P",i,"B",sep="")
    e <-  paste("$P",i,"E",sep="")
    
    keyval <- list()
    
    if(!is.na(desc[i]))
      keyval[[s]] <- desc[i]
    
    keyval[[n]] <- colnames(exprs.m)[i]
    
    if(data.range == "data")
      keyval[[r]] <- ceiling(max(exprs.m[,i], na.rm = TRUE))
    else if(is.numeric(data.range))
      keyval[[r]] <- data.range
    else
      stop("Invalid data.range parameter")
    
    keyval[[b]] <- 32
    keyval[[e]] <- "0,0"
    flowCore::keyword(flowFrame) <- keyval
    
    
    pdata[i,"minRange"] <- min(exprs.m[,i], na.rm = TRUE)
    pdata[i,"maxRange"] <- max(exprs.m[,i], na.rm = TRUE)
    
  }
  
  flowCore::pData(params) <- pdata
  flowCore::parameters(flowFrame) <- params
  
  # keyval[["$DATATYPE"]] <- "F"
  return(flowFrame)
}


#' Write a flowFrame as FCS file
#'
#' This function writes a flowFrame as an FCS file, taking care of updating the \code{$FILENAME} keyword
#'
#' @param flowFrame the \code{flowFrame} to write
#' @param path destination path
#' @export
write_flowFrame <- function(flowFrame, path) {
  
  f.name <- basename(path)
  flowCore::keyword(flowFrame)[["$FIL"]] <- f.name
  flowCore::write.FCS(flowFrame, path)
  return(invisible(NULL))
}



get_panel_table <- function(files.list) {
  
  message("Reading FCS parameters...")
  panel.table <- read_parameters(files.list)
  message("Done")
  common.names <- get_common_names(panel.table)
  problem.idx <- get_problem_idx(panel.table, common.names)
  panel.table <- panel.table[, order(colSums(problem.idx), decreasing = T), drop = F]
  
  panel.table <- data.frame(Parameter = row.names(panel.table), common.names, panel.table, check.names = F, stringsAsFactors = F)
  names(panel.table)[2] <- "Most common"
  return(panel.table)
}


read_parameters <- function(files.list) {
  
  ret <- lapply(files.list, function(f) {
    fcs <- flowCore::read.FCS(f, which.lines = 1)
    df <- data.frame(name = as.character(flowCore::parameters(fcs)$name),
                     desc = as.character(flowCore::parameters(fcs)$desc), check.names = F,
                     stringsAsFactors = F)
    df$desc[is.na(df$desc)] <- ""
    names(df) <- gsub("desc", basename(f), names(df))
    return(df)
    
  })
  
  ret <- Reduce(function(a, b) {
    merge(a, b, by = "name", all.x = T, all.y = T)
  }, ret)
  
  ret <- data.frame(ret, check.names = F, stringsAsFactors = F)
  
  row.names(ret) <- ret$name
  ret$name <- NULL
  
  return(ret)
}


#' Rename FCS parameter names
#'
#' This function renames FCS parameter names (i.e. the $PnN keyword)
#'
#' @param fcs FCS data as returned by \code{read_fcs}
#' @param names.map A named character vector mapping old names to new names.
#'   The \code{names} of this vector should be the old names, and the values
#'   the corresponding new names. New names need to be unique
#'
#' @return Returns the FCS data with the updated names
#'
#' @export
rename_fcs_parameters_name <- function(fcs, names.map) {
  
  ret <- fcs
  old.names <- colnames(ret$m)
  new.names <- as.character(names.map[old.names])
  stopifnot(!any(duplicated(new.names)))
  colnames(ret$m) <- new.names
  names(ret$desc) <- new.names
  
  spill.keyword <- grep("SPILL", names(ret$keywords), value = T)
  
  if(length(spill.keyword) > 0) {
    m <- ret$keywords[[spill.keyword]]
    if (is.matrix(m)){
      colnames(m) <- as.character(names.map[colnames(m)])
      ret$keywords[[spill.keyword]] <- m}
  }
  return(ret)
}

#' Rename FCS parameter descriptions
#'
#' This function renames FCS parameter descriptions (also called "short names", i.e.
#'   the $PnS keyword)
#'
#' @param fcs FCS data as returned by \code{read_fcs}
#' @param names.map A named character vector mapping parameter names (i.e. $PnN) to
#'   the desired parameter descriptions ($PnS). The \code{names} of this vector
#'   should be the parameter names, and the values the parameter descriptions. All
#'   the parameter names of the \code{fcs} data must be present in \code{names(names.map)}
#'
#' @return Returns the FCS data with the updated descriptions
#'
#' @export
rename_fcs_parameters_desc <- function(fcs, names.map) {
  
  ret <- fcs
  stopifnot(all(colnames(ret$m) %in% names(names.map)))
  
  ret$desc <- names.map[colnames(ret$m)]
  return(ret)
}


#' Renames FCS parameters in a set of files
#'
#' This function renames FCS parameters in a set of files, as specified by a template
#'
#' @param working.dir The directory in which the files are located
#' @param out.dir The output directory. This will be created as a subfolder of
#'   \code{working.dir}
#' @param tab The table specifying how the renaming should be performed. This
#'   is a structure similar as the \code{data.frame} returned from \code{read_parameters}. \code{row.names}
#'   correspond to the current parameter names in the files, and columns correspond to files (which
#'   should be located in the \code{working.dir}). The value at \code{[row, column]} should contain the
#'   value you want to set for the corresponding parameter description (i.e. $PnS). In addition the
#'   \code{data.frame} must contain two extra columns with the following names
#'   \itemize{
#'      \item{\code{Remove}}{ logical indicating whether the parameter should be removed from the files}
#'      \item{\code{Parameter}}{ character indicating the desired new name of the parameter (i.e. the new
#'      $PnN)}
#'  }
#' @export
rename_parameters_in_files <- function(working.dir, out.dir, tab) {
  
  out.dir <- file.path(working.dir, out.dir)
  esiste <- dir.exists(out.dir)
  
  if (!esiste){
    dir.create(path = out.dir, recursive = T)}
  
  file.cols <- grep("Remove|Parameter", names(tab), invert = T)
  for(i in file.cols) {
    f.name <- file.path(working.dir, names(tab)[i])
    print(sprintf("Processing %s ...", f.name))
    fcs <- read_fcs(f.name)
    to.remove <- row.names(tab)[tab$Remove]
    if(length(to.remove) > 0)
        {fcs <- remove_parameters(fcs, to.remove)}
    names.map <- tab$Parameter
    names(names.map) <- row.names(tab)
    fcs <- rename_fcs_parameters_name(fcs, names.map)
    names.map <- tab[, i]
    names(names.map) <- tab$Parameter
    names.map <- names.map[!is.na(tab[, i])]
    fcs <- rename_fcs_parameters_desc(fcs, names.map)
    out.name <- file.path(out.dir, names(tab)[i])
    write_fcs(fcs, out.name)}
}


get_common_names <- function(tab) {
  
  ret <- apply(tab, 1, function(x) {
    x <- as.character(x)
    return(names(sort(table(x), decreasing = T)[1]))
  })
  
  return(as.vector(ret))
}

get_problem_idx <- function(tab, common.names) {
  
  ret <- tab != common.names
  ret[is.na(ret)] <- 1
  return(ret)
}

#to.remove = names of parameters to remove
remove_parameters <- function(fcs, to.remove) {
  
  ret <- fcs
  ret$m <- ret$m[, !(colnames(ret$m) %in% to.remove)]
  ret$desc <- ret$desc[colnames(ret$m)]
  return(ret)
}


add_parameter <- function(fcs, v, name, desc) {
  
  ret <- fcs
  ret$m <- cbind(ret$m, v)
  colnames(ret$m)[ncol(ret$m)] <- name
  
  ret$desc <- c(ret$desc, desc)
  names(ret$desc)[length(ret$desc)] <- name
  return(ret)
}


panel_fcs <- function(flow_Frame, stringa) {
  
  result <- try(read.FCS(filename = flow_Frame), silent = FALSE)
  
  if (inherits(x = result,"flowFrame")){
    nome <- write.FCS(x = result, filename = stringa)
    options(warn = 0)
    return(nome)}
  else{ # Process any error messages
    # Ignore warnings while processing errors
    options(warn = -1)
    msg <- geterrmessage()
    if (grepl("missing value", msg)) {
      result <- paste("USER ERROR: Did you supply any argument?")} 
    else if (grepl("This does not seem to be a valid FCS2.0, FCS3.0 or FCS3.1 file",msg)) {
      # not valid entry
      result <- "the files selected cannot be parsed as regular FCS files: try with another flowSet"} 
    else
    {result <- paste0("Some errors found on parsing the selected files: system reports: ", msg)}
    
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}

add_dim <- function(flow_Frame, dim_name, dim_vect) {

  expr <- exprs(flow_Frame)
  ggdf <- data.frame(expr, dim_vect, check.names = F, stringsAsFactors = F)
  names(ggdf)[names(ggdf) == "dim_vect"] <- dim_name
  newdim <- as.matrix(dim_vect)
  colnames(newdim) <- dim_name
  new_flowframe <- fr_append_cols(fr = flow_Frame, cols = newdim)
  lung <- length(names(new_flowframe@parameters@data$name))
  new_flowframe@description[["$PAR"]] <- lung
  return(new_flowframe)
}


cell_Id_fS <- function(flow_Set, stringhe) {
  
  fSlength <- length(flow_Set)
  fS_Id <- new(Class = "flowFrame") 
  
  for (i in 1:fSlength){
    pDname <- as.character(pData(parameters(flow_Set[[i]]))$name)
    pDdesc <- as.character(pData(parameters(flow_Set[[i]]))$desc)
    name_mask <- pDdesc[!is.na(pDdesc)]
    mask <- !(pDdesc %in% name_mask)
    chan_undescribed <- which(mask) #gives the position of the un-descripted channel
    chan_fisici <- grep("FSC-|SSC-|TIME|DENSITY|SCORE|CELL_ID", pDname, value = FALSE, ignore.case = TRUE) #to find the lasers and the time
    
    mask <- !(chan_undescribed %in% chan_fisici)
    chan_to_describe <- chan_undescribed[which(mask)] # this gives the position of the undescribed channels only
    
    if (length(chan_to_describe)>0){
      flow_Set[[i]]@parameters@data$desc[chan_to_describe] <- flow_Set[[i]]@parameters@data$name[chan_to_describe]}
    #This part has been introduced to put a name in the channel without description
    #The description will then match the name of the marker
    
    nomi_proibiti <- grep("cell_Id|density|score", pDname, value = FALSE, ignore.case = TRUE) 
    
    if (length(nomi_proibiti)>0) {fF_expr <- exprs(flow_Set[[i]])[,-nomi_proibiti]}
    else {fF_expr <- exprs(flow_Set[[i]])}
    
    fF_expr <- as.data.frame(x = fF_expr, optional = TRUE, stringsAsFactors = FALSE)
    #this is to add the cell_Id
    fF_expr$ID <- seq.int(nrow(fF_expr))
    fF_Id <- add_dim(flow_Frame = flow_Set[[i]], dim_name = "cell_Id", dim_vect = fF_expr$ID)
    
    if (fSlength < 100){
      if (i < 10) 
      {fFname = paste0("./tmpdata/Id_0",i,"_",stringhe[i])}
      else
      {fFname = paste0("./tmpdata/Id_",i,"_", stringhe[i])}} 
    else # if fSlength >= 100
    {if (i < 10) 
    {fFname = paste0("./tmpdata/Id_00",i,"_", stringhe[i])}
      else { 
        if (i < 100)
        {fFname = paste0("./tmpdata/Id_0",i,"_", stringhe[i])} 
        else 
        {fFname = paste0("./tmpdata/Id_",i,"_", stringhe[i])}}}
    
    write.FCS(x = fF_Id, filename = fFname)}
  
  fFfiles <- list.files(path = "./tmpdata", pattern="^Id_*", full.names = TRUE) 
  options(warn = 2) # Turn warnings into errors so they can be trapped
  result <- try(expr = read.flowSet(files = fFfiles, truncate_max_range = FALSE))
  
  if (inherits(x = result,"flowSet")){
    options(warn = 0)
    return(result)}
  else {  # Process any error messages
    # Ignore warnings while processing errors
    options(warn = -1)
    msg <- geterrmessage()
    result <- paste("Something went wrong with the selected flowSet. cell_Id_fS internal routine returns: ", 
                    msg, sep = " ")
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


save_fS <- function(in_flow_Set, handled_flow_Set, stringhe) {
  # This function takes in the in_flow_Set which is the original flowSet from which we take the marker expression in their 
  #  orginal form, and the handled_flowSet which is the flowSet handled by the pre-clustering algorithms so with (in general) a
  # subset of event with the marker expression changed by the transformation and or the alignment. It then save a flowSet with
  # all the flowFrames in .fcs with the both the handled and the original expression but of course with the related subset of 
  # events
  
  length_in <- length(in_flow_Set)
  fFvectorName <- vector(mode = "character", length_in)
  pDname <- as.character(pData(parameters(in_flow_Set[[1]]))$name)
  pDdesc <- as.character(pData(parameters(in_flow_Set[[1]]))$desc)
  name_mask <- pDdesc[!is.na(pDdesc)]
  mask <- pDdesc %in% name_mask
  pDname_subset <- c(unname(pDname[mask]))
  
  fFvectorName <- paste0("./tmpdata/complete_",stringhe)
  
  fSfiltered <-  vector("list",length_in)
  outfS <-  new(Class = "flowSet")
  
  for (i in 1:length_in){
    fF_expr_in <- exprs(in_flow_Set[[i]])
    fF_expr_in <- as.data.frame(x = fF_expr_in, optional = TRUE, stringsAsFactors = FALSE)
    cell_Id_in <- fF_expr_in$cell_Id  
    fF_expr_pre <- exprs(handled_flow_Set[[i]])
    fF_expr_pre <- as.data.frame(x = fF_expr_pre, optional = TRUE, stringsAsFactors = FALSE)
    cell_Id_pre <- fF_expr_pre$cell_Id
    #fF_expr_filtered <- fF_expr_in %>% filter(cell_Id_in %in% cell_Id_pre) doesn't work on shiny environment
    bool_filter <- cell_Id_in %in% cell_Id_pre
    fF_expr_in$bool <- bool_filter
    fF_expr_filtered <- fF_expr_in[fF_expr_in$bool == TRUE,]
    fF_expr_filtered$bool <- NULL
    fSfiltered[[i]] <- as.matrix(fF_expr_filtered)} 
  #this has a subset of the matrix of the original (in_flowSet) filtered by the cell_Id of the handled (pre_flowSet) 

  for (i in 1:length_in){
    outfF <- handled_flow_Set[[i]]
    for(k in 1:length(pDname_subset)){
      name <- pDname_subset[k]
      new_name <- paste0("<",name,">")
      if (!(new_name %in% pDname_subset)){
      outfF <- add_dim(flow_Frame = outfF, dim_name = new_name, dim_vect = fSfiltered[[i]][,name])}}
    
    write.FCS(x = outfF, filename = fFvectorName[i])
    print(paste0("sample nr.",i," produced"))}
}


dim_comparison <- function(beforefS, afterfS) {
  
  summary.b4 <- fsApply(beforefS,function(frame){(dim(exprs(frame)))})
  summary.b4 <- as.data.frame(summary.b4, optional = TRUE)
  summary.b4$sample_id <- rownames(summary.b4)
  colnames(summary.b4)<- c("entry num of events","num of dim", "sample_id")
  summary.b4 <- summary.b4[, c(3,1,2)]
  
  if(inherits(x = afterfS,"flowSet")){
    summary <- fsApply(afterfS,function(frame){(dim(exprs(frame)))})
    colnames(summary)<- c("num of events after handling","num of dim")
    summary <- as.data.frame(summary, optional = TRUE)}
  else {
    summary <- dim(exprs(afterfS))
    summary <- data.frame("num of events after handling" = summary[1], "num of dim" = summary[2], check.names = FALSE)}
  
  diffLength <- nrow(summary.b4) - nrow(summary)
  numrowsum <- nrow(summary)
  
  if (diffLength > 0) {
    for (i in seq_along(1:diffLength)){
      summary <- summary %>% add_row() 
      summary$`num of events after handling`[i+ numrowsum] <- 0}}
  
  res <- summary.b4
  res$`num of dim` <- NULL
  res$`num of events after handling` <- summary$`num of events after handling`
  res$difference <- summary.b4$`entry num of events` - res$`num of events after handling`
  res$`removed percentage` <- round((res$difference / summary.b4$`entry num of events`)*100, digits = 2)
  
  file_name <- paste0("./tmpdata/fS_comparison_table.csv")
  write.csv(x = res, file = file_name)
  
  melt_res <- res
  res$sample_id = NULL
  melt_res$difference <- NULL
  melt_res$`removed percentage` <- NULL
  melt_res <- melt(data = melt_res)
  
  fS_summary <- data.frame("Nr of samples" = length(afterfS),
                           "Total nr of events of the original flowSet" = sum(summary.b4$`entry num of events`),
                           "Total nr of events of the handled flowSet" = sum(summary$`num of events after handling`),
                           "Percentage of the removed events" = paste0(
                             round((1 - sum(summary$`num of events after handling`) / sum(summary.b4$`entry num of events`))*100, digits = 2),"%"),
                           
                           "Min number of events through the original samples" = min(summary.b4$`entry num of events`),
                           "Original sample with the min number of events" = paste0("sample nr.",as.character(which.min(summary.b4$`entry num of events`))),
                           "Max number of events through the original samples" = max(summary.b4$`entry num of events`),
                           "Original sample with the max number of events" = paste0("sample nr.",as.character(which.max(summary.b4$`entry num of events`))),
                           "mean number of events of the original flowSet" = round(mean(summary.b4$`entry num of events`), digits = 1),
                           "standard deviation of events of the original flowSet" = round(sd(summary.b4$`entry num of events`),digits = 1),
                           
                           "Min number of events through the handled samples" = min(summary$`num of events after handling`),
                           "Handled sample with the min number of events" = paste0("sample nr.",as.character(which.min(summary$`num of events after handling`))),
                           "Max number of events through the handled samples" = max(summary$`num of events after handling`),
                           "Handled sample with the max number of events" = paste0("sample nr.", as.character(which.max(summary$`num of events after handling`))), 
                           "mean number of events of the handled flowSet" = round(mean(summary$`num of events after handling`),digits = 1),
                           "standard deviation of events of the handled flowSet" = round(sd(summary$`num of events after handling`), digits = 1),
                           check.names = FALSE, row.names = "flowSet statistic")
  fS_summary <- t(fS_summary)
  
  pl <- ggplot(data = melt_res, aes(sample_id, value)) + geom_bar(aes(fill = variable), 
                                                                  width = 0.4, position = position_dodge(width=0.5), stat="identity") +  
    theme(legend.position="top", legend.title = 
            element_blank(),axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1),
          axis.title.y=element_blank())
  
  lista <- list(res, pl, fS_summary)
  return(lista)
}


plot_comparison <- function(b4fS, afterfS) {
  
  fSlength <- length(b4fS)
  pDname <- as.character(pData(parameters(b4fS[[1]]))$name)
  pDdesc <- as.character(pData(parameters(b4fS[[1]]))$desc)
  name_mask <- pDdesc[!is.na(pDdesc)]
  mask <- pDdesc %in% name_mask
  pDname_subset <- pDname[mask]
  pDdesc_subset <- pDdesc[mask]
  b4_df <- new(Class = "data.frame")
  after_df <- new(Class = "data.frame")
  fFvectorName <- vector(mode = "character", fSlength)
  
  for (i in seq_along(1:fSlength)){
    if (fSlength < 100){
      if (i < 10) 
        {fFvectorName[i] = paste0("sample_0", as.character(i))}
      else
        {fFvectorName[i] = paste0("sample_", as.character(i))}} 
    else # if fSlength >= 100
    {if (i < 10) 
      {fFvectorName[i] = paste0("sample_00", as.character(i))}
      else { 
        if (i < 100)
          {fFvectorName[i] = paste0("sample_0", as.character(i))} 
        else 
          {fFvectorName[i] = paste0("sample_", as.character(i))}}}
    
    if (inherits(x = b4fS,"flowSet")) 
    {expr <- exprs(b4fS[[i]])}
    else 
    {expr <- exprs(b4fS)}
    
    if (inherits(x = afterfS,"flowSet"))
    {exprd <- exprs(afterfS[[i]])}
    else 
    {exprd <- exprs(afterfS)}
    
    b4  <- as_tibble(expr)
    after <- as_tibble(exprd)
    tag <- 1:nrow(b4)
    tag[] <- 0
    b4$tag <- tag 
    tag <- 1:nrow(after)
    tag[] <- 1
    after$tag <- tag
    nomi <- names(b4)
    b4_left <- b4 %>% left_join(x = b4 , y = after, by = c("cell_Id"))
    
    #in tag.y there are now three values: 
    #0 are the original events; 
    #1 are the events left after manipulation; 
    #NULL are the discarded events only present in b4_df but absent in after_df
    
    nomi <- names(b4_left)
    nomisub <- str_subset(nomi, ".x$")
    nomisub <- append(nomisub, values = "tag.y", after = length(nomisub))
    b4_left <- b4_left %>% dplyr::select(all_of(nomisub))
    nomi <- names(b4_left)
    nomi <- gsub("\\.x", "", nomi)
    names(b4_left) <- nomi
    
    b4_left$tag.y[b4_left$tag.y == 1] <- 0
    b4_left$tag.y[is.na(b4_left$tag.y)] <- 1
    
    #in tag.y there are now two values: 0 are the remaining events; 1 are the discarded (tagged) events 
    
    b4_left$tag <- NULL
    b4 <- b4_left
    
    nomi_proibiti <- grep(pattern = "FSC-|SSC-|TIME|DENSITY|SCORE|CELL_ID", x = pDname_subset, ignore.case = TRUE, value = FALSE)
    
    if (length(nomi_proibiti)!=0){
      pDname_subset <- pDname_subset[-c(nomi_proibiti)]
      pDdesc_subset <- pDdesc_subset[-c(nomi_proibiti)]}
    
    pDname_subset <- unname(append(pDname_subset, values = "tag", after = length(pDname_subset)))
    
    b4 <- b4 %>% rename(tag=tag.y) %>% dplyr::select(all_of(pDname_subset))
    sample_ids <- rep(fFvectorName[[i]], times = nrow(b4))
    b4$sample_id <- sample_ids
    sample_ids <- rep(fFvectorName[[i]], times = nrow(after))
    after$sample_id <- sample_ids
    #b4_df <- rbind(b4_df, b4)
    b4_df <- b4
    after_df <- rbind(after_df, after)}
  
  ggdfmelt <- melt(b4_df, id.var = c("sample_id","tag"), value.name = "expression", variable.name = "antigen")
  if ((nrow(ggdfmelt)) > RENDER_DATA) {
    ggdfmelt <- sample_n(tbl = ggdfmelt, size = RENDER_DATA, replace = F)}
  
  ggdfmelt$tag <- as.factor(ggdfmelt$tag)
  levels(ggdfmelt$tag)[levels(ggdfmelt$tag)==0] <- "remaining event"
  levels(ggdfmelt$tag)[levels(ggdfmelt$tag)==1] <- "discarded event"
  
  legend_title <- list(yref='paper',xref="paper",y=1.05,x=1.1, text="type of events",showarrow=F)
  
  p <- plot_ly(data = ggdfmelt, x = ~expression, y = ~antigen, color = ~tag, colors = c("green", "red"), type = "box") %>% 
    layout(yaxis = list(title = "markers"), xaxis = list(title = "expression ranges"), annotations = legend_title)
  return(p)
}


concatenating_fS <- function(flow_Set, stringa) {
  
  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", length = fSlength)
  nome <- vector(mode = "character", length = 1)
 
  for (i in seq_along(1:fSlength)){
    if (fSlength < 100){
      if (i < 10) 
      {fFvectorName[i] = paste0(stringa, "_0", as.character(i))}
      else
      {fFvectorName[i] = paste0(stringa, "_", as.character(i))}} 
    else # if fSlength >= 100
    {if (i < 10) 
    {fFvectorName[i] = paste0(stringa,"_00", as.character(i))}
      else { 
        if (i < 100)
        {fFvectorName[i] = paste0(stringa, "_0", as.character(i))} 
        else 
        {fFvectorName[i] = paste0(stringa, "_", as.character(i))}}}}
  
  for (i in seq_along(1:fSlength)) {
    fFvectorName[[i]] <- paste0("./tmpdata/", fFvectorName[[i]], ".fcs")
    flowframe <- flow_Set[[i]]
    SampleID <- as.matrix(vector(mode = "numeric", length = dim(flow_Set[[i]])[1]))
    colnames(SampleID) <- c("SampleID")
    #if (flow_Jo == TRUE){
    #  point <- dim(flow_Set[[i]])[1]
    #  fFrnorm <- rnorm(n = point, mean = i*10000, sd = 200)
    #  SampleID[ ,1] <- fFrnorm
    #} else SampleID[,1] <- rep(i, dim(flowframe)[1])
    
    SampleID[,1] <- rep(i, dim(flowframe)[1])
    options(warn = 2) # Turn warnings into errors so they can be trapped
    result <- try(fr_append_cols(fr = flow_Set[[i]], cols = SampleID))
    
    if (inherits(x = result,"flowFrame")){
      options(warn = 0) # reset warning status
      write.FCS(result, filename = fFvectorName[[i]])}
    else{# Process any error messages
      # Ignore warnings while processing errors
      options(warn = -1)
      # If this script were a function, warning() and stop() could be called to pass errors upstream
      msg <- geterrmessage()
      messaggio <- paste0("System reports: ", msg, 
                          "... means 'something wrong with the fr_append_cols function inside internal concatenating_fs'")
      
      # Restore default warning reporting
      options(warn=0)
      return(messaggio)
    }}
  nomi_files <- paste0("^",stringa, "_*")
  fFfiles <- list.files(path = "./tmpdata", pattern = nomi_files, full.names = TRUE) 
  #fFfiles <- str_sort(x = fFfiles, numeric = TRUE) #this is necessary to sort in the correct way the flowFrames in the flowSet
  fS <- read.flowSet(files = fFfiles, truncate_max_range = FALSE)
  expr <- fsApply(fS, exprs)
  conc_flowframe <- fS[[1]]
  ggdf <- as.matrix(expr)
  exprs(conc_flowframe) <- ggdf
  return(conc_flowframe)
}


cleaning_fS <- function(flow_Frame, fFname, par1, par2, seed = NULL) {
 
  physical <- vector(mode = "character", length = 4) #just in case we have to deal with 4 physical channels
  pDname <- c(flow_Frame@parameters@data$name)
  nomi <- grep("FSC-|SSC-", pDname, value = FALSE, ignore.case = TRUE)
  if (length(nomi)>0){
    for(i in seq_along(1:length(nomi))){
      physical[i] <- unname(pDname[i])}
    for(i in seq_along(1:length(physical))){
      if (physical[i] == "")
      {physical <- physical[-i]}}}
  else {physical <- NULL}
  
  if ('cell_Id' %in% pDname)
    {physical <- c(physical, "cell_Id")}

  set.seed(seed)
  options(warn = 1) # Error is just printed - This is to avoid the funtion stops with the following message: 
  #"Error in ord_fcs_time(set[[i]], timeCh) : (converted from warning) Expression data in the file sample_03.fcs were not originally ordered by time.
  result <- try(expr = flow_auto_qc(fcsfiles = flow_Frame, remove_from = "all", output = 1, timeCh = NULL, 
                                    second_fractionFR = "timestep", alphaFR = par1, decompFR = FALSE, outlier_binsFS = FALSE,
                                    pen_valueFS = par2*10, max_cptFS = 3, sideFM = "both", fcs_QC = "_QC", 
                                    mini_report = "_QCmini", ChExcludeFS = physical,
                                    fcs_highQ = FALSE, fcs_lowQ = FALSE, html_report = FALSE), silent = FALSE)
  
  if (inherits(x = result,"flowFrame")){ 
    options(warn = 0)
    nome <- write.FCS(x = result, filename = fFname)
    return(nome)}
  
  # Process any error messages
  if (inherits(x = result, "try-error")) {
    # Ignore warnings while processing errors
    options(warn = -1)
    msg <- geterrmessage()
    options(warn=0)
    return(msg)}
}


visualize_rangefS <- function(flow_Set) {

  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", fSlength)
  for (i in seq_along(1:fSlength)){
    if (fSlength < 100){
      if (i < 10) 
        {fFvectorName[i] = paste0("sample_0", as.character(i))}
      else
        {fFvectorName[i] = paste0("sample_", as.character(i))}} 
    else # if fS_length >= 100
      {if (i < 10) 
        {fFvectorName[i] = paste0("sample_00", as.character(i))}
      else { 
        if (i < 100)
          {fFvectorName[i] = paste0("sample_0", as.character(i))} 
        else 
          {fFvectorName[i] = paste0("sample_", as.character(i))}}}}
  
  sample_ids <- rep(fFvectorName, fsApply(flow_Set, nrow))
  res.expr <- exprs_sub(flow.Set = flow_Set)
  exprsub <- res.expr[[1]]
  
  ggdf <- data.frame(sample_id = sample_ids, exprsub, check.names = FALSE)
  ggdf <- melt(ggdf, id.var = "sample_id", value.name = "expression", variable.name = "antigen")
  
  if ((nrow(ggdf)) > RENDER_DATA) {
    ggdf <- sample_n(tbl = ggdf, size = RENDER_DATA, replace = F)}
  
  p <- plot_ly(ggdf, x = ~expression, y = ~antigen) %>%
    add_boxplot() %>%
    layout(yaxis = list(title = "sample's markers")) 
  #%>% toWebGL()
  return(p)
}


visualize_expr_fS <- function(flow_Set, colonne = NULL,  width = NULL, height = NULL) {

  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", fSlength)
  
  for (i in seq_along(1:fSlength)){
    if (fSlength < 100){
      if (i < 10) 
      {fFvectorName[i] = paste0("sample_0", as.character(i))}
      else
      {fFvectorName[i] = paste0("sample_", as.character(i))}} 
    else # if fSlength >= 100
    {if (i < 10) 
    {fFvectorName[i] = paste0("sample_00", as.character(i))}
      else { 
        if (i < 100)
        {fFvectorName[i] = paste0("sample_0", as.character(i))} 
        else 
        {fFvectorName[i] = paste0("sample_", as.character(i))}}}}
  
  sample_color <- color.sample
  sample_ids <- rep(fFvectorName, fsApply(flow_Set, nrow))
  res.expr <- exprs_sub(flow.Set = flow_Set)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
  posizione <- which(nomi == "density")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]
    exprsub <- exprsub[,-posizione]
    desc <- desc[-posizione]
    colnames(exprsub) <- desc}
  
  posizione <- which(nomi == "cell_Id")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]
    exprsub <- exprsub[,-posizione]
    desc <- desc[-posizione]
    colnames(exprsub) <- desc}
  
  ggdf <- data.frame(sample_id = sample_ids, exprsub, check.names = FALSE)
  ggdf <- melt(ggdf, id.var = "sample_id", value.name = "expression", variable.name = "antigen")
  
  if ((nrow(ggdf)) > RENDER_DATA) {
    ggdf <- sample_n(tbl = ggdf, size = RENDER_DATA, replace = F)}
  
  if (is.null(colonne)){
    colonne <- ceiling(length(res.expr[[3]])/2)-1
    if (colonne<=0) {colonne=1}
    righe <- ceiling(length(res.expr[[3]])/2)}
  else {righe = length(res.expr[[3]])
  colonne = colonne}
  
  p <- ggplot(ggdf, aes(x = expression, color = sample_id, group = sample_id)) +
    geom_freqpoly(stat = "density") + 
    facet_wrap(~ antigen, scales = "free", ncol = colonne, nrow = righe) +
    scale_color_manual(values = sample_color)
  if (is.null(width))
    {p <- ggplotly(p = p, originalData = TRUE)}
  else
    {p <- ggplotly(p = p, originalData = TRUE, width, height)}
  return(p)
}


visualize_expr_fF <- function(flow_Frame) {

  pDname <- flow_Frame@parameters@data$name
  pDdesc <- flow_Frame@parameters@data$desc
  name_mask <- pDdesc[!is.na(pDdesc)]
  mask <- pDdesc %in% name_mask
  pDname_subset <- pDname[mask]
  pDdesc_subset <- pDdesc[mask]
  
  sample_color <- "red"
  expr <- flow_Frame@exprs
  
  if ("SampleID" %in% pDname){
    pDname_subset <- append(x = pDname_subset, values = "SampleID", after = length(pDname_subset))
    pDdesc_subset <- append(x = pDdesc_subset, values = "SampleID", after = length(pDdesc_subset))}
  
  exprsub <- expr[, pDname_subset]
  colnames(exprsub) <- pDdesc_subset
  sample_ids <- rep("conc sample", dim(flow_Frame@exprs)[1])
  ggdf <- data.frame(sample_id = sample_ids, exprsub, check.names = FALSE)
  ggdf <- melt(ggdf, id.var = "sample_id", value.name = "expression", variable.name = "antigen")
  
  if ((nrow(ggdf)) > RENDER_DATA) {
    ggdf <- sample_n(tbl = ggdf, size = RENDER_DATA, replace = F)}
  
  p <- ggplot(ggdf, aes(x = expression, color = sample_id, 
                        group = sample_id)) +
    geom_density() +
    facet_wrap(~ antigen, nrow = 4, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
    guides(color = guide_legend(ncol = 1)) +
    scale_color_manual(values = sample_color)
  
  p <- ggplotly(p)
  return(p)
}


trans_fS <- function(flow_Set, type_of, param, min_quantile, max_quantile){

  fSlength <- length(flow_Set)
  fSnapshot <-  vector("list",fSlength)
  for (i in 1:fSlength){
    fSnapshot[[i]] <- as.data.frame(exprs(flow_Set[[i]]))}
  
  pDname <- as.character(pData(parameters(flow_Set[[1]]))$name)
  res.expr <- exprs_sub(flow.Set = flow_Set)
  #exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  #desc <- res.expr[[3]]
  
  posizione <- which(nomi == "density")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]}
  
  posizione <- which(nomi == "cell_Id")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]}
  
  pDname_subset <- nomi
  
  ###################### arcsinh ------
  
  if (type_of=="arcsinh"){
    tfS <- fsApply(x = flow_Set, FUN = function(x, cofactor = param) {
      
      m <- exprs(x)[, pDname_subset]
      m <- asinh(m / param)
      exprs(x) <- m
      x
    })}
  
  ###################### scale ------
  
  if (type_of=="scale"){
    normalize.function <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))}
    
    tfS <- fsApply(flow_Set, FUN = function(x){
      m <- exprs(x)[, pDname_subset]
      m <- normalize.function(m)
      exprs(x) <- m
      x})}  
  
  ###################### zero ------
  if (type_of == "zero"){
    zeroize.function <- function(x) {
      rng <- colQuantiles(x, probs = c(min_quantile, max_quantile), drop = TRUE)
      x01 <- t((t(x) - rng[, 1]) / (rng[, 2] - rng[, 1]))
      x01[x01 < 0] <- 0
      x01[x01 > 1] <- 1
      x01}
    
    tfS <- fsApply(flow_Set, FUN = function(x){
      m <- exprs(x)[, pDname_subset]
      m <- zeroize.function(m)
      exprs(x) <- m
      x})}
  
  ###################### arcsinh+ ------
  if (type_of == "arcsinhplus"){
    tfS <- fsApply(flow_Set, FUN = function(x, cofactor = param) {
      m <- exprs(x)[, pDname_subset]
      m <- m - min(m)
      m <- asinh(m / param)
      exprs(x) <- m
      x})}
  
  ###################### none ------
  if (type_of == "none"){tfS <- flow_Set}
  
  ###################### flowTrans ------
  if (type_of == "mclMultivArcSinh" || type_of == "mclMultivBiexp" || type_of == "mclMultivBoxCox" || type_of == "mclMultivLinLog"){
    options(warn = 2) 
    lista <- vector(mode = "list", length = fSlength)
    tfS <- flow_Set
    for (i in 1:fSlength){
      tfF <- try(expr = flowTrans::flowTrans(dat = flow_Set[[i]], fun = type_of, dims = nomi, n2f = FALSE, parameters.only = FALSE))
      
      if (!(inherits(x = tfF,"try-error"))){ 
        options(warn = 0)
        tfS[[i]] = tfF$result
        param <- flowTrans::extractParams(x = tfF)[[1]]
        lista[[i]] <- param}
    else{
      options(warn = -1) # Ignore warnings while processing errors
      msg <- geterrmessage()
      flow_Set <- paste("Something went wrong with the scaling process. flowTrans() function returns: ", 
                        msg, sep = " ")
      # Restore default warning reporting
      options(warn=0)}}}  
  
 ###################### flowVS ------
  if (type_of == "flowVS"){
    #options(warn = 2) #commented because it often happens "Error in optimize(afun, b, maximum = TRUE) : (converted from warning) NA/Inf replaced by maximum positive value"
    result <- try(expr = estParamFlowVS(flow_Set,channels=pDname_subset))
    if (!(inherits(x = result,"try-error"))){ 
      options(warn = 0)
      #unfortuantely the transFlowVS put a description also in the prohibited channels (see the cell_Id_fS function)
      #the pDname used below is to one in the input flowSet
      lista_descrizioni <-  vector("list",fSlength)
      description_list <- vector("list", fSlength)
      for (i in 1:fSlength){
        lista_descrizioni[[i]] <- flow_Set[[i]]@parameters@data$desc
        description_list[[i]] <- flow_Set[[1]]@description}
      
      tfS = transFlowVS(flow_Set, channels=pDname_subset, cofactors = result)
      
      #nomi_proibiti <- grep("FSC-A|FSC-H|SSC-A|SSC|TIME|CELL_ID", pDname, value = FALSE, ignore.case = TRUE) #these are the prohibited names
      #for (i in 1:length(flow_Set)){
      #  for (k in 1:length(nomi_proibiti)){
      #    tfS[[i]]@parameters@data$desc[nomi_proibiti[k]] <- NA}
      #    write.FCS(x = tfS[[i]], filename = fFvectorName[i])}}
      
      for (i in 1:fSlength){
        tfS[[i]]@parameters@data$desc <- lista_descrizioni[[i]]
        tfS[[i]]@description <- description_list[[i]]}}
    else{
      options(warn = -1) # Ignore warnings while processing errors
      msg <- geterrmessage()
      flow_Set <- paste("Something went wrong with the scaling process. transFlowVS() function returns: ", 
                      msg, sep = " ")
      # Restore default warning reporting
      options(warn=0)}}
  
  ##################### flow_Set output 
  
  if (inherits(x = flow_Set, "flowSet"))
    {for (i in 1:fSlength){
      m <- as.matrix(fSnapshot[[i]])
      tm <- exprs(tfS[[i]])
      m[, pDname_subset] <- tm[, pDname_subset]
      exprs(flow_Set[[i]])[, pDname_subset] <- m[, pDname_subset]}}
  
  return(flow_Set)
}
  

align_flowStats <- function(flow_Set, marker_list, trans) {
  
  pDname <- as.character(pData(parameters(flow_Set[[1]]))$name)
  pDdesc <- as.character(pData(parameters(flow_Set[[1]]))$desc)
  dimfS <- dim_fS(flow_Set)
  #dimfS <- tibble(`dimension name` = pDname, `marker description` = pDdesc) %>% rowid_to_column("Id") %>% drop_na()  
  dimfS.sel <- dimfS[marker_list,]
  marker_names <- dimfS.sel$`dimension name`
  
  name_mask <- pDdesc[!is.na(pDdesc)]
  mask <- pDdesc %in% name_mask
  pDname_subset <- c(unname(pDname[mask]))
  
  res.expr <- exprs_sub(flow.Set = flow_Set)
  nomi <- res.expr[[2]]
  
  posizione <- which(nomi == "density")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]}
  
  posizione <- which(nomi == "cell_Id")
  if (length(posizione)>0){
    nomi <- nomi[-posizione]}
  
  pDname_subset <- c(unname(nomi))
  
  fSlength <- length(flow_Set)
  fSnapshot <-  vector("list",fSlength)
  for (i in 1:fSlength){
    fSnapshot[[i]] <- as.data.frame(exprs(flow_Set[[i]]))}
  
  #wf <- flowCore::workFlow(flow_Set) #'flowCore::workFlow' is deprecated. Use 'flowWorkspace::GatingSet' instead.
  
  righe <- ceiling(sqrt(length(marker_names)))
  if (righe==0) {righe <- 1}
  colonne <- floor(sqrt(length(marker_names)))
  if ((righe*colonne)<length(marker_names)) {righe <- righe +1}
  #wf <- flowCore::workFlow(flow_Set) #'flowCore::workFlow' is deprecated. Use 'flowWorkspace::GatingSet' instead.
  #plotb4 <- densityplot(name~., Data(wf[["asinh"]]), channels=marker_names, scales=list(y=list(draw=F)), filter=lapply(marker_names, curv1Filter), layout=c(colonne,righe), warn.unused = FALSE)
  #tl <- transformList(from = pDname_subset, tfun = asinh, transformationId="asinh")
  #flowWorkspace::add(wf, tl)
  #norm <- normalization(normFun=function(x, parameters, ...)
  #  warpSet(x, parameters, ...), parameters=marker_names, normalizationId="Warping")
  #add(wf, norm, parent="asinh")
  
  gs <- GatingSet(flow_Set)
  #trans = FALSE # even is the transformation has been performed, it is better to re-execute it here because the warpSet procedure
  #fails when not performed with error "ror in chol.default(Asym) : the leading minor of order #num is not positive definite" 
  #The problem is that the densit plot will be squeezed and quite awful
  if (trans==FALSE){
    transformation <- transformerList(from = marker_names, asinhtGml2_trans())
    gs <- transform(gs, transformation)}
  
  fS <- gs_cyto_data(gs)
  fS <- cytoset_to_flowSet(fS)
  plotb4 <- flowViz::densityplot(name~., fS, channels=marker_names, scales=list(y=list(draw=F)),
                                 filter=lapply(marker_names, curv1Filter), layout=c(colonne,righe), warn.unused = FALSE)
  #the densityplot function does not work properly
  options(warn = 1) # Turn warnings into errors so they can be trapped
  fS_out <- try(expr = warpSet(x = fS, stains = marker_names, warpFuns = F), silent = FALSE)
  #fs_out <- warpSet(x = fs, stains = marker_names, warpFuns = F)
  #fS_out <- cytoset_to_flowSet(fS_out)
  
  if (inherits(x = fS_out,"flowSet")){ 
    options(warn = 0)
    fS_out <- cytoset_to_flowSet(fS_out)
    plotafter <- flowViz::densityplot(name~., fS_out, channels=marker_names, scales=list(y=list(draw=F)), 
                             filter=lapply(marker_names, curv1Filter), layout=c(colonne,righe), warn.unused	= FALSE)}
  else{
    plotafter <- NULL
    msg <- fS_out #it was msg <- geterrmessage()
    fS_out <- paste0("Something went wrong with the aligning process. System reports: ", msg)
    options(warn = 0)
    return(list(fS_out, plotb4, plotafter))}
  
  if (inherits(x = fS_out,"flowSet")){
    for (i in 1:fSlength){
      m <- as.matrix(fSnapshot[[i]])
      tm <- exprs(fS_out[[i]])
      m[, pDname_subset] <- tm[, pDname_subset]
      exprs(flow_Set[[i]])[, pDname_subset] <- m[, pDname_subset]}}
  
  png(filename = "./tmpdata/PlotAlignB4.png", width = 1400, height = 210*fSlength, units = "px", res = 125)
  print(plotb4)
  dev.off()
  
  png(filename = "./tmpdata/PlotAlignAfter.png", width = 1400, height = 210*fSlength, units = "px", res = 125)
  print(plotafter)
  dev.off()
  
  lista <- list(flow_Set, plotb4, plotafter)
  return(lista)
}

#___________________________________________________________________________________________________
# The following functions have been added because of a problem with the original spade function which is obsolete with R4.0
# They are the original SPADE.addDensityToFCS and the SPADE.downsampleFCS. This is the error with the original function
#  Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'object' in selecting a method for function 'exprs': (converted from warning) '.local' is deprecated.
#  Use 'keyword' instead.
#____________________________________________________________________________________________________


spade.addDensityToFCS <- function (infilename, outfilename, cols = NULL, arcsinh_cofactor = NULL, 
          transforms = flowCore::arcsinhTransform(a = 0, b = 0.2), 
          kernel_mult = 5, apprx_mult = 1.5, med_samples = 2000, comp = TRUE) 
{
  if (!is.null(arcsinh_cofactor)) {
    warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
    transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
  }
  
  #this was the original function: 
  #in_fcs <- SPADE.transform.FCS(SPADE.read.FCS(infilename, comp = comp), transforms)
  in_fcs <- read.FCS(filename = infilename)
  in_data <- exprs(in_fcs)
  params <- parameters(in_fcs)
  pd <- pData(params)
  if (is.null(cols)) {
    cols <- as.vector(pd$name)
  }
  idxs <- match(cols, pd$name)
  if (any(is.na(idxs))) {
    stop("Invalid column specifier")
  }
  density <- SPADE.density(in_data[, idxs], kernel_mult = kernel_mult, 
                           apprx_mult = apprx_mult, med_samples = med_samples)
  if (max(density) == 0) 
    warning(paste(infilename, "has degenerate densities, possibly due to many identical observations", 
                  sep = " "))
  #this was the original function: 
  #in_fcs <- SPADE.read.FCS(infilename, comp = FALSE, transformation = FALSE)
  in_fcs <- read.FCS(filename = infilename)
  in_data <- exprs(in_fcs)
  channel_number <- ncol(in_fcs) + 1
  channel_id <- paste("$P", channel_number, sep = "")
  channel_name <- "density"
  channel_range <- max(density) + 1
  plist <- matrix(c(channel_name, channel_name, channel_range, 
                    0, channel_range - 1))
  rownames(plist) <- c("name", "desc", "range", 
                       "minRange", "maxRange")
  colnames(plist) <- c(channel_id)
  pd <- rbind(pd, t(plist))
  pData(params) <- pd
  out_data <- cbind(in_data, density = density)
  #this was the original function: 
  #out_frame <- flowFrame(out_data, params, description = description(in_fcs))
  descr = keyword(in_fcs)
  out_frame <- new("flowFrame", exprs = out_data, parameters = params, description = descr)
  keyval <- list()
  keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
  keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
  keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
  keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
  keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
  keyword(out_frame) <- keyval
  write.FCS(out_frame, outfilename)
}


spade.downsampleFCS <- function (infilename, outfilename, exclude_pctile = 0.01, target_pctile = NULL, 
          target_number = NULL, target_percent = 0.1) 
{ 
  #this was the original function: 
  #in_fcs <- SPADE.read.FCS(infilename, comp = FALSE, transform = FALSE)
  in_fcs <- read.FCS(filename = infilename)
  in_data <- exprs(in_fcs)
  params <- parameters(in_fcs)
  pd <- pData(params)
  d_idx <- match("density", pd$name)
  if (is.na(d_idx)) {
    stop("No density parameter in FCS file")
  }
  boundary <- quantile(in_data[, d_idx], c(exclude_pctile, 
                                           target_pctile), names = FALSE)
  out_data <- subset(in_data, in_data[, d_idx] > boundary[1])
  if (!is.null(target_percent)) {
    target_number = round(target_percent * nrow(out_data))
    message("Targeting ", target_number, " events for ", 
            infilename)
  }
  density <- out_data[, d_idx]
  if (is.null(target_number)) {
    boundary <- boundary[2]
    out_data <- subset(out_data, boundary/density > runif(nrow(out_data)))
  }
  else if (target_number < nrow(out_data)) {
    density_s <- sort(density)
    cdf <- rev(cumsum(1/rev(density_s)))
    boundary <- target_number/cdf[1]
    if (boundary > density_s[1]) {
      targets <- (target_number - 1:length(density_s))/cdf
      boundary <- targets[which.min(targets - density_s > 
                                      0)]
    }
    out_data <- subset(out_data, boundary/density > runif(length(density)))
  }
  else if (target_number > nrow(out_data)) {
    stop("More events requested than present in file")
  }
  
  #this was the original function: 
  #out_frame <- flowFrame(out_data, params, description = description(in_fcs))
  descr = keyword(in_fcs)
  out_frame <- new("flowFrame", exprs = out_data, parameters = params, description = descr)
  write.FCS(out_frame, outfilename)
}

spade_fS <- function(flow_Set, stringhe, type, p1, p2, seed) {
  # the function returns also the flowSet with all the flowframes with the density only
  # it con be used to make some comparison before and after the downsampling
  
  fSlength <- length(flow_Set)
  pDname <- as.character(pData(parameters(flow_Set[[1]]))$name)
  pDdesc <- as.character(pData(parameters(flow_Set[[1]]))$desc)
  
  if ('density' %in% pDname)
  {return ("density channel is present in the sample's expression matrix. Please remove it before to downsampling")}
  
  name_mask <- pDdesc[!is.na(pDdesc)]
  mask <- pDdesc %in% name_mask
  pDname_subset <- pDname[mask]
  pDdesc_subset <- pDdesc[mask]
  nomi <- grep("FSC-|SSC-|TIME|DENSITY|SCORE|CELL_ID", name_mask, value = FALSE, ignore.case = TRUE)
  
  if (length(nomi)!=0){
    pDname_subset <- pDname_subset[-c(nomi)]
    pDdesc_subset <- pDdesc_subset[-c(nomi)]}
  
  fFvectorName <- vector(mode = "character", length = fSlength)
  
  fFvectorName <- stringhe
  fFvectorName_out <- fFvectorName
  if (type == "down_perc"){ #case of downsampling by perc
    suffix <- paste0("Dspade_perc",as.character(p1),"_",as.character(p2))}
  else
    suffix <- paste0("Dspade_value",as.character(p1),"_",as.character(p2))
  
  for (i in seq_along(1:fSlength)) {
    fFvectorName[[i]] <- paste0("./tmpdata/", fFvectorName[[i]]) # "./tmpdata/downsample_01.fcs"
    fFvectorName_out[[i]] <- paste0("./tmpdata/downdens_", fFvectorName_out[[i]]) #" ./tmpdata/downdens_downsample_01.fcs"
    set.seed(seed)
    options(warn = 1) # turn warnings into errors so they can be trapped
    result <- try(expr = spade.addDensityToFCS(infilename = fFvectorName[[i]], 
                                               outfilename = fFvectorName_out[[i]], cols = pDname_subset)) 
    
    if (!(inherits(x = result,"try-error"))){ 
      options(warn = 0)}
    else{
      # Ignore warnings while processing errors
      options(warn = -1)
      msg <- geterrmessage()
      result <- paste("Something went wrong with the downsampling process. SPADE.addDensity() function returns: ", 
                      msg, sep = " ")
      # Restore default warning reporting
      options(warn=0)
      return(result)}
    
    if (type == "down_perc"){ #case of downsampling by perc: the downsampling is made here and not in the donwsampling_fS
      if (result == fFvectorName_out[[i]]){
        options(warn = 0)
        nome = substr(x = fFvectorName[i], start = 11, stop=256)
        fFvectorName[i] <- paste0("./tmpdata/",suffix,nome)
        set.seed(seed)
        spade.downsampleFCS(infilename = fFvectorName_out[[i]], outfilename = fFvectorName[[i]], 
                            exclude_pctile = p1/100, target_percent = p2/100, target_pctile = NULL)
        nome = substr(x = fFvectorName[i], start = 11, stop = 25) #30
        nome = paste0("^", nome, "*")
        fFfiles <- list.files(path = "./tmpdata", pattern=nome, full.names = TRUE)}}
    else #case of downsampling by value
      {fFfiles <- list.files(path = "./tmpdata", pattern="^downdens_", full.names = TRUE)}}
    
    options(warn = 2) # Turn warnings into errors so they can be trapped
    result <- try(expr = read.flowSet(files = fFfiles, truncate_max_range = FALSE)) 
  
  if (!(inherits(x = result,"try-error"))){ 
    options(warn = 0)
    return(result)}
  else{
    # Ignore warnings while processing errors
    options(warn = -1)
    msg <- geterrmessage()
    result <- paste("Something went wrong with the downsampling process in the internal spade_fS function. It returns: ", 
                    msg, sep = " ")
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


outlier_score <- function(flow_Set, stringhe, algo, kappa, kappa_max = NULL, par = 1, seed){
  
  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", length = fSlength)
  pDname <- as.character(pData(parameters(flow_Set[[1]]))$name)
  pDdesc <- as.character(pData(parameters(flow_Set[[1]]))$desc)
  
  if ('score' %in% pDname)
  {return ("score channel is present in the sample's expression matrix. Please remove it before to downsampling")}
  
  name_mask <- pDdesc[!is.na(pDdesc)]
  mask <- pDdesc %in% name_mask
  pDname_subset <- pDname[mask]
  pDdesc_subset <- pDdesc[mask]
  nomi <- grep("FSC-|SSC-|TIME|DENSITY|SCORE|CELL_ID", name_mask, value = FALSE, ignore.case = TRUE)
  
  if (length(nomi)!=0){
    pDname_subset <- pDname_subset[-c(nomi)]
    pDdesc_subset <- pDdesc_subset[-c(nomi)]}

  for (i in seq_along(1:fSlength)){
    m <- exprs(flow_Set[[i]])
    m <- m[, pDname_subset]
  
    fFvectorName <- stringhe
    suffix <- paste0("_k",as.character(kappa),"_kmax",as.character(kappa_max))
    fFvectorName[[i]] <- paste0("./tmpdata/", suffix,fFvectorName[[i]]) #./tmpdata/dfF_samplei.fcs
    
    if (kappa_max > nrow(m))
      {kappa_max = nrow(m)}
    if (kappa_max < kappa)
      {kappa_max = kappa +1}
    
    set.seed(seed)
    options(warn = 2) # Turn warnings into errors so they can be trapped
    if (algo == "RKOF") {score <- try(DDoutlier::RKOF(dataset = m, k = kappa, C = 1, alpha = par))} # alpha = 1 and sigma2 = 1 deafult
    if (algo == "LOOP") {score <- try(DDoutlier::LOOP(dataset = m, k = kappa, lambda = par))} #lambda = 3 default
    if (algo == "LOF") {score <- try(DDoutlier::LOF(dataset = m, k = kappa))}
    if (algo == "LDF"){score <- try(DDoutlier::LDF(dataset = m, k = kappa, h = par, c = 1))} #h = 1 default
    #LDF A vector of Local Density Factor for observations. The greater the LDF, the greater the outlierness
    if (algo == "LDE"){score <- try(DDoutlier::LDF(dataset = m, k = kappa, h = par, c = 1))} #h = 1 default
    #LDE A vector of Local Density Estimate for observations. The greater the LDE, the greater centrality
    if (algo == "KNN_SUM") {score <- try(DDoutlier::KNN_SUM(dataset = m, k = kappa))}
    if (algo == "KNN_AGG") {score <- try(DDoutlier::KNN_AGG(dataset = m, k_min = kappa, k_max = kappa_max))}
    # Process any error messages
    
    if (!(inherits(x = score,"try-error")))
      {options(warn = 0)} # Restore default warning reporting
    else{
      # Ignore warnings while processing errors
      options(warn = -1)
      msg <- geterrmessage()
      result <- paste("Something went wrong with the downsampling process. The outlier_score internal routine returns: ", 
                      msg, sep = " ")
      # Restore default warning reporting
      options(warn=0)
      return(result)}
    
    if (algo == "LDF") {score <- score$LDF}
    if (algo == "LDE") {score <- score$LDE}
    
    # some algorithms like RKOF could compute some score values af infinite (Inf)
    max_sample_score <- score[which(!is.na(score))]
    max_sample_score <- max_sample_score[which(max_sample_score < Inf)]
    max_sample_score <- max(max_sample_score)
    score[which(is.infinite(score))] <- max_sample_score
    
    print(paste0("outlier's score of sample nr.",i," processed"))

    if ((length(score) == nrow(m))){
      newfF <- add_dim(flow_Frame = flow_Set[[i]], dim_name = "score", dim_vect = score)
      write.FCS(x = newfF, filename = fFvectorName[[i]])}}

  nome = substr(x = fFvectorName[i], start = 11, stop = 20)
  nome = paste0("^", nome, "*")
  fFfiles <- list.files(path = "./tmpdata", pattern=nome, full.names = TRUE) 
  #fFfiles <- str_sort(x = fFfiles, numeric = TRUE) #this is necessary to sort in the correct way the flowFrames in the flowSet
  newfS <- read.flowSet(files = fFfiles, truncate_max_range = FALSE)
  return(newfS)
}


downsampling_fS <- function(flow_Set, stringhe, algo, type, perc, value, seed, cut_off = NULL) {
  
  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", length = fSlength)
  
  if (type == "down_perc"){ #case of downsampling by perc
    suffix <- paste0("D",algo,"_perc",as.character(perc))}
  if (type == "down_val"){ #case of downsampling by val
    suffix <- paste0("D",algo,"_value",as.character(value))}
  if (type == "equalize"){ #case of downsampling by equilization
    suffix <- paste0("D",algo,"_equal",as.character(cut_off))}
  
    fFvectorName <- paste0("./tmpdata/", suffix, stringhe)
    
  for (i in seq_along(1:fSlength)){
    m <- exprs(flow_Set[[i]])
    m <- as.data.frame(x = m, optional = FALSE)
    nomi <- colnames(m)
    if (!(is.null(cut_off))){# This holds for equilize only since only using equilize the cut_off in not NULL
      if (nrow(m) >  cut_off + MIN_SAMPLE_LENGTH) 
      {eventi <- cut_off} 
      else {eventi <- nrow(m)}}
    else{eventi <- round(x = (100 - perc) * nrow(m) /100,digits = 0)}
    
    if (!(algo == "Random")){
      if ("score" %in% nomi){
        m_score <- m[,'score']
        max_sample_score <- m_score[which(!is.na(m_score))]
        max_sample_score <- max_sample_score[which(max_sample_score < Inf)]
        max_sample_score <- max(max_sample_score)
        m_score[which(is.infinite(m_score))] <- max_sample_score
        m[,'score'] <- m_score
        ordm <- m %>% arrange(desc(x = score)) %>% mutate(score = (score / max(score))*100)} 
      #case of DDoutlier algo in descending order: the outliers are in the head
     if ("density" %in% nomi){ #in ascending order for spade (smaller value bigger outlierness)
        ordm <- m %>% arrange(x = density) %>% mutate(density = (density / max(density))*100)} #case of spade in ascending order
      #case of spade algo in ascending order: the outliers are in the head 
      if ((type == "down_perc")||(type == "equalize")){
        #in both cases the events to choose are the ones with the higher densities (in the tail for spade) 
        #or the lower score (in the tail for the other algos)
        if (algo == "spade") {
          filtm <- tail(x = ordm, n = eventi)}
        else{
          filtm <- tail(x = ordm, n = eventi)}} 
      if (type == "down_val")
        {if ("density" %in% nomi) # case of spade
        {filtm <- ordm %>% dplyr::filter(density > value)} 
          else
         {filtm <- ordm %>% dplyr::filter(score < value)}} 
      if ("time" %in% nomi)
        {m <- filtm %>% arrange(time)}  
      if ("Time" %in% nomi)
        {m <- filtm %>% arrange(Time)}
      if ((!("Time" %in% nomi))&&(!("time" %in% nomi)))
        {m <- filtm}}
    
    else{ #this is in case of RANDOM downsampling
      set.seed(seed)
      m <- sample_n(tbl = m, size = eventi, replace = FALSE)
      if ("time" %in% nomi){
        m <- m %>% arrange(time)}  
      if ("Time" %in% nomi){
        m <- m %>% arrange(Time)}}
    
      m <- as.matrix(x = m)
      
    if (dim(m)[1] > MIN_SAMPLE_LENGTH_AFTER_DOWNSAMPLING){
      exprs(flow_Set[[i]]) <- m
      write.FCS(x = flow_Set[[i]], filename = fFvectorName[i])}
    else
    {return(paste0("System reports: The result of the downsampling process is invalid because sample nr. ", as.character(i),
                   " is too small ( less than ", as.character(MIN_SAMPLE_LENGTH_AFTER_DOWNSAMPLING), " events)"))}}
    
  nome = substr(x = fFvectorName[1], start = 11, stop = 24)
  nome = paste0("^", nome, "*")  
  fFfiles <- list.files(path = "./tmpdata", pattern=nome, full.names = TRUE)
  
  options(warn = 2) # Turn warnings into errors so they can be trapped
  result <- try(expr = read.flowSet(files = fFfiles, truncate_max_range = FALSE)) 
  
  # Process any error messages
  if (inherits(x = result,"try-error")) 
  {# Ignore warnings while processing errors
    options(warn = -1)
    # Ignore warnings while processing errors
    options(warn = -1)
    
    # If this script were a function, warning() and stop() could be called to pass errors upstream
    msg <- geterrmessage()
    messaggio <- paste0("System reports: ", msg, "... means 'something went wrong with the spade algorithm'")
    
    # Restore default warning reporting
    options(warn=0)
    return(messaggio)}
  else
    {options(warn=0)
      return(fFfiles)} #it do not return the flowSet
}


plot_score <- function(flow_Set, sample_color, type) {

  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", fSlength)
  vect_event <- data.frame()
  down_score <- data.frame()
  dens_score <- data.frame()
 
  for (i in seq_along(1:fSlength)){
    if (fSlength < 100){
        if (i < 10) 
        {fFvectorName[i] = paste0("sample_0", as.character(i))}
        else
        {fFvectorName[i] = paste0("sample_", as.character(i))}} 
      else # if fSlength >= 100
      {if (i < 10) 
        {fFvectorName[i] = paste0("sample_00", as.character(i))}
          else { 
            if (i < 100)
              {fFvectorName[i] = paste0("sample_0", as.character(i))} 
            else 
            {fFvectorName[i] = paste0("sample_", as.character(i))}}}
    vect_event <- rbind(vect_event, dim(exprs(flow_Set[[i]]))[1])}
    
  
  min_event <- min(vect_event)
  ###
  max_event <- max(vect_event)
  ###
  for (i in seq_along(1:fSlength)){
    nomi <- colnames(exprs(flow_Set[[i]]))
    if ("density" %in% nomi)
    {pos <- which(nomi == "density")
      sample_score <- exprs(flow_Set[[i]])[,pos]}
    if ("score" %in% nomi)
    {pos <- which(nomi == "score")
      sample_score <- exprs(flow_Set[[i]])[,pos]}
    if (type == "spade")
      {sample_score <- sort(x = sample_score, decreasing = FALSE)}
    else
      {sample_score <- sort(x = sample_score, decreasing = TRUE)}
    
    max_sample_score <- sample_score[which(!is.na(sample_score))]
    max_sample_score <- max_sample_score[which(max_sample_score < Inf)]
    max_sample_score <- max(max_sample_score)
    sample_score[which(is.infinite(sample_score))] <- max_sample_score
    
    if(max_sample_score > 0){
      sample_score <- round((sample_score / max_sample_score)*100, digits = 3)
      sample_score <- as.data.frame(sample_score)
      if (type == "spade")
        {names(sample_score) <- "density"}
      else
        {names(sample_score) <- "score"}
      sample_score$ID <- seq.int(nrow(sample_score))
      sample_score <- head(x = sample_score, n = round(min_event,digits = 0))
      min_event_txt <- as.character(min_event)
      sample_score$sample_id <- as.character(rep(fFvectorName[i], round(min_event, digits = 0)))
      down_score <- rbind(down_score, sample_score)}
    else
      return("There was a problem on the result of the downsampling algorithm. 
             Try another algorithm or another related parameter setting")}
  
  if (type == "spade")
  {pl_geom <- ggplot(down_score, aes(x = ID, y = density, color = sample_id, group = sample_id)) +
    geom_smooth() +
    labs(x = "number of events") +
    scale_color_manual(values = sample_color) +
    scale_x_continuous(breaks = pretty(down_score$ID, n = 10)) +
    labs(y="ordered densities", x = "number of events") +
    ggtitle(paste0("ordered densities for the first ", min_event_txt, " events"))}
  else
  {pl_geom <- ggplot(down_score, aes(x = ID, y = score, color = sample_id, group = sample_id)) +
    geom_smooth() +
    labs(x = "number of events") +
    scale_color_manual(values = sample_color) +
    scale_x_continuous(breaks = pretty(down_score$ID, n = 10)) +
    labs(y="ordered scores", x = "number of events") +
    ggtitle(paste0("ordered scores for the first ", min_event_txt, " events"))}
  
  for (i in seq_along(1:fSlength)){
    if ("density" %in% nomi)
    {pos <- which(nomi == "density")
    sample_score <- exprs(flow_Set[[i]])[,pos]}
    if ("score" %in% nomi)
    {pos <- which(nomi == "score")
    sample_score <- exprs(flow_Set[[i]])[,pos]}
   
    if (type == "spade")
      {sample_score <- sort(x = sample_score, decreasing = FALSE)}
    else
      {sample_score <- sort(x = sample_score, decreasing = TRUE)}
    
    max_sample_score <- sample_score[which(!is.na(sample_score))]
    max_sample_score <- max_sample_score[which(max_sample_score < Inf)]
    max_sample_score <- max(max_sample_score)
    sample_score[which(is.infinite(sample_score))] <- max_sample_score
    
    if(max_sample_score > 0){
    sample_score <- round((sample_score /  max_sample_score)*100, digits = 2)
    tot <- length(sample_score)
    
    sample_score <- as.data.frame(sample_score)
    if (type == "spade")
      {names(sample_score) <- "density"}
    else
      {names(sample_score) <- "score"}
    sample_score$ID <- seq.int(tot)
    sample_score$ID <- round((sample_score$ID/tot)*100, digits = 0)
    if (type == "spade")
      {sample_score <- sample_score %>% group_by(ID) %>% summarise(mean_density = mean(density, na.rm = TRUE)/100)}
    else
      {sample_score <- sample_score %>% group_by(ID) %>% summarise(mean_score = mean(score, na.rm = TRUE)/100)}
    sample_score$ID <- sample_score$ID/100
    
    righe_che_mancano <- 101 - nrow(sample_score)
    if (righe_che_mancano >0){
      righe <- tail(x = sample_score, righe_che_mancano)
      sample_score <- rbind(sample_score, righe)}
    sample_score$sample_id <- as.character(rep(fFvectorName[i], nrow(sample_score))) 
    #it was  sample_score$sample_id <- as.character(rep(fFvectorName[i], 101))
    dens_score <- rbind(dens_score, sample_score)}
    else
      return("problem on the downsampling algorithm. Try another algorithm or another related parameter setting")}
    
  if (type == "spade")
    {pl_dens <- ggplot(dens_score, aes(x = mean_density, color = sample_id, group = sample_id)) +
    #geom_freqpoly(stat = "density", size = 1) +
    geom_freqpoly(stat = "density", linewidth = 1) +#geom_freqpoly(stat = "density", size = 1) +
    scale_color_manual(values = sample_color) +
    scale_x_continuous(breaks = pretty(down_score$mean_density, n = 10)) +
    labs(y="normalized values", x = "normalized densities") +
    ggtitle(paste0("normalized densities along the samples"))}
  else
    {pl_dens <- ggplot(dens_score, aes(x = mean_score, color = sample_id, group = sample_id)) +
    geom_freqpoly(stat = "density", size = 1) +
    scale_color_manual(values = sample_color) +
    scale_x_continuous(breaks = pretty(down_score$mean_score, n = 10)) +
    labs(y="normalized values", x = "normalized scores") +
    ggtitle(paste0("normalized scores along the samples"))}
  
  p_geom <- ggplotly(p = pl_geom, originalData = TRUE)
  p_dens <- ggplotly(p = pl_dens, originalData = TRUE)
  lista <- list(down_score, p_geom, p_dens)
  return(lista)
}


edit_table <- function(flow_Set = NULL, file_name = NULL) {
  
  if (!is.null(flow_Set)){
    fSlength <- length(flow_Set)
    vectorFile <- vector(mode = "character", fSlength)
    vectorSampleid <- vector(mode = "character", fSlength)
    vectorTag1 <- vector(mode = "character", fSlength)
    vectorTag2 <- vector(mode = "character", fSlength)
    vectorTag3 <- vector(mode = "character", fSlength)
    vectorTag4 <- vector(mode = "character", fSlength)
    vectorTime <- seq(from = Sys.Date(), by = "days", length.out = fSlength)
    
    for (i in seq_along(1:fSlength)){
      if (!is.null(file_name)) {vectorFile[i] <- file_name[[i]]}
      else{
        vectorFile[i] = flow_Set@phenoData@data$name[i]}
      
      #vectorFile <- make.names(names = vectorFile, unique = TRUE) #I made unique even though Id_nr alwasy make the names different
      vectorFile <- forbidden_char(caratteri = vectorFile)
      
      if (fSlength < 100){
        if (i < 10) { 
          vectorSampleid[i] = paste0("sample_0", as.character(i))
          vectorTag1[i] = paste("patient_0", as.character(i), sep = "")
          vectorTag2[i] = paste("type_0", as.character(i), sep = "")}
        else {
          vectorSampleid[i] = paste0("sample_", as.character(i))
          vectorTag1[i] = paste("patient_", as.character(i), sep = "")
          vectorTag2[i] = paste("type_", as.character(i), sep = "")}}
      else # if fSlength >= 100
        {if (i < 10) {
          vectorSampleid[i] = paste0("sample_00", as.character(i))
          vectorTag1[i] = paste("patient_00", as.character(i), sep = "")
          vectorTag2[i] = paste("type_00", as.character(i), sep = "")}
        else { 
          if (i < 100) {
            vectorSampleid[i] = paste0("sample_0", as.character(i))
            vectorTag1[i] = paste("patient_0", as.character(i), sep = "")
            vectorTag2[i] = paste("type_0", as.character(i), sep = "")} 
          else {
            vectorSampleid[i] = paste0("sample_", as.character(i))
            vectorTag1[i] = paste("patient_", as.character(i), sep = "")
            vectorTag2[i] = paste("type_", as.character(i), sep = "")}}}
    
      if (i<27) {vectorTag3[i] = LETTERS[i]}
      if ((i>=27)&&(i<53)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[i-26])}
      if ((i>=53)&&(i<79)) {vectorTag3[i] = paste0(LETTERS[2],LETTERS[i-52])}
      if ((i>=79)&&(i<105)) {vectorTag3[i] = paste0(LETTERS[3],LETTERS[i-78])}
      if ((i>=105)&&(i<131)) {vectorTag3[i] = paste0(LETTERS[4],LETTERS[i-104])}
      if ((i>=131)&&(i<157)) {vectorTag3[i] = paste0(LETTERS[5],LETTERS[i-130])}
      if ((i>=157)&&(i<183)) {vectorTag3[i] = paste0(LETTERS[6],LETTERS[i-156])}
      if ((i>=193)&&(i<209)) {vectorTag3[i] = paste0(LETTERS[7],LETTERS[i-182])}
      if ((i>=209)&&(i<235)) {vectorTag3[i] = paste0(LETTERS[8],LETTERS[i-208])}
      if ((i>=235)&&(i<261)) {vectorTag3[i] = paste0(LETTERS[9],LETTERS[i-234])}
      if ((i>=261)&&(i<287)) {vectorTag3[i] = paste0(LETTERS[10],LETTERS[i-260])}
      if ((i>=287)&&(i<313)) {vectorTag3[i] = paste0(LETTERS[11],LETTERS[i-286])}
      if ((i>=313)&&(i<339)) {vectorTag3[i] = paste0(LETTERS[12],LETTERS[i-312])}
      if ((i>=339)&&(i<365)) {vectorTag3[i] = paste0(LETTERS[13],LETTERS[i-338])}
      if ((i>=365)&&(i<391)) {vectorTag3[i] = paste0(LETTERS[14],LETTERS[i-364])}
      if ((i>=391)&&(i<417)) {vectorTag3[i] = paste0(LETTERS[15],LETTERS[i-390])}
      if ((i>=417)&&(i<443)) {vectorTag3[i] = paste0(LETTERS[16],LETTERS[i-416])}
      if ((i>=443)&&(i<469)) {vectorTag3[i] = paste0(LETTERS[17],LETTERS[i-442])}
      if ((i>=469)&&(i<495)) {vectorTag3[i] = paste0(LETTERS[18],LETTERS[i-468])}
      if ((i>=495)&&(i<521)) {vectorTag3[i] = paste0(LETTERS[19],LETTERS[i-494])}
      if ((i>=521)&&(i<547)) {vectorTag3[i] = paste0(LETTERS[20],LETTERS[i-520])}
      if ((i>=547)&&(i<573)) {vectorTag3[i] = paste0(LETTERS[21],LETTERS[i-546])}
      if ((i>=573)&&(i<599)) {vectorTag3[i] = paste0(LETTERS[22],LETTERS[i-572])}
      if ((i>=599)&&(i<625)) {vectorTag3[i] = paste0(LETTERS[23],LETTERS[i-598])}
      if ((i>=625)&&(i<651)) {vectorTag3[i] = paste0(LETTERS[24],LETTERS[i-624])}
      if ((i>=651)&&(i<677)) {vectorTag3[i] = paste0(LETTERS[25],LETTERS[i-650])}
      if ((i>=677)&&(i<703)) {vectorTag3[i] = paste0(LETTERS[26],LETTERS[i-676])}
      if ((i>=703)&&(i<729)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[1],LETTERS[i-702])}
      if ((i>=729)&&(i<755)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[2],LETTERS[i-728])}
      if ((i>=755)&&(i<781)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[3],LETTERS[i-754])}
      if ((i>=781)&&(i<807)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[4],LETTERS[i-780])}
      if ((i>=807)&&(i<833)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[5],LETTERS[i-806])}
      if ((i>=833)&&(i<859)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[6],LETTERS[i-832])}
      if ((i>=859)&&(i<885)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[7],LETTERS[i-858])}
      if ((i>=885)&&(i<911)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[8],LETTERS[i-884])}
      if ((i>=911)&&(i<937)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[9],LETTERS[i-910])}
      if ((i>=937)&&(i<963)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[10],LETTERS[i-936])}
      if ((i>=963)&&(i<989)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[11],LETTERS[i-969])}
      if ((i>=989)&&(i<1015)) {vectorTag3[i] = paste0(LETTERS[1],LETTERS[12],LETTERS[i-988])}
      
      if (i < 11) {vectorTag4[i] <- paste0("time_0", as.character(i-1))}
      if ((i >=11)&&(i < 101)) {vectorTag4[i] <- paste0("time_", as.character(i-1))}
      if ((i >=101)&&(i < 1000)) {vectorTag4[i] <- paste0("time", as.character(i-1))}
  
      metadf <- data.frame("file_name" = vectorFile, 
                       "sample_id" = vectorSampleid, 
                       "patient_id" = factor(vectorTag1),
                       "type" = factor(vectorTag2), 
                       "condition" = factor(vectorTag3), 
                       "time_step" = factor(vectorTag4),
                       "date" = vectorTime,
                       check.names = FALSE,
                       stringsAsFactors = FALSE)}}
  else {
    metadf <- data.frame("file_name" = character(),
                        "sample_id" = character(),
                        "patient_id" = character(),
                        "type" = character(),
                        "condition" = character(),
                        "time_step" = character(),
                        "date" = character(),
                        check.names = FALSE,
                        stringsAsFactors = FALSE)}
  return(metadf)
}


edit_markers <- function(flow_Set = NULL) {
  if(!is.null(flow_Set)){

  res.expr <- exprs_sub(flow.Set = flow_Set)
  pDname_subset <- res.expr[[2]]
  pDdesc_subset <- res.expr[[3]]
  vectorSelected  <- vector(mode = "logical")
  vectorSelected <- rep(F, length(pDname_subset))
  Id <-  seq.int(length(pDname_subset))
  
  metadf <- data.frame("Id" = Id,
                       "dimension_name" = pDname_subset, 
                       "marker_description" = pDdesc_subset,
                       "selected" = vectorSelected,
                       check.names = FALSE,
                       stringsAsFactors = FALSE)}
  else {
    metadf <- data.frame("Id" = character(),
                         "dimension_name" = character(),
                         "marker_description" = character(),
                         "selected" = character(),
                         check.names = FALSE,
                         stringsAsFactors = FALSE)}
  return(metadf)
}


edit_pheno <- function(flow_Set = NULL, selected = NULL) {
  
  if(!is.null(flow_Set)){
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected)
  pDname_subset <- res.expr[[2]]
  pDdesc_subset <- res.expr[[3]]
  nr_righe <- 3
  livelli <- c("+","-","*") #this to define the possible values
  values = factor(levels = livelli, ordered = F, nmax = 3)
  
  m <- matrix(data = values, ncol = length(pDdesc_subset)+1, nrow = nr_righe)
  metadf = data.frame(m, stringsAsFactors = T, check.names = FALSE)  
  
  values <- c("*","+","+","+","-","-","-") # this to sample the values to the table
  
  for (j in 1:length(pDdesc_subset)){
    metadf[,j+1] <- factor(levels = livelli, nmax = 3)
    for (i in 1:nr_righe)
      metadf[i,j+1] <- sample(x = values, size = 1)
  }
  
  vectorName <- c("phenotype_name", pDdesc_subset)
  names(metadf) <-vectorName
  vectorName <- vector(mode = "character", nr_righe)
  
  for (i in seq_along(1:nr_righe)){
    vectorName[i] = paste0("pheno_", as.character(i))}
  
  metadf$phenotype_name<- vectorName}
  else
      metadf = data.frame("pheno matrix" = character(), stringsAsFactors = T, check.names = FALSE)  
  
  return(metadf)
}


MDSplot_new <- function(flow_Set, meta, type, selected = NULL, save = NULL, metodo) {
  
  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", fSlength)
  
  for (i in seq_along(1:fSlength)){
    if (fSlength < 100){
      if (i < 10) 
      {fFvectorName[i] = paste0("sample_0", as.character(i))}
      else
      {fFvectorName[i] = paste0("sample_", as.character(i))}} 
    else # if fSlength >= 100
    {if (i < 10) 
    {fFvectorName[i] = paste0("sample_00", as.character(i))}
      else { 
        if (i < 100)
        {fFvectorName[i] = paste0("sample_0", as.character(i))} 
        else 
        {fFvectorName[i] = paste0("sample_", as.character(i))}}}}
  
  sample_ids <- rep(fFvectorName, fsApply(flow_Set, nrow))
  
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
  expr_median_sample_tbl <- data.frame(sample_id = sample_ids, exprsub, check.names = FALSE) %>%
    group_by(sample_id) %>% 
    summarize_all(list(median = median))
  
  expr_median_sample <- t(expr_median_sample_tbl[, -1])
  colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id
  color_tag1 <- color.tag1
  color_tag2 <- color.tag2
  color_tag3 <- color.tag3
  color_tag4 <- color.tag4
  color_sample <- color.sample
  
  options(warn = 2) # Turn warnings into errors so they can be trapped
  #result <- try(plotMDS(expr_median_sample, plot = FALSE), silent = FALSE)
  
  distanza <- factoextra::get_dist(x = t(expr_median_sample), method = metodo)
  
  if (inherits(x = distanza,"dist"))
  {
    options(warn = 0)
    
    mds <- isoMDS(d = distanza, maxit = 50, trace = TRUE)
    result <- as.data.frame(mds$points)
    colnames(result) <- c("x", "y")
    plotto <- NULL
    
    ggdf <- data.frame(MDS1 = result$x, MDS2 = result$y, sample_id = colnames(expr_median_sample), check.names = FALSE)
    mm <- match(ggdf$sample_id, meta$sample_id)
    
    if (type == "tag1") {
      ggdf$tag1 <- meta[,3]
      ggdf$tag1 <- factor(x = ggdf$tag1, levels = levels(factor(x = ggdf$tag1, labels = as.character(unique(ggdf$tag1)))))
      names(ggdf)[4] <- names(meta)[3]
      
      foo <- paste0("ggplot(ggdf, aes(x = MDS1, y = MDS2, color = ", names(meta)[3], ")) +")
      foo1 <- paste0(foo, "geom_point(size = 2, alpha = 0.8) + theme_bw() + scale_color_manual(values = color_tag1) ")
      parsed <- parse(text = foo1)
      p <- eval(expr = parsed)
      if (!is.null(save)){
        ggsave(filename = "tag1_MDS.png", plot = p, device = "png", path = "./tmpdata/")}
    } else if (type == "tag2")
    {
      ggdf$tag2 <- meta[,4]
      ggdf$tag2 <- factor(x = ggdf$tag2, levels = levels(factor(x = ggdf$tag2, labels = as.character(unique(ggdf$tag2)))))
      names(ggdf)[4] <- names(meta)[4]
      
      foo <- paste0("ggplot(ggdf, aes(x = MDS1, y = MDS2, color = ", names(meta)[4], ")) +")
      foo1 <- paste0(foo, "geom_point(size = 2, alpha = 0.8) + theme_bw() + scale_color_manual(values = color_tag2) ")
      parsed <- parse(text = foo1)
      p <- eval(expr = parsed)
      if (!is.null(save)){
        ggsave(filename = "tag2_MDS.png", plot = p, device = "png", path = "./tmpdata/")}}
    else if (type == "tag3")
    {
      ggdf$tag3 <- meta[,5]
      ggdf$tag3 <- factor(x = ggdf$tag3, levels = levels(factor(x = ggdf$tag3, labels = as.character(unique(ggdf$tag3)))))
      names(ggdf)[4] <- names(meta)[5]
      
      foo <- paste0("ggplot(ggdf, aes(x = MDS1, y = MDS2, color = ", names(meta)[5], ")) +")
      foo1 <- paste0(foo, "geom_point(size = 2, alpha = 0.8) + theme_bw() + scale_color_manual(values = color_tag3) ")
      parsed <- parse(text = foo1)
      p <- eval(expr = parsed)
      if (!is.null(save)){
        ggsave(filename = "tag3_MDS.png", plot = p, device = "png", path = "./tmpdata/")}}
    else if (type == "tag4")
    {
      ggdf$tag4 <- meta[,6]
      ggdf$tag4 <- factor(x = ggdf$tag4, levels = levels(factor(x = ggdf$tag4, labels = as.character(unique(ggdf$tag4)))))
      names(ggdf)[4] <- names(meta)[6]
      
      foo <- paste0("ggplot(ggdf, aes(x = MDS1, y = MDS2, color = ", names(meta)[6], ")) +")
      foo1 <- paste0(foo, "geom_point(size = 2, alpha = 0.8) + theme_bw() + scale_color_manual(values = color_tag4) ")
      parsed <- parse(text = foo1)
      p <- eval(expr = parsed)
      if (!is.null(save)){
        ggsave(filename = "tag4_MDS.png", plot = p, device = "png", path = "./tmpdata/")}}
    else if (type == "sample_id")
    {
      ggdf$sample_id <- meta$sample_id[mm]
      p <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = sample_id)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_bw() +
        scale_color_manual(values = color_sample)
      
      plotto <- factoextra::fviz_dist(dist.obj = distanza, order = TRUE)
      
      if (!is.null(save)){
        ggsave(filename = "sample_MDS.png", plot = p, device = "png", path = "./tmpdata/")
        ggsave(filename = "HD_MDS.png", plot = plotto, device = "png", path = "./tmpdata")}}
    
    p <- ggplotly(p)
    lista <- list(p,plotto)
    return(lista)}
  else{ # Process any error messages
    # Ignore warnings while processing errors
    options(warn = -1)
    # If this script were a function, warning() and stop() could be called to pass errors upstream
    msg <- geterrmessage()
    #the complete error message is "Error: (converted from warning) Error in plotMDS.default: Only 2 columns of data: need at least 3"
    if (grepl("Only 2 columns of data: need at least 3", msg)) {
      result <- paste("Error: in plotMDS.default: Only 2 samples provided: need at least 3 to produce an MDS plot")
    } else if (grepl("error",msg)) {
      # not valid entry
      result <- "error from plotMDS"
    } else if (grepl("custom error message", msg)) {
      result <- "Could call stop() to pass an error upstream."}
    
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


diagnostic_plot <- function(flow_Set, meta, type, selected) {
  
  fSlength <- length(flow_Set)
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  sample_ids <- rep(meta$sample_id, fsApply(flow_Set, nrow))
  ggdf <- data.frame(sample_id = sample_ids, exprsub, check.names = FALSE)
  
  if ((nrow(ggdf)) > RENDER_DATA) {
    ggdf <- sample_n(tbl = ggdf, size = RENDER_DATA/5, replace = F)}
  
  ggdf <- melt(ggdf, id.var = "sample_id", value.name = "expression", variable.name = "antigen")
  mm <- match(ggdf$sample_id, meta$sample_id)
  
  color_tag1 <- color.tag1
  color_tag2 <- color.tag2
  color_tag3 <- color.tag3
  color_tag4 <- color.tag4
  
  if (type == "tag1"){
    ggdf$tag1 <- meta[,3][mm]
    ggdf$tag1 <- factor(x = ggdf$tag1, levels = levels(factor(x = ggdf$tag1, labels = as.character(unique(ggdf$tag1)))))
    
    p <- ggplot(ggdf, aes(x = expression, color = tag1, group = sample_id)) +
      geom_density() +
      facet_wrap(~ antigen, nrow = 4, scales = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = color_tag1)}
  else if (type == "tag2"){
    ggdf$tag2 <- meta[,4][mm]
    ggdf$tag2 <- factor(x = ggdf$tag2, levels = levels(factor(x = ggdf$tag2, labels = as.character(unique(ggdf$tag2)))))
    
    p <- ggplot(ggdf, aes(x = expression, color = tag2, group = sample_id)) +
      geom_density() +
      facet_wrap(~ antigen, nrow = 4, scales = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = color_tag2)}
  else if (type == "tag3"){
    ggdf$tag3 <- meta[,5][mm]
    ggdf$tag3 <- factor(x = ggdf$tag3, levels = levels(factor(x = ggdf$tag3, labels = as.character(unique(ggdf$tag3)))))
    
    p <- ggplot(ggdf, aes(x = expression, color = tag3, group = sample_id)) +
      geom_density() +
      facet_wrap(~ antigen, nrow = 4, scales = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = color_tag3)}
  else if (type == "tag4"){
    ggdf$tag4 <- meta[,6][mm]
    ggdf$tag4 <- factor(x = ggdf$tag4, levels = levels(factor(x = ggdf$tag4, labels = as.character(unique(ggdf$tag4)))))
    
    p <- ggplot(ggdf, aes(x = expression, color = tag4, group = sample_id)) +
      geom_density() +
      facet_wrap(~ antigen, nrow = 4, scales = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = color_tag4)}
    return(p)
}


evaluation_plotPCA <- function(flow_Set, seed, selected=NULL) {
 
  fSlength <- length(flow_Set)
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
  set.seed(seed)
  colnames(exprsub) <- desc
  pca_res <- prcomp(exprsub)
  ppca <- fviz_eig(pca_res)
  sum_pca <- summary(pca_res)
  
  gpc <- get_pca_var(pca_res)
  contribution <- (gpc$contrib)
  
  cos2_plot <- fviz_cos2(pca_res, choice = "var", axes = 1:2, cex.main=1.25, cex.lab=1.25, cex.axis=1.75)
  cos2_var <- fviz_pca_var(pca_res, col.var = "cos2", gradient.cols = c("#001EBB", "#E7B800", "#FC0707"), repel = TRUE, cex.main=1.25, cex.lab=1.25, cex.axis=1.75)
  PCA_importance_df <- data.frame("PC1" = c(sum_pca$importance["Cumulative Proportion", "PC1"]), 
                   "PC1+PC2" = c(sum_pca$importance["Cumulative Proportion", "PC2"]), check.names = FALSE)
  PCA_contribution_df <- data.frame(contribution, check.names = FALSE)
  PCA_contribution_df$var_names <- rownames(contribution)
  lista <- list(ppca, PCA_importance_df, PCA_contribution_df, cos2_plot, cos2_var)
  
  return(lista)
}


evaluation_plotKmeans <- function(data_Set, k_max, col_select=NULL, selected = NULL) {
  
  if (inherits(x = data_Set,"flowSet")){
    res.expr <- exprs_sub(flow.Set = data_Set, selected.markers = selected)
    exprsub <- res.expr[[1]]
    nomi <- res.expr[[2]]
    desc <- res.expr[[3]]}
  else {
    exprsub <- exprs(data_Set)
    exprsub <- subset(x = exprsub, select=col_select)}
  
  wss <- sapply(1:k_max, function(k){kmeans(exprsub, k, nstart=10,iter.max = 200, algorithm = "MacQueen")$tot.withinss})
  
#The usage of the selected algorithm is because the default one (Hartigan-Wong) gives the following error:
#Quick-TRANSfer stage steps exceeded maximum (= 4504200) probably due to the fact that some of the points (rows of ‘x’) are 
# extremely close, the algorithm may not converge in the “Quick-Transfer” stage. 
#See: https://stackoverflow.com/questions/21382681/kmeans-quick-transfer-stage-steps-exceeded-maximum
#Notice: wss means "within-cluster sum of square"
  
  wss <- as.data.frame(wss) 
  wss$ID <- seq(k_max)
  wss = wss[-1,]
  
  cplot <- ggplot(data = wss, aes(ID,wss)) + geom_point() + geom_line()+ #it was with geom_smooth()
    ggtitle("number of cluster evaluation for the concatenated fS") +
    labs(y="Within groups sum of squares", x = "Number of Clusters")
  
  lista <- list(cplot, wss)
  return(lista)
}


evaluation_Mclust <- function(flow_Set, g_max, selected = NULL) {
  
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
#d_clust <- Mclust(as.matrix(exprsub), G=1:g_max, modelNames = mclust.options("emModelNames"))
# for the models see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/
# complete list:  c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV")  
  
  models <- c("VVE", "VEV", "EVV", "VVV")
  mod_len <- length(models)
  m.bic <- mclustBIC(as.matrix(exprsub), G=1:g_max, modelNames = models)
  m.bic.df <- data.frame(matrix(ncol = length(models), nrow = g_max), check.names = FALSE)
  
  for(righe in 1:g_max) {
    for(colonne in 1:length(models)) 
    {m.bic.df[righe,colonne] <- m.bic[righe,colonne]}}
  
  m.bic.df <- m.bic.df[-1,]
  m.bic.df$n_cluster <- 2:g_max
  
  nomi <- c(models, "n_cluster")
  names(m.bic.df) <- nomi
  m.bic.df <- m.bic.df %>%  gather(`VVE`, `VEV`, `EVV`, `VVV`, key = "models", value = "BIC")
 
  pl <- ggplot(data = m.bic.df, aes(x = n_cluster, y = BIC, group = models)) +
    geom_line(aes(color=models))+
    geom_point(aes(color=models)) 
  
  m.bic.df <- m.bic.df %>% arrange(desc(n_cluster)) %>% arrange(desc(BIC)) 
  m.bic.df <- m.bic.df[1:5,]
  
  lista <- list(pl, m.bic.df)
  return(lista)
}

kmeans_clust <- function(flow_Set, selected_markers, seed, K = K, Algo) {
  
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected_markers)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
  set.seed(seed)
  
  options(warn = 2) # Turn warnings into errors so they can be trapped
  result <- try(expr = kmeans(x = exprsub, centers = K, iter.max = 100, nstart = K, algorithm = Algo))
  # add the algorithm and the iter.max and the center paramters
  #I put iter.max = 100 because of the error: Quick-TRANSfer stage steps exceeded maximum (= 4504200)
  if (inherits(x = result,"kmeans")){
    options(warn = 0)
    clust <- as.integer(result$cluster)
    codes <- data.frame(exprsub, mapping = clust, check.names = FALSE) %>% group_by(mapping) %>% 
      summarize_all(median) #[number of cluster, ]
    codes <- as.matrix(x = codes)
    
    lista <- list(clust, codes)
    return(lista)}
  else
  {# Process any error messages
   # Ignore warnings while processing errors
      options(warn = -1)
      msg <- geterrmessage()
      result <- paste("Something went wrong with the ReadInput function of the FlowSOM package. It reports: ",
                      msg, sep = " ")
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


flowSOM_clust <- function(flow_Set, selected_markers, seed) {
  
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected_markers)
  exprsub <- res.expr[[1]]
  pDname_subset <- c(unname(res.expr[[2]]))
  
  set.seed(seed)

  options(warn = 2) # Turn warnings into errors so they can be trapped
  result <- try(expr = ReadInput(flow_Set, transform = FALSE, scale = F))
  if (inherits(x = result,"FlowSOM")){
    options(warn = 0)
    fsom <- result 
    set.seed(seed)
    som <- BuildSOM(fsom, colsToUse = pDname_subset, silent = F)
    set.seed(seed)
    mst <- BuildMST(fsom = som, tSNE = TRUE, silent = F)
    lista <- list(fsom, som, mst)
    return(lista)}
  else
    {# Process any error messages
      if (inherits(result) == "try-error") {
        # Ignore warnings while processing errors
        options(warn = -1)
          
        # If this script were a function, warning() and stop()
        # could be called to pass errors upstream
        msg <- geterrmessage()
        result <- paste("Something went wrong with the ReadInput function of the FlowSOM package. It reports: ",
                          msg, sep = " ")} 
        
        # Restore default warning reporting
        options(warn=0)
        return(result)}
}


phenograph_clust <- function(flow_Set, selected_markers, seed, K = num_clust, implementation) {
  
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected_markers)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
  set.seed(seed)
 
  options(warn = 2) # It needs to not ignore warnings into because sometime it generates and warnings but the result is fine but length(clust)!=nrow(exprsub)
  if (implementation == "Rphenograph"){
    result <- try(expr = Rphenograph::Rphenograph(data = exprsub, k = K))} # case of Rphenograph
  if (implementation == "FastPG"){
    result <- try(expr = FastPG::fastCluster(data = exprsub, k = K, verbose = TRUE, progress = 'bar'))} # case of FastPG
  
  if (inherits(x = result,"list")){
    options(warn = 0)
    if (implementation == "Rphenograph"){
      clust <- as.integer(result[[2]]$membership)} #case Rphenograph
    if (implementation == "FastPG"){
      clust <- as.integer(result$communities + 1)} #case FastPG
    
    if (length(clust)==nrow(exprsub)){
      codes <- data.frame(exprsub, mapping = clust, check.names = FALSE) %>% group_by(mapping) %>% 
        summarize_all(median) #[number of cluster, ]
      codes <- as.matrix(codes)
      if (implementation == "Rphenograph"){modularity <- round(x = result[[2]]$modularity[[3]], digits = 3)} #for Rphenograph
      if (implementation == "FastPG"){modularity <- round(x = result$modularity, digits = 3)} #for FastPG
      lista <- list(clust, codes, modularity)
      return(lista)} 
    else
    {return(paste("Something went wrong with the selected phenograph implementation. Try another value for K"))}}
  else
  {# Process any error messages
    if (inherits(x = result, "try-error")) {
      # Ignore warnings while processing errors
      options(warn = -1)
      
      # If this script were a function, warning() and stop()
      # could be called to pass errors upstream
      msg <- geterrmessage()
      result <- paste("Something went wrong with the phenograph algorithm. It reports: ",
                      msg, sep = " ")} 
    
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


DepecheR_clust <- function(flow_Set, selected_markers, K, seed){
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected_markers)
  exprsub <- res.expr[[1]]
  nomi <- res.expr[[2]]
  desc <- res.expr[[3]]
  
  set.seed(seed)
  sSs <-  seq_len(nrow(exprsub))
  #penalties <- 2^seq(0, 3, by = 0.5)
  num_core <-  detectCores() 
  
  temp_dir = "./DepecheR"
  if (dir.exists(temp_dir)) {unlink(x = temp_dir, recursive = TRUE, force = TRUE)}
  options(warn = 0) # too many warnings (it was 2 to turn warnings into errors so they can be trapped)
  
  result <- try(expr = depeche(inDataFrame = exprsub, samplingSubset = sSs, sampleSize = "default",
                               selectionSampleSize = "default", k = K, minARIImprovement = 0.01, optimARI = 0.95, maxIter = 100, 
                               log2Off = FALSE, center = "peak", scale = TRUE, nCores = num_core, plotDir = temp_dir))
  
  # All the parameters are set to default except for iter.max (it was 50) and addNoise (was)
  
  if (inherits(x = result,"list")){
    options(warn = 0)
    clust <- as.integer(result$clusterVector)
    if (length(clust)==nrow(exprsub)){
      codes <- data.frame(exprsub, mapping = clust, check.names = FALSE) %>% group_by(mapping) %>% 
        summarize_all(median) #[number of cluster, ]
      
      #tSNE <- Rtsne(X = exprsub, pca=FALSE, num_threads = 0) #all the rest is set to default 
      # This is commented until I implement the latest steps of https://bioconductor.org/packages/release/bioc/vignettes/DepecheR/inst/doc/DepecheR_test.html
      #remember to add dI.map(mapres) in the mapping reactive procedure in the server.R
      tSNE <- NULL
      codes <- as.matrix(codes)
      lista <- list(clust, codes, tSNE)
      return(lista)} else
      {return(paste("Something went wrong with the ReadInput function of the Phenograph package. Try another value for K"))}}
  else
  {# Process any error messages
    if (inherits(x = result,"try-error")) {
      # Ignore warnings while processing errors
      options(warn = -1)
      
      # If this script were a function, warning() and stop()
      # could be called to pass errors upstream
      msg <- geterrmessage()
      result <- paste("Something went wrong with the phenograph algorithm. It reports: ",
                      msg, sep = " ")} 
    
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


flowSOM_plot <- function(mst, markers_df){

  markers_df <- transform(markers_df, dimension_name = as.character(dimension_name)) 
  #otherwise you have factors with $selected
  markers_df <- transform(markers_df, marker_description = as.character(marker_description))
  marker_name <- markers_df[markers_df$selected==TRUE,]$dimension_name
  marker_desc <- markers_df[markers_df$selected==TRUE,]$marker_description
  lunghezza <- length(marker_desc)
  vectorMarker <- vector(mode = "character",length = lunghezza)
  vectorTitle <- vector(mode = "character",length = lunghezza)
  marker_desc <- make.names(marker_desc)
  
  for (i in seq_along(1:lunghezza)){
    vectorMarker[i] = paste0("./tmpdata/PlotMarker", marker_desc[[i]],".png")
    vectorTitle[i] = paste0("flowSOM MST for ", marker_desc[[i]]," marker")}
  
  for (i in seq_along(1:lunghezza)){
    png(filename = vectorMarker[[i]])
    print(PlotMarker(fsom = mst, marker = marker_name[[i]], title = vectorTitle[[i]], refMarkers = mst$map$colsUsed[i]))
    dev.off()}
  
  listPlot <- lapply(1:lunghezza, function(i) readPNG(source = vectorMarker[[i]], native = FALSE))
  #listPlot <- lapply(1:lunghezza, function(i) PlotMarker(fsom = mst, marker = marker_name[[i]], title = vectorTitle[[i]], refMarkers = mst$map$colsUsed[i])) 
  
  png(filename = './tmpdata/PlotMarker.png', width = 1000, height = 2400, units = "px", res = 125)
  listGrob <- lapply(1:lunghezza, function(i) rasterGrob(listPlot[[i]]))
  exit_plot <- grid.arrange(grobs = listGrob,ncol=2)
  dev.off()
  
  png(filename = './tmpdata/PlotNumbers.png', width = 1200, height = 1200, units = "px", res = 125)
  print(PlotNumbers(fsom = mst, title = "Cluster Numbers in MST"))
  dev.off()
  
  png(filename = './tmpdata/PlotGrid.png', width = 1200, height = 1200, units = "px", res = 125)
  #UpdateNodeSize(fsom = mst, reset=TRUE) Warning message: 'UpdateNodeSize' is deprecated. Use 'the equalNodeSize and maxNodeSize parameters in the plot functions' instead.
  
  print(PlotStars(fsom = mst, markers = marker_name, title = "markers in grid layout", view = "grid"))
  dev.off()
  
  png(filename = './tmpdata/PlotStars.png', width = 1200, height = 1200, units = "px", res = 125)
  #UpdateNodeSize(fsom = mst, reset=TRUE)
  print(PlotStars(fsom = mst, markers = marker_name, title = "markers in MST layout"))
  dev.off()
  
  listPlot[[1]] <- readPNG(source = './tmpdata/PlotNumbers.png', native = FALSE)
  listPlot[[2]] <- readPNG(source = './tmpdata/PlotGrid.png', native = FALSE)
  listPlot[[3]] <- readPNG(source = './tmpdata/PlotStars.png', native = FALSE)
  
  listGrob <- lapply(1:3, function(i) rasterGrob(listPlot[[i]]))
  
  png(filename = './tmpdata/PlotComp.png', width = 1000, height = 1200, units = "px", res = 250)
  listGrob <- lapply(1:3, function(i) rasterGrob(listPlot[[i]]))
  exit_plot <- grid.arrange(grobs = listGrob,ncol=2)
  dev.off()
}


meta_clustering <- function(codes, maxK, reps, pItem, clusterAlg, distance, seed, savedir, save) {
## Metaclustering into a fixed amount of meta clusters (NMC) with ConsensusClusterPlus
  
  options(warn = 2) # Turn warnings into errors so they can be trapped
  if (save==TRUE){
  result <- try(expr = ConsensusClusterPlus(d = t(codes), maxK = maxK, reps = reps, pItem = pItem, clusterAlg = clusterAlg, 
                             distance = distance, seed = seed, title = savedir, plot = "png", writeTable = TRUE))}
  else {result <- try(expr = ConsensusClusterPlus(t(codes), maxK = maxK, reps = reps, pItem = pItem, clusterAlg = clusterAlg, 
                                              distance = distance, seed = seed, plot = "png", writeTable = FALSE))}
  if(length(result)==maxK){
    mc <- result
    options(warn = 0)
    return(mc)}
  
  if (inherits(x = result,"try-error")) {
    # Ignore warnings while processing errors
    options(warn = -1)
    
    # If this script were a function, warning() and stop() could be called to pass errors upstream
    msg <- geterrmessage()
    result <- paste("Something went wrong with the ConsensusClusterPlus function of the homonymous package. It reports: ",
                    msg, sep = " ")
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


cluster_eva <- function(fF ,clustering, color, distance, power=NULL, setseed, type = NULL){

  exprs_fF <- exprs(fF)
  exprs_fF <- cbind(exprs_fF, clustering)
  maxCell <- 25000#MAX_EVA #con 30000 si blocca
  set.seed(setseed)
  if (maxCell < nrow(exprs_fF)){
    inds_reduced <- sample.int(n = nrow(exprs_fF), size = maxCell)
    exprs_reduced <- exprs_fF[inds_reduced,]}
  else
    {exprs_reduced <- exprs_fF}
  clustering_reduced <- exprs_reduced[,"clustering"]
  
  nomi_proibiti <- grep("FSC-|SSC-|cell_Id|density|score|Time|SampleID|clustering", colnames(exprs_fF), value = FALSE, ignore.case = TRUE) 
  if (length(nomi_proibiti)>0) {exprs_reduced <- exprs_reduced[,-nomi_proibiti]}
  
  exprs_scaled <-  scale.default(exprs_reduced)
  
  # which is the same of:
  # for (i in 1:dim(exprs_reduced)[[2]]){
  #  exprs_scaled[,i] <- (exprs_reduced[,i] - mean(exprs_reduced[,i]))/sd(exprs_reduced[,i])}
  set.seed(setseed)
  
  distance_matrix <- dist(x = exprs_scaled, method = distance, p = power)
  
  exprs_scaled <- cbind(exprs_scaled, clustering_reduced)
  set.seed(setseed)
  sil <- cluster::silhouette(x = as.integer(exprs_scaled[,"clustering_reduced"]), dist = distance_matrix)
  colori <- color[1:length(unique(clustering_reduced))]
  
  sum_view <- summary(sil)
  width <- as.data.frame(sum_view[2]$clus.avg.widths)
  colnames(width) <- "silhouette average width"
  size <- as.data.frame(sum_view[3]$clus.sizes)
  size$Freq <- (size$Freq / sum(size$Freq))*100
  colnames(size) <- c("cluster","evaluated freq")
  sum_df <- data.frame(size, width)
  # in alternative to https://github.com/kassambara/factoextra
  #plotSil <- function (sil.obj, label = FALSE, print.summary = TRUE, colori = color, ...){
  plotSil <- function (sil.obj, label = FALSE, print.summary = TRUE, colori = color, ...){
    
    df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE)
    df <- df[order(df$cluster, -df$sil_width), ]
    if (!is.null(rownames(df))) 
      df$name <- factor(rownames(df), levels = rownames(df))
    else df$name <- as.factor(1:nrow(df))
    df$cluster <- as.factor(df$cluster)
    
    colori <- setNames(object = colori, nm = 1:length(colori))
    mapping <- aes_string(x = "name", y = "sil_width", color = "cluster", fill = "cluster")
    p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
           scale_color_manual(values = colori) + 
           #see https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
           #see https://ggplot2.tidyverse.org/reference/scale_manual.html
           scale_fill_manual(values = colori) +
           labs(y = "Silhouette width Si", x = "", title = paste0("Clusters silhouette plot ", "\n Average silhouette width: ", 
                                                                  round(mean(df$sil_width), 2)), " Distance criteria: ", distance) + 
           ggplot2::ylim(c(NA, 1)) + geom_hline(yintercept = mean(df$sil_width), linetype = "dashed", color = "red")
    
    p <- ggpubr::ggpar(p, ...)
    if (!label) {p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())}
    else if (label) {
      p <- p + theme(axis.text.x = element_text(angle = 45))}
    
    if (type == "clust"){  
    ggsave(filename = "silhouette_clust.png", plot = p, device = "png", path = "./tmpdata/")}else{
      ggsave(filename = "silhouette_meta.png", plot = p, device = "png", path = "./tmpdata/")}
    #p <- ggplotly(p) #was not possibile because it requires too much memory
    ave <- tapply(df$sil_width, df$cluster, mean)
    n <- tapply(df$cluster, df$cluster, length)
    sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 2), stringsAsFactors = TRUE)
    percentuali <- sum_df$evaluated.freq
    sil.sum$perc <- paste(round(percentuali, 2), "%", sep="")
    sil.sum <- sil.sum[,c(1,3,2,4)] #to change the column order
    lista <- list(p, sil.sum)
    return(lista)}
  
  gg_sil <- plotSil(sil.obj = sil, label = FALSE, print.summary = FALSE, colori = colori)
  pl <- gg_sil[[1]]
  df_summary <- gg_sil[[2]]
  if (type == "clust"){
  write.csv(x = df_summary, file = "./tmpdata/silhouette_clust.csv")} else
    {write.csv(x = df_summary, file = "./tmpdata/silhouette_meta.csv")}
  lista <- list(pl, df_summary)
  return(lista)
}
  

plot_clustering_heatmap_wrapper <- function(flow_Set, type_HM, clustering, selected_markers, cluster_labelling = NULL, 
                                            pheno_table = NULL, Nclust, central, min_quantile, max_quantile, quant,
                                            pos_threshold, neg_threshold, dist_method, hclust_method, color_clusters, 
                                            color_label = NULL, work = NULL, cluster.signature = NULL) {

  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected_markers)
  exprsub <- res.expr[[1]]
  pDname_subset <- res.expr[[2]]
  pDdesc_subset <- res.expr[[3]]
  
  if ('cell_Id' %in% colnames(exprsub)){
    posizione <- which(colnames(exprsub) == "cell_Id")
    exprsub <- exprsub[,-posizione]}
  if ('cell_Id' %in% pDdesc_subset){
    posizione <- which(pDdesc_subset == "cell_Id")
    pDdesc_subset <- pDdesc_subset[-posizione]}
  
  colnames(exprsub) <- pDdesc_subset
  #browser()
  normalize.function <- function(x, quantiles=quant) { #used to normalize from 0 to 1 the expression range
    if (!quant){
      rng <- colQuantiles(x = exprsub, probs = c(min_quantile, max_quantile), drop = TRUE) #[nr of markers x min_quantile [,1] , max quitile [,2]]
      expr01 <- t((t(exprsub) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    }else{
      expr01 <- (x - min(x)) / (max(x) - min(x))}
    expr01[expr01 < 0] <- 0
    expr01[expr01 > 1] <- 1
    return (expr01)}
  
  expr01 <- normalize.function(exprsub, quant)
  
  # Calculate the central tendency expression
  if ((central == "median")||(central == "mean")){
    expr01_median <- data.frame(expr01, cell_clustering = clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% 
    summarize_all(.funs = central)} #[number of cluster, ]
  
  if (central == "mode"){
    my_mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]}
    
    expr01_median <- data.frame(expr01, cell_clustering = clustering, check.names = FALSE) %>%
      group_by(cell_clustering) %>% 
      summarize_all(.funs = my_mode)} #[number of cluster, ]
  
      # Calculate cluster frequencies
  clustering_table <- as.numeric(table(clustering)) #[nr of cluster]
  
  cell_clustering <- expr01_median$cell_clustering
  expr01_median$cell_clustering <- NULL #before to calculate the distance with dist let' remove the cell_clustering column
  
#  if (!is.null(cluster.signature)) {
#    expr01_median <- as_tibble(cluster.signature)} 
  #in case input$signature_finding_method=="Densities" - Please notice than in this case the heatmap computation is not based to the central 
  # (typically the "median") values but on the peaks
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr01_median, method = dist_method) #[(n * (n-1))/ 2] where n is the number of clusters
  
  options(warn = 2) # Turn warnings into errors so they can be trapped
  hclust_result <- try(expr = hclust(d, method = hclust_method)) # a list of 7 objects -> for the pheatmap function
  if (inherits(x = hclust_result,"hclust")){
    options(warn = 0)
    expr_heat <- as.matrix(expr01_median[, c(colnames(expr01))])
    rownames(expr_heat) <- cell_clustering #[nr of cluster x markers]
    expr01_median$cell_clustering <- cell_clustering
   
    # Row annotation for the heatmap
    annotation_row <- data.frame(cluster = expr01_median$cell_clustering, stringsAsFactors = FALSE) # a factor [nr clusters x 1]
    if (Nclust > nrow(expr_heat)){
      nomi_righe <- rownames(x=expr_heat, do.NULL = FALSE)
      if (NA %in% nomi_righe){
        name_mask <- is.na(nomi_righe)
        posizione <- which(x = name_mask)
        nomi_righe[posizione] <- "UNASSIGNED"} #sometimes not all the clusters are assigned -> additional NA (not assigned)
      
      annotation_row <- data.frame(cluster = nomi_righe, stringsAsFactors = FALSE)} #dataframe [nr.clusters x 1] tipically is < Nclust
    #this because sometime the expr_heat could have a number of rows which is less than expected (the Nclust)
    
    mutate_func1 <- function(var) {ifelse(var >= pos_threshold, "+", (ifelse(var < neg_threshold, "-", NA)))}
    if(!is.null(cluster_labelling)){
      cluster_labelling$new_cluster <- factor(cluster_labelling$new_cluster) 
      #cluster_labelling is a dataframe of [nr.of clusters x 2 variables: original_cluster and new_cluster] 
      annotation_row$cluster_labelling <- cluster_labelling$new_cluster 
      
      df <- as.data.frame(expr_heat)
      df_log <- df %>% mutate_all(.funs = mutate_func1) 
      df_log <- df_log[!duplicated(df_log), ]
      numero_righe <- nrow(x = df_log)
      phenotype_name <- rep("phenotype_", numero_righe)
      for (i in seq_along(1:numero_righe)){if (i < 10) {phenotype_name[i] = paste0("phenotype_0", as.character(i))}
        else 
        {phenotype_name[i] = paste0("phenotype_", as.character(i))}}
      
      df_log <- cbind(phenotype_name, df_log)
      for(i in seq_along(colnames(expr_heat))){df_log[order(df_log[,length(colnames(expr_heat))+2-i]),]} #this is to order by all the columns
      df_log$cluster <- rownames(df_log) #new column
      if (is.null(work)){
        file_name <- paste0("./tmpdata/pheno_", type_HM ,".csv")}
      else {file_name <- paste0("./tmpdata/pheno_", type_HM ,"_wf",as.character(work),".csv")}
      write.csv(x = df_log, file = file_name)
      
      df$cluster <- annotation_row$cluster
      df$cluster_label <- annotation_row$cluster_labelling
      df$perc <- round(clustering_table / sum(clustering_table) * 100,2)
      df$value <- clustering_table #new column
      if (is.null(work)){
      file_name <- paste0("./tmpdata/HM_", type_HM ,".csv")}
      else {file_name <- paste0("./tmpdata/HM_", type_HM ,"_wf",as.character(work),".csv")}
      write.csv(x = df, file = file_name)}
    else
    {
      df <- as.data.frame(expr_heat)
      df_log <- df %>% mutate_all(.funs = mutate_func1) 
      df_log <- df_log[!duplicated(df_log), ]
      numero_righe <- nrow(x = df_log)
      phenotype_name <- rep("phenotype_", numero_righe)
      for (i in seq_along(1:numero_righe)){if (i < 10) {phenotype_name[i] = paste0("phenotype_0", as.character(i))}
        else 
        {phenotype_name[i] = paste0("phenotype_", as.character(i))}}
      df_log <- cbind(phenotype_name, df_log)
      for(i in seq_along(colnames(expr_heat))){df_log <- df_log[order(df_log[,length(colnames(expr_heat))+2-i]),]} 
      #this was to order by all the columns
      df_log$cluster <- rownames(df_log) #new column
      if (is.null(work)){
      file_name <- paste0("./tmpdata/pheno_", type_HM ,".csv")}
      else {file_name <- paste0("./tmpdata/pheno_", type_HM ,"_wf",as.character(work),".csv")}
      write.csv(x = df_log, file = file_name)
      #df$cluster <- annotation_row$cluster
      df$cluster <-  expr01_median$cell_clustering
      df$perc <- round(clustering_table / sum(clustering_table) * 100,2)
      df$value <- clustering_table #new column
      if (is.null(work)){
      file_name <- paste0("./tmpdata/HM_", type_HM ,".csv")}
      else {file_name <- paste0("./tmpdata/HM_", type_HM ,"_wf",as.character(work),".csv")}
      write.csv(x = df, file = file_name)}
    #This part is to prepare all the labels for the heatmap ###    
    if (is.null(cluster_labelling)){
      labels_row <- paste0(rownames(expr_heat), " (",round(clustering_table / sum(clustering_table) * 100, 2), "%)")} #[nr of cluster]
    else
    {labels_row <- paste0(rownames(expr_heat), " (",round(clustering_table / sum(clustering_table) * 100, 2), "%) ", 
                          cluster_labelling$new_cluster)} #[nr of cluster]
    
    labels_col <- colnames(expr_heat) #[selected markers]
    if ((!(type_HM=="pheno"))&&(!(type_HM=="full_label_plus"))){
      colori <- color_clusters[1:length(annotation_row$cluster)] 
      df <- as.data.frame(expr_heat, stringsAsFactors = FALSE)
      df$value <- clustering_table 
      df$colori <- colori 
      df <- df[order(df$value, decreasing = TRUE),]
      color_clusters <- df$colori}
    #this was inserted to arrange the color on the clusters by "value". In this way you can try to compare different outcomes but the colors 
    #are always based on the same arrangements based on the bigger clusters
    
    
    temp_expr_heat <- expr_heat
    if(!(type_HM == "pheno")){
    names(color_clusters) <- annotation_row$cluster
    color_clusters <- color_clusters[1:length(annotation_row$cluster)]}#[nr. of clusters]}
    else{
      bool_vect <- names(color_clusters) %in% annotation_row$cluster
      color_clusters <- color_clusters[bool_vect]}
   
    temp_color_clusters <- color_clusters
    temp_annotation_row <- annotation_row
    
    #this for cycle only repairs the numbering: no holes at the end of this procedure. This because the heatmap does not behave well when there
    # are holes in the exprs_heat, the annotation_row and the color_clusters: the real ones are the exprs_heat, the annotation_row and     
    # the color_clusters, but the heatmap is built with these fake temp_ entries
    #if (!all(is.character(names(temp_color_clusters)))){
      #for(i in seq_along(1:length(temp_color_clusters))){
      #  if(!(names(temp_color_clusters[i]) == as.character(i))){
      #    names(temp_color_clusters)[i] <- as.character(i)
      #    temp_annotation_row$cluster[i] <- i
      #    row.names(temp_expr_heat)[i] <- as.character(i)}}
  #}
    
    #if(type_HM == "pheno"){
    #  for(i in seq_along(1:length(temp_color_clusters))){
    #    if(!(names(temp_color_clusters[i]) == as.character(i))){
    #      pos <- which(temp_annotation_row$cluster == names(temp_color_clusters[i]))
    #      names(temp_color_clusters)[i] <- as.character(pos)}}
    #      temp_annotation_row$cluster <- names(temp_color_clusters)
    #      row.names(temp_expr_heat) <- names(temp_color_clusters)}
    
    #We use exrp_heat as the real data while temp_expr_heat is just a fake matrix with the row without holes
    if(!is.null(cluster_labelling)){
      if(type_HM == "meta"){
        annotation_colors <- list(cluster = temp_color_clusters) #a list with a [nr of cluster] vector
        annotation_legend <- FALSE}
      if(type_HM == "pheno"){
        annotation_colors <- list(cluster = temp_color_clusters) #a list with a [nr of cluster] vector
        annotation_colors$cluster_labelling <- color_clusters
        pos <- which(names(annotation_colors$cluster_labelling)=="UNASSIGNED")
        annotation_colors$cluster_labelling[pos] <- '#D3D3D3' #lightgrey
        pos <- which(names(annotation_colors$cluster_labelling)=="OVERLAP")
        annotation_colors$cluster_labelling[pos] <- '#000000' #black
        
        cl <- annotation_colors$cluster_labelling
        cl <- cl[sort(names(cl))]
        annotation_colors$cluster_labelling <- cl
        temp_annotation_row$cluster <- sort(temp_annotation_row$cluster)
        
        annotation_legend <- FALSE}
      if(type_HM == "label"){
        annotation_colors <- list(cluster = temp_color_clusters, cluster_labelling = color_label)
        pos <- which(names(annotation_colors$cluster_labelling)=="UNASSIGNED")
        annotation_colors$cluster_labelling[pos] <- '#D3D3D3' #lightgrey
        pos <- which(names(annotation_colors$cluster_labelling)=="OVERLAP")
        annotation_colors$cluster_labelling[pos] <- '#000000' #black
        annotation_legend <- FALSE}
      if(type_HM == "full_label"){
        annotation_colors <- list(cluster = temp_color_clusters) #a list with a [nr of cluster] vector
        annotation_colors$cluster_labelling <- color_clusters
        annotation_legend <- TRUE}
      if(type_HM == "full_label_plus"){
        annotation_colors <- list(cluster = temp_color_clusters, cluster_labelling = color_label)
        pos <- which(names(annotation_colors$cluster_labelling)=="UNASSIGNED")
        annotation_colors$cluster_labelling[pos] <- '#D3D3D3' #lightgrey
        pos <- which(names(annotation_colors$cluster_labelling)=="OVERLAP")
        annotation_colors$cluster_labelling[pos] <- '#000000' #black
        annotation_legend <- TRUE}}
    else{
      annotation_colors <- list(cluster = temp_color_clusters) #a list with a [nr of cluster] vector
      annotation_legend <- FALSE
      if(type_HM == "pheno"){
        annotation_colors <- list(cluster = temp_color_clusters) #a list with a [nr of cluster] vector
        annotation_colors$cluster_labelling <- color_clusters
        pos <- which(names(annotation_colors$cluster_labelling)=="UNASSIGNED")
        annotation_colors$cluster_labelling[pos] <- '#D3D3D3' #lightgrey
        pos <- which(names(annotation_colors$cluster_labelling)=="OVERLAP")
        annotation_colors$cluster_labelling[pos] <- '#000000' #black
        
        rownames(temp_expr_heat) <- rownames(annotation_row)
        annotation_legend <- FALSE}}
    
    color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100) #is the array of the 100 colors used in the heatmap
    #heatmaply(x = expr_heat, colors = color, plot_method = "plotly", 
    #xlab = "markers", ylab = "clusters", main = "meta-clustering heatmap", 
    #               hide_colorbar = F, k_row = Nclust, groupLabels=T, 
    #               row_side_colors = annotation_row, labRow = labels_row, width = 900, height = 1200)
    if (is.null(work)){file_png_name = paste0("./tmpdata/HM_", type_HM ,".png")}
    else {file_png_name = paste0("./tmpdata/HM_", type_HM ,"_wf",as.character(work),".png")}
    
    if(type_HM=="meta"){png(filename = file_png_name, width = 1400, height = 1600, units = "px", res = 125)}
    if(type_HM=="label"){png(filename = file_png_name, width = 1400, height = 1600, units = "px", res = 125)}
    if(type_HM=="full_label"){png(filename = file_png_name, width = 1600, height = 4800, units = "px", res = 125)}
    if(type_HM=="full_label_plus"){png(filename = file_png_name, width = 1600, height = 4800, units = "px", res = 125)}
    if(type_HM=="pheno"){png(filename = file_png_name, width = 1300, height = 1000, units = "px", res = 125)}
    
    if(type_HM=="full_label")
    {pheatmap(mat = expr_heat, color = color, cluster_cols = FALSE, cluster_rows = hclust_result, 
              labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", 
              fontsize = 11, fontsize_number = 9, treeheight_col = 150, treeheight_row = 150,
              annotation_row = NULL, annotation_colors = NULL, 
              annotation_legend = annotation_legend)
      dev.off()}
    
    if(type_HM=="full_label_plus")
    {pheatmap(mat = temp_expr_heat, color = color, cluster_cols = FALSE, cluster_rows = hclust_result, 
              labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", 
              fontsize = 11, fontsize_number = 9, treeheight_col = 150, treeheight_row = 150,
              annotation_row = temp_annotation_row, annotation_colors = annotation_colors, 
              annotation_legend = annotation_legend)
      dev.off()}
    
    if(type_HM=="label"||type_HM=="meta"||type_HM=="pheno")
    {pheatmap(mat = temp_expr_heat, color = color, cluster_cols = FALSE, cluster_rows = hclust_result, 
              labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", 
              fontsize = 12, fontsize_number = 8, treeheight_col = 150, treeheight_row = 150,
              annotation_row = temp_annotation_row, annotation_colors = annotation_colors, 
              annotation_legend = annotation_legend)
      dev.off()}
    
    expr_heat <- as.data.frame(expr_heat)
    
    if (!is.null(cluster_labelling)){
      lista <- list(expr_heat = expr_heat, labels_row = labels_row, color_clusters = color_clusters, 
                    color_label = annotation_colors$cluster_labelling)} 
    else
      {lista <- list(expr_heat = expr_heat, labels_row = labels_row, color_clusters = color_clusters)}
    
    return(lista)}
  # Process any error messages
  if (inherits(x = result,"try-error")) {
    options(warn = -1)  # Ignore warnings while processing errors
    msg <- geterrmessage()
    result <- paste("Something went wrong with the hclust function of the stats package. It reports: ",
                    msg, sep = " ")
    # Restore default warning reporting
    options(warn=0)
    return(result)}
}


map_gen <- function(flow_Set, type, two_three = NULL, selected_markers, seed_nr, limit, maxCell,
                    #tSNE parameters
                    Rtsne_pca = NULL, Rtsne_perp = NULL, Rtsne_theta = NULL, Rtsne_iter = NULL, Rtsne_eta = NULL,
                    #UMAP paramters
                    UMAP.n_neighbors = NULL, UMAP.metric = NULL, UMAP.min_dist = NULL, UMAP.spread = NULL, UMAP.random_state = NULL){

  fSlength <- length(flow_Set)
  fFvectorName <- vector(mode = "character", fSlength)
  
  for (i in seq_along(1:fSlength)){
    if (fSlength < 100){
      if (i < 10) 
      {fFvectorName[i] = paste0("sample_0", as.character(i))}
      else
      {fFvectorName[i] = paste0("sample_", as.character(i))}} 
    else # if fSlength >= 100
    {if (i < 10) 
      {fFvectorName[i] = paste0("sample_00", as.character(i))}
        else { 
          if (i < 100)
          {fFvectorName[i] = paste0("sample_0", as.character(i))} 
        else 
          {fFvectorName[i] = paste0("sample_", as.character(i))}}}}
  
  sample_ids <- rep(fFvectorName, fsApply(flow_Set, nrow)) # [#nr. of cells] vector of char "sample_x" with x = num of flowFrame
  expr <- fsApply(flow_Set, exprs) 
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected_markers) 
  exprsub <- res.expr[[1]] #exprs of flowSet with columns of selected markers only
  
  colnames(exprsub) <- res.expr[[3]]
  
  inds <- split(1:length(sample_ids), sample_ids)
  map_ncells <- table(sample_ids)
  map_inds <- unlist(inds)
  map_expr <- exprsub
  
  inds_reduced <- map_inds
  map_expr_reduced <- map_expr
  
  set.seed(seed_nr)
  
  if (limit){
    if (maxCell < nrow(map_expr)){
      inds_reduced <- sort(sample.int(n = nrow(map_expr), size = maxCell))
      map_expr_reduced <- map_expr[inds_reduced,]}}
  
  options(warn = 2) # Turn warnings into errors so they can be trapped
  if(type == "tSNE"){
    options(warn = 0)
    if ((two_three == FALSE)||is.null(two_three))
      {map_out <- try(expr = Rtsne(map_expr_reduced, dims = 2, check_duplicates = FALSE, pca = Rtsne_pca, perplexity = Rtsne_perp, 
                          theta = Rtsne_theta, max_iter = Rtsne_iter, eta = Rtsne_eta, verbose = TRUE, num_threads = 0),
                     silent = FALSE)}
    else
        {map_out <- try(expr = Rtsne(map_expr_reduced, dims = 3, check_duplicates = FALSE, pca = Rtsne_pca, perplexity = Rtsne_perp, 
                                                 theta = Rtsne_theta, max_iter = Rtsne_iter, eta = Rtsne_eta, verbose = TRUE, num_threads = 0), 
                                    silent = FALSE)}
      #if (!(class(map_out)=="try-error"))
    if (!(inherits(x = map_out,"try-error")))
    { 
      if ((two_three == FALSE)||is.null(two_three)){
        dr <- data.frame(exprsub, check.names = FALSE)
        dr$map_x <- NA
        dr$map_y <- NA
        dr$map_x[inds_reduced] <- map_out$Y[,1]
        dr$map_y[inds_reduced] <- map_out$Y[,2]}
      else{
        dr <- data.frame(exprsub, check.names = FALSE)
        dr$map_x <- NA
        dr$map_y <- NA
        dr$map_z <- NA
        dr$map_x[inds_reduced] <- map_out$Y[,1]
        dr$map_y[inds_reduced] <- map_out$Y[,2]
        dr$map_z[inds_reduced] <- map_out$Y[,3]}
      
      dr$sample_id <- sample_ids[map_inds]
      lista <- list(dr, map_inds, inds_reduced)
      return(lista)
    }
  }
  
  if (type == "UMAP"){
    options(warn = 0)
    if ((two_three == FALSE)||is.null(two_three))
    {map_out <- try(expr = umap(d = map_expr_reduced, verbose = TRUE, n_components = 2, n_neighbors = UMAP.n_neighbors, 
                                metric = UMAP.metric, min_dist = UMAP.min_dist, spread = UMAP.spread, random_state = UMAP.random_state), silent = FALSE)}
    else
      {map_out <- try(expr = umap(d = map_expr_reduced, verbose = TRUE, n_components = 3, n_neighbors = UMAP.n_neighbors, 
                                metric = UMAP.metric, min_dist = UMAP.min_dist, spread = UMAP.spread, random_state = UMAP.random_state), silent = FALSE)}

    if (!(inherits(x = map_out,"try-error"))){ 
      if ((two_three == FALSE)||is.null(two_three))
      {
        dr <- data.frame(exprsub, check.names = FALSE)
        dr$map_x <- NA
        dr$map_y <- NA
        for (i in seq_along(1:nrow(map_expr_reduced))){
          dr$map_x[inds_reduced[i]] <- map_out$layout[,1][i]
          dr$map_y[inds_reduced[i]] <- map_out$layout[,2][i]}}
      else{
        dr <- data.frame(exprsub, check.names = FALSE)
        dr$map_x <- NA
        dr$map_y <- NA
        for (i in seq_along(1:nrow(map_expr_reduced))){
          dr$map_x[inds_reduced[i]] <- map_out$layout[,1][i]
          dr$map_y[inds_reduced[i]] <- map_out$layout[,2][i]
          dr$map_z[inds_reduced[i]] <- map_out$layout[,3][i]}}
      
    dr$sample_id <- sample_ids[map_inds]
    lista <- list(dr, map_inds, inds_reduced)
    return(lista)}}
    
    # Process any error messages
  if (inherits(x = map_out,"try-error")) {
    # Ignore warnings while processing errors
    options(warn = -1)
      
    # If this script were a function, warning() and stop() could be called to pass errors upstream
    msg <- geterrmessage()

    # Restore default warning reporting
    options(warn=0)
      
    messaggio <- paste0("System reports: ", msg, "... means 'something went wrong with the maping algorithm'")
      # Restore default warning reporting
    options(warn=0)
    return(messaggio)}
}


map_plot_comp <- function(map_df, type = NULL, two_three = NULL, map_inds=NULL, map_reduced=NULL, metadata=NULL, clust=NULL, 
                          metaclust=NULL, NMC=NULL, color_clusters=NULL, save = NULL) 
{
  
  #map_df <- na.omit(map_df)
  if (!(type == "marker")){
    if (!is.null(metaclust))
    {
      
      code_clustering <- metaclust # [number of clusters]
      #clust <- clust[map_reduced] #introduced after the mapping limitation
      
      cell_clustering <- code_clustering[clust] # [number of cells] with #nr. of metacluster value
      map_df$cell_clustering <- factor(cell_clustering[map_inds], levels = 1:NMC) # add cell_clustering to map_df
      
      map_df <- na.omit(map_df) #introduced after the mapping limitation
      
      if (nrow(map_df) > RENDER_DATA) { 
        map_df <- sample_n(tbl = map_df, size = RENDER_DATA, replace = F)}
      
      ## Plot map colored by clusters
      if (two_three == FALSE){ #2D
        p <- ggplot(map_df,  aes(x = map_x, y = map_y, color = cell_clustering)) + geom_point(size = 0.6) +
          theme_bw() + #https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels/
          scale_color_manual(values = color_clusters) +
          guides(color = guide_legend(override.aes = list(size = 4), ncol = 2)) #was size = 20
        if (!is.null(save)){
          ggsave(filename = "map_metacluster.png", plot = p, device = "png", path = "./tmpdata/")}
        pl <- ggplotly(p)}
      else { #3D
        pl <- plot_ly(map_df, x = ~map_x, y = ~map_y, z = ~map_z, color = ~cell_clustering, colors = color_clusters) %>%
          add_markers(marker = list(size = 0.6)) %>%
          layout(scene = list(xaxis = list(title = 'x-map'),
                              yaxis = list(title = 'y-map'),
                              zaxis = list(title = 'z-map')),autosize = TRUE)}}
    else
    {
      if (!is.null(color_clusters))
      {
        cell_clustering <- clust # [number of cells] with #nr. of metacluster value
        map_df$cell_clustering <- factor(cell_clustering) # add cell_clustering to map_df
        map_df <- na.omit(map_df) #introduced after the mapping limitation
        
        pl <- ggplot(map_df, aes(x = map_x, y = map_y, color = cell_clustering)) + geom_point(size = 0.5) +
          scale_color_manual(values = color_clusters) +
          theme(
            # get rid of panel grids
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # Change plot and panel background
            plot.background=element_rect(fill = "gray"),
            panel.background = element_rect(fill = 'black'))
        if (!is.null(save)){
          ggsave(filename = "map_cluster.png", plot = pl, device = "png", path = "./tmpdata/")}}
      else {
        map_df <- na.omit(map_df) #introduced after the mapping limitation
        xr <- range(map_df$map_x)
        yr <- range(map_df$map_y)
        get_density <- function(x, y, ...) {
          dens <- MASS::kde2d(x, y, lims = c(xr,yr))
          ix <- findInterval(x, dens$x)
          iy <- findInterval(y, dens$y)
          ii <- cbind(ix, iy)
          return(dens$z[ii])}
        
        density_map <- tryCatch(get_density(map_df$map_x, map_df$map_y, n = 100), 
                            error=function(err) 
                            {print("Fail to produce density!") 
                              return(NA)})
        if (!(is.numeric(density_map))){
          pl <- ggplot(map_df, aes(x = map_x, y = map_y)) + geom_point(aes(x = map_x, y = map_y)) +
            theme(
              # get rid of panel grids
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              # Change plot and panel background
              plot.background=element_rect(fill = "gray"),
              panel.background = element_rect(fill = 'black'))
          if (!is.null(save)){
            ggsave(filename = "map.png", plot = pl, device = "png", path = "./tmpdata/")}
          } 
        else {
          pl <- ggplot(map_df, aes(x = map_x, y = map_y)) + geom_point(aes(x = map_x, y = map_y, color = density_map)) +
            scale_color_gradientn(colours = terrain.colors(20)) +
            theme(
              # get rid of panel grids
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              # Change plot and panel background
              plot.background=element_rect(fill = "gray"),
              panel.background = element_rect(fill = 'black'))
          if (!is.null(save)){
            ggsave(filename = "map.png", plot = pl, device = "png", path = "./tmpdata/")}}}}
    return(pl)}
  else{
    nomi <- make.names(names(map_df))
    colnames(map_df) <- nomi
    
    nomi <- nomi[!nomi %in% grep(paste0(c("map_x", "map_y", "sample_id"), collapse = "|"), nomi, value = T)] 
    #to get rid of the "map_x", "map_y" and "sample_id names" 
    lista_plt <- vector("list", length = length(nomi))
    #example: plot_cd4 <- ggplot(map_df,  aes(x = tSNEx, y = tSNEy, color = CD4)) +  geom_point(size = 0.9) +  theme_bw() +
    #  scale_color_gradientn("CD4", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
    
    stringa_lista <- ""
    stringa_nomi <- ""
    
    map_df <- map_df[complete.cases(map_df), ]  
    for (i in seq_along(1:length(nomi))){
      foo <- paste0("aes(x = map_x, y = map_y, color = ", nomi[i])
      foo2 <- paste0("ggplot(map_df, ", foo, "))")
      foo3 <- paste0(foo2, "+ geom_point(size = 0.9) + theme_bw() + scale_color_gradientn(")
      foo4 <- paste0(foo3,"'", nomi[i],"'")
      foo5 <- paste0(foo4, ", colours = colorRampPalette(rev(brewer.pal(n = 11, name = 'Spectral')))(50))")
      parsed <- parse(text = foo5)
      lista_plt[[i]] <- eval(expr = parsed)
      stringa_lista <- paste0(stringa_lista, "lista_plt[[",as.character(i),"]],")
      stringa_nomi <- paste0(stringa_nomi, "nomi[",as.character(i),"],")
    }
    stringa_nomi <- substr(stringa_nomi, 1, nchar(stringa_nomi)-1)
    numero_righe <- ceiling(length(nomi)/2)
    #example: #figure <- ggarrange(lista_plt[[1]], lista_plt[[2]], lista_plt[[3]],lista_plt[[4]],lista_plt[[5]], 
    #   labels = c(nomi[1], nomi[2], nomi[3], nomi[4], nomi[5]), ncol = 2, nrow = numero_righe)
    foo <- paste0("ggarrange(", stringa_lista,"labels = c(", stringa_nomi,"), ncol = 2, nrow =",as.character(numero_righe),")")
    parsed <- parse(text = foo)
    figura <- eval(expr = parsed)
    if (!is.null(save)){
      altezza <- (length(mk$selected[mk$selected==TRUE]) * 200)/100 # dimensions are inches
      larghezza <- 1350/100
      ggsave(filename = "marker_map.png", plot = figura, device = "png", path = "./tmpdata/", width = larghezza, height = altezza)}
    return(figura)}
}


phenocluster <- function(flow_Set, selected_markers, clustering, phenoquery, central, min_quantile, 
                         max_quantile, quant, pos_threshold, neg_threshold, cluster.signature = NULL){
  
  res.expr <- exprs_sub(flow.Set = flow_Set, selected.markers = selected_markers)
  exprsub <- res.expr[[1]]
  if ('cell_Id' %in% colnames(exprsub)){
    posizione <- which(colnames(exprsub) == "cell_Id")
    exprsub <- exprsub[,-posizione]}
  
  pDdesc_subset <- res.expr[[3]]
  if ('cell_Id' %in% pDdesc_subset){
    posizione <- which(pDdesc_subset == "cell_Id")
    pDdesc_subset <- pDdesc_subset[-posizione]}
  
  colnames(exprsub) <- pDdesc_subset
  
  ##### function description 
  #The goal is to produce a dataframe in which, for each of the calculated metacluster, you assign a phenotype  
  #described by the the phenoquery
  #
  # input description:
  #
  # mat_expr = matrix of marker expression. It is the whole expression matrix of the concatenated flowFrame 
  # collecting all the transformed flowFrame in a single flowSet: expr <- fsApply(tfS, exprs) 
  # cv is the integer vector of the cluster id: one cluster for each event
  # phenoquery = is a dataframe of dimension M*(N+1) where M is the number of distinct phenotypes expressed 
  # by each cell and N is the number of related markers (e.g: 	
  #   CD62L	CD45RA	CD95
  # Naive	+	+	-
  # CM	+	-	+
  # EM	-	-	+
  # TM	-	+	+
  # TSCM	+	+	-
  # each of the elements could have one of the three characters value: 
  #  "+" (positevely expressed), 
  #  "-" (non express) or 
  #  "*" indifferent, to take into account the case in which not all the markers are relevant for a specific
  # phenotype
  # the syntax of the Cell Ontology (flowCL package). Of course the following relationship must hold:
  #  pos_threshold (+) >= neg_threshold (-)
  #
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  # output description - res_log is a data.frame [number of clusters x 3] where the first column reports the original clusters 
  # having the number of levels equals to the ones found with the clustering algorithm (100 with flowSOM) the second column is 
  # the phenotype name as reported by the phenoquery table (ph) and the third reports the related number of the phenotype mapped 
  # on the various clusters: if you enter 12 phenotype you will have 1:12 levels that is each single row could be 
  # e.g: 81; "Naive cells", which means: 
  # the cluster 81 contains (is mapped) as a Naive cells (that is the 8th phenotype entered)
  # 
  
  normalize.function <- function(x, quantiles=quant) { #used to normalize from 0 to 1 the expression range
    
    if (!quant){
      rng <- colQuantiles(x = exprsub, probs = c(min_quantile, max_quantile), drop = TRUE) #[nr of markers x min_quantile [,1] , max quitile [,2]]
      expr01 <- t((t(exprsub) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    }else{ expr01 <- (x - min(x)) / (max(x) - min(x))}
    expr01[expr01 < 0] <- 0
    expr01[expr01 > 1] <- 1
    return (expr01)}
  
  expr01 <- normalize.function(x = exprsub, quantiles = quant)
  
  # Calculate the central tendency expression
  if ((central == "median")||(central == "mean")){
    expr01_median <- data.frame(expr01, cell_clustering = clustering, check.names = FALSE) %>%
      group_by(cell_clustering) %>% 
      summarize_all(central)} #[number of cluster, ]
  
  if (central == "mode"){
    my_mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]}
  
      expr01_median <- data.frame(expr01, cell_clustering = clustering, check.names = FALSE) %>%
      group_by(cell_clustering) %>% 
      summarize_all(my_mode)} #[number of cluster, ]
  
  if (!is.null(cluster.signature)) {
    expr01 <- cluster.signature
    n_clust <- nlevels(factor(clustering))
    expr01_median <- data.frame(expr01, cell_clustering = 1:n_clust, check.names = FALSE)} 
  #in case input$signature_finding_method=="Densities"
  
  
  expr_heat <- as.matrix(expr01_median[, c(colnames(exprsub))]) #[#meta-cluster x all markers]
  phenotype <- phenoquery[,1]
  phenoquery <- phenoquery[, -1]
  phenoquery <- as_tibble(phenoquery)
  expr_heat <- as_tibble(expr_heat[, pDdesc_subset])
  
  phenoquery<-phenoquery[,colnames(expr_heat)] 
  
  mutate_func1 <- function(var) {ifelse(var >= pos_threshold, 1L, (ifelse(var < neg_threshold, 0L, NA)))}
  
  expr_log <- expr_heat %>% mutate_all(.funs = mutate_func1) %>% mutate(cluster=expr01_median$cell_clustering) %>% 
    dplyr::select(-cluster,everything())
  expr_log_as_matrix <- as.matrix(x = expr_log)
  
  mutate_func2 <- function(var) {ifelse(var == "+", 1L, (ifelse(var == "-", 0L, 2)))}
  
  phenolog <- phenoquery %>% mutate_all(.funs = mutate_func2)
  phenolog <- phenolog %>% mutate(phenotype = phenotype)
  phenolog_as_matrix <- as.matrix(x = phenolog[,-(ncol(phenolog))])
  
  #left_log <- expr_log %>% left_join(y = phenolog)
  logCol <- logical(length = nrow(expr_log))
  logMat <- matrix(nrow = nrow(expr_log), ncol = ncol(expr_log)-1)
  phenotype <- rep(NA,nrow(expr_log))
  phenos <- matrix(nrow = nrow(expr_log), ncol = nrow(phenolog))
  up_log <- expr_log
  
  # LogMat, the nrow(expr_log)Xnrow(expr_log) logical matrix, collects every logCol column for every marker. It could have only TRUE or FALSE
  # It is TRUE when the expr_log is above the pos_threshold and the phenolog is "+"" or "*" or when expr_log is below the neg_threshold and the correspondent phenolog is "-" or "*"
  # Only the rows with all elements TRUE will be assiged to the signature in the phenolog
  # phenos collects every phenotypes collected after each row fo the phenolog. Only the one with one signature only will be assigned. Otherwise it will be "overlap"
  for (i in 1:nrow(phenolog_as_matrix)){        # phenolog is the translation in 0,1,2 of the phenotype table [i x j] = [number of entries x markers]
    for (j in 1:(ncol(phenolog_as_matrix))){
      ph <- phenolog_as_matrix[i,j]
      elCol <- expr_log_as_matrix[,j]           # expr_log is the 0 1 table derived from the expr_heat [i x j] = [clusters x markers]  
      for (k in 1:nrow(expr_log_as_matrix)){
        el <- elCol[k]
        if (ph == 2) 
          {logCol[k] <- TRUE}
        else 
          {logCol[k] <- ifelse(ph == el, TRUE, FALSE)}}
      logMat[,j] <- logCol}
    check <- FALSE
    for (l in 1:nrow(expr_log_as_matrix)){
      tmpRow <- logMat[l,]
      tmpRow[is.na(tmpRow)] <- FALSE
      if (all(tmpRow))
      {phenotype[l] <- as.character(phenolog$phenotype[i])}}
    phenos[,i] <- phenotype}
  
  for(i in 1:nrow(phenos)){
    if (!all(is.na(phenos[i,]))){   #if the row does not have all NA
      #if (sum(is.na(phenos[i,])) < (ncol(phenos))-1){ #if the row does have only one element different from the other
      if (length((unique(phenos[i,])) > 1)){
        row_clean <- na.omit(phenos[i,]) 
        if(length(unique(as.list(row_clean))) == 1){
          phenos[i,] <- row_clean[1]}
        else{phenos[i,] <- "OVERLAP"}}}}
  
  for(i in 1:nrow(phenos)){
    if (is.na(unique(phenos[i,]))) {
      phenotype[i] <- NA
      next}
    if (phenos[i,1] == "OVERLAP") {phenotype[i] <- "OVERLAP"}
    if (phenos[i,1] != "OVERLAP") {phenotype[i] <- phenos[i,1]}}
 
  up_log$phenotype <- phenotype 
  p <- up_log %>% dplyr::select(phenotype) %>% mutate_if(is.factor, as.character)
  
  for(i in 1:nrow(expr_log)){
    if(is.na(p$phenotype[i])){p$phenotype[i] <- "UNASSIGNED"}}
  res_log <- up_log %>% mutate(phenotype = p$phenotype) %>% dplyr::select(cluster, phenotype) %>% 
    rename(original_cluster = cluster, new_cluster = phenotype) 
  
  livelli <- factor(res_log$new_cluster)
  numerico <- as.numeric(livelli)
  res_log$numeric <- numerico
  return (res_log)
}


cluster_signature <- function(fS, selected.marker, clustId, color_marker, tipo = NULL, central, min_quantile, 
                              max_quantile, quant, pos_threshold, neg_threshold){
# The function is split in two: 
# - the first is to produce the density function for each single marker and to find the minimum of the sum of densities for each color 
# for each cluster. This minimum on the x position (normalized) provides the MFI reference for each single cluster to be used in place of the 0.5 (50%)
# used in the "traditional" way.
# - the second is to study the function of each single marker inside the cluster to find double peaks of densities: these double peaks
# are used to produce their values (x and f(x)) and to produce an heatmap of the marker/cluster matrix to suggest a possible increase of clusters (or 
# meta-cluster) to further split those events that potentially could be placed in more numerous clusters
#
# Inside the function there are two main functions: finite.differences to produce the derivative of the density's trend and the 
# studio_di_funzione to find max and mins
if (tipo == "cluster"){stringa <- character()}else{stringa <- "meta"} #this ti describe in the various strings, the correct naming
  
res.expr <- exprs_sub(flow.Set = fS, selected.markers = selected.marker)
exprsub <- res.expr[[1]]

pDdesc_subset <- res.expr[[3]]
pDname_subset <- res.expr[[2]]

colnames(exprsub) <- pDdesc_subset
marker_name <- pDdesc_subset
names(color_marker) <- marker_name

Nr_event <- nrow(exprsub) #total number of events
matrice_expr <- cbind(exprsub, clusterId = clustId)

# The following is to perform the derivative of a function
#References (https://rstudio-pubs-static.s3.amazonaws.com/295650_406b32cdd2d34ca3a150bbe95010d665.html)
#Burden, R. L., & Faires, J. D. (2011). Numerical analysis (9th ed.). Boston, MA: Brooks/Cole, Cengage Learning.
#Weisstein, Eric W. “Numerical Differentiation.” From MathWorld–A Wolfram Web Resource. 
#http://mathworld.wolfram.com/NumericalDifferentiation.html 
# In alternative see https://rviews.rstudio.com/2021/05/04/functional-data-analysis-in-r/ 
finite.differences <- function(x, y) {
  if (length(x) != length(y)) {stop('x and y vectors must have equal length')}
  n <- length(x)
  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = length(x))
  
  for (i in 2:(n-1)) {fdx[i] <- (y[i+1] - y[i-1]) / 2*(x[i+1] - x[i-1])}
  
  # Originally was fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i]) the forward differencing method, but this generates a too jagged function 
  
  # For the last value, since we are unable to estimate the x[0] and the x[n], since most of the times are equal to 0 let's fix it this way
  fdx[is.na(fdx)] <- 0 # to replace NA with 0
  
  fdx[n] <- 0
  fdx[1] <- 0
  return(fdx)}

normalize.function <- function(x, quantiles=quant) { #used to normalize from 0 to 1 the expression range
  
  if (!quant){
    rng <- colQuantiles(x = exprsub, probs = c(min_quantile, max_quantile), drop = TRUE) #[nr of markers x min_quantile [,1] , max quitile [,2]]
    expr01 <- t((t(exprsub) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  }else{ expr01 <- (x - min(x)) / (max(x) - min(x))}
  expr01[expr01 < 0] <- 0
  expr01[expr01 > 1] <- 1
  return (expr01)}

studio_di_funzione <- function(x_granularity, dens_func, derivate, ilcluster, ilmarker, central) {
  # Starting from the left to find the change in sign of df_dx. The massimo.minimo() function seeks the point of maximum of the
  # densities. It performs the maximum lookup from left to right. The entry values should be normalized (the y axis goes from 0 to 1)
  # It accept the following frame parameters: 
  # - peak_space : % minimum space between two peaks
  # - peak_gorge : % peak gorge between the first peak and the related minimum 
  
  # after the first maximum found it looks for a second maximum (if it exists) still proceeding from left to right until
  # it finds a second maximum peak which cannot be closer than the peak_space, the difference between the first peak and the 
  # gorge between the two peaks cannot be less than the minimum peak_gorge and above the density background 
  
  # and these are the density variables along the procedure below
  # - max1: the first maximum (relative value, in percentage)
  # - max2: the second maximum (relative value, in percentage)
  # - x_max1: the x of the first maximum (absolute value in pixel vector length that is 512. To get in percentage: x * 100 / vector_length)
  # - x_max2: the x of the second maximum (absolute value in pixel vector length that is 512)
  # - minimum:  minimum between the two peaks
  # - x_min:  minimum between the two peaks
  # - dens_bkg: the density background value that is the full-scale above the minimum densities (aka the most present "low" density above the 0)
  # - min_dens: the minimum normalized density % to be taken into account 
  
  # To be included in the UI #
  
  peak_space <- 0.05 # to be added on the general io mask to allow configuration. It sets also the minimum space between occurrences (max or min)
  peak_gorge <- 0.15 #to be added on the general io mask to allow configuration was 0.5)
  min_dens <- 0.04 #to be added on the general io mask to allow configuration - set to 0 because it happens that a change in the derivative sign happens 
  #around a very low level of density
  
  # all the values should be expressed in % (e.g. 0.05 means 5%)
  
  found_something <- FALSE
  max1_found <- FALSE
  max2_found <- FALSE
  min_found <- FALSE; min1_found <- FALSE
  max1 <- 0; x_max1 <- 0; max2 <- 0; x_max2 <- 0; min1 <- 0; 
  x_min1 <- 0; min2 <- 0; x_min2 <- 0; minimum <- 0; x_min <- 0; 
  min_derive <- 1.0e-5
  
  left_guard <- 0.01; right_guard <- 0.01
  #due to the possible peaks of the derivative it's is better to exclude the firsts values (0.01 = 1% = more than 5 values per side)
  
  ####################################### this is just for debugging 
  ######################################## use mad's samples with CD4CAR_handled 
  #if (ilcluster==1){browser()}
  #this entry part of the procedure is to select the most frequent density value to find the "background" full_scale. To do this let's split the density
  # value in finite number of chunk
  dens_low <- dens_func[dens_func > 0.005] # let's remove the lowest densities (less than 0.5% of density)
  splits <- floor(length(dens_low)/5) # this is quite arbitrary: it works with most of the functions
  if (splits < 2){splits <- 2}
  out_cut <- cut.default(x = dens_low, breaks = splits) #produce values splitting in ceiling(length(dens_low)/2) numbers of chunks
  chunks <- levels(out_cut)
  cut_table  <- data.frame(table(out_cut)) #counting the entries for each chunk
  cut_table <- cut_table[-c(ceiling(nrow(cut_table)/2):nrow(cut_table)),] #removing half of the table (leaving the lower density levels only)
  max_chunk <- cut_table[which.max (x = cut_table$Freq),] #select the maximum, that is the density value with the greatest number of entries (the most frequent)
  dens_bkg <- 1/(length(chunks))*(as.numeric(row.names(max_chunk))) #this selects the density background that is (max frequent chunk + 1)/(number of chunks)
  k <- 0
  skip <- TRUE
  spacing <- floor(peak_space*x_granularity) - 1 #these are the points which has to separates two distinguished events (min or max or changes in the derivative sign)
  
  check_rle<-rle(abs(derivate)<min_derive)
  check_rle_rev <- rle(abs(rev(derivate))<min_derive)
  k_max <- x_granularity - check_rle_rev$lengths[1]
  left_border <- floor(left_guard*x_granularity) 
  right_border <- x_granularity-floor(left_guard*x_granularity)
  rle_points <- rle(sign(derivate))
  relevant_points <- vector(mode = "integer", length = (length(rle_points$values)-1))
  ix <- 1; kx <- 0
  while (ix < (length(rle_points$values))) { #this is just to find the relevant points which are the once where the derivative changes in sign
    kx <-  rle_points[[1]][[ix]] + kx
    relevant_points[ix] <- kx
    ix <- ix + 1}
  #rbowser()
  relevant_points_reference <- rbind(relevant_points, round(relevant_points/x_granularity, 3)) #this is to help in the debugging
  for (i in left_border:right_border){ #to find the min/max points
    #if((i==382)||(i==455)||(i==460)||(i==758)||(i==854)||(i==1066)||(i==1129)||(i==1222)||(i==1228)||(i==1252)){rbowser()}  
    #the range is limited in the left and in the right of the two guards
    if ((k <= check_rle$lengths[1])&&(i<right_border)&&(skip)){k <- k+1; if (k==check_rle$lengths[1]){skip<-FALSE}; next} #to skip the lowest derivate's 
    #value typically at the beginning or the end of the normalized derivate vector. It should be done only once (skip)
    if ((i > k_max)&&(i < right_border)){next} #to skip the latest derivate values which are less than min_derive
    #if ((k < spacing - 1)&&(found_something)){k <- k+1; next} #this is used to skip the next steps based on spacing
    if ((k < floor(spacing/2))&&(found_something)){k <- k+1; next} #this is used to skip the next steps based on spacing
    # spacing/2 because sometimes the density function could be unstable for a long interval 
    found_something <- FALSE #in case the algorithm found something along the way in the next loop
    if ((i < (right_border - floor(spacing)))&&((derivate[i]*derivate[i+1]< 0)&&(!max1_found))) { #if the derivative changes the sign and it is the first max found
      if (dens_func[i-1] < dens_func[i]) { #case of maximum
        #if ((ilcluster==1)&&(ilmarker==3)){rbowser()}
        ####################################### this is just for debugging 
        #rbowser()#
        #the max1 is a maximum candidate: we should verify if:
        if ((dens_func[i] > (dens_bkg + min_dens))&&(var(sign(derivate[(i+1):(i+floor(spacing))]))==0)){ 
          # the max must be above the background
          #var(sign(derivate[(i+1):(i+floor(spacing/2))]))==0 is to filter too much rapid changes in the derivative (to avoid those show tooth in the 
          #derivative) before and after the current point. It was spacing/2 on both side but it turns out it is too much restricted: better just on the 
          #forward side
          max1 <- dens_func[i]
          x_max1 <- i # to express in percentage: x * 100 / x_granularity
          max1_found <- TRUE
          found_something <- TRUE
          k <- 0}}}
    if ((i < (right_border - floor(spacing)))&&((derivate[i]*derivate[i+1]< 0)&&(i > x_max1)&&(!max2_found))) { 
      #if the derivative changes the sign for the second time
      if (dens_func[i-1] > dens_func[i]){ #case of minimum - after a maximum there is always a minimum. 
        #if ((ilcluster==1)&&(ilmarker==3)){rbowser()}
        ####################################### this is just for debugging 
        #rbowser()#
        if ((abs(max1 - peak_gorge) > dens_func[i])&&(var(sign(derivate[(i+1):(i+floor(spacing*2/3))]))==0)){ 
          #This is to avoid a second max too much close to the min (e.g. the case you have a saw tooth derivative)
          # 2) case of minimum enough different from the density peak (otherwise skip the minimum) && to avoid saw tooth in the derivative
          minimum <- dens_func[i]
          x_min <- i
          min_found <- TRUE
          found_something <- TRUE
          k <- 0}}
      if (dens_func[i] > dens_func[i-1]){ #case of maximum: let's see if it is a max2 candidate  
        #rbowser()#
        if (((abs(i - x_max1) > spacing) > 0)&&(min_found)&&(var(sign(derivate[(i+1):(i+floor(spacing/3))]))==0)) { 
          #3) to verify the if the second maximum is far enough from the first one (peak_space should be in %, it's inside the spacing var) 
          # also to verify if there has been a minimum before (if not, it means that the two possible peaks are not enough different)
          #rbowser()#
          #if ((ilcluster==1)&&(ilmarker==3)){rbowser()}
          ####################################### this is just for debugging 
          if (dens_func[i] > dens_bkg + min_dens){ 
            # 4) to verify if the second peak stands out from the background and if ti is enough different from the other peak
            max2_found <- TRUE
            max2 <- dens_func[i]
            x_max2 <- i
            found_something <- TRUE
            k <- 0
            if (!max1_found) { # in case the second peak is found, but the first one is still missing (they swap)
              max1_found <- TRUE; max1 <- dens_func[i]; x_max1 <- i
              max2_found <- FALSE; max2 <- 0; x_max2 <- 0}} 
        }}}
    } # end of for cycle
  ############# two cycles: one for the maxs and one for the mins
  k <- 0; skip <- TRUE; min1 <- 0; x_min1 <- 0; min1_found <- FALSE; min2 <- 0; x_min2 <- 0; min2_found <- FALSE
  for (i in left_border:right_border){ #to find the min points (starting and ending before the limits)
    #if((i==505)||(i==758)||(i==782)||(i==782)||(i==902)||(i==1042)||(i==1252)||(i==1340)) #{rbowser()}
    if ((k <= check_rle$lengths[1])&&(i<right_border)&&(skip)){k <- k+1; if (k==check_rle$lengths[1]){skip<-FALSE}; next} #to skip the lowest 
    #derivate's value typically at the beginning or the end of the normalized derivate vector. It should be done only once (skip)
    if ((i > k_max)&&(i < right_border)){next} #to skip the latest derivate values which are less than min_derive
    if ((floor(spacing/2))&&(found_something)){k <- k+1; next}
    found_something <- FALSE
    if ((i < (right_border - floor(spacing)))&&((derivate[i]*derivate[i+1]< 0)&&(!min1_found)&&(i > x_max1))) { #if the derivative changes the sign and it is the first min found
      if (dens_func[i-1] > dens_func[i]) { #case of minimum
        #rbowser()#
        #if ((ilcluster==1)&&(ilmarker==3)){rbowser()}
        ####################################### this is just for debugging 
        if ((abs(max1 - peak_gorge) > dens_func[i])&&
            (abs(x_min1 - i) > spacing)&&(abs(x_min2 - i) > spacing)&&(abs(x_max1 - i) > spacing)&&(var(sign(derivate[(i+1):(i+floor(spacing/2))]))==0)){
          # This to avoid than the minimum points are closer than the "spacing" due to the peak space
          min1 <- dens_func[i]
          x_min1 <- i # to express in percentage: x * 100 / x_granularity
          minimum <- min1
          x_min <- x_min1
          min1_found <- TRUE
          found_something <- TRUE
          k <- 0}}}  
    if ((i < (right_border - floor(spacing)))&&(derivate[i]*derivate[i+1]< 0)&&(i > x_max1)&&(i < x_max2)) { #if the derivative changes the sign for the second time
      if (dens_func[i-1] > dens_func[i]) {#case of minimum
        #rbowser()#
        #if ((ilcluster==1)&&(ilmarker==3)){rbowser()}
        ####################################### this is just for debugging 
        if ((dens_func[i]>dens_bkg)&&(abs(x_min1 - i) < spacing)&&(abs(x_min2 - i) < spacing)&&(abs(x_max1 - i) < spacing)&&(abs(x_max2 - i) < spacing)&&
            (var(sign(derivate[(i+1):(i+floor(spacing/2))]))==0)){
          min2 <- dens_func[i]
          x_min2 <- i # to express in percentage: x * 100 / x_granularity
          min2_found <- TRUE
          found_something <- TRUE
          k <- 0}}}}# end for cycle
  
  #if ((ilcluster==1)&&(ilmarker==3)){rbowser()}
  ####################################### this is just for debugging 
    
if((min1_found)||(min2_found)){min_found <- TRUE}
  if ((min2<min1)&&(abs(x_min2 - 0.5) > abs(x_min1 - 0.5))){minimum <- min2; x_min <- x_min2} 
  if (max1_found) {max1 <- max(dens_func); x_max1 <- which.max(dens_func)} #the max1 is always the highest
  if (!((max1_found)&&(max2_found))){minimum <- 0; x_min<-0; min_found<- FALSE} #the minimum exists only if both max1 and max2 has been found
  if ((max1_found && max2_found)&&!((x_min > min(x_max1, x_max2))&&(x_min < max(x_max1,x_max2)))){
    minimum <- 0; x_min<-0; min_found<- FALSE
    max2 <- 0; x_max2 <- 0; max2_found <- FALSE} #the minimum exists only if it is in between max1 and max2 
  if (!max1_found){max1 <- 0; x_max1 <- 0; max1_found <- FALSE}
  
  x_max1 <- x_max1/x_granularity
  x_max2 <- x_max2/x_granularity
  x_min <- x_min/x_granularity
  
  lista <- list("max1 found" = max1_found, "max1" = max1,"x of max1" = x_max1, 
                "max2 found" = max2_found, "max2" = max2,"x of max2" = x_max2,
                "min found" = min_found, "min" = minimum, "x of minimum" = x_min,
                "density background" =dens_bkg) #fill the list with the max points and the list of index found2
  return (lista) 
} # end of function

################ Start analysis: cluster (i) x marker (j) ###############
clustId <- matrice_expr[,"clusterId"] #it could be meta-cluster in case of ...
n_clust <- nlevels(factor(clustId))
mat_list <- vector(mode = "list",length = n_clust) #These are the expression matrix with all the events for a cluster
new_MFI = double(n_clust)
result_marker_list <- vector(mode = "list", length = length(color_marker))
result_dens_sum_list <- vector(mode = "list", length = n_clust)
result_list <- vector(mode = "list", length = n_clust) #for (aName in marker_name){result_list[[aName]] <- list(aName)} #to name the markers inside the cluster

ld_list_MFI <- vector(mode = "list", length = length(color_marker)) 
#This is the list of the calculated ggplot layers reporting the sum of densities for the selected cluster
ld_list <- vector(mode = "list", length = length(color_marker)) 
#This is the list of the ggplot layers reporting the marker density for the selected cluster
expr01_peak <- data.frame(matrix(ncol = length(marker_name), nrow = n_clust)) 
expr01_no_delta <- data.frame(matrix(ncol = length(marker_name), nrow = n_clust)) 
expr_log <- data.frame(matrix(ncol = length(marker_name), nrow = n_clust)) 
warning_msg_1 <- ""; warning_msg_2 <- ""; warning_msg_3 <- ""; warning_msg_4 <- "";  
issue <- 0

if (tipo == "cluster"){warning_df <- data.frame (issue_nr = 0, cluster_nr = 0,marker_name = "-", warning_txt = "")}else
  {warning_df <- data.frame (issue_nr = 0, metacluster_nr = 0,marker_name = "-", warning_txt = "")}
warning_df <- warning_df[-1,] 

matrice <- matrice_expr[, setdiff(colnames(matrice_expr), "clusterId")]
exprsub <- normalize.function(matrice, quant) #this is to normalize the entire expression matrix
matrice_expr <- cbind(exprsub, clusterId = clustId)

for (i in 1:n_clust){ #from now on every analysis and function's study is related to the single cluster
  #if (i==1){browser()} #to debug selecting the analysis to study for debugging
  mat_list[[i]] <- matrice_expr[matrice_expr[, "clusterId"] == i,,drop=TRUE] 
  exprsub <- mat_list[[i]][, pDdesc_subset]
  #exprsub <- normalize.function(exprsub) #this is to normalize the entire expression matrix
  ggdf <- as.data.frame(exprsub)
  nr_event <- nrow(ggdf) #these r the number of events of the clusters: % = nr_event/Nr_event
  
  #----------------------------------------------------------------------------------------------------------------------#
  # WARNING GENERATOR: if the minimum of the sum of densities has not been found
  #warning_msg <- paste("The minimum of the sum of densities of cluster nr.",as.character(i), 
  #                     " has not been found: The related MFI cannot be correctly evaluated")
  #warning_msg_1 <- paste(warning_msg_2,warning_msg, sep = "\n")
  #----------------------------------------------------------------------------------------------------------------------#
  
  # The following is to isolate the density function for each marker in the selected cluster_id (i)
  # This will produce the plot of the marker expression densities. I'll extract with layer_data() the various densities
  ggdf <- reshape2::melt(data = ggdf, value.name = "expression_range", variable.name = "marker", id.vars=NULL) #nrow(ggdf) nr_event*length(marker_name)
  dens_plot <- ggplot2::ggplot(ggdf, ggplot2::aes(x = expression_range, color = marker, group = marker)) +
    ggplot2::geom_freqpoly(stat = "density", linewidth = 1) + 
    #facet_wrap(~ marker, scales = "free", ncol = colonne, nrow = righe) +
    ggplot2::scale_color_manual(values = color_marker)
  #pl_dens_plot <- plotly::ggplotly(dens_plot) #in case of plotly
  
  ld_dens <- ggplot2::layer_data(plot = dens_plot, i = 1) #pg <- ggplot2::ggplot_build(p) #pgdf <- pg[["data"]] 
  # The other way is to use ld_dens_data <- ggplot2::ggplot_build(plot = dens_plot). data are the same
  
  vector_length <- length(ld_dens$x)/length(marker_name) #see the layer_data vector length for each marker = 512
  sum_of_density <- vector(mode="numeric", length=vector_length)
  # The following is to sum all the densities functions 
  # This to find a center in the marker expression where to find the minimum level of density 
  #( namely, a splitting point between the negatively and the positively expressed markers )
 
  for (k in 1:length(color_marker)){
    marker_density <- dplyr::filter(.data = ld_dens, colour == color_marker[k])
    ld_list_MFI[[k]]$layer_data = dplyr::select(.data = marker_density, colour,y,x,density,scaled,ndensity,count)
    ld_list_MFI[[k]]$marker <- marker_name[k]
    ld_list_MFI[[k]]$number_of_event <- nr_event 
    sum_of_density <- sum_of_density + ld_list_MFI[[k]]$layer_data$density} 
  
  df <- data.frame(expression_range = ld_list_MFI[[1]]$layer_data$x, sum_of_densities = sum_of_density/max(sum_of_density))
  #using the ld_list_MFI[[1]]$layer_data$x which is common to every x axis' density
  p_dens_sum <- ggplot2::ggplot(data = df, ggplot2::aes(x = expression_range, y = sum_of_densities)) + 
    #the plot should be normalized for the y axis to be used in the "studio di funzione" function
    ggplot2::geom_area(mapping = ggplot2::aes(y = sum_of_densities), linewidth =1.2, color = "black", fill = "green") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
  
  #The following is to analyze the sum of densities - The routine should be similar to the one below used for each marker
  #In any case a new "studio_di_funzione" developed explicitly to find the best minimum density to divide in two the positive vs the negative marker zone 
  ld_dens_sum <- ggplot2::layer_data(plot = p_dens_sum, i = 1)
  
  df_dx <- finite.differences(x = ld_dens_sum$x, y = ld_dens_sum$y)
  scaled_df_dx <- df_dx/max(df_dx) #the same length of ld_dens_sum$x
  #scaled_density <-  ld_dens_sum$y / max(ld_dens_sum$y) # useless: max(ld_dens_sum$y) = 1
  
  f <- data.frame(expression_range = ld_dens_sum$x, marker_density = ld_dens_sum$y, scaled_df_dx)
  
  plot_dens_sum_con_la_derivata <- ggplot2::ggplot(data = f, mapping = ggplot2::aes(x = expression_range)) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = marker_density), linewidth=1.2, color = "red") +
    ggplot2::geom_line(mapping = ggplot2::aes(y = scaled_df_dx), linewidth=1, color = "blue") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  
  #if (i==1){rbowser()} #to debug selecting the analysis to study for debugging
  vector_length <- length(ld_dens_sum$x)
  result_MFI <- studio_di_funzione(x_granularity = vector_length, dens_func = ld_dens_sum$y, derivate = scaled_df_dx, ilcluster = i, ilmarker = 0)
# ilmarker = 0 because here we do not refer to any marker  
  ####################################################################### just for debugging purpose
  # In case the minimum has not been found: the new MFI will be the max1 or the maximum among the max1 and the max2
  
  if (result_MFI$`min found`){
    if ((result_MFI$`max1 found`)&&(!result_MFI$`max2 found`)){
      studio_MFI <- p_dens_sum  + 
        ggplot2::geom_vline(xintercept = result_MFI$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of max1`, label=paste("x of max density","=",as.character(round(result_MFI$`x of max1`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_vline(xintercept = result_MFI$`x of minimum`, linetype="twodash", color = "blue", linewidth=1.5) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of minimum`, label=paste("estimated MFI","=",as.character(round(result_MFI$`x of minimum`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_hline(yintercept = result_MFI$`density background`, linetype="dashed", color = "black", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(y=result_MFI$`density background`, label=paste("density background","=",as.character(round(result_MFI$`density background`, digits = 2))), x=0.5, hjust = -0.2)) +
        ggplot2::labs(title = paste(stringa,"cluster nr.", i, " : normalized sum of marker densities"))}
    
    if ((result_MFI$`max1 found`)&&(result_MFI$`max2 found`)){
      studio_MFI <- p_dens_sum + 
        ggplot2::geom_vline(xintercept = result_MFI$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of max1`, label=paste("x of max density","=",as.character(round(result_MFI$`x of max1`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_vline(xintercept = result_MFI$`x of minimum`, linetype="twodash", color = "blue", linewidth=1.5) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of minimum`, label=paste("estimated MFI","=",as.character(round(result_MFI$`x of minimum`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_vline(xintercept = result_MFI$`x of max2`, linetype="dashed", color = "red", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of max2`, label=paste("x of 2nd max density","=",as.character(round(result_MFI$`x of max2`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_hline(yintercept = result_MFI$`density background`, linetype="dashed", color = "black", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(y=result_MFI$`density background`, label=paste("density background","=",as.character(round(result_MFI$`density background`, digits = 2))), x=0.5, hjust = -0.2)) +
        ggplot2::labs(title = paste(stringa, "cluster nr.", i, " : normalized sum of marker densities"))}}
    #new_MFI[i] <- 0.5
    #if ((result_MFI$'x of minimum' < 0.75)&&(result_MFI$'x of minimum' > 0.25)){new_MFI[i] <- new_MFI[i] <- result_MFI$'x of minimum'}} 
    # was if (result_MFI$`min found`) only but I decided to use the minimum only when it is not too far from the center
    # I put 0.5 as MFI except in this case 
  
  else{ #if !(result_MFI$`min found`)
    if ((result_MFI$`max1 found`)&&(!result_MFI$`max2 found`)){
      studio_MFI <- p_dens_sum  + 
        ggplot2::geom_vline(xintercept = result_MFI$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of max1`, label=paste("x of max density","=",as.character(round(result_MFI$`x of max1`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_hline(yintercept = result_MFI$`density background`, linetype="dashed", color = "black", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(y=result_MFI$`density background`, label=paste("density background","=",as.character(round(result_MFI$`density background`, digits = 2))), x=0.5, hjust = -0.2)) +
        ggplot2::labs(title = paste(stringa, "cluster nr.", i, " : normalized sum of marker densities"))}
        #ggplot2::labs(title = paste(stringa, "cluster nr.", i, " : normalized sum of marker densities - NO Minimum found!"))}
    
    if ((result_MFI$`max1 found`)&&(result_MFI$`max2 found`)){
      studio_MFI <- p_dens_sum + 
        ggplot2::geom_vline(xintercept = result_MFI$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of max1`, label=paste("x of max density","=",as.character(round(result_MFI$`x of max1`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_vline(xintercept = result_MFI$`x of max2`, linetype="dashed", color = "red", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(x=result_MFI$`x of max2`, label=paste("x of 2nd max density","=",as.character(round(result_MFI$`x of max2`, digits = 2))), y=0.5, angle=90, vjust = 1.2)) +
        ggplot2::geom_hline(yintercept = result_MFI$`density background`, linetype="dashed", color = "black", linewidth=1) +
        ggplot2::geom_text(ggplot2::aes(y=result_MFI$`density background`, label=paste("density background","=",as.character(round(result_MFI$`density background`, digits = 2))), x=0.5, hjust = -0.2)) +
        ggplot2::labs(title = paste(stringa, "cluster nr.", i, " : normalized sum of marker densities"))}
        #ggplot2::labs(title = paste(stringa, "cluster nr.", i, " : normalized sum of marker densities - NO Minimum found!"))}
    result_MFI$`x of minimum` <- 0.5 #if not min found the MFI cannot be estimated
  } #else
  
  new_MFI[i] <- result_MFI$`x of minimum`
  result_dens_sum_list[[i]] <- result_MFI
  
  if (!result_MFI$`min found`){ #save the plot if result_MFI$`min found`==FALSE
    if (n_clust < 100){ ggplot2::ggsave(filename = paste0("_CS_",stringa,"cluster_","nr_",as.character(i),"_MFI_study.png"), 
                    plot = studio_MFI, device = "png", path = ("tmpdata/"))}
  #----------------------------------------------------------------------------------------------------------------------#
  # WARNING GENERATOR: if the minimum of the sum of densities has not been found
    warning_msg <- paste("The minimum of the sum of densities of ",stringa,"cluster nr.",as.character(i), 
                           " has not been found: The related MFI cannot be correctly evaluated")
    warning_msg_2 <- paste(warning_msg_2,warning_msg, sep = "\n")
    issue <- issue+1
    warning_df[nrow(warning_df) + 1,] <- c(as.character(issue), as.character(i), "-", "MFI cannot be estimated")
    #----------------------------------------------------------------------------------------------------------------------#
  } 
  
  if ((result_MFI$`x of minimum` < 0.25)||(result_MFI$`x of minimum` > 0.75)&&(n_clust < 100)){ #save the plot if the min is "extreme"
    if (n_clust < 100){
      ggplot2::ggsave(filename = paste("_CS_",stringa,"cluster_","nr_",as.character(i)," MFI_study.png"), 
                    plot = studio_MFI, device = "png", path = paste0("tmpdata/"))}
    #----------------------------------------------------------------------------------------------------------------------#
    # WARNING GENERATOR: the minimum is too close to the borders
    warning_msg <- paste("The MFI estimated for ", stringa, "cluster nr.",as.character(i), " is very close to the upper/lower
                         interval limit of the marker expressions")
    warning_msg_2 <- paste(warning_msg_2,warning_msg, sep = "\n")
    issue <- issue+1
    warning_df[nrow(warning_df) + 1,] <- c(as.character(issue), as.character(i), "-", "MFI close to the borders")
    #----------------------------------------------------------------------------------------------------------------------#
    }
  
  for (j in 1:length(color_marker)){
    #if ((i==67)&&(j==8)){rbowser()} #to debug selecting the analysis with issues
    marker_density <- dplyr::filter(.data = ld_dens, colour == color_marker[j])
    ld_list[[j]]$layer_data = dplyr::select(.data = marker_density, colour,y,x,density,scaled,ndensity,count)
    ld_list[[j]]$marker <- marker_name[j]
    ld_list[[j]]$number_of_event <- nr_event 
    ld_list[[j]]$median <- median(exprsub[,j])
    # for each marker the layer_data stores the index position of the x (which is the expression range) and the y (that is the expression value of the density)
    f <- data.frame(expression_range = ld_list[[j]]$layer_data$x, marker_density = ld_list[[j]]$layer_data$density)
    ld_list[[j]]$dens_plot <- ggplot2::ggplot(data = f, ggplot2::aes(x = expression_range, y = marker_density)) + 
      ggplot2::geom_area(linewidth=1.2, color = "blue", fill = "lightyellow") + 
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::labs(title = paste(names(color_marker[j]), "marker density",
                                  " of the ", stringa, "cluster nr.",as.character(i), "   ",stringa,"cluster weight =", 
                                  round(100*ld_list[[j]]$number_of_event/Nr_event, digits = 2),"%"), 
                    subtitle = paste(" Number of events: ", as.character(ld_list[[j]]$number_of_event))) 
    
    #df_dx <- finite.differences(x = ld_dens_plot_normalized$x, y = ld_dens_plot_normalized$y)
    df_dx <- finite.differences(x = ld_list[[j]]$layer_data$x, y = ld_list[[j]]$layer_data$density)
    #this is to scale the derivate. It does not matter the real values, but the zero of the derivate
    scaled_df_dx <- df_dx/max(df_dx)
    
    scaled_density <- ld_list[[j]]$layer_data$density / max(ld_list[[j]]$layer_data$density)
    
    f <- data.frame(expression_range = ld_list[[j]]$layer_data$x, marker_density = scaled_density, scaled_df_dx)
    ld_list[[j]]$scaled_df_dx <- scaled_df_dx
    
    ld_list[[j]]$dens_plot_con_la_derivata <- ggplot2::ggplot(data = f, mapping = ggplot2::aes(x = expression_range)) +
      ggplot2::geom_line(mapping = ggplot2::aes(y = marker_density), linewidth=1.2, color = "red") +
      ggplot2::geom_line(mapping = ggplot2::aes(y = scaled_df_dx), linewidth=1, color = "blue") +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
    
    vector_length <- length(ld_list[[j]]$layer_data$x) #see the layer_data vector length for each marker = 512 in case of geom_freqpoly (only) 
    result <- studio_di_funzione(x_granularity = vector_length, ld_list[[j]]$layer_data$density/max(ld_list[[j]]$layer_data$density), 
                                 derivate = ld_list[[j]]$scaled_df_dx, ilcluster = i, ilmarker = j) 
    ####################################################################### just for debugging purpose
    # This is the point of minimum in the sum of densities trend 
    if(result$`min found`){
      if ((result$`max1 found`)&&(!result$`max2 found`)){
        studio_plot <- ld_list[[j]]$dens_plot + 
          ggplot2::geom_vline(xintercept = result$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
          ggplot2::geom_vline(xintercept = result$`x of minimum`, linetype="dashed", color = "darkorchid", linewidth=1) +
          ggplot2::geom_vline(xintercept = new_MFI[i], linetype="twodash", color = "blue", linewidth=1.5) +
          ggplot2::geom_hline(yintercept = result$`density background`, linetype="dashed", color = "black", linewidth=1) +
          ggplot2::geom_hline(yintercept = max(ld_dens$density), linetype="dashed", color = "brown", linewidth=1) +
          ggplot2::geom_vline(xintercept = ld_list[[j]]$median, linetype="dashed", color = "green", linewidth=1) +
        
          ggplot2::geom_text(ggplot2::aes(x=result$`x of max1`, label=paste("x of max density","=",as.character(round(result$`x of max1`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(x=result$`x of minimum`, label=paste("x of min density","=",as.character(round(result$`x of minimum`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +      
          ggplot2::geom_text(ggplot2::aes(x=new_MFI[i], label=paste("estimated MFI","=",as.character(round(new_MFI[i], digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(y=result$`density background`, label=paste("density background","=",as.character(round(result$`density background`, digits = 2))), x=0.5))+
          ggplot2::geom_text(ggplot2::aes(y=max(ld_dens$density), label=paste("max marker density in the ",stringa,"cluster","=",as.character(round(max(ld_dens$density), digits = 2))), x=0.5))+
          ggplot2::geom_text(ggplot2::aes(x=ld_list[[j]]$median, label=paste("median of marker = ",as.character(round(ld_list[[j]]$median, digits = 2))), 
                                          y=max(ld_dens$density)/2), angle=90, vjust = 1.2)} 
      
      if ((result$`max1 found`)&&(result$`max2 found`)){
        studio_plot <- ld_list[[j]]$dens_plot + 
          ggplot2::geom_vline(xintercept = result$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
          ggplot2::geom_vline(xintercept = result$`x of minimum`, linetype="dashed", color = "darkorchid", linewidth=1) +
          ggplot2::geom_vline(xintercept = result$`x of max2`, linetype="dashed", color = "red", linewidth=1) +
          ggplot2::geom_vline(xintercept = new_MFI[i], linetype="twodash", color = "blue", linewidth=1.5) +
          ggplot2::geom_hline(yintercept = result$`density background`, linetype="dashed", color = "black", linewidth=1) +
          ggplot2::geom_hline(yintercept = max(ld_dens$density), linetype="dashed", color = "brown", linewidth=1) +
          ggplot2::geom_vline(xintercept = ld_list[[j]]$median, linetype="dashed", color = "green", linewidth=1) +
          
          ggplot2::geom_text(ggplot2::aes(x=result$`x of max1`, label=paste("x of max density","=",as.character(round(result$`x of max1`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(x=result$`x of minimum`, label=paste("x of min density","=",as.character(round(result$`x of minimum`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(x=result$`x of max2`, label=paste("x of 2nd max density","=",as.character(round(result$`x of max2`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1,2)) +
          ggplot2::geom_text(ggplot2::aes(x=new_MFI[i], label=paste("estimated MFI","=",as.character(round(new_MFI[i], digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(y=result$`density background`, label=paste("density background","=",as.character(round(result$`density background`, digits = 2))), x=0.5)) +
          ggplot2::geom_text(ggplot2::aes(y=max(ld_dens$density), label=paste("max marker density in the ", stringa, "cluster","=",as.character(round(max(ld_dens$density), digits = 2))), x=0.5))+
          ggplot2::geom_text(ggplot2::aes(x=ld_list[[j]]$median, label=paste("median of marker = ",as.character(round(ld_list[[j]]$median, digits = 2))), 
                                          y=max(ld_dens$density)/2), angle=90, vjust = 1.2)}}
    else{
      if ((result$`max1 found`)&&(!result$`max2 found`)){
        studio_plot <- ld_list[[j]]$dens_plot + 
          ggplot2::geom_vline(xintercept = result$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
          ggplot2::geom_vline(xintercept = new_MFI[i], linetype="twodash", color = "blue", linewidth=1.5) +
          ggplot2::geom_hline(yintercept = result$`density background`, linetype="dashed", color = "black", linewidth=1) +
          ggplot2::geom_hline(yintercept = max(ld_dens$density), linetype="dashed", color = "brown", linewidth=1) +
          ggplot2::geom_vline(xintercept = ld_list[[j]]$median, linetype="dashed", color = "green", linewidth=1) +
          
          ggplot2::geom_text(ggplot2::aes(x=result$`x of max1`, label=paste("x of max density","=",as.character(round(result$`x of max1`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(x=new_MFI[i], label=paste("estimated MFI","=",as.character(round(new_MFI[i], digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(y=result$`density background`, label=paste("density background","=",as.character(round(result$`density background`, digits = 2))), x=0.5))+
          ggplot2::geom_text(ggplot2::aes(y=max(ld_dens$density), label=paste("max marker density in the ", stringa, "cluster","=",as.character(round(max(ld_dens$density), digits = 2))), x=0.5))+
          ggplot2::geom_text(ggplot2::aes(x=ld_list[[j]]$median, label=paste("median of marker = ",as.character(round(ld_list[[j]]$median, digits = 2))), 
                                          y=max(ld_dens$density)/2), angle=90, vjust = 1.2)} 
      
      if ((result$`max1 found`)&&(result$`max2 found`)){
        studio_plot <- ld_list[[j]]$dens_plot + 
          ggplot2::geom_vline(xintercept = result$`x of max1`, linetype="dashed", color = "red", linewidth=1) +
          ggplot2::geom_vline(xintercept = result$`x of max2`, linetype="dashed", color = "red", linewidth=1) +
          ggplot2::geom_vline(xintercept = new_MFI[i], linetype="twodash", color = "blue", linewidth=1.5) +
          ggplot2::geom_hline(yintercept = result$`density background`, linetype="dashed", color = "black", linewidth=1) +
          ggplot2::geom_hline(yintercept = max(ld_dens$density), linetype="dashed", color = "brown", linewidth=1) +
          ggplot2::geom_vline(yintercept = ld_list[[j]]$median, linetype="dashed", color = "green", linewidth=1) +
          ggplot2::geom_text(ggplot2::aes(x=result$`x of max1`, label=paste("x of max density","=",as.character(round(result$`x of max1`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(x=result$`x of max2`, label=paste("x of 2nd max density","=",as.character(round(result$`x of max2`, digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1,2)) +
          ggplot2::geom_text(ggplot2::aes(x=new_MFI[i], label=paste("estimated MFI","=",as.character(round(new_MFI[i], digits = 2))), 
                                          y=max(ld_dens$density)/2, angle=90, vjust = 1.2)) +
          ggplot2::geom_text(ggplot2::aes(y=result$`density background`, label=paste("density background","=",as.character(round(result$`density background`, digits = 2))), x=0.5)) +
          ggplot2::geom_text(ggplot2::aes(y=max(ld_dens$density), label=paste("max marker density in the ", stringa, "cluster","=",as.character(round(max(ld_dens$density), digits = 2))), x=0.5))+
          ggplot2::geom_text(ggplot2::aes(x=ld_list[[j]]$median, label=paste("median of marker = ",as.character(round(ld_list[[j]]$median, digits = 2))), 
                                          y=max(ld_dens$density)/2), angle=90, vjust = 1.2)}
    }#else 
    
    result_marker_list[[j]] <- result
    
    if (nr_event < MIN_EVENTS_TO_TRUST_PLOT_REPORT){
      #----------------------------------------------------------------------------------------------------------------------#
      # WARNING GENERATOR: if the nr of event is less than MIN_EVENTS_TO_TRUST_PLOT_REPORT
      warning_msg <- paste("The number of event is less than", as.character(MIN_EVENTS_TO_TRUST_PLOT_REPORT), "in ", stringa,
                           "cluster nr.", as.character(i), "for marker ", marker_name[j],". Check using manual gating")
      warning_msg_1 <- paste(warning_msg_1,warning_msg, sep = "\n")
      issue <- issue+1
      warning_df[nrow(warning_df) + 1,] <- c(as.character(issue+1), as.character(i), "-", "too few events")
    }
    #----------------------------------------------------------------------------------------------------------------------#
    if (n_clust < 100){ #the plot is generated only if the number of cluster is less than 100 (flowSOM)
    ggplot2::ggsave(filename = paste0(stringa,"cluster_nr_",as.character(i),"_density_of_marker_",as.character(j),".png"), 
                    plot = studio_plot, device = "png", path = ("tmpdata/"))}
    # this was inside the condition below: this is the attempt to show the density of the marker anyway
    
    if ((result$`max1 found`)&&(result$`max2 found`)){
      #ggplot2::ggsave(filename = paste0(stringa,"cluster_nr_",as.character(i),"_density_of_marker_",as.character(j),".png"), 
      #                plot = studio_plot, device = "png", path = ("tmpdata/"))
      #----------------------------------------------------------------------------------------------------------------------#
      # WARNING GENERATOR: if there are more than one single peak in the marker density
      warning_msg <- paste("There are more than a single density peak in ", stringa, "cluster nr.",as.character(i), 
                             ". Check the related density plots for marker ", marker_name[j])
      warning_msg_3 <- paste(warning_msg_3, warning_msg, sep = "\n")
      issue <- issue+1
      warning_df[nrow(warning_df) + 1,] <- c(as.character(issue), as.character(i), marker_name[j], "more than a single density peak")
      }
    #----------------------------------------------------------------------------------------------------------------------#
    
    
    if (((abs(new_MFI[i]- result$`x of max1`) < MIN_PEAK_DISTANCE_TO_SHOW_PLOT_REPORT)&&
         (result$`max1 found`))||((abs(new_MFI[i]- result$`x of max2`)< MIN_PEAK_DISTANCE_TO_SHOW_PLOT_REPORT)&&
                                  (result$`max2 found`))){
      #ggplot2::ggsave(filename = paste(stringa,"cluster_nr_",as.character(i),"_density_of_marker",as.character(j),".png"), 
      #                plot = studio_plot, device = "png", path = ("tmpdata/"))
      
      #----------------------------------------------------------------------------------------------------------------------#
      # WARNING GENERATOR: if the estimated MFI is closer than 3% to one on the maximum 
      warning_msg <- paste("The estimated normalized MFI is closer than", as.character(MIN_PEAK_DISTANCE_TO_SHOW_PLOT_REPORT * 100), 
                             "% to the single density peak in ", stringa,"cluster nr.",as.character(i), "for marker ", marker_name[j], 
                             ". Check using manual gating")
      warning_msg_4 <- paste(warning_msg_4,warning_msg, sep = "\n")
      issue <- issue+1
      warning_df[nrow(warning_df) + 1,] <- c(as.character(issue), as.character(i), marker_name[j], "MFI too close to the maximum")
      }
    #----------------------------------------------------------------------------------------------------------------------#
    
    
  } # for each marker
  result_list[[i]] <- result_marker_list} #for each cluster

logic_matrix <- matrix(data = 0, nrow = n_clust, ncol = length(color_marker))
rownames(logic_matrix) <- 1:n_clust
colnames(logic_matrix) <- names(color_marker)
nr_of_plots <- 0

# To plot an heatmap with the "red" double peak cells - 
# see https://stackoverflow.com/questions/3789549/display-a-matrix-including-the-values-as-a-heatmap
for (i in 1:n_clust){
  for (j in 1:length(color_marker))
  {if (result_list[[i]][[j]]$`max2 found` == TRUE){
    nr_of_plots = nr_of_plots + 1
    logic_matrix[i,j] <- 1} else {logic_matrix[i,j] <- 0}
  } # end of marker
} #end of clust

##### this is to plot the double peaks heat map
#----------------------------------------------------------------------------------------------------------------------#
# put the clean_graph() function from the helpers.R

#clean_graph()
#file_png_name = paste0("./tmpdata/Non_Single_Peak.png")
#png(filename = file_png_name, width = 1300, height = 1000, units = "px", res = 125)
test <- unique(as.vector(logic_matrix))
if(length(test)!=1){png(filename = paste0("./tmpdata/_FS_not_unique_peak_per_marker_per_",stringa,"cluster.png"), width = 1200, height = 900, units = "px", res = 125)
  write.csv(x = logic_matrix, file = paste("./tmpdata/_FS_not_unique_peak_per_marker_per_",stringa,"cluster.csv"))
  pheatmap::pheatmap(mat = logic_matrix, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE, 
                   main = "Non single density peak for each cluster x marker", color = c("white", "red"))
  dev.off()}

#----------------------------------------------------------------------------------------------------------------------#

# Notice that the more difficult markers (e.g. not so sharp density, tetramer, ...) do not have any red cell 
#(tht is, no double peaks on them)

#-----------------------------------------------------------------------------------------------------------------------#
# Here is to build the expression matrix with the level of expression x clusters and for each marker
# In the original implementation (without this density study) is was the "expr01_median" table 
# The expr01_peak table in this case reports the peaks with the center value (0.5) "re-centered" with the calculated MFI.
# At the end: the same distance (on the normalized x axis) between the x_peak and the x_new_MFI should be the same that the
# the distance between the new calculated new_x_peak and 0.5 that is: delta = x_new_MFI - 0.5 so that the new_peak is 
# calculated with new_peak = peak - delta
# The expr_log table has the same meaning of the expr_log in the phenocluster function

colnames(expr01_peak) <- marker_name 
colnames(expr01_no_delta) <- marker_name 
colnames(expr_log) <- marker_name 


###### Expression matrix with the estimated MFI ###########################################
# This is to build a new expression matrix based on the estimated MFI. The estimated MFI is the new center of expression (0.5)
#
# Notice: the pos_threshold and neg_threshold should be taken into account
# 

for (i in 1:n_clust){
  delta <- new_MFI[i] - 0.5
  for (j in 1:length(color_marker)){
    if ((result_list[[i]][[j]]$`max1 found`)&&(!result_list[[i]][[j]]$`max2 found`)){ # case of only one max found
      new_x_peak <- result_list[[i]][[j]]$`x of max1` - delta #delta could be positive or negative
      new_x_no_delta <-  result_list[[i]][[j]]$`x of max1`}   
    if ((result_list[[i]][[j]]$`max1 found`)&&(result_list[[i]][[j]]$`max2 found`)){
      if ((result_list[[i]][[j]]$'max1')>(result_list[[i]][[j]]$'max2')){new_x_peak <- result_list[[i]][[j]]$`x of max1` - delta}
      else{new_x_peak <- result_list[[i]][[j]]$`x of max2` - delta
           new_x_peak_no_delta <- result_list[[i]][[j]]$`x of max2`}} # case of both max found the highest wins
    if (new_x_peak < 0) {new_x_peak <- 0}
    if (new_x_no_delta < 0) {new_x_no_delta <- 0}
    if (new_x_peak > 1) {new_x_peak <- 1}
    if (new_x_no_delta > 1) {new_x_no_delta <- 1}
    expr01_peak[[j]][[i]] <- new_x_peak
    expr01_no_delta[[j]][[i]] <- new_x_no_delta
    if (new_x_peak >= pos_threshold){expr_log[[j]][[i]] <- 1}
    if (new_x_peak < neg_threshold){expr_log[[j]][[i]] <- 0}
  }
}

expr_log[expr_log == 1] <- "+"
expr_log[expr_log == 0] <- "-"


expr_log <- cbind(cell_clustering = sequence(n_clust), expr_log)
warning_text <- paste(warning_msg_1,warning_msg_2, warning_msg_3, warning_msg_4, sep = "\n") 

expr01_peak_matrix <- as.matrix(expr01_peak)
expr01_no_delta_matrix <- as.matrix(expr01_no_delta)

rownames(expr01_peak_matrix) <- 1:n_clust
rownames(expr01_no_delta_matrix) <- 1:n_clust
#expr01_peak_matrix <- cbind(cell_clustering = sequence(n_clust), expr01_peak_matrix)

if(length(test)!=1){png(filename = paste0("./tmpdata/new_HM_",stringa,".png"), width = 1200, height = 900, units = "px", res = 125)
  pheatmap::pheatmap(mat = expr01_peak_matrix, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                   main = paste("New peak based",stringa,"cluster heatmap"), display_numbers = TRUE, number_color = "black", fontsize_number = 6, 
                   fontsize_row = 7, fontsize_col = 7)
  dev.off()}

if(length(test)!=1){png(filename = paste0("./tmpdata/new_HM_no_delta",stringa,".png"), width = 1200, height = 900, units = "px", res = 125)
  pheatmap::pheatmap(mat = expr01_no_delta_matrix, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                     main = paste("New peak based",stringa,"cluster heatmap"), display_numbers = TRUE, number_color = "black", fontsize_number = 6, 
                     fontsize_row = 7, fontsize_col = 7)
  dev.off()}

write.csv(as.data.frame(new_MFI),file = paste("./tmpdata/",stringa,"cluster_new_MFI.csv"), row.names = TRUE)
write.csv(as.data.frame(expr_log),file = paste("./tmpdata/",stringa,"cluster_expr_log.csv"), row.names = FALSE)
write.csv(as.data.frame(expr01_peak),file = paste("./tmpdata/",stringa,"cluster_expr01_peak.csv"), row.names = TRUE)
write.csv(as.data.frame(expr01_no_delta),file = paste("./tmpdata/",stringa,"cluster_expr01_no_delta.csv"), row.names = TRUE)
write.csv(warning_df,file = paste("./tmpdata/",stringa,"cluster_warning.csv"), row.names = FALSE)

# use cat(warning_text) to show with <CR>

############## Try to organize the warnings in a table. See http://www.sthda.com/english/wiki/writing-data-from-r-to-txt-csv-files-r-base-functions

lista_output <- list(expr01_peak, expr_log, warning_df)
return(lista_output)}
#-------------------------------------------------------------------------------------------------------------------------#
#This is to compose the warning strings

# passing variables using \n
#cat(string1,"\n",string2,"\n",string3) see: https://www.geeksforgeeks.org/r-program-to-print-a-new-line-in-string/

#lista_output <- list(result_list, result_plot, logic_matrix, heat_map)  
#return(lista_output)
#end of function


plot_labelling <- function(flow_Set, phenoclust, clustering, m_clustering, map_res, map_cell, 
                         selected_markers, metadata, color_clusters, save = NULL, work = NULL) {

  cluster_labelling <- as.data.frame(phenoclust)
  
  mm <- match(map_res$sample_id, metadata$sample_id)
  map_res$tag1 <- metadata[,3][mm]
  map_res$tag2 <- metadata[,4][mm]
  map_res$tag3 <- metadata[,5][mm]
  map_res$tag4 <- metadata[,6][mm]

  map_res$cell_label <- map_cell # add cell_label to map_df
  
  map_res <- na.omit(map_res) #introduced after the mapping limitation
  if ((nrow(map_res)) > RENDER_DATA) {
    map_res <- sample_n(tbl = map_res, size = RENDER_DATA, replace = F)}
  
  gg <- ggplot(map_res,  aes(x = map_x, y = map_y, color = cell_label)) +
    geom_point(size = 0.8) +
    theme_bw() +
    scale_color_manual(values = color_clusters) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  gg_map <- ggplotly(gg)
  
  if (!is.null(save)){
    if (work==1) {ggsave(filename = "map_cluster_wf1.png", plot = gg, device = "png", path = "./tmpdata/")}
    if (work==2) {ggsave(filename = "map_cluster_wf2.png", plot = gg, device = "png", path = "./tmpdata/")}
    if (work==3) {ggsave(filename = "map_cluster_wf3.png", plot = gg, device = "png", path = "./tmpdata/")}}
  
  gg_map_sample <- gg + facet_wrap(~ sample_id) ## Facet per sample
  if (!is.null(save)){
    if (work==1) {ggsave(filename = "map_sample_wf1.png", plot = gg_map_sample, device = "png", path = "./tmpdata/")}
    if (work==2) {ggsave(filename = "map_sample_wf2.png", plot = gg_map_sample, device = "png", path = "./tmpdata/")}
    if (work==3) {ggsave(filename = "map_sample_wf3.png", plot = gg_map_sample, device = "png", path = "./tmpdata/")}}
  
  gg_map_tag1 <- gg + facet_wrap(~ tag1) ## Facet per condition
  if ((!is.null(save))&&(uni.tag1 > 1)){
    nome_file <- paste0("map_", names(metadata)[3])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_map_tag1, device = "png", path = "./tmpdata/")}
  gg_map_tag2 <- gg + facet_wrap(~ tag2) ## Facet per patient
  if ((!is.null(save))&&(uni.tag2 > 1)){
    nome_file <- paste0("map_", names(metadata)[4])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_map_tag2, device = "png", path = "./tmpdata/")}
  gg_map_tag3 <- gg + facet_wrap(~ tag3) ## Facet per time_Step
  if ((!is.null(save))&&(uni.tag3 > 1)){
    nome_file <- paste0("map_", names(metadata)[5])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_map_tag3, device = "png", path = "./tmpdata/")}
  gg_map_tag4 <- gg + facet_wrap(~ tag4) ## Facet per time_Step
  if ((!is.null(save))&&(uni.tag4 > 1)){
    nome_file <- paste0("map_", names(metadata)[6])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_map_tag4, device = "png", path = "./tmpdata/")}
  lista <- list(gg_map, gg_map_sample, gg_map_tag1, gg_map_tag2, gg_map_tag3, gg_map_tag4)  
  return(lista)
}


join_tag <- function(frame, meta, type, tag1col=NULL, tag2col=NULL, tag3col=NULL, tag4col=NULL) {
  # This is to add properly to the concatenated flowFrame, the dimensions related to each of the tag defined in the meta_sample data (meta)
  # and to plot the quantity value with the related csv saved in the tmpdata directory
  
  #handling data entries
  if (type == "meta"){
    df <- as.data.frame(exprs(frame)[,c("SampleID", "clusterId", "meta_clusterId")])} # all numeric as a class
  else
    {df <- as.data.frame(exprs(frame)[,c("SampleID", "clusterId")])}

  meta$SampleID <- as.numeric(1:nrow(meta))
  rj <- right_join(meta,df)
  rj <- rj[,-(1:2)]
  rj$date <- NULL
  nomi_tag <- colnames(meta)
  nomi_tag <- nomi_tag[-(1:2)]
  ggplot_tag <- vector(mode = "list", length = length(nomi_tag) + 1)
  ggplot_meta_tag <- vector(mode = "list", length = length(nomi_tag) + 1)
  
  for (i in 1:4){rj[,i] <- as.factor(as.numeric(factor(rj[,i])))}
  rj$SampleID <- factor(rj$SampleID)
  rj$clusterId <- factor(rj$clusterId)
  if (type == "meta"){rj$meta_clusterId <- factor(rj$meta_clusterId)}
  meta_factor <- meta
  #for (i in 1:4){meta_factor[,i+2] <- (as.numeric(factor(meta_factor[,i+2])))}
  for (i in 1:4){meta[,i+2] <- (as.numeric(factor(meta[,i+2])))}
  
  #sample by clusterId
  sample_rj <- rj %>% group_by(clusterId,SampleID) %>% summarise(count = n())
  sample_rj <- sample_rj %>% ungroup %>% complete(clusterId, SampleID, fill = list(count=0))
  temp_rj <- sample_rj %>% group_by(clusterId) %>% mutate(perc=paste0(round(count/sum(count)*100, 2), "%")) %>% ungroup()
  write.csv(x = temp_rj, file = paste0("./tmpdata/clusterId_", "sample",".csv"))
  sample_rj <- sample_rj %>% group_by(clusterId) %>% mutate(perc=paste0(round(count/sum(count)*100, 2))) %>% ungroup()
  
  #the following shame is because I cannot manage to keep all the columns when a group_by
  for (k in 1:4){
    sample_rj[,nomi_tag[k]] <- as.numeric(NA)
    for (i in 1:nrow(meta)){
      sample_name <- meta$SampleID[i]
      tag_name <- meta[sample_name,nomi_tag[k]]
      (idx <- which(sample_rj[,"SampleID"] == sample_name))
      for (j in 1:length(idx)){
        sample_rj[idx[j],nomi_tag[k]]<-tag_name}}
    #cast_rj <- dcast(data = sample_rj, formula = clusterId ~ condition + SampleID, value.var = "perc")
    foo <- paste0("dcast(data = sample_rj, formula = clusterId ~ ",nomi_tag[k], " + SampleID, value.var = 'perc')")
    parsed <- parse(text = foo)
    cast_rj <- eval(expr = parsed)
    filename <- paste0("./tmpdata/clusterId_",nomi_tag[k],"_cast")
    write.csv(x = cast_rj, file = paste0(filename, ".csv"))}
  
  sample_rj$perc <- as.numeric(sample_rj$perc)
  plot_title <- paste0("cluster quantity distribution for each sample")
  gg_plot_sample <- ggplot(sample_rj, aes(x = clusterId, y = perc, fill = SampleID)) + geom_bar(stat = 'identity') + 
    scale_fill_manual(values = color.sample) + ggtitle(plot_title)
  ggplot_tag[[1]] <- gg_plot_sample
  ggsave(filename = paste0("clusterId_sample.png"), plot = ggplot_tag[[1]], device = "png", path = "./tmpdata/")
  
  #sample by meta_clusterId
  if (type == "meta"){ 
    sample_rj <- rj %>% group_by(meta_clusterId,SampleID) %>% summarise(count = n())
    sample_rj <- sample_rj %>% ungroup %>% complete(meta_clusterId, SampleID, fill = list(count=0))
    temp_rj <- sample_rj %>% group_by(meta_clusterId) %>% mutate(perc=paste0(round(count/sum(count)*100,2), "%")) %>% ungroup()
    write.csv(x = temp_rj, file = paste0("./tmpdata/meta_clusterId_", "sample",".csv"))
    sample_rj <- sample_rj %>% group_by(meta_clusterId) %>% mutate(perc=paste0(round(count/sum(count)*100, 2))) %>% ungroup()
    
    for (k in 1:4){
      sample_rj[,nomi_tag[k]] <- as.numeric(NA)
      for (i in 1:nrow(meta)){
        sample_name <- meta$SampleID[i]
        tag_name <- meta[sample_name,nomi_tag[k]]
        (idx <- which(sample_rj[,"SampleID"] == sample_name))
        for (j in 1:length(idx)){
          sample_rj[idx[j],nomi_tag[k]]<-tag_name}}
      #cast_rj <- dcast(data = sample_rj, formula = clusterId ~ condition + SampleID, value.var = "perc")
      foo <- paste0("dcast(data = sample_rj, formula = meta_clusterId ~ ",nomi_tag[k], " + SampleID, value.var = 'perc')")
      parsed <- parse(text = foo)
      cast_rj <- eval(expr = parsed)
      filename <- paste0("./tmpdata/meta_clusterId_",nomi_tag[k],"_cast")
      write.csv(x = cast_rj, file = paste0(filename, ".csv"))}
    
    foo <- paste0("dcast(data = sample_rj, formula = meta_clusterId ~ ",nomi_tag[1], " + SampleID, value.var = 'perc')")
    parsed <- parse(text = foo)
    cast_rj <- eval(expr = parsed)
    write.csv(x = cast_rj, file = paste0("./tmpdata/cast_meta_clusterId_", "sample",".csv"))
    
    sample_rj$perc <- as.numeric(sample_rj$perc)
    plot_title <- paste0("meta_cluster quantity distribution for each sample")
    gg_plot_sample <- ggplot(sample_rj, aes(x = meta_clusterId, y = perc, fill = SampleID)) + geom_bar(stat = 'identity') + 
      scale_fill_manual(values = color.sample) + ggtitle(plot_title)
    ggplot_meta_tag[[1]] <- gg_plot_sample
    ggsave(filename = paste0("meta_clusterId_sample.png"), plot = ggplot_meta_tag[[1]], device = "png", path = "./tmpdata/")}
  
  tag_pro <- function(tag_nr,tagcol){
    
    tag <-  nomi_tag[tag_nr]
    tagdata <- as.numeric(factor(rj[,tag_nr]))
    #frame <- add_dim(flow_Frame = frame, dim_name = tag, dim_vect = tagdata)
    #sample by cluster
    tag_rj <- rj %>% group_by(clusterId, rj[,tag_nr]) %>% summarise(count = n())
    names(tag_rj)[2] <- tag #it is always the second 
    tag_rj <- tag_rj %>% ungroup %>% complete(clusterId, tag_rj[,2], fill = list(count=0))
    tag_rj <- tag_rj %>% group_by(clusterId) %>% mutate(perc=paste0(round(count/sum(count)*100, 2), "%")) %>% ungroup
    write.csv(x = tag_rj, file = paste0("./tmpdata/clusterId_", tag,".csv"))
    
    tag_rj <- tag_rj %>% group_by(clusterId) %>% mutate(perc = paste0(round(count/sum(count)*100, 2))) %>% ungroup()
    tag_rj$perc <- as.numeric(tag_rj$perc)
    
    plot_title <- paste0("cluster quantity distribution for the ", tag, " tag")
    foo <- paste0("ggplot(tag_rj, aes(x = clusterId, y = perc, fill = ",tag, ")) + ")
    foo1 <- paste0(foo, "geom_bar(stat = 'identity', position = 'dodge') + ")
    foo2 <- paste0(foo1, "scale_fill_manual(values = tagcol) + ggtitle(plot_title)")
    parsed <- parse(text = foo2)
    ggplot_tag[[tag_nr+1]] <- eval(expr = parsed)
    ggsave(filename = paste0("clusterId_dodge_", tag,".png"), plot = ggplot_tag[[tag_nr+1]], device = "png", path = "./tmpdata/")
    
    foo1 <- paste0(foo, "geom_bar(stat = 'identity') + ")
    foo2 <- paste0(foo1, "scale_fill_manual(values = tagcol) + ggtitle(plot_title)")
    parsed <- parse(text = foo2)
    ggplot.tag <- eval(expr = parsed)
    ggsave(filename = paste0("clusterId_", tag,".png"), plot = ggplot.tag, device = "png", path = "./tmpdata/")
    
    #sample by meta_cluster
    if (type == "meta"){
      tag_rj <- rj %>% group_by(meta_clusterId, rj[,tag_nr]) %>% summarise(count = n())
      names(tag_rj)[2] <- tag
      tag_rj <- tag_rj %>% ungroup %>% complete(meta_clusterId, tag_rj[,2], fill = list(count=0))
      tag_rj <- tag_rj %>% group_by(meta_clusterId) %>% mutate(perc=paste0(round(count/sum(count)*100, 2), "%")) %>% ungroup
      write.csv(x = tag_rj, file = paste0("./tmpdata/meta_clusterId_", tag,".csv"))
      
      tag_rj <- tag_rj %>% group_by(meta_clusterId) %>% mutate(perc = paste0(round(count/sum(count)*100, 2))) %>% ungroup()
      tag_rj$perc <- as.numeric(tag_rj$perc)
      
      plot_title <- paste0("meta_cluster quantity distribution for the ", tag, " tag")
      foo <- paste0("ggplot(tag_rj, aes(x = meta_clusterId, y = perc, fill = ",tag, ")) + ")
      foo1 <- paste0(foo, "geom_bar(stat = 'identity') + ")
      foo2 <- paste0(foo1, "scale_fill_manual(values = tagcol) + ggtitle(plot_title)")
      parsed <- parse(text = foo2)
      ggplot_meta.tag <- eval(expr = parsed)
      ggsave(filename = paste0("meta_clusterId_", tag,".png"), plot = ggplot_meta.tag, device = "png", path = "./tmpdata/")
      
      foo1 <- paste0(foo, "geom_bar(stat = 'identity', position = 'dodge') + ")
      foo2 <- paste0(foo1, "scale_fill_manual(values = tagcol) + ggtitle(plot_title)")
      parsed <- parse(text = foo2)
      ggplot_meta.tag <- eval(expr = parsed)
      ggsave(filename = paste0("meta_clusterId_dodge_", tag,".png"), plot = ggplot_meta.tag, device = "png", path = "./tmpdata/")}
    if (type == "meta"){lista <- list(ggplot.tag, ggplot_meta.tag)}else{lista <- list(ggplot.tag, NULL)}
    return <- lista}
  
  if (uni.tag1 > 1){
    res <- tag_pro(tag_nr = 1, tagcol = tag1col); ggplot_tag[[1+1]] <- res[[1]]; 
    if (type == "meta"){ggplot_meta_tag[[1+1]] <- res[[2]]}}
  if (uni.tag2 > 1){
    res <- tag_pro(tag_nr = 2, tagcol = tag2col); ggplot_tag[[1+2]] <- res[[1]]; 
    if (type == "meta"){ggplot_meta_tag[[1+2]] <- res[[2]]}}
  if (uni.tag3 > 1){
    res <- tag_pro(tag_nr = 3, tagcol = tag3col); ggplot_tag[[1+3]] <- res[[1]]; 
    if (type == "meta"){ggplot_meta_tag[[1+3]] <- res[[2]]}}
  if (uni.tag4 > 1){
    res <- tag_pro(tag_nr = 4, tagcol = tag4col); ggplot_tag[[1+4]] <- res[[1]]; 
    if (type == "meta"){ggplot_meta_tag[[1+4]] <- res[[2]]}}
  if (type == "meta"){lista <- list(ggplot_tag, ggplot_meta_tag)}else{lista <- list(ggplot_tag)}
  return(lista)
}


plot_labellingbar <- function(flow_Set, flag, cell_label, NMC, metadata, color_pheno, save = NULL, work = NULL) {

  ################## sample_id
  
  sample_ids <- rep(metadata$sample_id, fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, sample_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  #For each sample, we plot its PBMC cell type composition represented with colored bars, 
  #where the size of a given stripe reflects the proportion of the corresponding cell type in a given sample

  if (flag)
    {cluster_list <- factor(rownames(props), levels = 1:NMC)}
  else
    {cluster_list <- rownames(props)}
  
  ggdf_sample <- melt(data.frame(cluster = cluster_list, props), 
                      id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
  
  gg_sample <- ggplot(ggdf_sample, aes(x = sample_id, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    scale_fill_manual(values = color_pheno) +
    ggtitle("Sample")
  
  gg_sample_dodge <- ggplot(ggdf_sample, aes(x = sample_id, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    scale_fill_manual(values = color_pheno) +
    ggtitle("Sample")
  
  if (!is.null(save)){
    if (work == 1){ggsave(filename = "bar_sample_wf1.png", plot = gg_sample, device = "png", path = "./tmpdata/")}
    if (work == 2){ggsave(filename = "bar_sample_wf2.png", plot = gg_sample, device = "png", path = "./tmpdata/")}
    if (work == 3){ggsave(filename = "bar_sample_wf3.png", plot = gg_sample, device = "png", path = "./tmpdata/")}
    if (work == 1){ggsave(filename = "dodge_sample_wf1.png", plot = gg_sample_dodge, device = "png", path = "./tmpdata/")}
    if (work == 2){ggsave(filename = "dodge_sample_wf2.png", plot = gg_sample_dodge, device = "png", path = "./tmpdata/")}
    if (work == 3){ggsave(filename = "dodge_sample_wf3.png", plot = gg_sample_dodge, device = "png", path = "./tmpdata/")}}
  bar_plot_sample <- ggplotly(gg_sample)
  dodge_plot_sample <- ggplotly(gg_sample_dodge)
  bar_plot_tag1 <- NULL
  bar_plot_tag2 <- NULL
  bar_plot_tag3 <- NULL
  bar_plot_tag4 <- NULL
  dodge_plot_tag1 <- NULL
  dodge_plot_tag2 <- NULL
  dodge_plot_tag3 <- NULL
  dodge_plot_tag4 <- NULL
  
  ################## tag1
  
  if(uni.tag1>1){
  tag1_ids <- rep(metadata[,3], fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, tag1_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  #For each sample, we plot its PBMC cell type composition represented with colored bars, 
  #where the size of a given stripe reflects the proportion of the corresponding cell type in a given sample
  
  if (flag)
    {cluster_list <- factor(rownames(props), levels = 1:NMC)}
  else
    {cluster_list <- rownames(props)}
  
  ggdf_tag1 <- melt(data.frame(cluster = cluster_list, props), 
                      id.vars = "cluster", value.name = "proportion", variable.name = names(metadata)[3])
  if (is.null(work))
  {filename <- paste0("./tmpdata/df_",names(metadata)[3],".csv")}
  else {filename <- paste0("./tmpdata/df_",names(metadata)[3],"_wf",as.character(work),".csv")}
  write.csv(x = ggdf_tag1, file = filename)
  
  plot_title <- paste0(names(metadata)[3])
  tag1 <-  names(ggdf_tag1)[2]
  foo <- paste0("ggplot(ggdf_tag1, aes(x = ", tag1,", y = proportion, fill = cluster)) + ")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity') + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1, "scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag1 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag1 > 1)){
    nome_file <- paste0("bar_", names(metadata)[3])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag1, device = "png", path = "./tmpdata/")}
  bar_plot_tag1 <- ggplotly(gg_tag1)
  
  plot_title <- paste0(names(metadata)[3])
  tag1 <-  names(ggdf_tag1)[2]
  foo <- paste0("ggplot(ggdf_tag1, aes(x = ", tag1,", y = proportion, fill = cluster)) + ")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity', position = 'dodge')  + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1, "scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag1 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag1 > 1)){
    nome_file <- paste0("dodge_", names(metadata)[3])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag1, device = "png", path = "./tmpdata/")}
  dodge_plot_tag1 <- ggplotly(gg_tag1)}
  
  ################## tag2
  
  if(uni.tag2>1){
  tag2_ids <- rep(metadata[,4], fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, tag2_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  #For each sample, we plot its PBMC cell type composition represented with colored bars, 
  #where the size of a given stripe reflects the proportion of the corresponding cell type in a given sample
  
  if (flag)
    {cluster_list <- factor(rownames(props), levels = 1:NMC)}
  else
    {cluster_list <- rownames(props)}
  
  ggdf_tag2 <- melt(data.frame(cluster = cluster_list, props), 
                      id.vars = "cluster", value.name = "proportion", variable.name = names(metadata)[4])
  
  if (is.null(work))
  {filename <- paste0("./tmpdata/df_",names(metadata)[4],".csv")}
  else {filename <- paste0("./tmpdata/df_",names(metadata)[4],"_wf",as.character(work),".csv")}
  write.csv(x = ggdf_tag2, file = filename)
  
  plot_title <- paste0(names(metadata)[4])
  tag2 <-  names(ggdf_tag2)[2]
  foo <- paste0("ggplot(ggdf_tag2, aes(x = ", tag2,", y = proportion, fill = cluster)) +")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity') + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1,"scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag2 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag2 > 1)){
    nome_file <- paste0("bar_", names(metadata)[4])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag2, device = "png", path = "./tmpdata/")}
  bar_plot_tag2 <- ggplotly(gg_tag2)
  
  plot_title <- paste0(names(metadata)[4])
  tag2 <-  names(ggdf_tag2)[2]
  foo <- paste0("ggplot(ggdf_tag2, aes(x = ", tag2,", y = proportion, fill = cluster)) +")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity', position = 'dodge') + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1,"scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag2 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag2 > 1)){
    nome_file <- paste0("dodge_", names(metadata)[4])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag2, device = "png", path = "./tmpdata/")}
  dodge_plot_tag2 <- ggplotly(gg_tag2)}

  ################## tag_3
  
  if(uni.tag3>1){
  tag3_ids <- rep(metadata[,5], fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, tag3_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  #For each sample, we plot its PBMC cell type composition represented with colored bars, 
  #where the size of a given stripe reflects the proportion of the corresponding cell type in a given sample
  
  if (flag)
    {cluster_list <- factor(rownames(props), levels = 1:NMC)}
  else
    {cluster_list <- rownames(props)}
  
  ggdf_tag3 <- melt(data.frame(cluster = cluster_list, props), 
                      id.vars = "cluster", value.name = "proportion", variable.name = names(metadata)[5])
  
  if (is.null(work))
  {filename <- paste0("./tmpdata/df_",names(metadata)[5],".csv")}
  else {filename <- paste0("./tmpdata/df_",names(metadata)[5],"_wf",as.character(work),".csv")}
  write.csv(x = ggdf_tag3, file = filename)
  
  plot_title <- paste0(names(metadata)[5])
  tag3 <-  names(ggdf_tag3)[2]
  foo <- paste0("ggplot(ggdf_tag3, aes(x = ", tag3,", y = proportion, fill = cluster)) +")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity') + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1, "scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag3 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag3 > 1)){
    nome_file <- paste0("bar_", names(metadata)[5])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag3, device = "png", path = "./tmpdata/")}
  bar_plot_tag3 <- ggplotly(gg_tag3)
  
  plot_title <- paste0(names(metadata)[5])
  tag3 <-  names(ggdf_tag3)[2]
  foo <- paste0("ggplot(ggdf_tag3, aes(x = ", tag3,", y = proportion, fill = cluster)) +")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity', position = 'dodge') + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1, "scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag3 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag3 > 1)){
    nome_file <- paste0("dodge_", names(metadata)[5])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag3, device = "png", path = "./tmpdata/")}
  dodge_plot_tag3 <- ggplotly(gg_tag3)}
  
  ################## tag4 - time_step
  
  if(uni.tag4>1){
  tag4_ids <- rep(metadata[,6], fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, tag4_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  #For each sample, we plot its PBMC cell type composition represented with colored bars, 
  #where the size of a given stripe reflects the proportion of the corresponding cell type in a given sample
  
  if (flag)
  {cluster_list <- factor(rownames(props), levels = 1:NMC)}
  else
  {cluster_list <- rownames(props)}
  
  ggdf_tag4 <- melt(data.frame(cluster = cluster_list, props), 
                    id.vars = "cluster", value.name = "proportion", variable.name = names(metadata)[6])
  
  if (is.null(work))
  {filename <- paste0("./tmpdata/df_",names(metadata)[6],".csv")}
  else {filename <- paste0("./tmpdata/df_",names(metadata)[6],"_wf",as.character(work),".csv")}
  write.csv(x = ggdf_tag4, file = filename)
  
  plot_title <- paste0(names(metadata)[6])
  tag4 <-  names(ggdf_tag4)[2]
  foo <- paste0("ggplot(ggdf_tag4, aes(x = ", tag4,", y = proportion, fill = cluster)) +")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity') + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1, "scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag4 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag4 > 1)){
    nome_file <- paste0("bar_", names(metadata)[6])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag4, device = "png", path = "./tmpdata/")}
  bar_plot_tag4 <- ggplotly(gg_tag4)
  
  plot_title <- paste0(names(metadata)[6])
  tag4 <-  names(ggdf_tag4)[2]
  foo <- paste0("ggplot(ggdf_tag4, aes(x = ", tag4,", y = proportion, fill = cluster)) +")
  foo1 <- paste0(foo, "geom_bar(stat = 'identity', position = 'dodge') + theme_bw() + guides(colour = guide_legend(override.aes = list(size=20))) + ")
  foo2 <- paste0(foo1, "scale_fill_manual(values = color_pheno) + ggtitle(plot_title)")
  parsed <- parse(text = foo2)
  gg_tag4 <- eval(expr = parsed)
  if ((!is.null(save))&&(uni.tag4 > 1)){
    nome_file <- paste0("dodge_", names(metadata)[6])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag4, device = "png", path = "./tmpdata/")}
  dodge_plot_tag4 <- ggplotly(gg_tag4)}

  lista <- list(bar_plot_sample, bar_plot_tag1, bar_plot_tag2, bar_plot_tag3, bar_plot_tag4,
                dodge_plot_sample, dodge_plot_tag1, dodge_plot_tag2, dodge_plot_tag3, dodge_plot_tag4)
  return(lista)
}


plot_median_expr <- function(flow_Set, phenoclust, cell_clustering, metadata, NMC, color_tag1, color_tag2,  color_tag3, 
                             color_tag4, save = NULL, work = NULL) {
 
  mm <- match(cell_clustering, phenoclust$original_cluster)
  cell_label <- phenoclust$new_cluster[mm]
  sample_ids <- rep(metadata$sample_id, fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, sample_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  #For each sample, we plot its PBMC cell type composition represented with colored bars, 
  #where the size of a given stripe reflects the proportion of the corresponding cell type in a given sample
  
  ggdf <- melt(data.frame(cluster = rownames(props), props), 
               id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
  
  ## Add condition info
  mm <- match(ggdf$sample_id, metadata$sample_id)
  ggdf$tag1 <- factor(metadata[,3][mm])
  ## Add patient info
  ggdf$tag2 <- factor(metadata[,4][mm])
  ## Add time step info
  ggdf$tag3 <- factor(metadata[,5][mm])
  ## Add time step info
  ggdf$tag4 <- factor(metadata[,6][mm])
  
  vect_shape1 <- rep(c(1:NMC), times = nlevels(factor(metadata[,3])))
  vect_shape2 <- rep(c(1:NMC), times = nlevels(factor(metadata[,4])))
  vect_shape3 <- rep(c(1:NMC), times = nlevels(factor(metadata[,5])))
  vect_shape4 <- rep(c(1:NMC), times = nlevels(factor(metadata[,6])))
  
  if (length(vect_shape1)<length(flow_Set)){vect_shape1 <- 1:length(flow_Set)}
  if (length(vect_shape2)<length(flow_Set)){vect_shape2 <- 1:length(flow_Set)}
  if (length(vect_shape3)<length(flow_Set)){vect_shape3 <- 1:length(flow_Set)}    
  if (length(vect_shape4)<length(flow_Set)){vect_shape4 <- 1:length(flow_Set)}    
  
  gg_tag1 <- ggplot(ggdf) +
    geom_boxplot(aes(x = tag1, y = proportion, color = tag1, fill = tag1), position = position_dodge(), alpha = 0.5, outlier.color = NA) +
    geom_point(aes(x = tag1, y = proportion, color = tag1, shape = sample_id), alpha = 0.8, position = position_jitterdodge()) +
    facet_wrap(~ cluster, scales = "free", nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 6)) +
    scale_color_manual(values = color_tag1) +
    scale_fill_manual(values = color_tag1) +
    scale_shape_manual(values=vect_shape1) 
  if ((!is.null(save))&&(uni.tag1 > 1)){
    nome_file <- paste0("median_", names(metadata)[3])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag1, device = "png", path = "./tmpdata/")}
  
  gg_tag2 <- ggplot(ggdf) +
    geom_boxplot(aes(x = tag2, y = proportion, color = tag2, fill = tag2), position = position_dodge(), alpha = 0.5, outlier.color = NA) +
    geom_point(aes(x = tag2, y = proportion, color = tag2, shape = sample_id), alpha = 0.8, position = position_jitterdodge()) +
    facet_wrap(~ cluster, scales = "free", nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 6)) +
    scale_color_manual(values = color_tag2) +
    scale_fill_manual(values = color_tag2) +
    scale_shape_manual(values=vect_shape2)
  if ((!is.null(save))&&(uni.tag2 > 1)){
    nome_file <- paste0("median_", names(metadata)[4])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag2, device = "png", path = "./tmpdata/")}
  
  gg_tag3 <- ggplot(ggdf) +
    geom_boxplot(aes(x = tag3, y = proportion, color = tag3, fill = tag3), position = position_dodge(), alpha = 0.5, outlier.color = NA) +
    geom_point(aes(x = tag3, y = proportion, color = tag3, shape = sample_id), alpha = 0.8, position = position_jitterdodge()) +
    facet_wrap(~ cluster, scales = "free", nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 6)) +
    scale_color_manual(values = color_tag3) +
    scale_fill_manual(values = color_tag3) +
    scale_shape_manual(values=vect_shape3)
  if ((!is.null(save))&&(uni.tag3 > 1)){
    nome_file <- paste0("median_", names(metadata)[5])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag3, device = "png", path = "./tmpdata/")}
  
  gg_tag4 <- ggplot(ggdf) +
    geom_boxplot(aes(x = tag4, y = proportion, color = tag4, fill = tag4), position = position_dodge(), alpha = 0.5, outlier.color = NA) +
    geom_point(aes(x = tag4, y = proportion, color = tag4, shape = sample_id), alpha = 0.8, position = position_jitterdodge()) +
    facet_wrap(~ cluster, scales = "free", nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 6)) +
    scale_color_manual(values = color_tag4) +
    scale_fill_manual(values = color_tag4) +
    scale_shape_manual(values=vect_shape4)
  if ((!is.null(save))&&(uni.tag4 > 1)){
    nome_file <- paste0("median_", names(metadata)[6])
    nome_file <- paste0(nome_file, "_wf", as.character(work),".png")
    ggsave(filename = nome_file, plot = gg_tag4, device = "png", path = "./tmpdata/")}
  
  med_expr_tag1 <- ggplotly(gg_tag1)
  med_expr_tag2 <- ggplotly(gg_tag2)
  med_expr_tag3 <- ggplotly(gg_tag3)
  med_expr_tag4 <- ggplotly(gg_tag4)
  
  lista <- list(med_expr_tag1, med_expr_tag2, med_expr_tag3, med_expr_tag4)
  return(lista)
}


plot_stream <- function(flow_Set, flag, phenoclust, cell_label, metadata, NMC, color_pheno) {
  
  tag4_ids <- rep(metadata[,6], fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, tag4_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  
  num_steps <- dim(props)[2]
  nomi_col <- 1000:(1000+(num_steps-1))
  colnames(props) <- as.Date.character(x = nomi_col, format = "%Y")
  colnames(counts) <- as.Date.character(x = nomi_col, format = "%Y")
  
  #props$cluster <- as.numeric(row.names(props))
  props$cluster <- row.names(props)
  #props <- props[order(props$cluster),]
  
  #counts$cluster <- as.numeric(row.names(counts))
  counts$cluster <- row.names(counts)
  
  #counts <- counts[order(counts$cluster),]
  
  #cluster_list <- as.numeric(rownames(props))
  cluster_list <- rownames(props)
  
  names(color_pheno) <- 1:length(cluster_list)
  
  if (!(length(x = cluster_list)==length(color_pheno))){
    if (!any(is.na(as.numeric(names(color_pheno[1:(length(x = cluster_list))]))))){ 
      #to check whether the first color_cluster names are all numerics
      for (i in seq_along(1:length(color_pheno))){
        if(is.na(cluster_list[i])){break}
        if  (!(cluster_list[i] == as.character(i))){color_pheno[i] <- color_pheno[i+1]}}}
  }#this only shift colors in correspondence of jumps in the numbering. The numbers remainains the same
  
  color_pheno <- color_pheno[1:length(cluster_list)] #[nr. of clusters]
  nomi_clust <- names(color_pheno)
  new_color_pheno <- data.frame(color_pheno,nomi_clust, stringsAsFactors = FALSE)
  new_color_pheno <- new_color_pheno[order(new_color_pheno$nomi_clust), ]
  new_color_pheno$nomi_clust <- NULL
  new_color_pheno <- new_color_pheno$color_pheno
  
  props_df <- data.frame(props, stringsAsFactors = FALSE, check.names = FALSE)
  props_ggdf <- melt(data = props_df, id.vars = "cluster")
  props_ggdf <- props_ggdf[order(props_ggdf$cluster, props_ggdf$variable),]
  
  stream_plot_prop <- props_ggdf%>% streamgraph(key = "cluster", value = "value", date = "variable", interpolate="cardinal") %>%
    sg_axis_x(tick_interval = num_steps, tick_format = "|") %>%
    sg_axis_y(tick_count = 10) %>% 
    sg_fill_manual(values = c(new_color_pheno)) %>% #doesn't work
    sg_legend(show=TRUE, label="cluster: ") #%>%
    #sg_title(title = "timestep proportional streamplot")
  
  
  counts_df <- data.frame(counts, stringsAsFactors = FALSE, check.names = FALSE)
  counts_ggdf <- melt(data = counts_df, id.vars = "cluster")
  counts_ggdf <- counts_ggdf[order(counts_ggdf$cluster, counts_ggdf$variable),]
  
  stream_plot_counts <- counts_ggdf%>% streamgraph(key = "cluster", value = "value", date = "variable", interpolate="cardinal") %>%
    sg_axis_x(tick_interval = num_steps, tick_format = "|") %>%
    sg_axis_y(tick_count = 10) %>% 
    sg_fill_manual(values = c(new_color_pheno))%>%
    sg_legend(show=TRUE, label="cluster: ") #%>%
    #sg_title(title = "time_step counts streamplot")
  
  lista <- list(stream_plot_prop, stream_plot_counts)
  return(lista)
}


sample_cluster_df <- function(flow_Set, flag, cell_label, NMC, metadata, work = NULL) {
  
  ################## sample_id
  
  sample_ids <- rep(metadata$sample_id, fsApply(flow_Set, nrow))
  counts_table <- table(cell_label, sample_ids)
  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table)
  #For each sample, we plot its PBMC cell type composition represented with colored bars, 
  #where the size of a given stripe reflects the proportion of the corresponding cell type in a given sample

  if (flag)
  {cluster_list <- factor(rownames(props), levels = 1:NMC)}
  else
  {cluster_list <- rownames(props)}
  
  ggdf_sample <- melt(data.frame(cluster = cluster_list, props), 
                      id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
  
  if (is.null(work)){file_name <- "./tmpdata/df_sample.csv"}
  else {file_name <- paste0("./tmpdata/df_sample_wf", as.character(work),".csv")}
  write.csv(x = ggdf_sample, file = file_name)
  return(ggdf_sample)
}


cluster_t_test <- function(df_sample_wf, NMC, metadata, work = NULL) {
#
# this procedure is to build a t_test table from the df_condition table in which the rows are related to 
# clusters and the columns are the sample and the related characteristic of the sample metadata 
#
  
  metadata <- mt #passing the global dataframe to keep the original values (for some reason the dataframe changes inside the function)
  sample_ids <- metadata$sample_id 
  nr_samples <- length(sample_ids)
  nr_clusters <- length(levels(factor(df_sample_wf$cluster)))
  nr_tags1 <- length(levels(factor(metadata[,3])))
  nr_tags2 <- length(levels(factor(metadata[,4])))
  nr_tags3 <- length(levels(factor(metadata[,5])))
  nr_tags4 <- length(levels(factor(metadata[,6])))
  
  v_tags1 <- levels(factor(metadata[,3]))
  v_tags2 <- levels(factor(metadata[,4]))
  v_tags3 <- levels(factor(metadata[,5]))
  v_tags4 <- levels(factor(metadata[,6]))
  
  p_val <- vector(mode = "numeric",length = NMC)
  eff_size <- vector(mode = "numeric",length = NMC)
  
  colonne <- colnames(metadata)[-c(1,2)]
  somma <- 0
  #the first loop i for producing the necessary columns to add to the tt_in dataframe. In general there will be 
  # n*(n-1)/2 numbers of columns to add, where n is the length of each tag
  
  tt_in <- reshape2::dcast(data = df_sample_wf, formula = cluster ~ sample_id, fun.aggregate = sum, value.var = "proportion")
  if(nr_tags1>1){
    combo <- vector(length = nr_tags1*(nr_tags1-1)/2)
    combinazioni <- t(combn(v_tags1, 2)) 
    for(i in seq_along(1:(nrow(combinazioni)))){ 
      combo[i] <- paste0(combinazioni[i,1],"_vs_",combinazioni[i,2])
    }
    for (j in 1:(nr_tags1*(nr_tags1-1)/2)){
      subs_tag<- subset(metadata, metadata[,3] == combinazioni[j,1])
      first_group <- subs_tag[,"sample_id"]
      subs_tag<- subset(metadata, metadata[,3] == combinazioni[j,2])
      second_group <- subs_tag[,"sample_id"]
      tt_in <- tt_in %>% mutate(p_value=p_val)
      tt_in <- tt_in %>% mutate(eff_size=eff_size)
      colnames(tt_in)[colnames(tt_in) == "p_value"] <- paste0("p_",combo[j])
      colnames(tt_in)[colnames(tt_in) == "eff_size"] <- paste0("ES_",combo[j])
    }
    tt_out <- tt_in
    
    for (i in 1:NMC){
      for (j in 1:(nr_tags1*(nr_tags1-1)/2)){
        subs_tag<- subset(metadata, metadata[,3] == combinazioni[j,1])
        first_group <- subs_tag[,"sample_id"]
        subs_tag<- subset(metadata, metadata[,3] == combinazioni[j,2])
        second_group <- subs_tag[,"sample_id"]
        
        primo <- t(tt_out[i, first_group])
        secondo <- t(tt_out[i, second_group])
        primo <- cbind(primo, rep(combinazioni[j, 1], length(primo)))
        secondo <- cbind(secondo, rep(combinazioni[j, 2], length(secondo)))
        df <- as.data.frame(rbind(primo, secondo))
        colnames(df) <- c("values", "group")
        df$group <- factor(df$group)
        df$values <- as.numeric(df$values)
        
        options(warn = 1) # Turn warnings into errors so they can be trapped
        #old_t_res <- try(t.test(x = tt_out[i,first_group], y =  tt_out[i,second_group], alternative = "two.sided", var.equal = TRUE), silent = FALSE)
        t_res <- try(rstatix::t_test(data = df, formula = values ~ group, alternative = "two.sided", var.equal = TRUE, paired = FALSE), silent = FALSE)
        cohen <- try(rstatix::cohens_d(data = df, formula = values ~ group, var.equal = TRUE, paired = FALSE), silent = FALSE)
        
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,1+nr_samples+j] <- t_res$p
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,1+nr_samples+j] <- NaN
          options(warn = 0)}
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,2+nr_samples+j] <- abs(cohen$effsize)
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,2+nr_samples+j] <- NaN
          options(warn = 0)}
      }
    }
    if (is.null(work)){file_name <- "./tmpdata/df_t_test_tag1.csv"} else 
    {file_name <- paste0("./tmpdata/df_t_test_tag1_wf", as.character(work),".csv")}
    write.csv(x = tt_out, file = file_name)
    print("t_test for tag1 produced")
  }
  
  ####### tag2
  tt_in <- reshape2::dcast(data = df_sample_wf, formula = cluster ~ sample_id, fun.aggregate = sum, value.var = "proportion")
  if(nr_tags2>1){
    combo <- vector(length = nr_tags2*(nr_tags2-1)/2)
    combinazioni <- t(combn(v_tags2, 2)) 
    for(i in seq_along(1:(nrow(combinazioni)))){ 
      combo[i] <- paste0(combinazioni[i,1],"_vs_",combinazioni[i,2])
    }
    for (j in 1:(nr_tags2*(nr_tags2-1)/2)){
      subs_tag<- subset(metadata, metadata[,4] == combinazioni[j,1])
      first_group <- subs_tag[,"sample_id"]
      subs_tag<- subset(metadata, metadata[,4] == combinazioni[j,2])
      second_group <- subs_tag[,"sample_id"]
      tt_in <- tt_in %>% mutate(p_value=p_val)
      tt_in <- tt_in %>% mutate(eff_size=eff_size)
      colnames(tt_in)[colnames(tt_in) == "p_value"] <- paste0("p_",combo[j])
      colnames(tt_in)[colnames(tt_in) == "eff_size"] <- paste0("ES_",combo[j])
    }
    tt_out <- tt_in
    
    for (i in 1:NMC){
      for (j in 1:(nr_tags2*(nr_tags2-1)/2)){
        subs_tag<- subset(metadata, metadata[,4] == combinazioni[j,1])
        first_group <- subs_tag[,"sample_id"]
        subs_tag<- subset(metadata, metadata[,4] == combinazioni[j,2])
        second_group <- subs_tag[,"sample_id"]
        
        primo <- t(tt_out[i, first_group])
        secondo <- t(tt_out[i, second_group])
        primo <- cbind(primo, rep(combinazioni[j, 1], length(primo)))
        secondo <- cbind(secondo, rep(combinazioni[j, 2], length(secondo)))
        df <- as.data.frame(rbind(primo, secondo))
        colnames(df) <- c("values", "group")
        df$group <- factor(df$group)
        df$values <- as.numeric(df$values)
        
        options(warn = 1) # Turn warnings into errors so they can be trapped
        old_t_res <- try(t.test(x = tt_out[i,first_group], y =  tt_out[i,second_group], alternative = "two.sided", var.equal = TRUE), silent = FALSE)
        t_res <- try(rstatix::t_test(data = df, formula = values ~ group, alternative = "two.sided", var.equal = TRUE, paired = FALSE), silent = FALSE)
        cohen <- try(rstatix::cohens_d(data = df, formula = values ~ group, var.equal = TRUE, paired = FALSE), silent = FALSE)
        
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,1+nr_samples+j] <- t_res$p
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,1+nr_samples+j] <- NaN
          options(warn = 0)}
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,2+nr_samples+j] <- abs(cohen$effsize)
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,2+nr_samples+j] <- NaN
          options(warn = 0)}
      }
    }
    if (is.null(work)){file_name <- "./tmpdata/df_t_test_tag2.csv"} else 
    {file_name <- paste0("./tmpdata/df_t_test_tag2_wf", as.character(work),".csv")}
    write.csv(x = tt_out, file = file_name)
    print("t_test for tag2 produced")
  }
  
  ####### tag3
  tt_in <- reshape2::dcast(data = df_sample_wf, formula = cluster ~ sample_id, fun.aggregate = sum, value.var = "proportion")
  if(nr_tags3>1){
    combo <- vector(length = nr_tags3*(nr_tags3-1)/2)
    combinazioni <- t(combn(v_tags3, 2)) 
    for(i in seq_along(1:(nrow(combinazioni)))){ 
      combo[i] <- paste0(combinazioni[i,1],"_vs_",combinazioni[i,2])
    }
    for (j in 1:(nr_tags3*(nr_tags3-1)/2)){
      subs_tag<- subset(metadata, metadata[,5] == combinazioni[j,1])
      first_group <- subs_tag[,"sample_id"]
      subs_tag<- subset(metadata, metadata[,5] == combinazioni[j,2])
      second_group <- subs_tag[,"sample_id"]
      tt_in <- tt_in %>% mutate(p_value=p_val)
      tt_in <- tt_in %>% mutate(eff_size=eff_size)
      colnames(tt_in)[colnames(tt_in) == "p_value"] <- paste0("p_",combo[j])
      colnames(tt_in)[colnames(tt_in) == "eff_size"] <- paste0("ES_",combo[j])
    }
    tt_out <- tt_in
    
    for (i in 1:NMC){
      for (j in 1:(nr_tags3*(nr_tags3-1)/2)){
        subs_tag<- subset(metadata, metadata[,5] == combinazioni[j,1])
        first_group <- subs_tag[,"sample_id"]
        subs_tag<- subset(metadata, metadata[,5] == combinazioni[j,2])
        second_group <- subs_tag[,"sample_id"]
        
        primo <- t(tt_out[i, first_group])
        secondo <- t(tt_out[i, second_group])
        primo <- cbind(primo, rep(combinazioni[j, 1], length(primo)))
        secondo <- cbind(secondo, rep(combinazioni[j, 2], length(secondo)))
        df <- as.data.frame(rbind(primo, secondo))
        colnames(df) <- c("values", "group")
        df$group <- factor(df$group)
        df$values <- as.numeric(df$values)
        
        options(warn = 1) # Turn warnings into errors so they can be trapped
        old_t_res <- try(t.test(x = tt_out[i,first_group], y =  tt_out[i,second_group], alternative = "two.sided", var.equal = TRUE), silent = FALSE)
        t_res <- try(rstatix::t_test(data = df, formula = values ~ group, alternative = "two.sided", var.equal = TRUE, paired = FALSE), silent = FALSE)
        cohen <- try(rstatix::cohens_d(data = df, formula = values ~ group, var.equal = TRUE, paired = FALSE), silent = FALSE)
        
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,1+nr_samples+j] <- t_res$p
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,1+nr_samples+j] <- NaN
          options(warn = 0)}
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,2+nr_samples+j] <- abs(cohen$effsize)
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,2+nr_samples+j] <- NaN
          options(warn = 0)}
      }
    }
    if (is.null(work)){file_name <- "./tmpdata/df_t_test_tag3.csv"} else 
    {file_name <- paste0("./tmpdata/df_t_test_tag3_wf", as.character(work),".csv")}
    write.csv(x = tt_out, file = file_name)
    print("t_test for tag3 produced")
  }
  
  ####### tag4
  tt_in <- reshape2::dcast(data = df_sample_wf, formula = cluster ~ sample_id, fun.aggregate = sum, value.var = "proportion")
  if(nr_tags4>1){
    combo <- vector(length = nr_tags4*(nr_tags4-1)/2)
    combinazioni <- t(combn(v_tags4, 2)) 
    for(i in seq_along(1:(nrow(combinazioni)))){ 
      combo[i] <- paste0(combinazioni[i,1],"_vs_",combinazioni[i,2])
    }
    for (j in 1:(nr_tags2*(nr_tags2-1)/2)){
      subs_tag<- subset(metadata, metadata[,6] == combinazioni[j,1])
      first_group <- subs_tag[,"sample_id"]
      subs_tag<- subset(metadata, metadata[,6] == combinazioni[j,2])
      second_group <- subs_tag[,"sample_id"]
      tt_in <- tt_in %>% mutate(p_value=p_val)
      tt_in <- tt_in %>% mutate(eff_size=eff_size)
      colnames(tt_in)[colnames(tt_in) == "p_value"] <- paste0("p_",combo[j])
      colnames(tt_in)[colnames(tt_in) == "eff_size"] <- paste0("ES_",combo[j])
    }
    tt_out <- tt_in
    
    for (i in 1:NMC){
      for (j in 1:(nr_tags4*(nr_tags4-1)/2)){
        subs_tag<- subset(metadata, metadata[,6] == combinazioni[j,1])
        first_group <- subs_tag[,"sample_id"]
        subs_tag<- subset(metadata, metadata[,6] == combinazioni[j,2])
        second_group <- subs_tag[,"sample_id"]
        
        primo <- t(tt_out[i, first_group])
        secondo <- t(tt_out[i, second_group])
        primo <- cbind(primo, rep(combinazioni[j, 1], length(primo)))
        secondo <- cbind(secondo, rep(combinazioni[j, 2], length(secondo)))
        df <- as.data.frame(rbind(primo, secondo))
        colnames(df) <- c("values", "group")
        df$group <- factor(df$group)
        df$values <- as.numeric(df$values)
        
        options(warn = 1) # Turn warnings into errors so they can be trapped
        old_t_res <- try(t.test(x = tt_out[i,first_group], y =  tt_out[i,second_group], alternative = "two.sided", var.equal = TRUE), silent = FALSE)
        t_res <- try(rstatix::t_test(data = df, formula = values ~ group, alternative = "two.sided", var.equal = TRUE, paired = FALSE), silent = FALSE)
        cohen <- try(rstatix::cohens_d(data = df, formula = values ~ group, var.equal = TRUE, paired = FALSE), silent = FALSE)
        
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,1+nr_samples+j] <- t_res$p
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,1+nr_samples+j] <- NaN
          options(warn = 0)}
        if (inherits(x = t_res,"rstatix_test")){ 
          tt_out[i,2+nr_samples+j] <- abs(cohen$effsize)
          options(warn = 0)}
        else{ # Ignore warnings while processing errors
          options(warn = -1)
          tt_out[i,2+nr_samples+j] <- NaN
          options(warn = 0)}
      }
    }
    if (is.null(work)){file_name <- "./tmpdata/df_t_test_tag2.csv"} else 
    {file_name <- paste0("./tmpdata/df_t_test_tag4_wf", as.character(work),".csv")}
    write.csv(x = tt_out, file = file_name)
    print("t_test for tag4 produced")
  }
}

