######################## Libraries -----------------------------
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
#library(shinyFiles)
#library(V8) # I loaded this to try to play a sound at the end of each long procedure. invain

library(parallel)
library(DT)
library(tidyverse)
library(lubridate)
library(plotly)
library(RColorBrewer)
library(reshape2) 
library(rhandsontable)
library(ggrepel)
#library(heatmaply)
library(pheatmap)
library(cluster) #This is for the cluster performance evaluation
library(png) # these are for plotting on a grid the png files
library(grid)
library(gridExtra)
library(streamgraph) #devtools::install_github("hrbrmstr/streamgraph")
library(ggpubr)
library(fs)
library(audio)
library(rstatix)

#library(openCyto)
#library(flowViz)
library(flowCore)
library(flowWorkspace) #This must be introduced to perform the alignment process
#library(ggcyto)
library(flowAI)
library(flowTrans)
#library(scales) #maybe to be added for the flowStats transformation 
#see https://bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/flowWorkspace-Introduction.html
library(flowStats)
library(DepecheR) #This must be added to the list
library(spade)  #https://github.com/nolanlab/spade the installation is based on github and devtools (install_github("nolanlab/spade@HEAD"))
#library(limma)
library(MASS) #it replaces limma for a more general MDSplot (with factoextra)
library(matrixStats)
library(factoextra)
library(mclust)
library(kohonen) # its presence to be checked since we are trying:
#library(meanShiftR) #to clusterize the tsne map. See https://towardsdatascience.com/the-5-clustering-algorithms-data-scientists-need-to-know-a36d136ef68 
#library(dbscan) ##to clusterize the tsne map. See http://www.sthda.com/english/wiki/wiki.php?id_contents=7940 

library(FlowSOM)
library(DDoutlier)
library(Rphenograph) #install_github("JinmiaoChenLab/Rphenograph@HEAD")
library(FastPG) #BiocManager::install("sararselitsky/FastPG")
library(ConsensusClusterPlus)
library(Rtsne)
library(umap)
library(flowVS)

######################## Header -----------------------------
dbHeader <- dashboardHeader(title = "cytoChain")


######################## Sidebar -----------------------------
dbSidebar <- dashboardSidebar(
    
    sidebarMenu(
        
        menuItem("Intro & Purposes", tabName = "intro", icon = icon("home", lib = "glyphicon")),
        menuItem("FlowSet optimization Workflow", tabName = "preClustering", icon = icon("vials")),
        menuItem("Metadata & Assays Workflow", tabName = "metadata", icon = icon("file-alt")),
        menuItem("High dim analysis Workflow", tabName = "clusteringworkflow", icon = icon("microscope"), startExpanded = TRUE, 
                 menuSubItem("classical workflow", tabName = "classical"),
                 menuSubItem("alternative workflow", tabName = "alternative")),
        menuItem("Panel Editor", tabName = "Panel_editor", icon = icon("pencil-alt"), badgeLabel = "new", badgeColor = "green"),
        menuItem("Credits & References", tabName = "credits", icon = icon("r-project"), badgeLabel = "new", badgeColor = "green")))


######################## Purpose TabPanel -----------------------------
tabPurpose <- tabPanel("Loading & parsing samples", icon = icon("home", lib = "glyphicon"),
                       column(width=12,
                       valueBoxOutput(outputId = "entryEvents", width = 3)),
                       column(width=12,
                       h4("This is a", strong("workflow"), "to prepare and to analyse your high dimensional flow cytometry data. It is organized 
                         as a", strong("storyboard"), "because the main purpose is to setup your", shiny::em("flowSet"), "through a step by step 
                         process, but you could jump directly to a specific", strong("frame"), "(each individual step corresponds to a tab in this 
                         dashboard) having the proper input for the selected step"),
                       # em is masked by nclust package, that's why I have to specify shiny::em
                       h4("The flowSet optimization workflow steps are:", 
                          tags$ol(tags$li(strong("Loading "), ", parsing and vizualizing your data: to display the generic properties of the loaded 
                                         flowSet"),
                                  tags$li(strong("Cleaning "), "your samples: to proceed in your analysis with a flowSet free of \"dirty\" events"), 
                                  tags$li(strong("Transforming "), "(scaling): to scale your dimensions"),
                                  tags$li(strong("Aligning"), ": applying an algorythm to align your dimensions along you samples"), 
                                  tags$li(strong("Downsampling"), ": an effective methods based on each event density"), 
                                  #tags$li(Subsetting: to produce a flowSet/flowFrame with a reduced number of events based on the selected 
                                  #        subpopulation"),
                                  tags$li("Concatenating: to aggregate in a single flowFrame your flowSet collecting the whole bunch of your 
                                         events"))),
                       #h4(strong("Attention! Follow these tips when submitting your samples:")),
                       h4("All samples must have a common panel (an identical marker set with both name and description for each stain). If some of 
                       your samples contains a different marker set, please consider to further split your analysis accordingly or try to fix them 
                          using the ", strong("Panel Editor")),
                       h4("All samples must be ", strong("already compensated. "), "The whole analysis is based on the expression matrix values: 
                       this means that all the channel's values have already applied the proper compensation process."),
                       h4("Before supplying inputs to cytoChain be sure to accomplish the ", strong("Forward scatter gating"), "for selecting the 
                       cell size and filtering out debris (typically by the FSC-A & FSC-H channel) and ", strong("Side scattering gating"), "related 
                       to the cell’s granularity to remove the doublets (typically by using SSC-A & FSC-A physical channel)"),
                       h4("This gating process is strongly related to the reasercher experience and it can hardly automated for any kind of events")),

                       
                       tags$hr(style="border-color: green;"),
                       h3("Select the proper process modules to define your flowSet optimization procedure"),
                       checkboxGroupButtons(inputId = "process", label = "Make your choice:", 
                                            #choiceValues = list("clean", "scale", "align", "downsample", "subset", "concatenate"),
                                            choiceValues = list("clean", "scale", "align", "downsample", "concatenate"),
                                            choiceNames = list(icon("shower"), 
                                                               icon("ruler-combined"), 
                                                               icon("align-center"), 
                                                               icon("vial"), 
                                                               #icon("cut"), 
                                                               icon("link")),
                                            #selected = list(F,F,F,F,F,F),
                                            selected = list(F,F,F,F,F),
                                            checkIcon = list(T = icon("ok", lib = "glyphicon"), F = icon("remove", lib = "glyphicon"), 
                                                             justified=TRUE),
                                            #checkIcon = list(yes = tags$i(style = "color: steelblue"), no = tags$i(style = "color: steelblue")),
                                            status = "success", size = "lg", individual = TRUE, justified = TRUE),
                       # theme styling ##???
                       tags$head(tags$style(HTML('#checkbox :after, #checkbox :before{background-color:#bff442;}'))),
                       
                       h3(textOutput("checktxt")),
                       h3("your selected items are the following:"),
                       textOutput("txt"),
                       tags$hr(style="border-color: green;"),
                       fluidRow(
                           column(3,
                                  h3("Load your samples as collection of FCS files"),
                                  fileInput(inputId = "fSinput", label = h4("File input"), multiple = T,
                                    buttonLabel = "FCS files", placeholder = "upload your flowFrames"),
                                  p("The samples must be valid FCS files of version 3.0+")), 
                           column(1, offset = 1,
                                  h3("...or")),
                           column(4,
                                  h3("Load your samples collected in a single zip file"),
                                  fileInput(inputId = "zipinput", label = h4("File input"), multiple = FALSE, 
                                            accept = c("application/zip", "ZIP archive", ".zip"),
                                     buttonLabel = "zip file", placeholder = "upload your single zip file"),
                                  p("The samples inside the zip file must be valid FCS files of version 3.0+. Only the zip format is supported"))),
                       
                       tags$hr(style="border-color: green;"),
                       h3("Perform loading & parsing process"),
                       actionButton("runLoad", "Load and parse your samples", icon("paper-plane"), 
                                      style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                       h3(textOutput("checktxt_parsing_fS")), 
                       conditionalPanel(
                         condition = ("output.load_panelStatus"),  
                         #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                         h4("Push to get some summaries of the loaded samples. Please notice that, all along this web-tool, to produce this kind of 
                            visualization you have to press this gray buttons after the main operation (usually performed after pushing red 
                            buttons): the tables and the plots won’t update automatically"), 
                         actionButton('runLoadShow', "Showing samples' details"),
                         h3(textOutput("checktxt_loading_fS")), br(),
                         conditionalPanel(
                           condition = ("input.runLoadShow !== 0"),
                           h4("This is a summary of the number of samples entered and the events involved"),
                           DTOutput(outputId = "fSsummary")%>% withSpinner(), br(),
                           h4("The following is a table with all the samples each one with its amount of event - Notice than the sample name is 
                              a temporary id assigned by this web app environement"),
                           DTOutput(outputId = "samples")%>% withSpinner(),
                           downloadButton('downloadfS_table_csv', "Download table in csv format"),br(),
                           h4("...In a graphic format:"),
                           plotOutput(outputId = "sample_plot")%>% withSpinner(), br(),
                           h4("This table reports the original sample name (the .FCS file) and its size in byte"),
                           DTOutput(outputId = "contents")%>% withSpinner(),br(),
                           h4("The following is the the list of the markers used with the correspending dye name and the description. In particular
                           this is the list of the markers having a description in the loaded FCS files. From this point on, to the markers without 
                           a description will be assigned the corresponding name"),
                           h4("cytoChain takes for granted that you follow some conventions. First of all, each of the dimension name and marker 
                           description should be distinct and unique inside each flowFrame. There are also some", 
                              strong("restrictions to be followed"), "otherwise some of the features of cytoChain may not work properly:"),
                           h4("   1) The physical dimensions' names should contain: ",strong("'FSC' or 'SSC'"), ". That's the offical cytometry 
                           standard names reserved for the forward and side scattering laser beams"),
                           h4("   2) The markers' names must not contain: ",strong("'FSC', 'SSC', 'TIME', 'score', 'density' or 'cell_Id'"), 
                           "whatever case used. That's because these names are reserved for the physical lasers, the 'time' dimensions, 
                           while the 'cell_Id' is an additional service channel introduced by cytoChain to tag the single event, and the 
                           'score' and ' 'density' are names used by the downsampling process to measure the 'outlierness' of an event"),
                           DTOutput(outputId = "markers")%>% withSpinner(),
                           downloadButton('downloadfS_dim_csv', "Download table in csv format"),br(),
                           h4("This plot shows the range of the markers. A correctly scaled marker range should spread, more or less, along 
                              a sigle order of magnitude"),
                           plotlyOutput(outputId = "marker_range")%>% withSpinner(),br(),
                           h4("This pie is just a visualization of the sample colors which will be used along all the analysis process"),
                           plotOutput(outputId = "samples_color"),
                           p("Since the sample colors are chosen randomly on the basis of the number of samples, sometime it could be that
                           the color set is unsatisfactory. In this case go back up again in this page, change the seed and push the 
                           'showing samples details' button again"),
                           
                           tags$hr(style="border-color: green;"),br(),
                           
                           tags$footer(p("A", strong("workflow"), "or a", shiny::em("pipeline"), "is a ordered set of procedures which can be 
                           perfomed in order to achieve some results. A", strong("stroyboard"), "and its frame set is the way a workflow can be 
                           organized within this webtool"), 
                           p("A ", strong("flowSet"), " is the conventional name used in R environment in order to handle and to parse the", 
                             shiny::em("fcs"), "files which are the basic file format of your flow cytometry output. A ",strong("flowFrame")," is the 
                             single sample, namely the single experiment (the fcs file is, in digital format, the content of the 'tube' full of your 
                             cells loaded to the cytometer); the flowSet is a collection of flowFrames. Each flowFrame contains many relevant data 
                             about the experiment, but the main content is the ", strong("expression matrix"), ": an ", strong("n x m"), " real 
                             number matrix values in which each row correspond to a single event measured by the cytometer. The measure of each 
                             event (typically a single cell) is the coordinate of a point in the m-dimensional space of the marker expression: 
                             every point has m level of expression, one for each marker (the stain)."),
                           p("For every detail about the fcs format see ", 
                             tags$a(href="https://bioconductor.org/packages/release/bioc/vignettes/flowCore/inst/doc/fcs3.html", "fcs3 format."),
                                            "All about flowFrame parsing and content handling in this tool is based on the ", shiny::em("flowCore"), 
                             "library see", tags$a(href = "https://www.nature.com/articles/srep20686", "openCyto"))),
                           tags$hr(style="border-color: green;"),
                           h4("Additional notes on the loaded samples:"),
                           p("As you can see on the table content shown above, the dimensions (aka the markers, the stains, ...) shown are only 
                           the ones which, in the original flowFrame, holds both names and descriptions. Please notice than, even for the rest of 
                           the operation in the following workflows, only those markers which are present on both side of the flowSet table 
                           description, could be taken in consideration.")
                           )#conditional panel input.runLoadShow
                         )#conditional panel input.runLoad
                       ) # tabPanel "Loading & parsing samples"

######################## Clean TabPanel -----------------------------
tabClean <- tabPanel(title = "Clean", icon =icon("shower"),
                     column(width=12,
                            valueBoxOutput(outputId = "c.entryEvents"),
                            valueBoxOutput(outputId = "c.handledEvents"),
                            valueBoxOutput(outputId = "c.percEvents")),
                     column(width=12,
                     p("In this section, it is possible to perform a 'cleaning' operation. The aim is to remove those events in each flowFrame, 
                       which could be the outcome of a bad data handling in the acquisition session of the cytofluorimeter."),
                     p("This part of the flowSet optimization workflow, is very much important to filter out those data which may results in kind of 
                     disturbing noise in the following analysis. It is based on the ", strong("flowAI"), " package of Gianni Monaco, Hao Chen (see 
                     the 'credit tab'),  a package for automatic and interactive quality control for flow cytometry data. In this implementation 
                     only the automatic quality control on the FCS data acquired is performed. This is done by evaluating three different 
                       properties (performed in a sequence of steps):"),
                     p("1) ", strong("flow rate"), " :This first step is based on the sample (the .fcs file) keyword $TIMESTEP and the 'time' 
                     dimension. The purpose is to clean the cell outliers of the acquisition step. For example it happens often that, especially in 
                     the starting phase and the very last phase, when you start or you stop the cell acquisition from the tube, the fluid is not in 
                     the laminar flow conditions: this always cause bad measurements affecting the whole stain panel"), 
                     p("2) ", strong("signal acquisition"), " :it is the control which involves each single signal (aka stain, dimension,...) of your 
                     panel. Of course the signal changes in time and this is the aim in order to build the expression matrix for the sample, but 
                       the mean and standard deviation of the median should remain constant over the course of the analysis: the signal acquisition 
                       is anomalous when there are changes which are too marked. The analysis try to leave out these 'out of range' events"), 
                     p("3) ", strong("dynamic range"), " :This is the last step in which the events from the upper and lower limits of the dynamic 
                       range are checked. The limit are related to the maximum value of the dynamic range related to the instrument itself, limit 
                       which are pre-set by the manufacturer"), 
                     h4("The statistical analysis of the flowSet is based on the time channel values: Please, verify that the 'time' channel is 
                        present inside the expression matrix of your flowFrames."),
                     h4("The second step of the quality check (signal acquisition check) is based on expression matrix values. Remember: the 
                        channels should be compensated!"),
                     p("Please go through the relevant parameters listed on the left side for the cleaning operation. You can control both the 
                       flow rate check (with the 'alphaFR' parameter) and the signal acquisition check (with the 'pen_valueFS' parameter). The 
                       dynamic range check is set to the default setting. ")),
                     h4("Perform cleaning process"),
                     p("Be patient! Depending on the sample's size (e.g. with FCS files bigger than 2 or 3 Mbyte), this process could take a long 
                       time!"),
                     actionButton("runClean", "Run cleaning algorithm", icon("paper-plane"), 
                                  style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                     h4(textOutput("checktxt_cleaning_fS")),
                     conditionalPanel(
                       #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                       condition = ("output.clean_panelStatus"),   
                       h3("Push to get some summaries of the samples loaded"), 
                       actionButton('runCleanShow', "Showing samples' details"), br(),
                       conditionalPanel(
                         condition = ("input.runCleanShow !== 0"),
                         DTOutput(outputId = "cleanedSummary")%>% withSpinner(), br(),
                         h4("The following is the table reporting the number of events in each flowFrames after the cleaning operation"),
                         DTOutput(outputId = "cleaned")%>% withSpinner(),
                         downloadButton('downloadCleanedfS_table_csv', "Download table in csv format"),br(),br(),
                         plotOutput(outputId = "clean_comparison_plot", width = "100%", height = "1200px"), 
                         h4("This table reports a comparison before and after the cleaning..."),
                         DTOutput(outputId = "clean_comparison")%>% withSpinner(),
                         downloadButton('downloadComparison_table_csv', "Download table in csv format"),br(),br(),
                         h4("...and this plots shows the differences between the two flowSet"),
                         p('Notice the red points for each marker expression: the farther are apart, the stronger is the chance that we are 
                           filtering outliers'),
                         plotlyOutput(outputId = "plot_clean_comparison")%>% withSpinner(),br(),
                         h4("The figure below shows the expression density of the various markers"),
                         #plotlyOutput(outputId = "marker_density_clean",  width = "100%", height = "1200px")%>% withSpinner(), br(),
                         plotlyOutput(outputId = "marker_density_clean")%>% withSpinner(), br(),
                         h4("The table below shows a report with a detailed statistic of the cleaning operation"),
                         DTOutput(outputId = "QC_clean")%>% withSpinner(),
                         downloadButton('downloadQC_table_csv', "Download table in csv format"),br(),br(),
                         downloadButton(outputId = "downloadCleanfS",label = "download your cleaned flowFrames")
                         ) #conditionalPanel input.runCleanShow !== 0
                       ) #conditionalPanel input.runClean !== 0
                     )

######################## Scale TabPanel -----------------------------
tabScale <- tabPanel(title = "Scale", icon =icon("ruler-combined"),
                     p("This part of the flowSet optimization workflow is related to the scaling process. It is well known that the range of the 
                     measured expression variable for each single marker could vary a lot and it is a common action to transform these range in a 
                     scaling process. The most used scale transformation for cytometry is the inverse hyperbolic sine ('arcsinh'), but other type of 
                     transformation could also be implemented. For example the simple normalization scaling based on the (x - min(x)) / 
                     (max(x) - min(x)) function, or the 'arcsinh+' in which the x - min(x) shifting is performed in order to deal with the classical 
                     'arcsinh' transfromation but returning only positive values (useful for certain clustering algorithms which does not allow 
                       negative values)"),
                     p("An interesting transformation is the 'zero' scaling in which, together with a normalization (setting the range between 
                       0 and 1) the values below a certain thrsholds could be set to 0 (and above a certain other threshold could be set to 1). This 
                       type of transformation is used along the high dimensional analysis workflow in this tool, in order to emphatyze the color 
                       differences in a the various heatmap plots, but it is usually not reccomanded in the normal analysis"),
                     p("A more enlightened method is to generalize the arcsinh transformation (or other well know transformations) with some 
                     parameters, from the ", strong("flowTrans")," package (see ", 
                       tags$a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-546"), ")"),
                     p("Another approach is the one allowed by the ",strong('flowVS')," package. flowVS (see ", 
                       tags$a(href="https://www.bioconductor.org/packages/release/bioc/html/flowVS.html"," )"),  
                       "is able to detect changes in populations across biological conditions, stabilizing the variance separately on each 
                       fluorescence channel. The same channel in all samples will be transformed with the same parameter (cofactor). In this way the 
                       arcsinh transformation is still the function used but with a tuned parameter, one for each channel (and not a common one for 
                       the whole bunch of channels - in this case the arcsinh parameter on the sidebar is not used). For this reason it is very much
                       inmportant to cross check the channel's densities along the samples in the plot produced at the end of the transformation.
                       The 'flowVS' scaling may take a quite long time depending on the number of channels"), br(),
                     h3("Perform transforming process"),
                     actionButton("runScale", "Run transforming algorithm", icon("paper-plane"), 
                                  style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                     h3(textOutput("checktxt_transforming_fS")), br(),
                     conditionalPanel(
                       condition = ("output.scale_panelStatus"), 
                       #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                       h3("Push to get some summaries of the transformed samples"), 
                       actionButton('runScaleShow', "Showing sample's details"), br(),
                       conditionalPanel(
                         condition = ("input.runScaleShow !== 0"),
                         h4("This plot shows the range of the markers. Check the spread of each marker's range"),
                         plotlyOutput(outputId = "marker_range_trans")%>% withSpinner(), br(),
                         h4("The figure below shows the expression density of the various markers"),
                         plotlyOutput(outputId = "marker_density_trans", width = "100%", height = "2400px")%>% withSpinner(),
                         tags$hr(style="border-color: green;"),
                         h4("Additional notes on the transformed samples:"),
                         p("With the transformed values it is actually possible to double check and to compare the marker trends, namely each
                           single dimensions. It is important to stress that, in case of a marker expression trend which is basically flat along its 
                           whole expression range, it is better to neglect the dimension itself for the rest of the analysis. In fact,  a virtually 
                           constant behaviour does not allow any inference with the rest of the markers: it is like a constant adding no information 
                           to correlate to the rest of the variables, the other cell related expressions. So then, the advice here is not to 
                           consider this 'flat' marker for the rest of the analysis and most of all in the metadata workflow"), br(),
                         actionButton('saveScale', "Producing your transformed and complete samples"), br(),
                         h4("Remember to produce your brand new samples each time you perform a data handling!"),br(),br(),
                         conditionalPanel(
                             condition = ("output.saveScale_panelStatus"),
                             downloadButton(outputId = "downloadTransfS",label = "download your transformed flowFrames"),br(),br(),
                             downloadButton(outputId = "downloadTransfS_complete", label = "download your .fcs samples with both the handled 
                             expression matrix and the original one")))))

######################## Align TabPanel -----------------------------
tabAlign <- tabPanel(title = "Align", icon =icon("align-center"),
                     h4("This tab of the flowSet optimization workflow is dedicated to the event alignment process. Especially for those experiments 
                     which involves a consistent number of samples, it could happen that the instrument itself may change, even for a slight 
                     fluctuation, the overall condition during the sample acquisition: this could spread along the whole set of samples to become a 
                     little drift on the marker expression values"),
                     h4("This alignment process tends to adjust along all the sample set, the different population of events, in order to align them 
                     following the density distributions of the selected markers to fix this drifting phenomenon. This is performed through the 
                     warpSet function of the ", strong("flowStat"), " package. The process is based on four steps for each single selected marker:", 
                       tags$ol(tags$li("Searching for the high density peaks (landmarks): each of these landmark should be related to a certain 
                                       population of events"),
                               tags$li("Number estimation of those landmarks: as it can be seen from the densities plots, some markers could have 
                                       multiple peaks"), 
                               tags$li("K-means clustering: to select and classify the relevant population areund the landmarks"),
                               tags$li("Landmarks alignment estimation: to compute the best fitting coordinates to align the population over the 
                                       different samples"))),
                     h4(" ",strong("Warning:")," The whole process works fine if the high density areas truly represent particular sub-types of 
                     cells. In order to align along the different samples , a certain population should be uniquely determined, which is true when 
                     the markers are ", strong("binary"), " that is the cells are either positive or negative for a particular marker (namely, 
                     the typical case of lineage markers). Avoid to align with markers which smear all along their expression!"),br(),
                     p("The event alignemnt is performed through the ", strong("flowStat"), " package developed by F.Hahne, N. Gopalakrishnan, 
                       A.H. Khodabakhshi, C. Wong and K. Lee"),
                     tags$hr(style="border-color: green;"),
                     h4("The following table reports the marker list. Please select the rows in the table in order to choose on which 
                        markers the alignment should be performed"),
                     
                     DT::dataTableOutput(outputId = 'marker2align')%>% withSpinner(),
                     verbatimTextOutput('mark2align.res'), br(),
                     h3("Perform marker's alignment"),
                     actionButton("runAlign", "Run the aligning algorithm", icon("paper-plane"), 
                                  style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                     h3(textOutput("checktxt_aligning_fS")),
                     conditionalPanel(
                       condition = ("output.align_panelStatus"), 
                       h3("Push to get some plots showing the density shapes of the samples"), 
                       actionButton('runAlignShow', "Show aligned flowSet details'"),br(),
                       #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                       conditionalPanel(
                         condition = ("input.runAlignShow !== 0"),
                         h4("this plots shows the selected marker density before the alignment process"),
                         plotOutput(outputId = "plot_align_b4", width = "100%", height = "2000px")%>% withSpinner(),br(),
                         h4("... while this shows the expression density after the alignment"),
                         plotOutput(outputId = "plot_align_after", width = "100%", height = "2000px")%>% withSpinner(),
                         h4("... and the classical densities view"),
                         plotlyOutput(outputId = "plot_align_after2", width = "100%", height = "2000px")%>% withSpinner(),
                         tags$hr(style="border-color: green;"),
                         h4("Additional notes on the aligned samples:"),
                         p("A safe and meaningful choice of the markers to perform the alignment, should be carefully selected. In particular, try 
                         to align those dimensions in which the expression of the scaled value trend, shows little differences between the sample's 
                         peaks. If the differences are unduly large, especially in those markers which have more than one peak, the alignment could 
                         be wrongly performed, trying to align peaks which are unrelated. So then, compare the markers trends before and after the 
                         alignment to verify the results and to avoid possible misalignments"),
                         p("Finally, please also notice the sensitive nature of this type of process: together with the scaling, this is the data 
                         transformation which fundamentally alters the data input and indeed the experiment values. But while scaling is mandatory 
                         to produce something visible and comparable in scale, performing an easy and reversible function transformation, this 
                         alignment process perform a complex data manipulation, so please the advice is 'handle the alignment  with care' motivating 
                         your choice. The best would be to perform the complete analysis with and without the alignment to cross-check whether 
                         there is a significant improvement in the maps production with clearer clusters and phenotypes"),br(),
                         actionButton('saveAlign', "Producing your transformed and complete samples"), br(),
                         h4("Remember to produce your brand new samples each time you perform a data handling!"),br(),br(),
                         conditionalPanel(
                             condition = ("output.saveAlign_panelStatus"),
                             downloadButton(outputId = "downloadAlignedfS",label = "download your aligned flowSet"), br(), br(),
                             downloadButton(outputId = "downloadAlignedfS_complete",label = "download your .fcs samples with both the handled 
                             expression matrix and the original one"))))) 


######################## Downsample TabPanel -----------------------------
tabDowns <- tabPanel(title = "Downsample", icon =icon("vial"),
                              column(width=12,
                                     valueBoxOutput(outputId = "d.entryEvents"),
                                     valueBoxOutput(outputId = "d.handledEvents"),
                                     valueBoxOutput(outputId = "d.percEvents")),
                              column(width=12,
                     p("A relevant part of the flowSet optimization data preparation is related to the down-sampling process. 
                     The main scope here is to remove the outliers, namely in this context, those events (observed cells) in the multi-dimensional 
                     space of the expression matrix, that lies in an abnormal distance from other points of the events population. Since it is a 
                     statistical concept, analyst should define what will be considered abnormal: before abnormal observations can be singled out, 
                     it is necessary to characterize normal observations. With this, having already cleaned the event space within the cleaning 
                     process, which are the outliers then?"),
                     p("This point is essential. Most of the time, in fact, the down-sampling is performed in order to limit a redundant number of 
                     events, to a more concise set, in order to deal with a reduced and less heavy amount of data without loosing the relevant 
                     information. This is not the case: the aim here is remove those events which lie alone or isolated from the rest of the points. 
                     The outcome is a dataset with only the 'dense' part of the event space. This should help the clustering algorithm to better 
                       define the clusters especially for those critical 'border points' which could limit the convergence of the clustering 
                       process."),  
                     h3("Perform downsampling process"),br(),
                     actionButton("runDown", "Run the downsampling algorithm", icon("paper-plane"), 
                                  style="color: #fff; background-color: #ff0000; border-color: #2e6da4"), br(),
                     h3(textOutput("checktxt_downsampling_fS")),br(),
                     conditionalPanel(
                         condition = ("output.down_panelStatus"),
                         #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                         h3("Push to get some summaries of the downloaded samples"), 
                         actionButton('runDownShow', "Show downsampled data details'"),br(),
                         conditionalPanel(
                             condition = ("input.runDownShow !== 0"),
                             DTOutput(outputId = "downsampledSummary")%>% withSpinner(),br(),
                             h4("The following is the table reporting the number of events in each flowFrames after the downsampling"),
                             DTOutput(outputId = "downsampled")%>% withSpinner(),
                             downloadButton('downloadDownsampled', "Download table in csv format"),br(),br(),
                             h4("The following barplot and the table reports a comparison before and after the downsampling..."),
                             plotOutput(outputId = "down_comparison_plot", width = "100%", height = "1200px"), br(),
                             DTOutput(outputId = "down_comparison")%>% withSpinner(), br(),
                             downloadButton('downloadComparisonTable', "Download table in csv format"),br(),
                             h4("...and this plots shows the marker's ranges differences between the two flowSet"),
                             p('Notice the red points for each marker expression: the farther are apart, the stronger is the chance that we are 
                               filtering outliers. The comparison is with the entry flowSet (before any handling)'),
                             plotlyOutput(outputId = "plot_down_comparison")%>% withSpinner(), br(),
                             h4("The figure below shows the expression density of the various markers"),
                             plotlyOutput(outputId = "marker_density_down", width = "100%", height = "1200px")%>% withSpinner(), br()),
                         conditionalPanel(
                             condition = ("output.down_algoStatus"),
                             h4('The following two couples of plots are related to the downsampling algorithms evaluation. The first couple is 
                             related to the original flowSet, before the down-sampling process, while the latter refers to the downsampled 
                                flowSet'),
                             p('The figure shows the outliers', strong("score "), 'or ', strong("density "), 'for each of the samples depending on
                             the type of algorithm selected (see below the additional notes on the downsampling mechanism). Each down-sampling 
                             algorithm assignes a score|density to every single event: in general, the higher|lower is the score|density, the higher 
                             is the probability that the event is isolated from the rest, namely there is a good chance that the related event is an 
                             outlier: the greater the score, the greater the outlierness. With this, the way to evaluate this score|densities values 
                             assignment is to order them starting from the higher|lower scores|density'),
                             p('Notice than the ordered scores|densities plot restricts the view taking into consideration the amount of events of 
                             the smallest sample. That is to better focus on the most interesting part of the trend: the higher|lower score|density 
                             for every single sample'),
                             p('While the ordered scores|densities plot is mostly useful to find a possible valuable cut-off point to set your
                               percentage amount or the cut-off value, the purpose of the normalized scores|densities plot is just to have a global 
                               view of the trends of the scores|densities calculated. In fact, in this case all the values are normalized in order 
                               to let you compare the trends along the different samples. Take this into account when you deal with experiment with 
                               samples collecting a very different number of events.'),
                             
                             plotlyOutput(outputId = "downsampled_b4Score")%>% withSpinner(),br(),
                             p('The following plot shows the densities of the score values. Since all the values are normalized both in terms of 
                             densities|scores - on the abscissa - and in terms of event amount - on the ordinate - the graph allows a comparison 
                             of the score trends of the for each single sample. This graph could be useful because the plot above focuses only on 
                             the first ordered number of events, while with this one you can realize how the values are distributed along the whole 
                             range of every samples. Notice that, a good distribution is the one in which all the samples have an evident tail with 
                             a pronunced elbow towards the "1" (the maximum normalized value). For ths spade algorithm, on the other hand, 
                               a good distribution is when all the samples have this kind of trend towards the "0"'),
                             plotlyOutput(outputId = "downsampled_b4Dens")%>% withSpinner(),br(),
                             h4('The following plots are related to the flowSet scores after the down-sampling process'),
                             plotlyOutput(outputId = "downsampled_afterScore")%>% withSpinner(),br(),
                             plotlyOutput(outputId = "downsampled_afterDens")%>% withSpinner(),br(),
                             h4('Additional notes on the downsampling mechanism:'),
                             p('Warning: In case the number of events, in one of the samples in your entry flowSet, is less than MIN_SAMPLE_LENGTH 
                             =50, the downsampling is not possible'),
                             p("You can filter the outliers in three different ways: ", strong("by percentage "), "or ", strong("by value "), " or ",
                             strong ("equalize"), ". When you select by percentage, the ", strong("exclude_pctile "), "is the threshold to filter out 
                             the events (e.g: with a threshold of 3% you filter out the 3% of the events with the higher score). In this case every 
                             single sample will be reduced by an amount which is related to its total number of events. When you choose by value, you 
                             cut out those events which are lower than the set threshold (set by the ", strong("density|score value"),". Of course in 
                             this case you are surely selecting only those events above the threshold but without knowing exactly their amount (the 
                             advice is to perform iteratively the selected algorithm with different thresholds. Keep in mind that the score values 
                             expressed in the plot are normalized (they are not the ones directly calculated by the selected algorithm - this to be 
                             able to compare the values on different samples). In any case, the downsampling by value could be a powerful option to 
                             filter out different amounts of events per sample. Especially when the sample's are very unbalanced, if the selected 
                             algorithm works well, you have the possibility to cut more the bigger samples, than the smaller ones, harmonizing with 
                               this the sample's size"),
                             p("When it needs to equilize the samples in terms of amount of events, select the ", strong("equalize "),"option. With
                               this you can filter out all the events which exceed the ", strong("equalize value"), "starting of course from those 
                               events which are considered outliers by the selected algorithm. Of course, in case of heavy cut-off, most of the 
                               filtered events come from the denser (bulky) zone of the hyperspace: the risk is to get rid of significative events, 
                               but that's the only way if you want to flat your samples"),
                             p('There are some exceptions. The first is obvious: a Random algorithm does not compute any score: it is possible 
                             to filter randomly only by percentage and not by value. It must be noticed that, in order to avoid any kind of biased 
                               outcome, the best is to select the Random algorithm with the Equalize type of down-sampling: select the amount of 
                               event you want to preserve, download your new samples and then, if you would like to remove the outliers, just 
                               perform a new down-sample by percentage or by value with your preferred algorithm. This is useful when the total 
                               amount of events exceeds of 1 million. In this case the system could run out of memory and to avoid this, the best is 
                               to select the down-sampling in the equalize way as a first step in order to reduce the overall sample’s weight'), 
                             p('The other exception is the spade algorithm (see ', 
                             tags$a(href="https://www.bioconductor.org/packages//2.12/bioc/html/spade.html", "spade"), ') that, as a score, it 
                             computes the event density: in this case the lower is the density the higher is the outlierness'),
                             p('The last exception is related to the LDE algorithm. The greater is the Local Density Estimate score, the greater 
                             the centrality. This last algorithm is used to downsample the most dense zone of the expression hyperspace and it 
                             can be used to perform a downsampling which allows the reduction the total amount of events. This can be useful for very 
                             big samples set, to reduce the amount of resources needed for the clustering analysis. The drawback is that you do not 
                             filter out the outliers but those part of events which are in the most dense zone of the hyperspace. Usually these are 
                             the most representative events, but in case you would like to focus on the rare population, downsampling the denser zones 
                             means to focus more on those events which are more isolated, namely by selectively favouring certain events, the quantities 
                             are distorted'),
                             p('If you can find any elbow, or a marked discontinuity along the ordered scores (or densities) plot, this could be a 
                                point of considerable interest to set your cut-off percentage for your downsampling, especially if your goal is to
                               remove the number of outliers, rather than simply downsample the amount of events. In this case, just compute the 
                               related percentage (recall that the ordered scores/densities plot focus on the first number of events, taking as 
                               a reference the smallest sample) to set your cut-off percentage'),
                             p('The same applies to a cut-off criteria based on values instead of percentage. In this case you set the value around  
                               the elbow point as reported by the order scores/densities. The only drawback is that you cannot infer a priori the 
                               amount of events you would filter out: iterate you down-sampling process until you find satisfying result'),
                             p('Not every algorithm is able to provide a sharp discontinuity, namely a clearly shaped elbow in order to be capable
                               of guide you through this analysis: the big question ', shiny::em("which are the outliers in my dataset?"), 'is one 
                               of the most complex and controversial you can find in statistical literature because (like for the clustering 
                               algorithm) the results strongly depends of the data to analyze, its complexity and even its magnitude in terms of 
                               number of elements, namely an algorithm that suit for every occasion does not exist. For this reason cytoChain 
                               provide a list of algorithm to choose from'),
                             p('If you do not notice any relevant difference between the plots reporting the trend before and after the 
                               down-sampling process, this means that you do not cut a sensible amount of events, which could be perfect choice if
                               you just do not want to strip-off a considerable amount of events to preserve the original relative percentage
                               when you proceed with the analysis. The problem is that, most of times you have to consider the downsampling process
                               to simplify and improve your dataset in order to safely perform any clustering process. More properly, to get 
                               rid of the outliers could be a must in order to deal with a qualitative treatable dataset. A straightforward way to 
                               evaluate this treatability and the related rendering process ability is to go through the Metadata & Assays workflow 
                               and in particular, to select the flowSet evaluation tab. The first analysis there, is to produce the scree plot of 
                               the eigenvalues from the PCA algorithm: this is the easiest way to evaluate the sparsity and the complexity of your 
                               dataset. A treatable flowSet is the one having the two principal components capable of handling at least the 70% of 
                               the relevant information. Again, this is not apodictic but it provides you a rule of thumb to know what is what and 
                               most of all, a quantitative index to compare different dataset')),br(), #condition = ("output.down_algoStatus")
                         actionButton('saveDown', "Produce your downsampled flowSet"), br(),
                         h4("Remember to produce your brand new samples each time you perform a data handling!"),br(),br(),
                         conditionalPanel(
                             condition = ("output.saveDown_panelStatus"),
                             downloadButton(outputId = "downloadDownfS",label = "download your downsampled flowFrames"),br(), br(),
                             downloadButton(outputId = "downloadDownfS_complete", label = "download your .fcs samples with both the handled  
                                    expression matrix and the original one")))))#tabDowns

######################## Concatenate TabPanel -----------------------------
tabConc <- tabPanel(title = "Concatenate", icon =icon("link"),
                    p("The concatenation in the cytometry jargon is to merge each single event of the entire flowSet in a single flowFrame"),
                    p("This is an essential step which is kind of hidden inside the rest of the analysis, and here is implemented just with the 
                      purpose to let the user to download the fcs file with the embedded 'sample_id' column"),
                    p("The concatenating is essential for the mapping phase of the high dimensional analysis workflow: only a concatenated sample set can be 
                    analyzed and compared along the various samples. In fact, a 'dimensionalty-reduced' map like tSNE or UMAP is a representation 
                    of the whole event of the entire data set within a single map: it is the representation of the concatenated flowSet"),
                    h3("Perform concatenating process"),
                    h3(textOutput("checktxt_concatenating_fS")),
                    actionButton("runConc", "Run the concatenating algorithm", icon("paper-plane"), 
                                 style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),br(),
                    conditionalPanel(
                        condition = ("output.conc_panelStatus"),
                        #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                        h3("Push to get some summaries of the concatenated samples"),
                        actionButton('runConcShow', 'Visualize the concatenating reports'),
                        conditionalPanel(
                            condition = ("input.runConcShow !==0"),
                            h4("These are the details of the concatenated flowFrame"),
                            tableOutput(outputId = "concatenated")%>% withSpinner(),
                            h4("This figure shows the densitys of the concatenated flowFrame"),
                            plotlyOutput(outputId = "marker_density_conc")%>% withSpinner(),br(),
                            downloadButton(outputId = "downloadConcfS",label = "download your concatenated flowFrame")))
                    )

######################## Meta TabPanel -----------------------------
tabMeta <- tabPanel("edit metadata", icon =icon("pen"), 
                    h4("The aim of this workflow is to set additional variables to combine with the sample_id. In the subsequent analysis these new 
                    entries could help in distinguish the various population types and to compare the different sub-population based on these 
                    additional parameters. With this is possible to enter some advanced characteristic of the sample and to associate them to the 
                    subsequent cluster analysis. You can appreciate the benefit of these meta-data entries in the plots generated of the high 
                    dimensional analysis workflow. The scope of the rest of the tabs ('", strong("flowSet analysis"), ", '"
                      , strong("flowSet evaluation"), " and '", strong("map evaluation"), ") is to provide some high level information about the 
                      samples, without entering into the cell's details."), br(),
                    h4("There are three tables to edit in order to set this meta-data: the", strong("sample meta-data"), " , the ", 
                       strong("markers meta-data")," and the", strong("phenotype meta-data")),
                    h4("Load the sample files inherited from the flowSet optimization process (if you followed the that workflow process)"), 
                    h4("If you prefer, you can load a new flowSet from the filebox 'File input' on the left sidebar panel"),
                    h4("In the same left sidebar panel, you have the possibility to load the meta-data by three .csv files: one for the sample 
                    meta-data, one for the markers and one for the phenotype meta-data. These tables are kind of labour-intensive to edit. 
                    Especially when you deal with the same type of analysis, it could be useful to load them once you already edited, and so you 
                    can use them without the need to re-edit them all the times. The advice is to edit once, save it and, in case, re-load them 
                      before to pass to the analysis step. In any case, "),
                    h4("The files loaded by the side-bar panel take priority over the ones loaded within the flowSet optimization workflow."),
                    h4("Please also notice that this is the only place where you can load your samples' description for the subsequent 
                    analysis. Once you loaded and costimized your table entries, push the red button below to actually load both table and 
                       samples"),
                    h4("Be sure your meta-data table is congruent with the samples you loaded!"),br(),
                    
                    actionButton("runMeta", "Load the meta-table of the selected samples", icon("paper-plane"), 
                                 style="color: #fff; background-color: #008000; border-color: #2e6da4"),br(),
                    rHandsontableOutput(outputId = "hot_meta")%>% withSpinner(),
                    h4("Please edit the table in accordance to your flowSet description and enter the modified table 
                       pushing the botton below"),
                    h4("Once you edited, the table could be downloaded to be further handled in your favourite way: the 'condition' and 'patient_id' 
                       are just tags suggested to fill the row entries the easiest way on the run, but you are free to define your own tags when you
                       manipulate your .csv file to be loaded for the meta_data definition"),
                    h4("Attention! The purpose of the file loding widget on the left is to let you load your samples to analize without passing by 
                    the ", strong('flowSet optimization workflow'), ". For this reason these files have priority. ", strong('Skip it'), "if you want 
                     to keep using the ones handled in the dataset optimization workflow"),
                    h4("Warning! The samples should be edited as ", strong("day/month/year"), " and in time order: the oldest dates must come first and selected as 'time_00'"),
                    p("Since each single sample is loaded holding the order in your directory, the best would be to name your files according to 
                      the order of the various steps. On the other hand, if you are not interested in the time levels, you can skip this, holding
                      the list automatically entered or setting just one date/step for the whole sample's set"),
                    p("The dates in the related column should be in the english format, namely: '%month/%day/%year'"),
                    h3(textOutput("checktxt_loading_meta_fS")), 
                    conditionalPanel(
                        condition = ("output.meta_panelStatus"),
                        actionButton("runMetaload", "Push to use this edited meta-table to be used 4 furter analysis", icon("paper-plane"), 
                                 style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                        hr(),
                        p("Following are the barplots depicting the colors used in the maps and in the histograms along your samples in the 
                          analysis. Just the tag with at least two different entries are shown. Keep a common symbolic entry if you are not 
                          interested in one tag"),
                        
                        conditionalPanel(
                            condition = "output.condition_tag1",
                            plotOutput(outputId = "tag1.color"),
                            plotOutput(outputId = "tag1.bar"),br()),
                        conditionalPanel(
                            condition = "output.condition_tag2",
                            plotOutput(outputId = "tag2.color"),
                            plotOutput(outputId = "tag2.bar"),br()),
                        conditionalPanel(
                            condition = "output.condition_tag3",
                            plotOutput(outputId = "tag3.color"),
                            plotOutput(outputId = "tag3.bar"),br()),
                        conditionalPanel(
                            condition = "output.condition_tag4",
                            plotOutput(outputId = "tag4.color"),
                            plotOutput(outputId = "tag4.bar"),br()),br(),
                    
                        conditionalPanel(
                            #tags$audio(src = "dingdong.wav", type = "audio/wav", autoplay = NA),
                            condition = ("input.runMetaload !== 0"),br(),
                            downloadButton('downloadmetaResume', "Download the resume of the event distributions in csv format"),br(),br(),
                            
                            h4("The scope of the table below is to select some markers. This selection is relevant for the rest of the analysis, 
                            because only the marker of your choice will be considered. As a matter of fact, even if the entire panel may include 
                            even more than 30 different markers (aka stains, aka dimensions, ...), it is certainly not recommended to use more 
                            than 10÷15 elements in the first steps of the analysis. In general it is not advisable to select all the markers, 
                            especially to find a specific population with a certain phenotype. You could select a subset of markers with the 
                            following editable table and the to extend on a step by step  basis, the total number of markers"), br(),
                            h4("Attention! All samples must have a common panel (an identical marker set) to be analized. If some of your samples
                               contains a different marker set, please consider to further split your analysis accordingly'"),br(),
                            h4("Load the markers details inherited from the flowSet optimization process (if you followed the process)..."),
                            actionButton("runMetaMarkers", "Load the markers' meta-table of the selected samples", icon("paper-plane"), 
                                 style="color: #fff; background-color: #008000; border-color: #2e6da4"),br(),
                            hr(),
                            rHandsontableOutput(outputId = "hot_markers")%>% withSpinner(),
                            h4("Edit the table in accordance to your flowSet description and enter the modified table pushing the botton below"),
                            h4("Warning! At least two markers must be selected"),
                            actionButton("runMetaMarkersLoad", "Push to use this edited markers' meta-table to be used 4 further analysis", 
                                     icon("paper-plane"), style="color: #fff; background-color: #ff0000; border-color: #2e6da4"), br(),
                            conditionalPanel(
                                condition = ("input.runMetaMarkersLoad !== 0"),
                                h3(textOutput("checktxt_loading_meta_marker")), 
                                #tags$audio(src = "dingdong.wav", type = "audio/wav", autoplay = NA), br(),
                                #tags$hr(style="border-color: green;"),br(),
                                p("With this last table, it is possible to set the phenodata table for the cell phenotype definition.
                                  This is the only table which is not mandatory, since you can perform every analysis in cytoChain even without 
                                  this entry"),
                                p("The user could define each single population of interest describing the type of marker expression that could be 
                                  either: ", strong("positive ('+' means 'well expressed')"), "or ", strong("negative: ('-' means 'thinly 
                                expressed')"), "or ", strong("non-relevant (with the wildcard symbol '*')"), ". This last choice is meaningful for 
                                those cases in which, for a particular phenotype is indifferent whether the related marker is expressed or not; in 
                                all, not every marker selected for the differentiation are relevant for the definition of the single phenotype."), 
                                  p("Usually, when multiple combination of marker is defined, not every phenotype will be mapped in a cluster set: 
                                  in this case some of the cluster (or meta-cluster) will be ", strong("'unassigned'"), ". In short, only some of 
                                  the meta-clusters could be assigned to a defined phenotypes, while some of them will remains unassigned."),
                                h4("Attention: for the phenotypes use some meaningful and, most of all, syntactically correct names, that is, do 
                                   not use any special characters (the only one allowed is the under",strong("_"),"score) and do not use just a
                                   number"),
                                h4("Edit the table entering the phenotype you want to analyze"), 
                                
                                actionButton("runPhenoTable", "Load the phenotype meta-table 4 the selected samples", icon("paper-plane"), 
                                         style="color: #fff; background-color: #008000; border-color: #2e6da4"),br(),
                                rHandsontableOutput(outputId = "hot_pheno")%>% withSpinner(),
                                p("Modify the table adding more phenotypes to analize (the entries showed are just examples)"),
                                p("Each of the elements could have one of the three characters value: '+', '-', or '*', with the following 
                                relationship: ", strong("pos_threshold (+) >= neg_threshold (-)")),
                                p("The negative and positive thresholds are set by the related labeling widgets on the left"),
                                actionButton("runPhenoTableload", "Push to use this edited pheno-table to be used 4 further analysis", 
                                             icon("paper-plane"), style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                                p("click with the mouse's right button if you want to save each table as .csv file"),
                                h4("Attention: the markers and the phenotype names are case sensitive! Please use the same marker's names in your 
                                phenotype tables when you edit the table and you load the related the .csv table from the left sidebar panel."),
                                h4("Warning: if you are going to load a phenotype table generated by cytoChain during an analyisis session 
                                   (tipically the one you may download from the 'Clustering' and the 'MetaClustering' tab of 'classical 
                                   workflow') please notice that the first column reports the clusterId. You should strip-off this column before 
                                   you load here'"),br(),
                                h3(textOutput("checktxt_loading_meta_pheno")), 
                                
                                conditionalPanel(
                                    condition = ("output.meta_phenoStatus"),
                                    #tags$audio(src = "dingdong.wav", type = "audio/wav", autoplay = NA), br(),
                                    p("Following are the colors used in the maps and in the histograms for the defined phenotypes"),
                                    plotOutput(outputId = "pheno.color"),
                                    
                                    p("Just because you may set up the analysis with or without a pheno-table input, this is buttton let you get 
                                        rid off a previously loaded table"),
                                    actionButton("resetMetaTable", "Push to reset the pheno-data table", 
                                             icon("paper-plane"), style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),br(),
                                    verbatimTextOutput("pheno_file_status"),
                                    
                                    h4("Additional notes on the phenotype table:"),
                                    p("It could also happen that, the same phenotype could span along more than one cluster, but also it 
                                    could be the case in which more than one phenotype should be assigned to the same set. This is the typical case 
                                    in which the phenotype is defined with one or more “*” (the wildcard). In this case the signature is undefined 
                                    because it is shown with the 'OVERLAP' tag (because it is an overlapping of two or more signatures), in all maps 
                                    and plots. An analysis of the heatmap in the 'meta-cluster' page tab should help to reconstruct the correct list."),
                                    p("Especially for the study of the various heatmaps produced by the tool along the analysis, it could be useful 
                                    to neatly lay out your set of stains, e.g. placing on left the most meanigful markers. To get this kind of 
                                    result, just dispose your favourite marker order in your phenotype table and you'll get the same order in the 
                                    various heatmaps"),
                                    p("There is no limit on the possible choices and this could be an extremely powerful tool to analyze the 
                                    phenotype of each single event following how the related proportions change overtime  or correlating the 
                                    population by conditions entered as tag1, tag2, tag3 and tag4 in the first meta-table. The chance to define all the 
                                    possible phenotype in terms of marker expression could lead the researcher to better picture the marker’s 
                                    interactions and to frame possible new cell lineages and their relation with certain subpopulation (even the 
                                    smallest one)."))))))


######################## Analysis TabPanel -----------------------------
tabAna <- tabPanel("flowSet analysis", icon =icon("balance-scale"), 
                   p("Here it is possible to analize the select flowSet in terms of differences or similarities between the samples"),
                   p("This is done through the", strong("Multidimensional scaling - MDS"), " which is a visual representation of distances or 
                   dissimilarities between sets of objects. The objects in this case are the selected samples (loaded and treated though the 
                     pre-processing workflow, or directly in the meta-data workflow)."),
                   p("the samples that are more similar (or have shorter distances) are closer together on the MDS graph than objects that are less 
                   similar (or have longer distances). MDS can also serve as a dimension reduction technique for high-dimensional data, like in our 
                     samples."),
                   p("If the sample metadata are entered, the MDS plots will let you figure out the calculated distances in term of similarities 
                     between the samples itself or depicted in terms of the conditions entered as 'tag1','tag2' or 'tag3' in the meta_data table"),
                   p("In brief, the MDS is calculated through a sequence of steps:"),
                   p(" - 1) Assign a number of points to coordinates in an n-dimensional space (in our case the dimensions are the selected markers).
                   The axes orientation and the coordinates themself are arbitrary and they do not have a physical meaning."), 
                   p(" - 2) Calculate the distances for all pairs of points (namely, in case of euclidean distance, this is the straight-line 
                   between two points x and y in Euclidean space (Pythagorean theorem) for n-dimensional space. This results in the 
                     similarity matrix."),
                   p(" - 3) Compare the similarity matrix with the original input matrix. This is done to fit the data to be represented  in a 'two 
                     dimension' space (which may be represented by a figure, like the plot that will be shown in this case)"),
                   p(strong('Notes: ')," The right distance is not always the euclidean! If you are pretty confident about the differences in 
                     your samples through different categories (tags), try to perform iteratively the MDS computation with different distances"),
                   p("If you want to know more about the distances and the so called ", shiny::em('curse of dimensionality'), " see: ", 
                     tags$a(href="https://www.datanovia.com/en/lessons/clustering-distance-measures/","clustering distances"), " , ", 
                     tags$a(href="https://pages.mtu.edu/~shanem/psy5220/daily/Day16/MDS.html","MDS"), " , ", 
                     tags$a(href="https://en.wikipedia.org/wiki/Curse_of_dimensionality","The curse of dimensionality")), br(),
                   
                   h4("Push the button to show several MDS (Multidimensional Scaling) and diagnostic plot to show the 
                       distances between the different samples"),
                   actionButton("runflowSetAna", 'Show plots'),
                   conditionalPanel(
                       condition = ("input.runflowSetAna !== 0"),
                       h3(textOutput("checktxt_plotting_mds")),
                       conditionalPanel(
                           condition = ("output.flowSetAna_panelStatus"), 
                           plotOutput(outputId = "MDS_HD",width = "1230px", height = "1200px"), br(),
                           #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                           plotlyOutput(outputId = "MDSsample_id",width = "1300px", height = "650px") %>%  withSpinner(), br(),
                           conditionalPanel(
                               tags$hr(style="border-color: green;"),
                               condition = "output.condition_tag1",
                               plotlyOutput(outputId = "MDStag1",width = "1300px", height = "650px")  %>% withSpinner(), br()),
                           conditionalPanel(
                               tags$hr(style="border-color: green;"),
                               condition = "output.condition_tag2",
                               plotlyOutput(outputId = "MDStag2",width = "1300px", height = "650px")  %>% withSpinner(), br()),
                           conditionalPanel(
                               tags$hr(style="border-color: green;"),
                               condition = "output.condition_tag3",
                               plotlyOutput(outputId = "MDStag3",width = "1300px", height = "650px")  %>% withSpinner(), br()),
                           conditionalPanel(
                               tags$hr(style="border-color: green;"),
                               condition = "output.condition_tag4",
                               plotlyOutput(outputId = "MDStag4",width = "1300px", height = "700px")  %>% withSpinner(), br()),
                           p("As you can see from the plot below, you can easily respond to the following questions: "),
                           p("The closer points represent similar samples as your expectation?"),
                           p("Which of the metadata could better describe the distance between the samples?"), 
                           p("Are the metadata reflecting the point distribution in the multidimensional space?"),br(),
                           tags$hr(style="border-color: green;"),br(), br(),
                           p("The plots below are for diagnostic, but it could be helpful to see the various stains trends, along the expression 
                           matrix, all grouped with the main metadata level information (tag1, tag2 and time_step)"),
                           plotOutput(outputId = "dia_tag1") %>% withSpinner(), br(),
                           tags$hr(style="border-color: green;") %>% withSpinner(),
                           plotOutput(outputId = "dia_tag2") %>% withSpinner(), br(),
                           tags$hr(style="border-color: green;"),
                           plotOutput(outputId = "dia_tag3") %>% withSpinner(), br(),
                           tags$hr(style="border-color: green;"),
                           plotOutput(outputId = "dia_tag4") %>% withSpinner() )))

######################## Evaluation TabPanel -----------------------------
tabEva <- tabPanel("flowSet evaluation", icon =icon("binoculars"), 
                   p("This page is devoted to the evaluation of the selected flowSet in terms data characteristic, like sparsity of the points in 
                   a multidimensional space, the grouping of the points in possible clusters, and even the shape of the group of points"),
                   p("In our case, the flowSet is composed of several flowFrames (the samples) and each flowFrame is basically a collection of 
                   events which are depicted as points in a multidimensional space. These evaluation algorithms are based on a single concatenated 
                   flowFrame with the whole set of events present in each of the original sample"),
                   p("First we propose a very easy method to evaluate the quality of the set of events. Definitely the Principal Component Analysis 
                   (PCA) is not suited for the flow cytometry data, because of the non-linearity of the event distribution in a typical sample. For 
                   this reason it is always preferred to utilize a tSNE-like representation, instead of ISOMAP or PCA. In any case, PCA is the most 
                   straightforward and light algorithm for dimensionality reduction (the concept behind PCA is very much similar to the MDS) and it 
                   can provide a good indication of how much spread are the events in the sample, with a equi-spaced set of events being the most 
                   sparse and difficult to reduce set, while in a event space with a perfectly linear point distribution, the point can be 
                   trivially represent in one axe only. On the other hand, in an extremely distributed set of events, the complex tSNE-like map 
                   could also not converge and this is the reason why we introduced this direct evaluation with PCA. In quantitative terms, the 
                   relative variance for each variable (aka channel or marker in this context) and the covariance between them can be a useful way 
                   to characterize the event space, the variable's behaviour and thier mutual interaction."),
                   p("After this first evalution the following steps are to evaluate the number of possible clusters to group the sample's events: 
                   the first evaluation is performed by the Kmeans clustering, while the last one is by the several other statistical methods by the 
                     mclust package"),br(),
                   h4("Push the button to show methods to evaluate the event sparsity inside the loaded flowSet"),
                   actionButton('runPCAEva', 'Show plots'), br(),
                   conditionalPanel(
                       condition = ("input.runPCAEva !== 0"),
                       h3(textOutput("checktxt_runPCAEva")),
                       conditionalPanel(
                           condition = ("output.flowSetEva_panelStatus"),
                           #tags$audio(src = "shipsbell.wav", type = "audio/wav", autoplay = NA),
                           p("Visualize eigenvalues (scree plot). Show the percentage of variances expressed by each principal component"),
                           plotOutput(outputId = "pca_res") %>% withSpinner(), br(),
                           p("When the trend shows a drastic reduction after the two first principal components and the cumulative proportion of 
                           PC1+PC2 is well above than 50%, other more complex dimensionality reductions algorithm like t-SNE could be performed
                           with no problem with a good clustering performance. On the other hand, if you are far from this safety situation, it is 
                           common to face some difficulties because the events are too sparse or too much equally distributed. In these cases, at 
                           worst, the dimensionality reduction algorithms may also not converge and fail to process and, at best, the maps could 
                           be not very efficient in differentiating the various clusters (clusters could spread along the map, there are many 
                           overlapping clusters and so on...). PCA is a linear transformation: the error we introduce, squeezing the", shiny::em("N"), 
                           " variables to ", shiny::em("2")," is relevant: ", shiny::em("1 - (the variance of the two first principal components)")),
                           p("the first two values are shown below"),
                           tableOutput(outputId = "pca_df"), br(),
                           p("the table below shows the contributions (in %) of the variables to the principal components. The variable (var) 
                             contribution to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component)
                             where the cos2 is related to variable covariance"),
                           tableOutput(outputId = "pca2_df"), br(),
                           p("The following is a barplot of the variable contribution. A high cos2 indicates the importance of the variable on the 
                             pricipal components which, in the correlation circle, is emphatize with longest strokes"),
                           
                           plotOutput(outputId = "cos2plot", width = "1300px", height = "1300px") %>% withSpinner(), br(),
                           p("This is the correlation circle, a representation of the covariance for the various variables"), br(), 
                           plotOutput(outputId = "cos2var", width = "1320px", height = "1320px") %>% withSpinner(), br(),
                           h4("Push the button to show the elbow method to estimate the number of cluster with kmeans algorithm"),
                           actionButton('runKmeansEva', 'Show plots'),
                           h4("Warning! Depending on the entry dataset volume and on the max number of cluster set, the calculation may take a long 
                              time!"), br(),
                           conditionalPanel(
                               condition = ("input.runKmeansEva !== 0"),
                               plotOutput(outputId = "kmeans_out", height = "1200px") %>% withSpinner(),
                               tableOutput(outputId = "kmeansTable_out") %>% withSpinner(),
                               p('The plot allows a certain evalution of the number of possible clusters. In particular it shows the wss, namely the 
                               within-cluster sum of square, against the number of clusters. The first clusters will add much information, but at 
                               some point the marginal gain will drop, giving an angle in the graph. The number of clusters is chosen at this point, 
                               hence the “elbow criterion”. Thanks to Sunny Anand (see ', 
                                 tags$a(href="https://www.r-bloggers.com/finding-optimal-number-of-clusters/", "finding number of clusters)")),
                                 p("This “elbow” point cannot always be unambiguously identified, and this is the case of most of the common 
                                 flow cytometry data. So then, for our application, it is more useful to infer the number of clusters like a 
                                 projection of a function limit: is the curve tangent with 0 slope? which will be the real number if you extrapolate 
                                 it from the curves? Try to increase the number of cluster evaluation with the related widget on the left sidebar if 
                                 you feel like being too far from the elbow to infer a meaningful point"),
                               h4("Push the button to show the plot to estimate the number of cluster with mclust algorithm"),
                               actionButton('runMclustEva', 'Show plots'),
                               h3("Attention! Depending on the entry dataset volume and on the max number of cluster set, the calculation may take a 
                              very long time!"), br(),
                               conditionalPanel(
                                   condition = ("input.runMclustEva !== 0"),
                                   plotOutput(outputId = "mclust_out") %>% withSpinner(),
                                   tableOutput(outputId = "mclustTable_out")%>% withSpinner(),
                                   p("The plot shows several different models. Choose the one with the best estimation. Also in this case the 
                                   optimal number cannot always be unambiguously identified"),
                                   p("To better understand the output you can study the mclust manual, (see ", 
                                     tags$a(href="https://cran.r-project.org/web/packages/mclust/index.html", "mclust references"), " or the ",
                                     tags$a(href="https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html", "mclust vignette"), ")."),
                                     p("In short, the mclust packaege comprise a set of algorithms which are based on the Gaussian finite mixture 
                                     model fitted by EM. While the kmeans algorithm is agnostic in terms of knowledge of the typology of events, 
                                     the algorithms used by mclust are statistically based: the main assumption is that the events have Gaussian 
                                       distribution, which is the the case of the flow cytometry data"),
                                       p("But there could be different ways to classify a Gaussian distribution. In particular mclust categorize 
                                       three ways to classify the distribution: shape, volume and orientation. The three capital letters referes 
                                       to those distrbution typology"),
                                       p("E.g.'VEV' means ellipsoidal, with equal shape, 'EEE' means ellipsoidal, equal volume, shape, and 
                                       orientation, 'VVV' which means ellipsoidal, varying volume, shape, and orientation. This is the most generic 
                                       type of gaussian distribution for mclust: each cluster could have different volume, different shape 
                                       (spherical, diagonal or ellipsoidal) and different orientation, which is the case for most of our 
                                         experimental samples"))))))


######################## Map evaluation TabPanel -----------------------------

tabtSNEeva <- tabPanel("Map evaluation", icon =icon("schlix"),
                   p("In this tab you can produce and download the tSNE maps related to the flowSet handled in the 'flowSet optimization workflow' tab"),
                   p("In case you entered the flowSet directly in the 'Metadata and Assays' workflow, this process works only for the 
                   initial flowSet aka the 'entry flowSet'. Select at least one of the flowSet optimization action if you want to analyze the 
                     tSNE maps and compare them: the initial unhandled flowSet versus the one manipulated, namely the outcome of the last performed 
                     process in the chain: cleaning, scaling, aligning and down-sampling"),
                   p("The aim here is to compare the two maps, before and after the flowSet handling. The two  different tSNE plots are 
                      related to the original flowSet and the one related to the flowSet handled after the latest relevant action you performed 
                      within the flowSet optimization workflow"),
                   p("You can also perform the clustering of the tSNE map dimensions only. This could be useful especially for complicated flowSet, 
                     when you want to distinguish and differentiate and assign a name the various 'islands' of the tSNE map"),
                   h4("Please also notice than the tSNE map is built using the linage markers selected in the 'edit metadata' tab, all the 
                      rest of dimensions are neglected here"),
                   p('thanks to Kamil Slowikowski	for the density color, see ', 
                     tags$a(href="https://slowkow.com/notes/ggplot2-color-by-density/", "color by density")),br(),
                   h4("Warning! tSNE map calculation is one of the longest task: so then, be patient..."), br(),
                   actionButton('runfStSNEinit', 'Generate the tSNE map of the entry flowSet'),br(),
                   h3(textOutput("checktxt_eva_init")),
                   conditionalPanel(
                       condition = ("output.tSNEinit_panelStatus"),
                       h4("This is the tSNE plot of the entry flowSet..."),
                       plotOutput(outputId = "mapinit_out",  width = "1450px", height = "1350px") %>% withSpinner(),br(),
                       actionButton('runKtSne', 'Perform the tSNE map clusterization'),br(),
                       conditionalPanel(
                           condition = ("output.runKtSne_panelStatus"),
                           h4("This is the clusterized tSNE map of the flowSet"),
                           plotOutput(outputId = "tSNEkmeans_clust_init", click = "plot_clik", width = "1450px", height = "1350px"),
                           verbatimTextOutput("tSNEkmeans_clust_info"),
                           downloadButton(outputId = "downloadKmeansTsne",label = "download your concatenated flowFrame with the tSNE's and the 
                                              kmeans dimensions"))),

                    actionButton('runfStSNEfinal', 'Generate the tSNE map of the handled flowSet'), br(),
                    h3(textOutput("checktxt_eva_final")),
                    conditionalPanel(
                        condition = ("output.tSNEfinal_panelStatus"),
                        h4("This is the tSNE map of the handled flowSet"),
                        plotOutput(outputId = "mapfinal_out",width = "1450px", height = "1350px") %>% withSpinner(), br(),
                        actionButton('runKtSne_unique', 'Perform the tSNE map clusterization'),br(),
                        conditionalPanel(
                            condition = ("output.runKtSneunique_panelStatus"),
                            h4("This is the clusterized tSNE map of the flowSet"),
                            plotOutput(outputId = "tSNEkmeans_clust_unique", click = "unique_plot_clik", width = "1450px", height = "1350px"),br(),
                            verbatimTextOutput("tSNEkmeans_clust_unique_info"),
                            downloadButton(outputId = "downloadKmeansunique",label = "download your concatenated flowFrame with the tSNE's and the 
                                              kmeans dimensions"))))

######################## Clustering TabPanel -----------------------------

tabClust <- tabPanel("Clustering", icon =icon("object-ungroup"), 
                     p("In this section you start the clustering analysis. Currently, one of the most used clustering algorithm used in cytometry is 
                       (", tags$a(
                           href="https://bioconductor.org/packages/release/bioc/html/FlowSOM.html", "flowSOM"), ") ", "because is by far, the best in 
                         terms of accuracy and in terms of resource-demanding"),
                     p("In brief, flowSOM is based on the ", strong("Self Organized Map - (SOM)"), "algorithm, a type of artificial neural network 
                     (ANN) an unsupervised learning approach which uses as the input training samples, the space of the events. The product is 
                     a two-dimensional discretized representation of the input space of the training samples, called a 'Kohonen map'. You can find 
                       some good introduction to the algorithm in the ", tags$a(href="https://en.wikipedia.org/wiki/Self-organizing_map", 
                                                                                "related wikipedia page"), ", and also in the following pages: ",
                        tags$a(href="https://annalyzin.wordpress.com/2017/11/02/self-organizing-map/", "SELF-ORGANIZING MAPS TUTORIAL"), 
                        "and in the ", tags$a(href="https://www.superdatascience.com/blogs/the-ultimate-guide-to-self-organizing-maps-soms/", 
                                              "The Ultimate Guide to Self Organizing Maps (SOM's)"), ". Like all the AI based algorithms, 
                                              the production of SOM's works better with a lot of events and it is tipically set with a large number 
                                              of clusters - For this reason it is essential to reduce the amount of clusters with some meta-clustering algorithm 
                     which will group the intial amount to a less numerous set. That's why the meta-clustering is an essential step in this wokflow, 
                        while using other clustering approaches, this additional step is often left out."),
                     p("Another valuable clustering algorithm is ", strong("Phenograph"), ", in the two different implementations ", 
                     strong("Rphenograph"), " and the faster ", strong("FastPG"), "which could perform better in terms of accuracy for 
                     little experiments (with less than half a million events). You cannot choose the amount of clusters you want, but is the 
                     Phenograph itself which select the best number of clusters based on the number of nearest neighbours events. The modularity,
                        which is one of Phenograph outcome is a measure of the connectedness of a clustered network. When comparing different 
                        clusterings of the same network (same 'K'), the one with the higher modularity is better"),
                     p(strong("Attention: The execution of the clustering process resets all the data related to previous meta-clustering and labeling 
                     processes (if already perfomed)")),
                     p("Start from here when you set up a brand new analysis among those available, namely:",
                        tags$ol(tags$li("Clustering & Meta-clustering - follow in this order: ", strong("Clustering "),">> ", 
                                        strong("Meta-Clustering "), ">> ", strong("Mapping "), ">> ..."),
                                tags$li("Clustering & Labeling: follow in this order: ", strong("Clustering "), ">> ", strong("Labeling I "), ">> ", 
                                        strong("Mapping "), ">> ", strong("Labeling II "), ">> ..."), 
                                tags$li("Clustering & Labeling & Meta-clustering - follow in this order: ", strong("Clustering "), ">> ", 
                                        strong("Labeling I "), ">> ", strong("Meta-clustering "), ">> ", strong("Mapping "), ">> ", strong("Labeling II "),
                                        ">>..."))),br(),
                     
                     p(strong("Additional notes on the high dimensional analysis workflow")),
                     p("The high dimensional analysis workflow provides a clustering operation (using one of the selected clustering 
                             algorithm) which may be followed by a meta-clustering. Reducing the amount of groups of similar events to a reasonable 
                             number of homogeneous sets could help in simplifying the analysis the data results, the map plots and so on. We named 
                             the operation of mapping a set of edited phenotypes (by the phenotype table entered in the 'Edit Metadata' tab) to the 
                             various calculated clusters, the ", strong("labeling"), "operation. It is possible to map the various entries expressed 
                             by the phenotype matrix to the clusters (within the 'Labeling I' tab) and to further map the phenotypes in the 
                             meta-clusters (within the 'Labeling II' tab)."),
                     p("The main high dimensional analysis strategies in this workflow are:", 
                       tags$ol(tags$li("1.	Clustering + Meta-clustering "),
                               tags$li("2.	Clustering + Label "), 
                               tags$li("3.	Clustering + Label + Meta-clustering"),)), 
                     p(strong("- 1: Clustering + Meta-clustering: "), "Let’s suppose that you do not have any idea about the possible 
                             phenotypes which could be present inside your experiment. Of course you select a certain panel to get some results 
                             from your experiment but as often happens, apart from some set of lineage markers, it could be difficult to compose a 
                             proper data frame to enter in the 'Metadata and Assays'. Another possibility is that you just do not care about any 
                             specific phenotype and this is a common case especially in high dimensional flow cytometry analysis when you just set 
                             up a big panel with more than 20 markers without dealing with any specific indication about the outcome of your 
                             experiment. In this case the clustering algorithm will split up your total set of event (collected in a single 
                             concatenated flowFrame) assigning each of them to a specific cluster of a set of 'Num_cluster': so then, each event 
                             belonging to a cluster will be assigned an additional dimension (a clusterId). Since most of times, especially with 
                             very complex experiments, you will deal with a lot of clusters, it is always better to reduce this number of clusters 
                             with the meta-clustering operation. The meta-clustering operation will just re-group each single cluster within a 
                             meta-cluster, namely a “cluster of clusters. "),  
                     p(strong("2.	Clustering + Label: "), "The meta-clustering could reduce the number of entries and thus simplifying 
                             all the analysis but, as a drawback, it can reduce the detailed sub-setting introduced by the selected clustering 
                             algorithm. In fact, since a cluster could have a certain median expression for each single channel, a meta-cluster will 
                             have a median of a set of clusters collecting altogether more groups of events. Furthermore, while a panel of N 
                             channels could have 2^N combinations corresponding to the number of different phenotypes, a clustering algorithm 
                             calculating the expression of K clusters while discriminate K > M number of phenotypes (in M is the number of 
                             meta-clusters). This strategy is suitable in those cases in which it is important to find a very specific phenotype 
                             which can be formed with very few amount of events, namely a case of rare events. In this case it can be meaningful to 
                             skip the meta-clustering and concentrate just to the clustering output with K (> M) different groups of events. On top 
                             of this advantage you can also put a label to your phenotype in the related table entry (through the 'Edit Metadata' 
                             tab). Unlike the previous strategy, you have to enter the phenotype description table in order to perform the labeling 
                               action"),
                     p(strong("3.	Clustering + Label + Meta-clustering: "), "This is the most complete analysis in which it is possible to 
                             join both clustering and meta-clustering. Like in the previous analysis strategy you have to enter the phenotype table 
                             first. In this way both clusters and meta-clusters phenotype will be named as per your edited table."), 
                     p(strong("Additional notes on some important parameters for this workflow")), 
                     p("Along with the number of cluster or meta-cluster, there are two pairs of important parameter to take into account, 
                               collected here on the left sidebar under the name of ",strong("Labeling parameters"), "because they tune how the 
                               various phenotypes are classified. The heatmaps produced under the ", strong("Labeling I "), "and the ", 
                       strong("Labeling II "), "tabs are based on the normalization of the expression matrix for the selected markers of the 
                               events grouped by the clustering (or meta-clustering) algorithm. The first pair of parameters are the ", 
                       strong("minQuantile/maxQuantile "), "that is the lower/upper limits for the zero transformation. These two thresholds 
                               (which we recommend to set them symmetrically like for the default setting: minQuantile = 0.05 and maxQuantile = 0.95) 
                               affect the normalization (like for the zero transfomation - see the manual), setting the threshold for the 
                               positive/negative MFI (Mean Fluorescent Intensity) values. As the ranges of marker intensities can vary substantially, 
                               they set the low and the high percentiles as a boundary to force the marker MFI to 0 or to 1. This normalization of the 
                               transformed data can give better representation of relative differences in marker expression between annotated cell 
                               populations in the heatmaps. This effect can be skipped checking the related checkbox (no Quantiles handling)."),
                     p("The other two parameters can be very sensitive as they constrained the normalized MFI which can be 
                               greater/less than the ", strong("MFI positive/negative threshold. "), "By default a phenotype is considered “+” (aka 
                               positively expressed) for certain marker if the normalized MFI is greater than 0.5 and considered “-” (aka negatively 
                               expressed) if the normalized MFI is inferior than 0.5. With a more stringent threshold (e.g: 0.6 and 0.4) you will 
                               identify less cluster (or meta-cluster) as belonging to a certain phenotype: you will have more “undetermined” type of 
                               cells, but the ones which can be assigned are more likely to be truly positive (or truly negative)"),
                     
                     h4("Push the button to run the selected clustering algorithm"),br(),
                     actionButton("runClust", label = "Run the selected clustering algorithm", icon("paper-plane"), 
                                  style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),br(),
                     h3(textOutput("checktxt_clust")), br(),   
                     conditionalPanel(
                         condition = ("output.clust_panelStatus"),
                         h4("Push the button to plot clustering result"),br(),
                         actionButton("plotClust", label = "Show the results 4 the selected clustering algorithm"),br(),
                         conditionalPanel(
                             condition = ("output.plotClust_panelStatus"),
                             conditionalPanel(
                                 condition = "input.ClustAlgo == 'flowSOM'",
                             imageOutput(outputId = "markerPlot", width = "1000px", height = "2400px") %>% withSpinner(),
                             tags$hr(style="border-color: green;"), br(),
                             imageOutput(outputId = "compPlot", width = "1000px", height = "1200px") %>% withSpinner(), br(),
                             downloadButton(outputId = "downloadClust",label = "download the plots"),br(),
                             
                             h4(strong("Additional notes on the flowSOM results:")),
                             p("Sofie Van Gassen, the author of the flowSOM package, implemented an additional feature along with the SOM's: the
                             organization of the resulting clustering in a ", strong("Minimal Spanning Tree (MST)"), ". A MST is an acyclic graph
                             which connects the nodes of a graph in such a way that the sum of the weights of the branches is minimal. By doing 
                             this, nodes will get connected to the ones they are the most similar. This feature put 'flowSOM' (like the 'spade' 
                             package) between those trajectory inference algorithms which aim to automatically reconstruct the developmental path 
                             of the various cell phenotypes. For this see 
                               (", tags$a(href="https://onlinelibrary.wiley.com/doi/pdf/10.1002/eji.201646347", "European Journal of Immunology - 
                               REVIEW: Computational methods for trajectory inference from single-cell transcriptomics - Robrecht Cannoodt, 
                                          Wouter Saelens and Yvan Saeys."), ") "), 
                             p("It is important to interpret the flowSOM's MST, but not just the marker's related content but even its shape, to
                             infer the algorithm outcomes. In general, a more twisted and inconsistent (namely 'a not so acyclic') graph tree is 
                             bad signal for your analysis: a good clustering result is a linear MST in which each branch can groups several clusters 
                               in a meta-cluster. Also notice the heatmap's dendograms and their relationship with the MST's branches")),br(),
                             
                             conditionalPanel(
                                 condition = "input.ClustAlgo == 'Rphenograph'||input.ClustAlgo == 'FastPG'",
                                 tags$hr(style="border-color: green;"), br(),
                                 htmlOutput("phenograph_text")),
                             conditionalPanel(
                               condition = "input.ClustAlgo == 'DepecheR'",
                               tags$hr(style="border-color: green;"), br(),
                               htmlOutput("depeche_text")),
                             conditionalPanel(
                                 condition = "input.ClustAlgo == 'KMeans'",
                                 tags$hr(style="border-color: green;"), br(),
                                 htmlOutput("kmeans_text"))), #plotClust_panelStatus
                             actionButton("saveClust", label = "Produce a concatenated sample with the additional 'clusterId' dimension"),br(),br(),
                             conditionalPanel(
                                 condition = ("output.clust_save_panelStatus"),
                                 downloadButton(outputId = "downloadFCSClust",label = "download the concatenated sample with cluster dimension")),
                             
                             conditionalPanel(
                                 condition = "output.clust_eva_panelStatus",
                                 tags$hr(style="border-color: green;"), br(),
                                 p("The following is the outcome of the clustering performance evaluation method. It is called Silhouette and it is 
                                 a method of interpretation and validation of the consistency of the clustering data. Please, do not be afraid if 
                                 the report dos not show a truly amazing performance quality, with a lot of events negatively evaluated inside 
                                 most of the clusters. Clustering is unsupervised and most of your experiments are complex to be uniquely clustered: 
                                 it's hard to know a priori which the best clustering is. In any case the method can provide a good test for 
                                 comparison of different clustering performances, in order to better select parameters, algorithms, 
                                   cluster amount... "),br(),
                                 plotOutput(outputId = "sil_plot", width = "1100px", height = "1000px") %>% withSpinner(),
                                 tableOutput(outputId = "sil_summary"),
                                 downloadButton(outputId = "downloadPerf_csv",label = "download the silhouette clustering table details")
                                 )#Clustering_Perf
                             
                         )#plotClust_panelStatus
                     )#tab panel


######################## Label I TabPanel -----------------------------

tabLabelClust <- tabPanel("Labeling I", icon =icon("blender"), value = "Clustering_p",
                     p("The clustering algorithm grouped the various events in different homogeneous groups. Each of this group, in turn, expresses a 
                     median value (for each of the selected markers) to form an expression profile. Each of this this expression profile could match 
                       the phenotype entered in the 'meta-data' related workflow."),
                     p("Of course, not every cluster belongs to a certain entered phenotype, but the ones having the median expression of
                       the corresponding phenotype description could be assigned ('labeled') as they belongs to a certain cell tipology"),
                     h4("Warning: Don't forget to edit the phenotype related meta-data in the 'Meta-data and Assays' workflow if you want to use 
                        this workflow step"),br(),
                     h4("Push the button to label the clustered cells"),
                     
                     actionButton("labelClust", label = "Perform the labeling algorithm", icon("paper-plane"), 
                                  style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),br(),
                     h3(textOutput("checktxt_labelClust")), br(),    
                     conditionalPanel(
                         condition = ("output.labelClust_panelStatus"),
                         
                     h4("At this point, the cell events are grouped and labeled according to the entered metadata. Every single cluster has a global  
                    expression which is calculated on the median of each single cell expression grouped in the selected cluster. If this median 
                    values matches with the profile of the phenotype metadata table, then the events which belong to that cluster expressed that 
                        particular phenotype."), br(),
                     
                     p("With this, it is possible to plot the heatmap of all the cell events as they are grouped in distinct clusters by the selected 
                       clustering algorithm. You could have less than the entered number of clusters (flowSOM for example tend to reduce the number 
                       of effective clusters"),
                     
                     h4("Push the button to show the heatmap plot 4 the selected clustering algorithm"),
                     actionButton("plot_labelClust", label = "Show the clustering heatmap"),br(),
                     conditionalPanel(
                         condition = ("output.plot_labelClust_panelStatus"),
                         imageOutput(outputId = "HM_full_label", width = "1200px", height = "3000px") %>% withSpinner(), br(),br(),br(),
                         downloadButton(outputId = "downloadClust_p",label = "download the plots"),br(),br(),
                         downloadButton(outputId = "downloadClust_csv",label = "download the .csv table of the heatmap"),br(),br(),
                         downloadButton(outputId = "downloadClustPheno_csv",label = "download the .csv pheno table"),br(),br(),
                         h4("Additional notes on the heatmap:"),
                         p("Each single row of the heatmap represent a cluster, namely a group of events with the same type of expression 
                           calculated on the median of the expression value for each single cell belonging to that cluster."),
                         p("This could help to visually evaluate the 'real' amount of clusters: in fact, if you have several cluster with 
                           the same profile (see the heatmap colors), probably you can 'meta-cluster', that is to group on the same cluster, more 
                           than a single cluster. This could happen for eaxmple if you select a number of clusters bigger than the effective 
                           ones"),
                         p("Another view of the labeled results can be seen within the 'labeling II' tab. In case of meta-cluster execution, the 
                           labeling is performed grouping the clusters in meta-clusters, while skipping the meta-clustering, the labeling is 
                           executed just evaluating the median expression of the various population."),
                         downloadButton(outputId = "downloadClustAll",label = "download all the produced analysis data"),br(),br()
                         )))


######################## MetaClust TabPanel -----------------------------

tabMetaClust <- tabPanel("Meta-Clustering", icon =icon("object-group"), 
                         p("This is the section related to the meta-clustering. This is executed by the ", strong("ConsensusClusterPlus")," 
                           package, an implementation of a kind of 'standard' algorithm, which is capable to group the various clusters obtained by 
                           the selected clustering algorithm. Producing these ", shiny::em("clusters of clusters"), " can be very much useful in 
                           case of a big amount of clusters, which is a typical situation with, for example, 'flowSOM' which by default works with 
                           a ", shiny::em("clusterization"), " of 100 clusters"),
                         p("Please, notice than, the values in the heatmaps showing the marker expression of the 
                         various clusters are very much sensitive to their expression scales: please select a properly transformed samples to 
                         obtain accurate values"),
                         h4("Push the button to run the meta-clustering algorithm"),
                         actionButton("runMetaClust", "Run the selected MetaClustering algorithm", icon("paper-plane"), 
                                      style="color: #fff; background-color: #ff0000; border-color: #2e6da4"), br(),
                         h3(textOutput("checktxt_metaClust")), br(),    
                         conditionalPanel(
                             condition = ("output.metaClust_panelStatus"),
                             h4("Push the button to plot the meta-clustering result"),
                             actionButton('plotMeta', 'Show the MetaClustering related plot'), br(),
                             conditionalPanel(
                                 condition = ("output.plotMeta_panelStatus"),
                                 imageOutput(outputId = "metaHeatmap_cluster", width = "1200px", height = "1600px") %>% withSpinner(),br(),
                                 downloadButton(outputId = "downloadMeta",label = "download the plots"),br(),
                                 h4("Additional notes on the heatmap:"),
                                 p("The meta-clustering procedure ,", strong("could always be skipped."), " It could be not important to reduce the 
                                 number of clusters, especially using clustering algorithm or flowSet which does not require so many clusters."),
                                 p("In any case, it should be noticed than the meta-clustering procedure works on the clusters themself and not on 
                                 group of events which are labeled together just because they present a similar phenotype like it is shown within 
                                 the 'Labeling I' tab. Using the ", strong("ConsensusClusterPlus "), "algorithm, means to group together clusters 
                                 which are closer inside the hyper-space described by the expression matrix of the working samples, while with the 
                                 criteria used within the 'Label I' tab, we label events which are just similar but do not necessarely lay together 
                                 in the hyperspace of events. Therefore the mate-clustering results are more consistent with the clustering 
                                   approach"),br(),
                                 downloadButton(outputId = "downloadMeta_csv", label = "download the csv table of the meta-cluster heatmap"),
                                 br(),br(),
                                 downloadButton(outputId = "downloadPhenoMeta_csv", label = "download the csv of the pheno-table"),
                                 br(),br(),
                                 downloadButton(outputId = "downloadConsensus_csv", 
                                                label = "download the csv files of the consensusClusterPlus package"),
                                 br(),br(),
                                 downloadButton(outputId = "downloadConsensus_png", 
                                                label = "download the plots of the consensusClusterPlus package"),
                                 br(),br(),
                                 actionButton("saveMetaClust", 
                                              label = "Produce a concatenated sample with the additional meta-cluster dimension"),br(),
                                 h4("Notice that, in case you handled the original flowSet within the 'flowSet optimization workflow', the 
                                 produced flowFrame (a unique '.fcs' file) will collect, along with the new computed dimensions, also the 
                                    original channels"),br(),#),br(),
                                 downloadButton(outputId = "downloadMetaAll", 
                                                label = "download all the produced analysis data"),br(),
                                 
                                 h3(textOutput("checktxt_saveMetaClust")), br(),    
                                    conditionalPanel(
                                    condition = ("output.saveMetaClust_panelStatus"),
                                    downloadButton(outputId = "downloadFCSMeta_Clust",
                                                label = "download the concatenated sample with cluster dimension"),br(),br()),
                             conditionalPanel(
                                 condition = ("output.metaClust_eva_panelStatus"),
                                 tags$hr(style="border-color: green;"), br(),
                                 p("The following is the outcome of the meta-clustering performance evaluation method"),br(),
                                 plotOutput(outputId = "meta_sil_plot", width = "1100px", height = "1000px") %>% withSpinner(),
                                 tableOutput(outputId = "meta_sil_summary"),
                                 downloadButton(outputId = "downloadMetaPerf_csv",label = "download the silhouette clustering table details"))
                             )))
                        #??? show the traking plots of the ConsensusClusterPlus

######################## MapComp TabPanel -----------------------------
tabMapComp <- tabPanel("Mapping", icon =icon("atlas"), 
                        p("How to show in an unique plot, all the cell events separated as effectively and efficiently as possible? You sould notice
                          than the original cell event space is an hyperspace of ", shiny::em("M"), " dimensions where ", shiny::em("M"), " is the 
                          amount of selected stains ('aka' markers) used for the analysis: with more than just three dimensions there is no way to 
                          represent such a complex plot. For this reason the dimensionality reduction algorithms help us out in this aim: the all 
                          amount of the events are coated in a two dimensions plot (but also the three dimensions versions are available). This is 
                          possible with the classical PCA (Principal Component Analysis) or the Isomap algorithms, but with complex non-linear data
                          like the flow cytometry expression matrix, ", strong("tSNE"), " (t-distributed Stochastic Neighbor Embedding) and ", 
                          strong("UMAP"), " (Uniform Manifold Approximation and Projection), with their three dimensions variants, are the most 
                          efficient dimensionality reduction algoritms suitable for our type of dataset"),
                        p("The biggest advantage of these algorithms is their ability to keep the various events separated and to group them in 
                          homogeneous subgroups. In this sense, tSNE and UMAP are kind of clustering algorithms themself, but their output is only 
                          a two (or three) dimensions map for the projection of the whole bunch of events: they do not assign every single event 
                          to a specific cluster, but they merely distribute each event in their two (or three) axes space keeping them separated in
                          homogenous groups. In particular, tSNE provides a consistent map with the various homogenous event grouped together, but 
                          these groups are distributed quite randomly while UMAP (if properly configured) should converge in a map with the much 
                          similar groups kept closer together."),
                        p("There are tons of web pages and articles about these maps and their variations (here we are using  a Barnes-Hut 
                          implementation of t-SNE) but it worth mentioning: ", tags$a(href="https://lvdmaaten.github.io/tsne/", "tSNE - Laurens van 
                          der Maaten"), " ,  ", tags$a(href="https://distill.pub/2016/misread-tsne/", "How to Use t-SNE Effectively"), " and ",
                          tags$a(href="https://lvdmaaten.github.io/tsne/", "Paper Dissected: “Visualizing Data using t-SNE” Explained"), "."),br(),
                        h4("Push the button to perform the selected dimensionality reduction algorithm"),
                        actionButton("runMapsComp", "Run the selected dimensionality reduction", icon("paper-plane"), 
                                     style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                        h3(textOutput("checktxt_runMapsComp")), br(),    
                        conditionalPanel(
                            condition = ("output.runMapsComp_panelStatus"),
                            h4("Push the button to plot the selected map"),
                            actionButton('plotmapcomp', "Show the maps and the cluster's related quantity plot"),br(),
                            conditionalPanel(
                                condition = ("output.runMapShow_panelStatus"),
                                plotlyOutput(outputId = "showmapcomp", width = "1100px", height = "1000px") %>% withSpinner(), br(),
                                
                                p("The plot below is the ratio of each sample for each single cluster"),
                                plotOutput(outputId = "showJtSample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                downloadButton(outputId = "downloadJtSample_csv",label = "download the .csv of this plot"),br(),br(),
                                tags$hr(style="border-color: green;"),br(),
                                
                                conditionalPanel(
                                    condition = "output.condition_meta",
                                    p("The plot below is the ratio of each sample for each single meta_cluster"),
                                    plotOutput(outputId = "showJtMetaSample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                    downloadButton(outputId = "downloadJtMetaSample_csv",label = "download the .csv of this plot"),br(),br(),
                                    tags$hr(style="border-color: green;")),br(),
                       
                                conditionalPanel(
                                    condition = "output.conditionJt_tag1",
                                    plotOutput(outputId = "showJt1_Sample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                    downloadButton(outputId = "downloadSampleTag1_csv",label = "download the .csv of this plot"),br(),br(),
                                    tags$hr(style="border-color: green;"),br(),
                                    conditionalPanel(
                                        condition = "output.condition_meta",
                                        plotOutput(outputId = "showJt1_MetaSample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                        downloadButton(outputId = "downloadMetaSampleTag1_csv",label = "download the .csv of this plot"),br(),br(),
                                        tags$hr(style="border-color: green;")),br()),
                                
                                conditionalPanel(
                                    condition = "output.conditionJt_tag2",
                                    plotOutput(outputId = "showJt2_Sample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                    downloadButton(outputId = "downloadSampleTag2_csv",label = "download the .csv of this plot"),br(),br(),
                                    tags$hr(style="border-color: green;"),br(),
                                    conditionalPanel(
                                        condition = "output.condition_meta",
                                        plotOutput(outputId = "showJt2_MetaSample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                        downloadButton(outputId = "downloadMetaSampleTag2_csv",label = "download the .csv of this plot"),br(),br(),
                                        tags$hr(style="border-color: green;")),br()),
                                
                                conditionalPanel(
                                    condition = "output.conditionJt_tag3",
                                    plotOutput(outputId = "showJt3_Sample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                    downloadButton(outputId = "downloadSampleTag3_csv",label = "download the .csv of this plot"),br(),br(),
                                    tags$hr(style="border-color: green;"),br(),
                                    conditionalPanel(
                                        condition = "output.condition_meta",
                                        plotOutput(outputId = "showJt3_MetaSample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                        downloadButton(outputId = "downloadMetaSampleTag3_csv",label = "download the .csv of this plot"),br(),br(),
                                        tags$hr(style="border-color: green;")),br()),
                                
                                conditionalPanel(
                                    condition = "output.conditionJt_tag4",
                                    plotOutput(outputId = "showJt4_Sample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                    downloadButton(outputId = "downloadSampleTag4_csv",label = "download the .csv of this plot"),br(),br(),
                                    tags$hr(style="border-color: green;"),br(),
                                    conditionalPanel(
                                        condition = "output.condition_meta",
                                        plotOutput(outputId = "showJt4_MetaSample", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                                        downloadButton(outputId = "downloadMetaSampleTag4_csv",label = "download the .csv of this plot"),br(),br(),
                                        tags$hr(style="border-color: green;")),br()),
                                
                                actionButton("saveMap", "Produce a concatenated flowSet with all the additional dimensions"),br(),
                                h4("With this you can download a concatenated flowframe with all the additional dimension used for the analysis: the 
                                clusterId, the meta-clusterId and all the map dimensions"),br(),br(),
                                actionButton("plot_marker_map", "Plot all the markers distributions within the map"),
                                h4("You can plot the various marker intensities within the map"),br(),br(),
                                conditionalPanel(
                                    condition = ("output.runMapMarker_panelStatus"),
                                    plotOutput(outputId = "showmarkermap", width = "auto", height = "auto") %>% withSpinner(),br(),
                                ),
                                conditionalPanel(
                                    condition = ("output.saveMap_panelStatus"),
                                    downloadButton(outputId = "downloadFCSMap", label = "download the concatenated sample with all the additional 
                                                   dimensions")))))
                        

######################## Label II TabPanel -----------------------------
tabLabel <- tabPanel("Labeling II", icon =icon("blender"),  value = "labeling",
                     h4("Push the button to perform the labeling algorithm"),br(),
                     p("remember: for the labeling operation, please fill the related phenotype table in the metadata menu"),
                     actionButton("runLabel", "Run the labeling algorithm", icon("paper-plane"), 
                                  style="color: #fff; background-color: #ff0000; border-color: #2e6da4"), br(),
                     h3(textOutput("checktxt_runLabel")), br(),    
                     conditionalPanel(
                         condition = ("output.label_panelStatus"),
                         h4("Push the button to plot the labeling results"),br(),
                         actionButton('runLabelPlot', 'Show the heatmaps'),
                         h3(textOutput("checktxt_runLabelPlot")), br(),    
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatus"),
                             h4("The following plots show the heatmaps with the labeling when the meta-clustering is performed"),br()),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatus"),
                             imageOutput(outputId = "metaHeatmap_label",  width = "1200px", height = "1200px") %>% withSpinner()),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatus"),
                             downloadButton(outputId = "downloadLabel_csv",label = paste0(("download the csv table of the heatmaps "),
                                            ("related to the phenotypes calculated with clustering")))),br(),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatus"),
                             imageOutput(outputId = "metaHeatmap_pheno",  width = "1200px", height = "900px") %>% withSpinner()),br(),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatus"),
                             downloadButton(outputId = "downloadPheno_csv",label = paste0(("download the csv table of the heatmaps "),
                                                                    ("related to the phenotypes calculated with clustering")))),br(),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatus"),
                             downloadButton(outputId = "downloadLabel",label = "download the plots")),
                         
                         tags$hr(style="border-color: green;"), 
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatusPlus"),
                             h4("The following plots show the heatmaps with the labeling without the meta-clustering process"),br()),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatusPlus"),
                             imageOutput(outputId = "metaHeatmap_labelPlus",  width = "1300px", height = "3000px") %>% withSpinner()),br(),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatusPlus"),
                             downloadButton(outputId = "downloadLabel_csvPlus",label = "download the csv table of the heatmap")),br(),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatusPlus"),
                             imageOutput(outputId = "metaHeatmap_phenoPlus",  width = "900px", height = "900px") %>% withSpinner()),br(),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatusPlus"),
                             downloadButton(outputId = "downloadPheno_csvPlus",label = "download the csv table of the heatmap")), br(),
                         conditionalPanel(
                             condition = ("output.labelPlot_panelStatusPlus"),
                             downloadButton(outputId = "downloadLabelPlus",label = "download the plots")),br(),
                         ))

######################## Pheno maps I TabPanel -----------------------------
tabLabel_p <- tabPanel("Pheno maps I",  icon =icon("search-plus"), value = "Label_p",
                     h4("Push the button to show labeling plots over the selected map"),br(),
                     p("The following plots are related to the meta_data tables"), br(),
                     actionButton('runLabelPlot_p', 'Show the plots'), br(),
                     h3(textOutput("checktxt_runPhenomapsI")), br(),    
                     conditionalPanel(
                         condition = ("output.labelPlotp_panelStatus"),
                         h4("The following plots show the maps with the labeling operations over the selected map. The latter take into account 
                            samplesIds"),br(),
                         plotlyOutput(outputId = "showLabeldMap", width = "1200px", height = "1200px") %>% withSpinner(),
                         tags$hr(style="border-color: green;"),br(),
                         conditionalPanel(
                             condition = "output.condition_tag",
                            plotlyOutput(outputId = "showSampleMap", width = "100%", height = "1000px") %>% withSpinner(),
                            tags$hr(style="border-color: green;"))))

######################## Pheno maps II TabPanel -----------------------------
tabLabel_pp <- tabPanel("Pheno maps II",  icon =icon("search-plus"), value = "Label_pp",
                       h4("Push the button to show other labeling plots over the selected map"),br(),
                       p("The following plots are related to the meta_data tables"), br(),
                       actionButton('runLabelPlot_pp', 'Show the plots'), br(),
                       h3(textOutput("checktxt_runPhenomapsII")), br(),    
                       conditionalPanel(
                           condition = ("output.labelPlotpp_panelStatus"),
                           h4("The following plots show the maps with the labeling operations phocusing on the conditions entered as 
                              tag1, tag2, tag3 and tag4 (if any)"),br(),
                           
                           conditionalPanel(
                               condition = "output.condition_tag1",
                               plotlyOutput(outputId = "showTag1Map", width = "100%", height = "1000px") %>% withSpinner(), br()),
                           conditionalPanel(
                               condition = "output.condition_tag2",
                               plotlyOutput(outputId = "showTag2Map", width = "100%", height = "1000px") %>% withSpinner(), br()),
                           conditionalPanel(
                               condition = "output.condition_tag3",
                               plotlyOutput(outputId = "showTag3Map", width = "100%", height = "1000px") %>% withSpinner(), br()),
                           conditionalPanel(
                               condition = "output.condition_tag4",
                               plotlyOutput(outputId = "showTag4Map", width = "100%", height = "1000px") %>% withSpinner(), br())))

######################## Quantity TabPanel -----------------------------
tabQuant <- tabPanel("Quantity plots", icon =icon("chart-bar"), 
                     h4("Push the button to display some quantity plots"),
                     p("The following plots are related to the meta_data tables"),br(),
                     actionButton('runQuantPlot', 'Show the quantity plots'),br(),
                     h3(textOutput("checktxt_runQuantPlot")), br(),    
                     conditionalPanel(
                         condition = ("output.QuantPlot_panelStatus"),
                         
                         plotlyOutput(outputId = "showBarSample", width = "100%", height = "1000px")  %>% withSpinner() ,br(),
                         plotlyOutput(outputId = "showDodgeSample", width = "100%", height = "1000px")  %>% withSpinner(),
                         downloadButton(outputId = "downloadQuantSample_csv",label = "download the .csv of this plot"),br(),br(),
                         tags$hr(style="border-color: green;"),br(),
                         h4("The barplots show the quantity of the labeled samples, taking into account the metadata elements sample_id and the conditions
                         edited in the meta_data table as tag1, tag2 and tag3 (if any)"),
                         
                         conditionalPanel(
                             condition = "output.condition_tag1",
                             plotlyOutput(outputId = "showBarTag1", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                             plotlyOutput(outputId = "showDodgeTag1", width = "100%", height = "1000px")  %>% withSpinner(),
                             downloadButton(outputId = "downloadQuantTag1_csv",label = "download the .csv of this plot"),br(),br(),
                             tags$hr(style="border-color: green;"),br()),
                         conditionalPanel(
                             condition = "output.condition_tag2",
                             plotlyOutput(outputId = "showBarTag2", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                             plotlyOutput(outputId = "showDodgeTag2", width = "100%", height = "1000px")  %>% withSpinner(),
                             downloadButton(outputId = "downloadQuantTag2_csv",label = "download the .csv of this plot"),br(),br(),
                             tags$hr(style="border-color: green;"),br()),
                         conditionalPanel(
                             condition = "output.condition_tag3",
                             plotlyOutput(outputId = "showBarTag3", width = "100%", height = "1000px")  %>% withSpinner(), br(),
                             plotlyOutput(outputId = "showDodgeTag3", width = "100%", height = "1000px")  %>% withSpinner(),
                             downloadButton(outputId = "downloadQuantTag3_csv",label = "download the .csv of this plot"),br(),br(),
                             tags$hr(style="border-color: green;"),br()),
                         
                         h4("The boxplots show, with jittered points, of the sample-level cluster proportions for each tag (if any)"),
                         conditionalPanel(
                             condition = "output.condition_tag1",
                             plotlyOutput(outputId = "med_expr_tag1", width = "100%", height = "1000px")  %>% withSpinner(),
                             tags$hr(style="border-color: green;"),br()),
                         conditionalPanel(
                             condition = "output.condition_tag2",
                             plotlyOutput(outputId = "med_expr_tag2", width = "100%", height = "1000px")  %>% withSpinner(),
                             tags$hr(style="border-color: green;"),br()), 
                         conditionalPanel(
                             condition = "output.condition_tag3",
                             plotlyOutput(outputId = "med_expr_tag3", width = "100%", height = "1000px")  %>% withSpinner(),
                             tags$hr(style="border-color: green;"),br())))

######################## Time-step TabPanel -----------------------------
tabQuant_p <- tabPanel("Time-step plots", icon =icon("shoe-prints"), value = "Quant_p",
                     h4("Push the button to disply the quantity plots"),
                     p("The following plots are related to the meta_data tables. If the time_step are correctly entered, you can see the barplot and 
                       the streamPlot as the evolution of your data_set"),br(),
                     actionButton('runQuantPlot_p', 'Show the time-step plots'),br(),
                     h3(textOutput("checktxt_runQuantPlot_p")), br(),    
                     conditionalPanel(
                         condition = ("output.QuantPlotp_panelStatus"),
                         h4("The barplots show the quantity of the labeled samples, taking into account the time_step metadata element (if any)"),br(),
                         plotlyOutput(outputId = "showBarTag4", width = "100%", height = "800px")  %>% withSpinner(), br(),
                         plotlyOutput(outputId = "showDodgeTag4", width = "100%", height = "800px")  %>% withSpinner(),
                         downloadButton(outputId = "downloadQuantTag4_csv",label = "download the .csv of this plot"),br(),br(),
                         tags$hr(style="border-color: green;"),br(),
                         streamgraphOutput(outputId = "stream_plot_prop", width = "100%", height = "800px")  %>% withSpinner(),
                         tags$hr(style="border-color: green;"),br(),
                         streamgraphOutput(outputId = "stream_plot_count", width = "100%", height = "800px")  %>% withSpinner(),
                         tags$hr(style="border-color: green;"),br(),
                         plotlyOutput(outputId = "med_expr_tag4", width = "100%", height = "800px")  %>% withSpinner()))

######################### sidebar panel flowSet optimization-----------------------------------------------------
sidebar_panel_pre <- sidebarPanel(width = 3,br(),
                                  tags$head(
                                      tags$style(HTML("hr {border-top: 1px solid #008000;}"))
                                  ),
                                  h4("console messages"),
                                  verbatimTextOutput("console_output_pre", placeholder = T),
                                  bsPopover(id = "console_output_pre", title = "console messages", 
                                            paste0('This collects the report from the various processes executed along the several steps ',
                                                           'of the selected workflow. Before to move along with the next procedure, please wait ', 
                                                           'until you get some feedback from the previous commands (namely, wait before pushing ',
                                                           'the next button!)')),
                                  
                                  tags$hr(style="border-color: green;"),
                                  
                                  # Input: Select the cleaning parameters ---
                                  numericInput(inputId = "set_seed", label = "set.seed() function", value = 0608, min = 1, max = 10000),
                                  bsPopover(id = "set_seed", title = "set seed", content = paste0('This number fixes the seed for the pseudo-random number generation, which is ', 
                                                           'important for most of the procedures within cytoChain. Based on this you can ', 
                                                           'also choose the color combination of your samples in the first "Loading and parsing" ',
                                                           ' tab')),
                                  tags$hr(style="border-color: green;"),
                                  h4("cleaning parameters"),
                                  
                                  # Input: Select the cleaning parameters ---
                                  sliderInput(inputId = "clean_param1", "cleaning: alphaFR",
                                              value = 5, min = 0.1, max = 50.0, step = 0.01),
                                  bsPopover(id = "clean_param1", title = "cleaning parameters alphaFR", 
                                            content = paste0('The level of statistical significance used to accept anomalies detected by the ESD',
                                                           '(Extreme Studentized Deviate) method used by "flowAI". It affects the flow rate check: ', 
                                                           'the lower, the lower is the impact')),
                                  
                                  sliderInput(inputId = "clean_param2",
                                              "cleaning: pen_valueFS",
                                              value = 40, min = 1, max = 200, step = 1),# it was 80
                                  bsPopover(id = "clean_param2", title = "cleaning parameters pen_valueFS",
                                            content = paste0('The value of the penalty for the changepoint detection algorithm for the signal ', 
                                            'acquisition check. The higher is the value, the less strict is the detection of the ',
                                            'anomalies. It can be very much sensitive, especially for small (< 10000 events) samples')),
                                  
                                  tags$hr(style="border-color: green;"),
                                  h4("scaling parameters"),
                                  # Input: Select the type of transformation---
                                  radioGroupButtons(
                                      inputId = "trans", label = "Choose your transformation function:", 
                                      choices = c("Arcsinh" = "arcsinh",
                                                  "flowVS" = "flowVS",
                                                  "Generalized Arcsinh" = "mclMultivArcSinh", 
                                                  "Biexponential" = "mclMultivBiexp", 
                                                  "Generalized Box-Cox" = "mclMultivBoxCox", 
                                                  "LinLog"= "mclMultivLinLog", 
                                                  "Arcsinh+" = "arcsinhplus", 
                                                  "Scale" = "scale", 
                                                  "Zero" = "zero",  
                                                  "None" = "none"), 
                                      justified = TRUE, status = "primary", direction = "vertical", selected = "arcsinh", 
                                      checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                                  bsPopover(id = "trans", title = "transformation functions", content = paste0("Select your favourite scaling criteria")),
                                  
                                  # Input: Slider for the transformation parameter ---
                                  sliderInput(inputId = "trans_param_arcsinh",label = "arcsinh parameter (cofactor)",
                                              value = 150, min = 1, max = 200),
                                  bsPopover(id = "trans_param_arcsinh", title = "arcsinh cofactor", 
                                            content = paste0('The slope of the linear zone within the arcsinh transformation. The default is the recommended', 
                                    'value mentioned in bunches of papers')),
                                  
                                  sliderInput(inputId = "trans_param_min", label = "zero parameter: min_quantile",
                                              value = 0.05, min = 0, max = 0.5),
                                  bsPopover(id = "trans_param_min", title = "trans parameter min",
                                    content = paste0('for the "Zero" transformation only. It is the lower percentile threshold to set the expression', 
                                    'value to 0 inside the expression matrix')),
                                  
                                  sliderInput(inputId = "trans_param_max", label = "zero parameter: max_quantile",
                                              value = 0.95, min = 0.5, max = 1),
                                  bsPopover(id = "trans_param_max", title = "trans parameter max", 
                                            content = paste0('for the "Zero" transformation only. It is the upper percentile threshold to set the expression', 
                                    'value to 1 inside the expression matrix')),
                                  
                                  tags$hr(style="border-color: green;"),
                                  h4("downsampling parameters"),
                                  
                                  radioGroupButtons(
                                      inputId = "down_algo", label = "Choose your down-sampling algorithm:", 
                                      choices = c("spade", "RKOF", "KNN_AGG", "KNN_SUM", "LOF", "LOOP", "LDF", "LDE", "Random"), 
                                      justified = TRUE, status = "primary", direction = "vertical", selected = "spade", 
                                      checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                                  bsPopover(id = "down_algo", title = "down-sampling algorithm",
                                            content = paste0('Select your preferred algorithm. The available selection is a subset of the algorithm ', 
                                                           'collected in the "DDoutlier" R package along with "spade". See the related manuals ',
                                                           'Please keep in mind that "spade" computes densities: smaller densities-> outlierness, ',
                                                           'while the other algorithms computes scores: higher score -> the outlierness')),
                                  
                                  radioGroupButtons(
                                      inputId = "down_type", label = "Type of down-sample filtering", 
                                      choices = c("by percentage" = "down_perc", "by value" = "down_val", "equalize" = "equalize"),
                                      justified = TRUE, status = "primary", direction = "horizontal", selected = "down_perc", 
                                      checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                                  bsPopover(id = "down_type", content = "down-sample filtering",
                                            title = paste0('You can choose if you want to filter out events selecting them by percentage, by value ',
                                                           'or equalize by a specific number')),
                                  
                                  sliderInput("down_percentile",
                                              "downsampling parameters: exclude_pctile",
                                              value = 3, min = 0.1, max = 90.0, step = 0.1 , post = " %"), 
                                  bsPopover(id = "down_percentile", title = "exclude percentile", 
                                            content = paste0('This parameter affects the downsampling quantities for every available algorithms.', 
                                                           'Densities (for "spade") , or Scores (for all the other algorithms) are properly ',
                                                           'ordered: a percentage of events below this % will be excluded, choosing among ', 
                                                           'the lower densities (in case of "spade") or the upper scores (for the rest)', 
                                                           '. E.g. a setting of 5% means that 5% of events will be left out from the analysis). ',
                                                           'Please notice that you filter the same amount of event from every single sample: ', 
                                                           'pay attention when you deal with very different sample weight and consider to evaluate ', 
                                                           'the downsampling operation, sample by sample')),
                                  
                                  sliderInput("down_value",
                                              "downsampling parameters: score|density value",
                                              value = 60,min = 2, max = 98, step = 0.1),
                                  bsPopover(id = "down_value", title = "score/density value", 
                                            content = paste0('It is possible to choose a density|score value (density for spade and score for every ', 
                                                           'other algorithm) as a threshold to filter out all the evens with score below/above this ',
                                                           'set value. The value itself has been normilized from 0 to 100 for every density/score ',
                                                           'calculated by the various down-sampling algorithms. In this case you can cut a different ', 
                                                           'amount of event for each of the sample. The drawback is that you cannot infer the number ', 
                                                           'of event you filter')),
                                  
                                  numericInput(inputId = "down_number", 
                                               label = "downsampling parameter: common equalize value", value = 10000, min = 100, max = 200000, step = 5),
                                  bsPopover(id = "down_number", title = "equalizing value", 
                                            content = paste0('It is possible to equalize all of your samples to a specific number of events')),
                            
                                  sliderInput("down_kappa",
                                              "downsampling parameter: K",
                                              value = 10,min = 0, max = 100, step = 1),
                                  bsPopover(id = "down_kappa", title = "K downsampling parameter", 
                                            content = paste0('This is a parameter related to every single algorithm except "spade". It is the K_min ',
                                            'for the "KNN_AGG", while it is simply the "K" for the rest of the algorithms')),
                                  
                                  sliderInput("down_kappa_max",
                                              "downsampling parameters: K_max",
                                              value = 30,min = 0, max = 200, step = 1),
                                  bsPopover(id = "down_kappa_max", title = "K_max downsampling parameter", 
                                            content = paste0('This is the "K_max" for the "KNN_AGG" algorithm')),
                                  
                                  sliderInput("down_par_lambda",
                                              "LOOP parameter lambda",
                                              value = 3,min = 0, max = 10, step = 1),
                                  bsPopover(id = "down_par_lambda", title = "LOOP downsampling parameter",
                                            content = paste0('This is the "lambda" for the LOOP algorithm')),
                                  
                                  sliderInput("down_par",
                                              "downsampling parameter",
                                              value = 1,min = 0, max = 10, step = 1),
                                  bsPopover(id = "down_par", title = "alpha downsampling parameter",
                                            content = paste0('This is the "alpha" for the RKOF algorithm. It is the "h" for the "LDE" and the "LDF" ',
                                                           'algorithms')))            

######################### sidebar panel meta-data --------------------------------------------
sidebar_panel_meta <- sidebarPanel(width = 3,
                                   
                                   h4("console messages"),
                                   verbatimTextOutput("console_output_meta", placeholder = T),
                                   
                                   tags$hr(style="border-color: green;"),
                                   
                                   p("Select your '.fcs' samples. The purpose of this file loding widgets is to let you load your samples to analize 
                                   without passing from the flowSet optimization workflow. For this reason these files have priority. Skip this 
                                     if you want to keep using the ones handled in the flowSet optimization workflow"),
                                   fileInput(inputId = "fSmeta_in", label = h3("Load your samples as collection of FCS files"), multiple = T,
                                             buttonLabel = "FCS files", placeholder = "upload your flowFrames"),
                                   h4("or..."),
                                   fileInput(inputId = "zipmeta_in", label = h3("Load your samples collected in a single zip file"), multiple = F,
                                             buttonLabel = "zip file", placeholder = "upload your flowFrames", 
                                             accept = c("application/zip", "ZIP archive", ".zip")),
                                   
                                   tags$hr(style="border-color: green;"),
                                   
                                   fileInput(inputId = "meta_file", label = h4("meta-sample table"), multiple = F,
                                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"),
                                             buttonLabel = "load your .csv file", placeholder = "upload your meta-data table"),
                                   
                                   h4("time/step entries in the meta-sample table"),
                                   # Input: Select the existence of a time/step entries in the meta-sample table
                                   checkboxInput(inputId = "time_step", label = "time/step entries in the sample meta-table", value = TRUE),
                                   bsPopover(id = "time_step", title = "time step entries existence", 
                                             content = paste0('This box must be un-checked in case you have a sample meta-table without any ',
                                                            'time/step oredered entry column. This will avoid an automatic check to control the ',
                                                            'correct consistency of the table')),
                                   
                                   tags$hr(style="border-color: green;"),
                                   
                                   fileInput(inputId = "marker_file", label = h4("marker-data table"), multiple = F,
                                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"),
                                             buttonLabel = "load your .csv file", placeholder = "upload your marker-data table"),
                                   
                                   tags$hr(style="border-color: green;"),
                                   
                                   fileInput(inputId = "pheno_file", label = h4("pheno-data table"), multiple = F,
                                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"),
                                             buttonLabel = "load your .csv file", placeholder = "upload your pheno-data table"),
                                   
                                   tags$hr(style="border-color: green;"),
                                   
                                   h4("MDS parameters"),
                                   # Input: Select the distance method ---
                                   p("The following parameter is to set distance measures for the MDS evaluation plots"),
                                   radioGroupButtons(
                                       inputId = "mdsMethod", label = "Select your method for computing the distance matrix for the MDS evalution:", 
                                       choices = c( "euclidean", "maximum", "minkowski", "pearson", "spearman", "kendall"), 
                                       justified = TRUE, status = "primary", direction = "vertical", selected = "euclidean", 
                                       checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                                   bsPopover(id = "mdsMethod", title = "MDS method", content = paste0('Select your favourite distance method:' ,
                                                                              'euclidean: is the root sum-of-squares of differences',
                                                                              'maximum is the maximum difference',
                                                                              'minkowski can be considered as a generalization of the Euclidean distance')),
                                   
                                   tags$hr(style="border-color: green;"),
                                   
                                   h4("number of cluster estimation parameters"),
                                   # Input: Select the cluster evaluation paramters ---
                                   p("The following parameter is to set the max number of cluster evaluation with kmeans"),
                                   p("This is for the elbow method evaluation. Please notice than the elbow cannot always be unambiguously 
                                     identified"),
                                   
                                   
                                   numericInput(inputId = "kmeans_cluster", label = "max number of cluster estimation with kmeans", 
                                                value = 10, min = 2, max = 50),br(),
                                   p("The following paramter is to set the max number of cluster evaluation with mclust"),
                                   p("The mclust is based on several different models. Choose the best model within the plot."),
                                   
                                   numericInput(inputId = "max_mclust", label = "max number of cluster estimation with mclust", 
                                                value = 10, min = 3, max = 30),br(),
                                   
                                   tags$hr(style="border-color: green;"),
                                   # Input: Select the type of clustering algorithm---
                                   h4("map evalution with clusterization"),
                                   radioGroupButtons(
                                       inputId = "tsne_algo", label = "Choose your clustering algorithm to clusterize yor tsne map:", 
                                       choices = c("k_means", "SOM"), 
                                       justified = TRUE, status = "primary", direction = "horizontal", selected = "k_means", 
                                       checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                                   bsPopover(id = "tsne_algo", title = "clustering method", 
                                             content = paste0("Select your favourite clustering algorithm to clusterize the tSNE map")),
                                   
                                   p("the following parameter is to set the number of cluster to use for the tSNE map evaluation"),
                                   numericInput(inputId = "NumClust", label = "number of clusters", 
                                                value = 20, min = 2, max = 200))

######################### sidebar panel clust --------------------------------------------
sidebar_panel_clust <- sidebarPanel(width = 3,
                                    tags$head(
                                        tags$style(HTML("hr {border-top: 1px solid #008000;}"))
                                    ),
                                    
                                    h4("console messages"),
                                    verbatimTextOutput("console_output_clust", placeholder = T),

                                    tags$hr(style="border-color: green;"),
                                    actionButton("runAuto", label = "Run the whole analysis with the current parameters", icon("paper-plane"), 
                                                 style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),br(),
                                    
                                    conditionalPanel(
                                        condition = ("output.runAuto_panelStatus"),
                                        h4("Your analysis report is ready"),
                                        downloadButton(outputId = "downloadRunAuto",label = "download your analysis material")),
                                    
                                    tags$hr(style="border-color: green;"),
                                    h4("Clustering parameters"),
                                    sliderInput(inputId = "k", label = paste0("select the number of clusters for K_means algorithm ", 
                                                                              "or the number of number of nearest neighbours for Phenograph"), 
                                                min = 2, max = 200, value = 30, ticks = TRUE, dragRange = FALSE),
                                    bsPopover(id = "k", title = "cluster quantity", 
                                              content = paste0("Select the number of clusters ", "Please notice than flowSOM currently works only with the default ",
                                                                              "value of 100 clusters")),
                                    
                                    radioGroupButtons(
                                        inputId = "ClustAlgo", label = "Choose your clustering algorithm:", 
                                        choices = c("FlowSOM" = "flowSOM", "Rphenograph" = "Rphenograph", "FastPG" = "FastPG", "DEPECHE" = "DepecheR",  
                                                    "K-means" = "KMeans"),
                                        justified = TRUE, status = "primary", direction = "vertical", selected = "flowSOM", 
                                        checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                                    bsPopover(id = "ClustAlgo", title = "Clustering algorithm", 
                                              content = paste0("Select your favourite clustering algorithm")),
                                    
                                    radioGroupButtons(
                                      inputId = "KmeansAlgo", label = "Kmeans algorithm", 
                                      choices = c("Hartigan-Wong" = "Hartigan-Wong", "Forgy" = "Forgy", "MacQueen" = "MacQueen"),
                                      selected = "Hartigan-Wong"), 
                                    bsPopover(id = "KmeansAlgo", title = "kmeans related method", 
                                              content = paste0("Select your favourite Kmeans algorithm")),
                                   
                                    radioGroupButtons(
                                      inputId = "ClusteringPdistance", label = "Clustering validation distance", 
                                      choices = c("Euclidean" = "euclidean", "Maximum" = "maximum", "Minkowski" = "minkowski"), 
                                      selected = "euclidean" ),
                                    bsPopover(id = "ClusteringPdistance", title = paste0("The clustering (and meta-clustering) performance evaluation 
                                                                                         is based on the way you measure ", "the event's distance: you 
                                                                                         can choose from 'euclidean' ", "'maximum', or 'minkowski",
                                                                                         "The power of the minkowsky distance is the amount of selected channels")),
                                    tags$hr(style="border-color: green;"), 
                                    h4("Meta-clustering parameters"),
                                    
                                    radioGroupButtons(
                                      inputId = "MetaClustAlgo", label = "MetaClustering algorithm", 
                                      choices = c("ConsensusClusterPlus" = "ConsensusClusterPlus", "Other" = "Other"), 
                                      selected = "ConsensusClusterPlus" ),
                                    bsPopover(id = "MetaClustAlgo", title = paste0("Select your favourite meta-clustering algorithm ", 
                                                                                   "...others will be soon implemented!")),
                                    
                                    sliderInput(inputId = "MetaClustNum", label = "select the number of Meta-Clusters", 
                                                min = 2, max = 50, value = 20, ticks = TRUE, dragRange = FALSE),
                                    bsPopover(id = "MetaClustNum", title = paste0("Set the number of meta-clusters - ", 
                                                                                   "Warning: it goes without saying that the number of meta-clusters ",
                                                                                   "should be less than the number of clusters")),
                                    
                                    radioGroupButtons(
                                      inputId = "MetaClustCriteria", label = "Meta-Clustering criteria", 
                                      choices = c("hierarchical (hclust)" = "hc", "paritioning around medoids (pam)" = "pam", "k-means" = "km"),
                                      selected = "hc"),
                                    bsPopover(id = "MetaClustCriteria", title = paste0("The meta-clustering algorithm is based on a criteria for ", 
                                                                                       "grouping the clusters: (meta-clustering is just a way to ", 
                                                                                       "'clusters' the data from clustering). You can choose from: ", 
                                                                                       "hierarchical clustering, paritioning around medoids or ",
                                                                                       "k-means")),
                                    
                                    radioGroupButtons(inputId = "MetaClustDist", 
                                                 label = "Meta-Clustering distance", 
                                                 choices = c("Pearson correlation" = "pearson",
                                                             "Euclidean" = "euclidean",
                                                             "Maximum" = "maximum",
                                                             "Minkowski" = "minkowski"), 
                                                 selected = "euclidean"),
                                    bsPopover(id = "MetaClustDist", title = paste0("Each meta-clustering alogithm is based on the way you measure ", 
                                                                                    "the event's distance")),
                                    
                                    tags$hr(style="border-color: green;"),
                                    h4("Labeling parameters"),
                                    
                                    radioGroupButtons(
                                        inputId = "signature_finding_method", label = "Choose your  signature finding method:", 
                                        choices = c("Central tendency (original)" = "Central_tendency", "Densitie's study" = "Densities"),
                                        justified = TRUE, status = "primary", direction = "horizontal", selected = "Central_tendency", 
                                        checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                                    bsPopover(id = "signature_finding_method", title = paste0("Select your favourite signature finding method")),
                                    
                                    checkboxInput(inputId = "no_Quant", label = "no Quantiles handling", value = FALSE),
                                    bsPopover(id = "no_Quant", 
                                              title = paste0("Select this to exclude the effect of quantiles on the heatmaps' expression values")),
                                    
                                    sliderInput(inputId = "minQuantile", label = "minQuantile: max effect when set to 0.5", 
                                                min = 0.01, max = 0.5, value = 0.05, ticks = TRUE, dragRange = FALSE),
                                    sliderInput(inputId = "maxQuantile", label = "maxQuantile: max effect when set to 0.5", 
                                                min = 0.5, max = 1.0, value = 0.95, ticks = TRUE, dragRange = FALSE),
                                    bsPopover(id = "maxQuantile", title = "min and max quantile", 
                                    content = paste0("The minQuantile and the maxQuantile parameter will act upon the heatmap values calculation")),
                                    
                                    sliderInput(inputId = "pos_threshold", label = "select the threshold for the positive MFI value", 
                                                min = 0.5, max = 0.8, value = 0.5, ticks = TRUE, dragRange = FALSE),
                                    
                                    sliderInput(inputId = "neg_threshold", label = "select the threshold for the negative MFI limit", 
                                                min = 0.2, max = 0.5, value = 0.5, ticks = TRUE, dragRange = FALSE),
                                    bsPopover(id = "neg_threshold", title = "min and max threshold", 
                                              content = paste0("The min and the max threshold will act upon the heatmap values calculation")),
                                    
                                    radioGroupButtons(
                                      inputId = "HeatMapCentral", label = "HeatMap Central tendency measure", 
                                      choices = c("Mean" = "mean", "Median" = "median", "Mode" = "mode"), selected = "median"),
                                    
                                    radioGroupButtons(
                                      inputId = "HeatMapMethod", label = "Heat-Map distance (method)", 
                                      choices = c("euclidean" = "euclidean", "maximum" = "maximum", "binary" = "binary", "minkowski" = "minkowski"), 
                                      selected = "euclidean"),
                                    #manhattan has been discarded since it is related to discrete distances
                                    
                                    radioGroupButtons(inputId = "HeatmapHCCriteria", 
                                                 label = "Hierchical-Clustering criteria for the Heat-map", 
                                                 choices = c("average" = "average", 
                                                             "mcquitty" = "mcquitty", 
                                                             "median" = "median",
                                                             "centroid" = "centroid"), 
                                                 selected = "median"),
                                    
                                    tags$hr(style="border-color: green;"),
                                    h4("Internal cluster validation"),
                                    checkboxInput(inputId = "clusteringPerf", label = "Cluster validation", value = FALSE),
                                    bsPopover(id = "clusteringPerf", title = "Clustering perfomance evaluation", 
                                              content = paste0("This box must be un-checked in case you do not want to evaluate the clustering 
                                                             (and meta-clustering) performance: This checkbox is introduced in order to save computatonal time")),
                                    
                                    tags$hr(style="border-color: green;"),
                                    h4("Map general parameters"),
                                    # Input: To tune the maximum number limitation for the map processing
                                    checkboxInput(inputId = "map_limit", label = "map processing limitation", value = TRUE),
                                    bsPopover(id = "map_limit", 
                                              title = paste0('This box must be un-checked in case you do not want any limitation in processing ',
                                                             'the selected dimensionality reduction algorithm. Removing the limit could mean to force ', 
                                                             'a too demanding process and the app could stop in the middle of the map computation')),
                                    
                                    numericInput(inputId = "mapMax", 
                                                 label = "Max number of point for the selected map", value = 100000, min = 100, max = 1000000, step = 5),
                                    bsPopover(id = "mapMax", 
                                              title = paste0('Non linear dimensionality reduction algorithms are among the most resource-intensive. ',
                                                             'Here is the maximum amount of processed event which applies if the "map limitation" is checked')),
                                    
                                    radioButtons(inputId = "mapType", label = "Select the type of map", 
                                                 choices = c("tSNE" = "tSNE", "UMAP" = "UMAP"), selected = "tSNE"),
                                    
                                    switchInput(inputId = "two_three", label = "2D or 3D", onLabel = "3D", offLabel = "2D", value = FALSE, 
                                                handleWidth = 200), ####?????
                                    
                                    tags$hr(style="border-color: green;"),
                                    h4("tSNE configuration"),
                                    
                                    switchInput(inputId = "tSNE_PCA", label = "pre-PCA step to be performed", value = TRUE),
                                    
                                    p("For the perplexity parameter it should hold the following: (3 * perplexity < nrow(X) - 1)"),
                                    sliderInput(inputId = "perplexity", label = "select the perplexity", 
                                                min = 1, max = 200, value = 30, ticks = TRUE, dragRange = FALSE),
                                    
                                    sliderInput(inputId = "theta", label = "select theta value", 
                                                min = 0.1, max = 1.0, value = 0.5, ticks = TRUE, dragRange = FALSE),
                                    
                                    sliderInput(inputId = "tSNEiter", label = "select number of iteration", 
                                                min = 10, max = 5000, value = 1000, ticks = TRUE, dragRange = FALSE),
                                    
                                    sliderInput(inputId = "tSNEeta", label = "select eta value", 
                                                min = 10, max = 500, value = 200, ticks = TRUE, dragRange = FALSE),
                                    
                                    tags$hr(style="border-color: green;"),
                                    h4("UMAP configuration"),
                                    sliderInput(inputId = "UMAP.n_neighbors", label = "nearest neighbors", 
                                                min = 2, max = 200, value = 20, ticks = TRUE, dragRange = FALSE),
                                    
                                    radioButtons(inputId = "UMAP.metric", label = "Select the type of metric", 
                                                 choices = c("euclidean" = "euclidean", "chebyshev" = "chebyshev", "minkowski" = "minkowski"), 
                                                 selected = "euclidean"),

                                    sliderInput(inputId = "UMAP.min_dist", label = "minimum distance", 
                                                min = 0.01, max = 0.99, value = 0.1, step = 0.01, ticks = TRUE, dragRange = FALSE),
                                    
                                    sliderInput(inputId = "UMAP.spread", label = "spread", 
                                                min = 0.1, max = 5, value = 2, step = 0.1, ticks = TRUE, dragRange = FALSE))

######################## body ---------------------------------------------------------------------------
dbBody <- dashboardBody(
    useShinyjs(),  # Include shinyjs
    #extendShinyjs(functions = "shinyjs.beep", script = "./www/beep.js"), trying to use javascript to produce sound: invain
    tabItems(
        ### Intro tabItem ----------------------------------------------
        tabItem(tabName = "intro",
                
                #extendShinyjs(text = 'shinyjs.ding = function() {var snd = new Audio("./www/ding.ogg"); snd.play(); }'),  
                # I tried this to play a sound at the end of each long procedure... invain
                h1("cytoChain: a web-app to analyze your high dimensional flow cytometry experiments"),
                br(),
                h4("You can download the manual here:"),
                downloadButton(outputId = "download_manual",label = "download cytoChain's manual"),br(),
                
                h4("The whole environment is conceived in a modular way, in workflows, to setup your analysis, performing the proper actions"), 
                h4("Three workflows are currently implemented:", 
                  tags$ol(tags$li("FlowSet optimization: to prepare your flowFrames for the further analysis"),
                          tags$li("Metadata & flowSet assays: for adding additional data to characterize and to evaluate your samples"), 
                          tags$li("High dimensional analysis: to perform the whole flowSet analysis"))),
                img(src = "Workflows.png", height = 500, width = 1000),
                h4("The flowSet optimization workflow is made of the following modules:"),
                img(src = "Puzzle-preII.png", height = 400, width = 1100),
                p("Please, notice that the 'Align' module is still under deplyment"),br(),
                h4("The metadata & flowSet assays workflow is made of the following modules:"),
                img(src = "Puzzle-metaII.png", height = 500, width = 900),br(),br(),
                h4("...and finally, the High dimensional workflow is made of the following modules:"),
                img(src = "Puzzle-clustII.png", height = 400, width = 850),
                h3("Attention! The app is still in a very draft version. 
                   The author declines all responsibility resulting from any kind of use")
        ),
        
        ### flowSet optimization tabItem ----------------------------------------------
        tabItem(tabName = "preClustering",
                titlePanel("flowSet optimization"),
                # Sidebar layout with input and output definitions ---
                sidebarLayout(
                    # Sidebar panel for inputs ---
                    sidebar_panel_pre,
                    # Main panel for displaying outputs ---
                    mainPanel(
                        navbarPage(title="Workflow", id = "preClusteringWorkflow", inverse = T, position = "static-top", 
                               #tabsetPanel(
                                   tabPurpose,  
                                   tabClean,
                                   tabScale,
                                   tabAlign,
                                   tabDowns,
                                 #  tabSubss,
                                   tabConc,
                                   #id = "preClusteringWorkflow", type = "tabs"), 
                                   fluid = T)))), #end of preClustering tabItem
        
        ### Metadata tabItem ----------------------------------------------
        tabItem(tabName = "metadata",
                
                titlePanel("MetaData editing"),
                # Sidebar layout with input and output definitions ---
                sidebarLayout(
                    # Sidebar panel for inputs ---
                    sidebar_panel_meta,
                    # Main panel for displaying outputs ---
                    mainPanel(
                        navbarPage(title="edit Meta", id = "metaDataWorkflow", inverse = T, position = "static-top", 
                                   tabMeta,
                                   tabAna,
                                   tabEva,
                                   tabtSNEeva), fluid = T))), #end of metadata tabItem
        
        ### Classical Workflow tabItem ----------------------------------------------
        tabItem(tabName = "classical",
                titlePanel("High dimensional analysis workflow"),
                sidebarLayout(
                    sidebar_panel_clust,
                    mainPanel(
                        navbarPage(title="Clustering pipeline", inverse = T, position = "static-top",
                                   tabClust,
                                   tabLabelClust,
                                   tabMetaClust,
                                   tabMapComp,
                                   tabLabel,
                                   tabLabel_p,
                                   tabLabel_pp,
                                   tabQuant,
                                   tabQuant_p), fluid = T))), #end of classical tabItem
        
        ### Alternative Workflow tabItem ----------------------------------------------
        tabItem(tabName = "alternative",
                h2("t.b.d")),
        
        ### Panel editor tabItem ----------------------------------------------
        tabItem(tabName = "Panel_editor",
                titlePanel(""),
                mainPanel(
                    navbarPage("Panel Editor",
                               div(
                                   tabPanel("edit",
                                        h4("Push if you need to edit the panel's names"),
                                        p("When you collect samples from different experiments, or when you import FCS files from other 
                                          environments, it could be some disalignments in the panel's name. Here you can try to fix them"),
                                        p("Thanks to Pier Federico Gherardini for this idea from '*premessa*' package, see: ",
                                            tags$a(href="https://github.com/nolanlab/scaffold", "SCAFFoLD")), 
                                        br(),
                                        
                                        #fileInput(inputId = "FCS_input", label = h4("File input"), multiple = T, 
                                        #          buttonLabel = "FCS files", placeholder = "upload your flowFrames"),br(),
                                        
                                        fluidRow(
                                            column(3,
                                                   h3("Load your samples as collection of FCS files"),
                                                   fileInput(inputId = "FCS_input", label = h4("File input"), multiple = T,
                                                             buttonLabel = "FCS files", placeholder = "upload your flowFrames"),
                                                   p("The samples must be valid FCS files of version 3.0 or FCS files of a more recent version")), 
                                            column(1, offset = 1,
                                                   h3("...or")),
                                            column(4,
                                                   h3("Load your samples collected in a single zip file"),
                                                   fileInput(inputId = "zip_input", label = h4("File input"), multiple = FALSE, 
                                                             accept = c("application/zip", "ZIP archive", ".zip"),
                                                             buttonLabel = "zip file", placeholder = "upload your single zip file"),
                                                   p("The samples inside the zip file must be valid FCS files of version 3.0+. Only the zip format is supported"))),
                                        
                                        
                                        h3("Perform loading & parsing process"),
                                        actionButton("panelEdit", "Load and parse your samples", icon("paper-plane")),
                                                     #style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                                        h3(textOutput("checktxt_parsing_panel")), 
                                        #verbatimTextOutput("filepaths"),
                                        br(),
                                        #textInput("paneleditorui_output_folder", 
                                        #          label = "Output folder name (it will be a sub-folder of your selected one)", value = "fixedFCS"),
                                        br(),
                                        rHandsontableOutput(outputId = "paneleditorui_panel_table"), style = 'width:2000px;'), #div
                                        actionButton("paneleditorui_process_files", "Process files", icon("paper-plane"),
                                                     style="color: #fff; background-color: #ff0000; border-color: #2e6da4"),
                                        conditionalPanel(
                                            condition = ("output.PE_panelStatus"), 
                                            downloadButton(outputId = "download_FCS",label = "download your new flowFrames"))
                                   , style = 'width:2000px;'), #div
                               #tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 12600px; /* or 950px */}"))
                               tags$head(tags$style(type="text/css", ".container-fluid {  max-width:12600px */}"))
                               ), fluid = T)),

        ### credits tabItem ----------------------------------------------
        
        tabItem(tabName = "credits",
                
                p("Here you can find the list of packages which are used in the various parts of the", strong("citoChain"), "code"),br(),
                p("Some of them constitute the sources of inspiration for this work. I should like to offer my heartfelt thanks to
                  all of the entire community of reasearcher and in particular to"),
                shiny::em("Malgorzata Nowicka, Carsten Krieg, Lukas M. Weber, Felix J. Hartmann, Silvia Guglietta, Burkhard Becher, 
                   Mitchell P. Levesque, Mark D. Robinson"), 
                p("for their workflow (see: ", 
tags$a(href="http://bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/CyTOF/doc/cytofWorkflow_BioC2017workshop.html#data-import",
       "BioC2017/Day2/Workshops"),
                ") and all the developers and bioinformatics who makes available some of their precious packages, in particular:"), 

                p("for ", strong('flowCore')),
                shiny::em("B Ellis, Perry Haaland, Florian Hahne, Nathan Le Meur, Nishant Gopalakrishnan, Josef Spidlen, Mike Jiang and Greg Finak 
                (2019). flowCore: flowCore: Basic structures for flow cytometry data"), br(),  br(),

                p("for ",  strong('flowAI')),
                shiny::em("Monaco,G. et al. (2016) flowAI: automatic and interactive anomaly discerning tools for flow 
                          cytometry data. Bioinformatics. 2016 Aug 15;32(16):2473-80"),br(),  br(),
                p("for ",  strong('flowStats')),
                shiny::em("Florian Hahne, Nishant Gopalakrishnan, Alireza Hadj Khodabakhshi, Chao-Jen Wong and Kyongryun Lee (2019). 
                          flowStats: Statistical methods for the analysis of flow cytometry data. http://www.github.com/RGLab/flowStats"),br(),  br(),
                p("for ",  strong('spade')),
                shiny::em("M. Linderman, P. Qiu, E. Simonds and Z. Bjornson (2018). spade: SPADE -- An analysis and visualization tool for Flow 
                          Cytometry. http://cytospade.org"), br(),  br(),
                p("for ",  strong('limma')),
                shiny::em(" Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential 
                expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47."),br(),  br(),
                p("for ", strong('matrixStats')),
                shiny::em("Henrik Bengtsson (2018). matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors). 
                          https://CRAN.R-project.org/package=matrixStats"), br(),  br(),
                p("for ",  strong('factoextra')),
                shiny::em("Alboukadel Kassambara and Fabian Mundt (2017). factoextra: Extract and Visualize the Results of Multivariate Data Analyses. 
                          https://CRAN.R-project.org/package=factoextra"), br(),  br(),
                p("for ",  strong('DDoutlier')),
                shiny::em("Jacob H. Madsen (2018). DDoutlier: Distance & Density-Based Outlier Detection. 
                          https://CRAN.R-project.org/package=DDoutlier"), br(),  br(),
                p("for ",  strong('mclust')),
                shiny::em(" Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation
                          using Gaussian finite mixture models The R Journal 8, pp. 205-233"),br(),  br(),
                p("for ",  strong('kohonen')),
                shiny::em(" Wehrens R, Kruisselbrink J (2018). Flexible Self-Organizing Maps in kohonen 3.0. _Journal of Statistical Software_, *87*(7), 1-18. 
                doi:10.18637/jss.v087.i07 (URL: https://doi.org/10.18637/jss.v087.i07)"),br(),  br(),
                p("for ",  strong('FlowSOM')),
                shiny::em("Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2019). FlowSOM: Using self-organizing maps for visualization and 
                          interpretation of cytometry data. http://www.r-project.org, http://dambi.ugent.be."), br(),  br(),
                p("for ",  strong('Rphenograph')),
                shiny::em(" Hao Chen (2015). Rphenograph: R implementation of the phenograph algorithm"), br(),  br(),
                p("for ",  strong('ConsensusClusterPlus')),
                shiny::em("Wilkerson, M.D., Hayes, D.N. (2010). ConsensusClusterPlus: a class discovery tool with confidence assessments and item 
                          tracking. Bioinformatics, 2010 Jun 15;26(12):1572-3."), br(),  br(),
                p("for ",  strong('flowVS')),
                shiny::em("Ariful Azad (2019). flowVS: Variance stabilization in flow cytometry (and microarrays)"), br(),  br(),
                p("for ",  strong('Rtsne')),
                shiny::em("L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning 
                Research 9(Nov):2579-2605, 2008. L.J.P. van der Maaten. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning 
                Research 15(Oct):3221-3245, 2014. Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation, 
                          URL: https://github.com/jkrijthe/Rtsne"), br(),  br(),
                p("for ",  strong('umap')),
                shiny::em("Tomasz Konopka (2019). umap: Uniform Manifold Approximation and Projection. https://CRAN.R-project.org/package=umap"),
                tags$hr(style="border-color: green;"),br(),  br(),
                p("for ",  strong('rstatix')),
                shiny::em("Alboukadel Kassambara. rstatix: Pipe-Friendly Framework for Basic Statistical Tests. https://cran.r-project.org/web/packages/rstatix/"),
                tags$hr(style="border-color: green;"),br(),  br(),

                p("I cannot forget to mention all the effort made by members of the R (and Python) developers community for all the features of 
                  such an indispensable set of accessory libraries:"),
                p("for ",  strong('DT')),
                shiny::em("Yihui Xie, Joe Cheng and Xianying Tan (2019). DT: A Wrapper of the JavaScript Library 'DataTables'. 
                          https://CRAN.R-project.org/package=DT"),br(),  br(),
                p("for ",  strong('tidyverse')),
                shiny::em("Hadley Wickham (2017). tidyverse: Easily Install and Load the 'Tidyverse'. https://CRAN.R-project.org/package=tidyverse"),br(),  br(),
                p("for ",  strong('lubridate')),
                shiny::em("Garrett Grolemund, Hadley Wickham (2011). Dates and Times Made Easy with lubridate. Journal of Statistical Software, 40(3), 1-25. 
                          URL http://www.jstatsoft.org/v40/i03/."),br(),  br(),
                p("for ", strong('plotly')),
                shiny::em("Carson Sievert (2018) plotly for R. https://plotly-r.com"),br(),  br(),
                p("for ",  strong('RColorBrewer')),
                shiny::em("Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. https://CRAN.R-project.org/package=RColorBrewer"),br(),  br(),
                p("for ",  strong('reshape2')),
                shiny::em("Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. 
                          URL http://www.jstatsoft.org/v21/i12/."), br(),  br(),
                p("for ",  strong('rhandsontable')),
                shiny::em("Jonathan Owen (2019). rhandsontable: Interface to the 'Handsontable.js' Library. http://jrowen.github.io/rhandsontable/"),br(),  br(),
                p("for ",  strong('pheatmap')),
                shiny::em("Raivo Kolde (2018). pheatmap: Pretty Heatmaps"), br(),  br(),
                p("for ",  strong('png')),
                shiny::em("Simon Urbanek (2013). png: Read and write PNG images. R package - https://CRAN.R-project.org/package=png"),br(),  br(),
                p("for ",  strong('gridExtra')),
                shiny::em("Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for 'Grid' Graphics. 
                          https://CRAN.R-project.org/package=gridExtra"),br(),  br(),
                p("for ",  strong('streamgraph')),
                shiny::em("Bob Rudis (2015). streamgraph: streamgraph is an htmlwidget for building streamgraph visualizations. 
                          http://github.com/hrbrmstr/streamgraph"), br(),  br(),
                tags$hr(style="border-color: green;"),
                
                p("and finally the shiny developers:"),
                p("for ",  strong('shiny')),
                shiny::em("Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2019). shiny: Web Application Framework for R. 
                          https://CRAN.R-project.org/package=shiny"),  br(),  br(),
                p("for ",  strong('shinythemes')),
                shiny::em("Winston Chang (2018). shinythemes: Themes for Shiny. https://CRAN.R-project.org/package=shinythemes"),br(),  br(),
                p("for ",  strong('shinydashboard')),
                shiny::em("  Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. 
                          https://CRAN.R-project.org/package=shinydashboard"),br(),  br(),
                p("for ",  strong('shinycssloaders')),
                shiny::em("Andras Sali (2017). shinycssloaders: Add CSS Loading Animations to 'shiny' Outputs. 
                          https://CRAN.R-project.org/package=shinycssloaders"), br(),  br(),
                p("for ",  strong('shinyWidgets')),
                shiny::em("Victor Perrier, Fanny Meyer and David Granjon (2019). shinyWidgets: Custom Inputs Widgets for Shiny.
                          https://CRAN.R-project.org/package=shinyWidgets"), br(),  br(),
                p("for ",  strong('shinyBS')),
                shiny::em(" Eric Bailey (2015). shinyBS: Twitter Bootstrap Components for Shiny. https://CRAN.R-project.org/package=shinyBS"), br(), br(),
                p("for ",  strong('shinyjs')),
                shiny::em("Dean Attali (2020). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. 
                          https://CRAN.R-project.org/package=shinyjs"), br(),

p("This is cytoChain_4.36 - the fcs samples are still handled in memory as a flowSet"),
p("In this version, the new implementation, along with the t-student test for each cluster/metacluster and each valuable tag, a new Effect Size 
  (https://www.rdocumentation.org/packages/rstatix/versions/0.7.2) evaluation is performed"), 
p("Saving of the two metadata files (mt = meta_sample.csv and mk =  meta_marker.csv)"),
p("The new implementation of the density study, to dig inside each of the cluster and to evaluate the possibility to increase the number of clusters"), 


))#end of tabItem
)#end of dashBoardBody

### main UI ----
UI <- dashboardPage(header = dbHeader, sidebar = dbSidebar, body = dbBody, skin = "green")

#??? http://htmlpreview.github.io/?https://github.com/KrishnaswamyLab/phateR/blob/master/inst/examples/bonemarrow_tutorial.html
#??? http://cole-trapnell-lab.github.io/monocle-release/ 
#??? https://stat.ethz.ch/R-manual/R-devel/library/base/html/Version.html 
#??? http://htmlpreview.github.io/?https://github.com/KrishnaswamyLab/phateR/blob/master/inst/examples/bonemarrow_tutorial.html 