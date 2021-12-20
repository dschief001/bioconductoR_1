pacman::p_load(tidyverse,here)
source(file = here("Functions","Functions_0_Project_file_management.R"))
# Filemanagement in Projects! ---------------------------------------------



# RAW files ---------------------------------------------------------------
# File-change dependent selection:
View(tibble(Filepath=list.files(
  path = str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),pattern = ".raw",full.names = T)
) %>%mutate(
  Filename = (str_extract(string = Filepath,pattern = "[^/]+$") %>% str_split(pattern = "\\.",simplify = T))[,1],
  Dateofchange=file.mtime(Filepath),
  RAWsize_MB= file.size(Filepath)*1e-6
) %>% arrange(Dateofchange)
)

# Select the framing Samples ----------------------------------------------
Selector_1 <- "Blank_1_1"
Selector_2 <- "Dilution_2_V4"

# Fileplotting timewise

tibble(Filepath=list.files(
  path = str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),pattern = ".raw",full.names = T)
) %>%mutate(
  Filename = (str_extract(string = Filepath,pattern = "[^/]+$") %>% str_split(pattern = "\\.",simplify = T))[,1],
  Dateofchange=file.mtime(Filepath),
  RAWsize_MB= file.size(Filepath)*1e-6
  # Type = function(Filename) {
  #   if (str_detect(string = Filename,pattern = "ution")) {
  #     # insert color filtering for different sample types here:
  #     1
  #   }else {NA}
  # }
) %>% arrange(Dateofchange) %>% slice(which(Filename == Selector_1):which(Filename ==Selector_2)) %>%
                                        ggplot(mapping = aes(x = Dateofchange,label=Filename,y=factor(1)))+
                                        geom_point(mapping = aes(x = Dateofchange,y=factor(0.5)))+
                                        geom_text(angle=90,size=5)+
                                        # scale_x_datetime(date_breaks = "day",date_labels = "%m") +
                                        scale_x_datetime(date_breaks = "hour",date_labels = "%a-%H")+
                                        scale_y_discrete()+
                                        theme_bw()


View(
  tibble(Filepath=list.files(
    path = str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),pattern = ".raw",full.names = T)
  ) %>%mutate(
    Filename = (str_extract(string = Filepath,pattern = "[^/]+$") %>% str_split(pattern = "\\.",simplify = T))[,1],
    Dateofchange=file.mtime(Filepath),
    RAWsize_MB= file.size(Filepath)*1e-6
  ) %>% arrange(Dateofchange) %>% slice(which(Filename == Selector_1):which(Filename ==Selector_2))
)

{
# library(shiny)
# library(miniUI)
# myGadgetFunc <- function() {
#   
#   ui <- miniPage(
#     gadgetTitleBar("My Gadget"),
#     miniContentPanel(
#       # Define layout, inputs, outputs
#       tableOutput(outputId = "Table1")
#     )
#   )
# 
#   server <- function(input, output, session) {
#     # Define reactive expressions, outputs, etc.
#     pacman::p_load(tidyverse,here)
#     output$Table1 <- renderTable(
#           tibble(Filepath=list.files(
#           path = str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),pattern = ".raw",full.names = T)
#         ) %>%mutate(
#           Filename = (str_extract(string = Filepath,pattern = "[^/]+$") %>% str_split(pattern = "\\.",simplify = T))[,1],
#           Dateofchange=file.mtime(Filepath),
#           RAWsize_MB= file.size(Filepath)*1e-6
#         ) %>% arrange(Dateofchange)
#     )
#     # When the Done button is clicked, return a value
#     observeEvent(input$done, {
#       returnValue <- "hi"
#       stopApp(returnValue)
#     })
#     observeEvent(input$cancel, {
#       stopApp(stop("No selection.", call. = FALSE))
#     })
#   }
# 
#   runGadget(ui, server, viewer = dialogViewer(dialogName = "1.",width = 1000,height = 1000))
# }
# myGadgetFunc()
} # old code stuff

Compare_vendor_vs_mzML_files(Selector_1 = Selector_1,Selector_2 = Selector_2)


# Check filenames and files.
# check for missing mzML files: 
if(nrow(Project_mzML.files %>% filter(is.na(Filepath))) !=0) {
  warning("Missing mzML Files!")
  print(Project_mzML.files %>% filter(is.na(Filepath)))
}
# RAW conversion tool: ----------------------------------------------------
# create the missing mzML files text file...
rstudioapi::jobRunScript(path = here("Functions","Vendor_specific_file_converter.R"),
                         name = "Files are Converted to mzML",
                         workingDir = here(),
                         importEnv = T)

# Snippet for future in callr: might be necessary, doesn't work wi --------
  # callr::r(function(Project_mzML.files) {
  #   write(x = Project_mzML.files %>% filter(is.na(Filepath)) %>% select(Filepath_RAW) %>% unlist(),file = file.path(stringr::str_remove(string = here::here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"files.txt"))
  #   setwd(stringr::str_remove(string = here::here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"))
  #   system('powershell -command "msconvert -f "files.txt" -o "mzML_Data""',intern = T)
  #   file.remove(file.path(stringr::str_remove(string = here::here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"files.txt"))
  #   getwd()
  # },args = list(Project_mzML.files))


##################################################
#If no raw - mzML check is required use:
##################################################
# mzML_info <- tibble(Filepath = list.files(file.path(str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"mzML_Data"),pattern = ".mzML",full.names = T)) %>%
# mutate(Filename = (str_extract(string = Filepath,pattern = "[^/]+$") %>% str_split(pattern = "\\.",simplify = T))[,1],mzMLsize_MB= file.size(Filepath)*1e-6)
# Project_files <- mzML_info %>% select(-mzMLsize_MB)
# Project_mzML.files <<- left_join(x = Project_files %>% rename(Filepath_RAW="Filepath"),
#                                  y = mzML_info,
#                                  by= "Filename")
# writexl::write_xlsx(x = Project_mzML.files,path = here(paste0("Proj_mzMLs.xlsx"))
##################################################


# Recheck the files - Selected Vendor vs mzML --------------------------------

Compare_vendor_vs_mzML_files(Selector_1 = Selector_1,Selector_2 = Selector_2)
if(nrow(Project_mzML.files %>% filter(is.na(Filepath))) !=0) {
  warning("Missing mzML Files!")
  print(Project_mzML.files %>% filter(is.na(Filepath)))
} else {message("Everything fine continue with script!")}

if (file.exists(here(paste0("Proj_mzMLs.xlsx")))) {
stop("Proj_mzMLs.xlsx File already exists!")
} else {
  writexl::write_xlsx(x = Project_mzML.files,path = here(paste0("Proj_mzMLs.xlsx")))
  # writexl::write_xlsx(x = Project_mzML.files,path = here(paste0("Proj_mzMLs_",Sys.Date(),".xlsx")))
  }

# Code to use in other scripts --------------------------------------------
# 
# Project_mzML.files <- readxl::read_xlsx(path = here(paste0("Proj_mzMLs.xlsx")))
# 
# Project_mzML.files %>% filter(str_detect(string = Filename,pattern = "Blood|Blut"))
# Project_mzML.files$Filepath