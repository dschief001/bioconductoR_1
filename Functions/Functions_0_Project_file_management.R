Compare_vendor_vs_mzML_files <- function(Selector_1,Selector_2) {
  rm(Project_mzML.files)
  # Filtering and saving
  Project_files <-   tibble(Filepath=list.files(
    path = str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),pattern = ".raw",full.names = T)
  ) %>%mutate(
    Filename = (str_extract(string = Filepath,pattern = "[^/]+$") %>% str_split(pattern = "\\.",simplify = T))[,1],
    Dateofchange=file.mtime(Filepath),
    RAWsize_MB= file.size(Filepath)*1e-6
  ) %>% arrange(Dateofchange) %>% slice(which(Filename == Selector_1):which(Filename ==Selector_2))# remove this for save of whole RAW data file
  dir.create(path = file.path(str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"mzML_Data"))
  
  tryCatch(expr = {
    mzML_info <- tibble(Filepath = list.files(file.path(str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"mzML_Data"),pattern = ".mzML",full.names = T)) %>%
      mutate(Filename = (str_extract(string = Filepath,pattern = "[^/]+$") %>% str_split(pattern = "\\.",simplify = T))[,1],mzMLsize_MB= file.size(Filepath)*1e-6)
  # Project_files <- mzML_info %>% select(-mzMLsize_MB)
    Project_mzML.files <<- left_join(x = Project_files %>% rename(Filepath_RAW="Filepath"),
                                   y = mzML_info,
                                   by= "Filename"
    ) %>% arrange(Dateofchange)},error = function(e) {
    dir.create(path = file.path(str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"mzML_Data"))
    Project_mzML.files <<- Project_files%>% rename(Filepath_RAW="Filepath")%>% mutate(Filepath=NA)})
  View(Project_mzML.files)
}
