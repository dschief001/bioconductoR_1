# RAW conversion tool: ----------------------------------------------------
# create the missing mzML files text file...
# Thu Dec 12 11:42:46 2019 ------------------------------
pacman::p_load(tidyverse,here)
Vendor_specific_file_converter <- function(Project_mzML.files){
  write(x = Project_mzML.files %>% filter(is.na(Filepath)) %>% select(Filepath_RAW) %>% unlist(),file = file.path(str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"files.txt"))
  setwd(stringr::str_remove(string = here::here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"))
  system('powershell -command "msconvert -f "files.txt" --zlib -o "mzML_Data""',intern = T)
  file.remove(file.path(str_remove(string = here(),pattern = "/[0-9,A-z,_ü]*(?=[^/]*$)"),"files.txt"))
  getwd()
}
Vendor_specific_file_converter(Project_mzML.files = Project_mzML.files)