
data_path = file.path("01_Data", "Input")
if (!dir.exists(data_path)) dir.create(data_path, recursive = TRUE)

egnyte_path = file.path("Z:", "Shared", "Projects", "2019", "D201901453.04 - DISE-ICF FIOB-1", 
                        "03 Working Docs_Analysis", "InputData")

input_files = list.files(egnyte_path, pattern = "csv|xlsx")

for (i in input_files){
  file.copy(file.path(egnyte_path, i),
            file.path(data_path, i),
            overwrite = TRUE)
}

