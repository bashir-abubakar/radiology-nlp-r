# load entity extractor function
source("D:/fall_risk_project/nlp_entity_extract.R")

# root folder
root_path <- "D:/fall_risk_project/mimic-cxr-reports/files/"

# list all txt files recursively
all_txt_files <- list.files(
  path = root_path,
  pattern = "\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

length(all_txt_files) 

# random sample 100
set.seed(123)
sample_files <- sample(all_txt_files, 100)

# read text
sample_texts <- purrr::map_chr(sample_files, ~ readr::read_file(.x))

# extract
entities_100 <- extract_entities(sample_texts, doc_ids = basename(sample_files))

View(entities_100) 

print(entities_100, n = Inf)

write.csv(entities_100,
          "D:/fall_risk_project/entities_sample100.csv",
          row.names = FALSE)
