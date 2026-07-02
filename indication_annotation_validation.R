# code to organize and process prematurity project annotations

library(bigrquery)
library(tidyverse)
library(glue)

out_folder <- "/mnt/arcus/lab/users/kafadare/indication_files/"

bq_auth()
proj_id = "scit605-healthymri-2d73d69f"

# Specify the directory with the annotations of interest
ann_dir <- "/mnt/arcus/lab/apps/annotator/premature_notes/*"

# Get a list of all the annotation files
allanots <- Sys.glob(glue("{ann_dir}/*"))
allanots.proc <- allanots %>% str_remove_all("\\.ann") %>% basename()

listofannotspersubject <- list()
length(listofannotspersubject) <- length(allanots.proc)
### Loop through each annotation file and extract the proc_ord_ids, annotations, and snomed_ids
listofannotspersubject <- list()
length(listofannotspersubject) <- length(allanots.proc)
for(ann in allanots.clean){
  #read line by line because characters like character ' " ' throws the new line encoding off.
  lines <- readLines(ann)
  lines <- gsub('"', '', lines)
  annot.i <- read_delim(I(lines), delim = "\t", col_names = FALSE, show_col_types = FALSE)
  #annot.i <- suppressWarnings(read_delim(ann, col_names=F,
  #                                       delim = '\t', show_col_types = FALSE))
  if (nrow(annot.i) > 0){
    annot.i <- annot.i %>%
      mutate(proc_ord_id = str_remove_all(ann, ".ann"))
    #get numerical range for indication
    indication.line <- annot.i[grep("Indication", annot.i[,2]),2]
    extracted_numbers <- str_extract_all(gsub(";", " ", indication.line), "\\d+\\.?\\d*") # Extracts digits, optionally with a decimal point and more digits
    ind_numbers_range <- range(as.numeric(unlist(extracted_numbers))) %>% { seq(.[1], .[2]) }
    #separate snomed IDs
    snomed.i <- annot.i %>%
      filter(grepl("SNOMED",X2)) %>%
      mutate(X3 = str_remove_all(X3,"\\t")) %>%
      mutate (text_block = str_extract(X2, "T\\d+")) %>%
      separate_wider_delim(X2, delim = ":", names = c("X2","snomed_term")) %>%
      select(snomed_term, snomed_name = X3, text_block) %>%
      distinct()
    #merge snomed IDs with annotation blocks based on text block ID
    annot_clean.i <- merge(annot.i[grepl("T", annot.i$X1),], snomed.i, by.x = "X1",  by.y = "text_block", all.x = TRUE)
    annot_clean.i <- annot_clean.i %>% mutate(numbers = map(str_extract_all(gsub(";", " ", annot_clean.i$X2), "\\d+\\.?\\d*"), as.numeric)) %>%
      mutate(X2 = gsub("[0-9;]", "", annot_clean.i$X2)) %>%
      filter(!str_detect(X2, "Indication"))
    #annot_clean.i$numbers <- map(annot_clean.i$numbers, as.numeric)
    #Create TRUE/FALSE column for within indication range
    annot_clean.i <- annot_clean.i %>%
      mutate(ind_range = map_lgl(numbers, ~ all(.x %in% ind_numbers_range))) %>%
      mutate(type = ifelse(ind_range, "Indication", "Finding")) %>%
      select(-c("numbers", "ind_range"))
    #Create Label column to combine snomed labels and non-snomed finding labels
    annot_clean.i$label <- ifelse(!is.na(annot_clean.i$snomed_name), annot_clean.i$snomed_name, annot_clean.i$X2)
    #Rename columns
    annot_clean.i <- annot_clean.i %>% rename(text_block = X1, main_tag = X2, tagged_text = X3)
  } else {
    annot_clean.i <- NULL
  }
  if (ann == allanots.clean[1]) {
    annots_bind = annot_clean.i
  } else {
    if (!is.null(annot_clean.i)) {
      if (nrow(annot_clean.i > 0)) {
        annots_bind = bind_rows(annots_bind, annot_clean.i)
      }
    }
  }
}
annots_bind$proc_ord_id <- basename(annots_bind$proc_ord_id)
annots <- annots_bind %>% 
  rename(snomed_code = snomed_term)

snomed_code_list <- paste0("'", sapply(strsplit(annots$snomed_code, ","),trimws), "'", collapse = ", ")

# Get snomed tables from arcus for snomed codes in the annotations table
snomed_icd10_query <- glue('SELECT
  icd10.icd10_code,
  icd10.snomed_code,
  icd10.snomed_name
  FROM arcus.icd10_to_snomed AS icd10
  WHERE snomed_code IN ({snomed_code_list});')

snomed_icd9_query <- glue('SELECT
  icd9.icd9_code,
  icd9.snomed_code,
  icd9.snomed_name
  FROM arcus.icd9_to_snomed AS icd9
  WHERE snomed_code IN ({snomed_code_list});')

snomed_icd10 <- bq_project_query(proj_id, snomed_icd10_query) %>% bq_table_download(., page_size=3500)
print(paste0("arcus.icd10_to_snomed table with 25.11 pat_ids, dimensions: ", dim(snomed_icd10)))

snomed_icd9 <- bq_project_query(proj_id, snomed_icd9_query) %>% bq_table_download(., page_size=3500)
print(paste0("arcus.icd9_to_snomed table with 25.11 pat_ids, dimensions: ", dim(snomed_icd9)))

indication <- annots %>% filter(type == "Indication")

length(unique(indication$proc_ord_id)) #712 reports


# match up with ICD codes (and dx_id from arcus)

indication_icd10 <- merge(indication, snomed_icd10 %>%
                select(-c("snomed_name")), by = "snomed_code", all.x = TRUE, all.y = FALSE)
indication_icd9 <- merge(indication, snomed_icd9 %>%
                           select(-c("snomed_name")),
                    by = "snomed_code", all.x = TRUE, all.y = FALSE)

write.csv(indication_icd10, file = glue("{out_folder}indication_icd10.csv"))
write.csv(indication_icd9, file = glue("{out_folder}indication_icd9.csv"))

