#!/usr/bin/env Rscript
##This is a template script to import SQL tables into R and perform 
##simple data operations such as merge/filter/summary/frequency table.

#Load packages. Install these packages with install.packages("pkg_name") if they are not installed.
library(bigrquery) #package to import SQL tables into R
library(tidyverse) #package suite for data manipulation, reading data into R, and plotting (ggplot2)
library(glue) #package to use expressions with strings

#Set up bigrquery to import tables using SQL
bq_auth()
proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1] #assign project ID

#Examples to import SQL tables

#Download SLIP project table
query <- glue("SELECT * FROM lab.proc_ord_projects") #the string here is the SQL query
projects <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500) #assign the name of the dataframe in R here
print(paste0("Projects table, dimensions: ", dim(projects))) #pastes the dimensions of the table downloaded

#Download SLIP grades table
query <- glue("SELECT * FROM lab.bak_grader_table_with_metadata_project_independent")
grades <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500)
print(paste0("Grades table, dimensions: ", dim(grades)))

#inspect names of variables in the tables
names(projects)
names(grades)

#see frequency table for a variable
table(projects$project)
#see unique levels of a variable
unique(grades$grader_name)
#get number of unique grader names
length(unique(grades$grader_name))

#merge based on common variable
merged <- merge(projects, grades, by = "proc_ord_id", all.x = TRUE)

#filter table for specific values of a variable (in this case for names of prematurity projects)

prematurity <- merged %>% filter(project %in% c("SLIP Infants Very Premature", 
                                                  "SLIP Infants Mod Premature", 
                                                  "SLIP Infants At Term"))



#filter based on grading criteria, graders (not coarse text or NLP models) and no grades that are 999 or 503
filtered <- prematurity
  filter(grade_criteria == "SLIP") %>% 
  filter(!(grader_name %in% c("NLP Models 2024-11-12", "NLP Models 2025-03-03", "Coarse Text Search 2023-07-24", 
                              "Coarse Text Search 2024-07-31", "Coarse Text Search 2024-01-29", "Coarse Text Search 2023-11-02",
                              "Coarse Text Search 2023-07-21"))) %>% filter(!(grade %in% c(999, 503)))
  
  
  
##You can also import arcus.patient table and merge by pat_id to include/inspect information about patients such as gestational age etc.

