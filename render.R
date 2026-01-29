#!/usr/bin/env Rscript

library(rmarkdown)
library(dplyr)
library(knitr)

# --------------------------------------------------
# Rscript render.R input.Rmd output_directory save_markdown_flag save_markdown_name [params in order]
# markdown_name var only relevant if saving markdown output, but needs to be supplied regardless
# save_markdown has to be 0/1 (true/false for saving output of markdown render as an html)
# --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript render.R input.Rmd output_directory save_markdown_flag(0/1) save_markdown_name(without extension) [params in order]")
}

rmd_script <- args[1]

print(paste0("Rendering: ", rmd_script))

# check output directory
output_dir <- args[2]
if(!dir.exists(output_dir)) {
stop("Output directory ", output_dir, " does not exist. Check path and spelling.")
} else{print(paste0("Saving output in subfolders under directory: ", output_dir))}

# assign save_markdown flag & save_markdown_name args
save_markdown_flag = args[3]
save_markdown_name = args[4]

if(length(args) > 4){
  print("You provided custom parameters for the R script. \nCurrently this script only works if you provided ALL available parameters. \nIt will throw an error if not. \nMight try to make it more flexible in the future, but for now this is a safety feature so that defaults in the Rmd do not get used when not desired.")
  # read in params for input.Rmd
  param_args <- args[-c(1,2,3,4)]
  # read YAML front matter
  yaml <- rmarkdown::yaml_front_matter(rmd_script)
  rmd_params <- ifelse(!is.null(yaml$params), yaml["params"], list())
  rmd_params <- rmd_params[[1]]
  print(paste0("Param names in ", rmd_script, ": ", names(rmd_params)))
  
  # parse key=value args
  input_params <- lapply(param_args, function(x) {
    kv <- strsplit(x, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Params must be key=value")
    kv
  })
  names(input_params) <- unlist(lapply(input_params, function(x) x[1]))
  input_params <- lapply(input_params, function(x) x[2])
  
  # keep only params defined in YAML
  if(sum(names(input_params) %in% names(rmd_params)) != length(names(rmd_params))){
    stop(paste0("Names of input params do not match names of rmd params. \nInput params: ", 
                paste(names(input_params), collapse = ", "), "\nRmd params:", paste(names(rmd_params), collapse = ", ")))
  } else{params_pass <- input_params[names(input_params) %in% names(rmd_params)]}
  
  # basic type coercion based on YAML defaults
  for (nm in names(params_pass)) {
    default <- rmd_params[[nm]]
    if (is.numeric(default)) params_pass[[nm]] <- as.numeric(params_pass[[nm]])
    if (is.logical(default)) params_pass[[nm]] <- as.logical(params_pass[[nm]])
  }
  print(paste0("Params you are running ", rmd_script, " with: ", params_pass))
} else {
  print("You provided no parameters for the Rmd. /nThe Rmd script will be rendered using the defaults defined in the script.")
  # read YAML front matter
  yaml <- rmarkdown::yaml_front_matter(rmd_script)
  rmd_params <- ifelse(!is.null(yaml$params), yaml["params"], list())
  rmd_params <- rmd_params[[1]]
  print(paste0("Params you are running", rmd_script, " with (defaults) : ", rmd_params))
  params_pass <- NULL
}

# add date to output filename
output_file <- paste0(tools::file_path_sans_ext(rmd_script), "_", save_markdown_name, "_", Sys.Date(), ".html")

# handle markdown saving based on flag
if(save_markdown_flag == 1){
  # make subfolder in output directory to save markdown output
  markdown_output_dir <- paste0(output_dir, "markdown_output/")
  if (!dir.exists(markdown_output_dir)) dir.create(markdown_output_dir)
  
  # Render R Markdown
  rmarkdown::render(
    input       = rmd_script,
    output_format = "html_document",
    output_file = output_file,
    output_dir  = markdown_output_dir,
    clean       = TRUE,
    params = params_pass,
    quiet = FALSE
  )
  
  print("Markdown rendered, output file saved to: ", paste0(markdown_output_dir, output_file))
  
} else if (save_markdown_flag == 0){
  # Render R Markdown
  rmarkdown::render(
    input       = rmd_script,
    output_format = "html_document",
    output_file = output_file,
    clean = TRUE,
    params = params_pass,
    quiet = FALSE
  )
  file.remove(paste0(markdown_output_dir, output_file))
  print("Markdown rendered, output file deleted as per save_markdown_flag == 0.")
} else {
  stop(paste0("Save_markdown_flag (argument order three) must be 0/1 flag. You provided: ", save_markdown_flag))
}
