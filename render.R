#!/usr/bin/env Rscript

library(rmarkdown)

# --------------------------------------------------
# Rscript render.R input.Rmd output_dir centile_path [output_file_base]
#
# --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript render.R input.Rmd output_dir centile_path [output_file_base]")
}

input_rmd <- args[1]

message("Rendering: ", input_rmd)

# check output directory
output_dir <- args[2]
if(!dir.exists(output_dir)) {
message("Output directory ", output_dir, " does not exist. Check path and spelling.")
} else{message("Saving all output in subfolders under directory: ", output_dir)}

# make subfolder in output directory to save markdown output
markdown_output_dir <- paste0(output_dir, "markdown_output/")
if (!dir.exists(markdown_output_dir)) dir.create(markdown_output_dir)

# assign variable to centiles if available
centiles_data_path <- args[3]

message("Using centile file: ", centiles_data_path)

# default file base = input filename without extension, if no file name supplied
default_name <- tools::file_path_sans_ext(basename(centiles_data_path))
output_filename <- if (length(args) > 3) args[4] else default_name

# add date to output filename
output_file <- paste0(output_filename, "_", Sys.Date(), ".html")

# Render R Markdown
  message("Params supplied: ", output_dir, centiles_data_path)
  rmarkdown::render(
    input       = input_rmd,
    output_format = "html_document",
    output_file = output_file,
    output_dir  = markdown_output_dir,
    clean       = TRUE,
    params = list(output_folder = output_dir,
                  centiles_path = centiles_data_path)
  )

message("Rendered file saved to: ", file.path(markdown_output_dir, output_file))
