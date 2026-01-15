#!/usr/bin/env Rscript

library(rmarkdown)

# --------------------------------------------------
# Rscript render.R input.Rmd [output_dir] [output_file_base]
#
# Defaults if only input.Rmd is supplied:
#   output_dir = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/markdown_output/"
#   output_file_base = defaults to name of input Rmd
# --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript render.R input.Rmd [output_dir] [output_file_base]")
}

input_rmd <- args[1]

message("Rendering: ", input_rmd)

# default output directory if no output directory arg supplied
output_dir <- if (length(args) >= 2) args[2] else "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/markdown_output/"

# default file base = input filename without extension, if no file name supplied
default_name <- tools::file_path_sans_ext(basename(input_rmd))
output_filename <- if (length(args) >= 3) args[3] else default_name

# add date to output filename
output_file <- paste0(output_filename, "_", Sys.Date(), ".html")

# check output directory existence, if not create
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Render R Markdown
rmarkdown::render(
  input       = input_rmd,
  output_format = "html document"
  output_file = output_file,
  output_dir  = output_dir,
  clean       = TRUE
)

message("Rendered file saved to: ", file.path(output_dir, output_file))
