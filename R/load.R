# Common setup for every script: run from the project root with
#   source('R/load.R')
# Loads the core libraries, the paths (config/paths.R) and every function file in R/.

if (!file.exists('config/paths.R')) {
    stop('Run the scripts from the project root (the folder that contains config/ and R/).')
}

suppressPackageStartupMessages({
    library(tidyverse)
})

source('config/paths.R')

for (f in setdiff(list.files('R', pattern = '\\.R$', full.names = TRUE), 'R/load.R')) {
    source(f)
}
rm(f)

# Save a ggplot as pdf + png with the same name in results/figures
# (pdf_device = cairo_pdf where the original script used it)
save_figure <- function(plot, name, width, height, units = 'in', dpi = 300, pdf_device = NULL) {
    out <- p('figures', name)
    ggsave(str_c(out, '.pdf'), plot, width = width, height = height, units = units, device = pdf_device)
    ggsave(str_c(out, '.png'), plot, width = width, height = height, units = units, dpi = dpi, bg = 'white')
    invisible(out)
}
