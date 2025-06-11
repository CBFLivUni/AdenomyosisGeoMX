# Load the quarto library
library(quarto)

# Render all the documents in order
purrr::map(list.files(pattern = "\\.qmd$")[1:12], ~quarto::quarto_render(.x))


