library(rmarkdown)
local({
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd("vignettes")

  render("SimSurvNMarker.Rmd", output_format = "html_document")
  render("SimSurvNMarker.Rmd", output_format = github_document(
    pandoc_args = "--webtex=https://latex.codecogs.com/svg.latex?"),
    output_file = "README.md")
  file.copy("README.md", file.path("..", "README.md"), overwrite = TRUE)
  unlink("README.md")
  unlink("README.html")
  unlink("man", recursive = TRUE)
})
