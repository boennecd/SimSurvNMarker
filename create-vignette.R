library(rmarkdown)
(function(){
  setwd("vignettes")

  render(input = "SimSurvNMarker.Rmd", output_format = "html_document")
  render(input = "SimSurvNMarker.Rmd", output_format = github_document(
    math_method = "default", html_preview = FALSE),
    output_file = "README.md")
  file.copy("README.md", file.path("..", "README.md"), overwrite = TRUE)
  unlink("README.md")
  unlink("README.html")
  fs <- list.files(file.path("man", "figures"))
  for(f in fs)
    file.copy(file.path("man", "figures", f),
              file.path("..", "man", "figures", f), overwrite = TRUE)
  unlink("man", recursive = TRUE)
})()
