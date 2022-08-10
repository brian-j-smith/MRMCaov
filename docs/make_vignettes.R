dest <- "vignettes"

out <- file(file.path(dest, "UserGuide.Rmd"), "w")

sourceLines <- function(file) {
  c(readLines(file.path("docs", "src", file)), "\n")
}

files <- c(
  "UserGuide.Rmd",
  "overview.Rmd",
  "using_intro.Rmd",
  "using_model.Rmd",
  "using_example.Rmd",
  "using_mrmc.Rmd",
  "using_srmc.Rmd",
  "using_stmc.Rmd",
  "using_curves.Rmd",
  "using_metrics.Rmd"
)

for (file in files) {
  writeLines(sourceLines(file), out)
}

writeLines("# References", out)

close(out)

file.copy("docs/src/bibliography.bib",
          file.path(dest, "bibliography.bib"),
          overwrite = TRUE)
