## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE,include=FALSE------------------------------------------

## ---- echo=FALSE,results="asis"-----------------------------------------------
files <- list.files("web", "\\.Rmd$")
for (file in files) {
  lines <- readLines(file.path("web", file), n = 10)
  title_idx <- grep("^title: ", lines)
  if (length(title_idx) > 0) {
    title <- sub("^title: ", "", lines[title_idx[1]])
    title <- sub('^"', "", title)
    title <- sub('"$', "", title)
    cat("* [",
      title,
      "]( https://eliaskrainski.github.io/INLAspacetime/vignettes/",
      sub("\\.Rmd", ".html", 
          sub("/web", "", file)),
      ")\n",
      sep = ""
    )
  }
}

