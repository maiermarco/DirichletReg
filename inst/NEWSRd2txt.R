library(tools)

old_op <- options(
  useFancyQuotes = "UTF-8",
  encoding       = "UTF-8"
)

opts_rd2txt <- Rd2txt_options(
  width            = 80L,
  minIndent        = 2L,
  extraIndent      = 2L,
  sectionIndent    = 2L,
  sectionExtra     = 2L,
  underline_titles = FALSE,
  unicode_symbols  = TRUE,
  itemBullet       = "• "
)

Rd2txt("NEWS.Rd", out="../NEWS", outputEncoding = "UTF-8", options = opts_rd2txt)

system("R CMD Rd2pdf --no-preview --encoding=UTF-8 --force NEWS.Rd")

compactPDF("NEWS.pdf", gs_quality="screen", verbose = TRUE)

options(old_op)

