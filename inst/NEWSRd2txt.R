library(tools)

Rd2txt_options(
  width            = 80L,
  minIndent        = 2L,
  extraIndent      = 2L,
  sectionIndent    = 2L,
  sectionExtra     = 2L,
  underline_titles = FALSE,
  unicode_symbols  = TRUE,
  itemBullet       = "• "
)

Rd2txt("NEWS.Rd", out="../NEWS", package = "DirichletReg", outputEncoding = "UTF-8")

system("R CMD Rd2pdf --no-preview --encoding=UTF-8 --force NEWS.Rd")
compactPDF("NEWS.pdf", gs_quality="printer", verbose = TRUE)
