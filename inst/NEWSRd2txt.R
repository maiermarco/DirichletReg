library(tools)

# Spell-checking NEWS.Rd
spelling::spell_check_files(path = "NEWS.Rd", ignore = spelling:::read_wordfile(wordfile = "WORDLIST"), lang = "en_US")

# NEWS.Rd --> ../NEWS.md
parsed_rd <- rd2markdown::get_rd(file = "NEWS.Rd")
parsed_rd <- parsed_rd[!sapply(parsed_rd, function(x){ x[[1]] == "UTF-8" })]
cat(rd2markdown::rd2markdown(parsed_rd, fragments = c()), file = "../NEWS.md")

# NEWS.Rd --> ../NEWS
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

# NEWS.Rd --> NEWS.pdf
system("R CMD Rd2pdf --no-preview --encoding=UTF-8 --force NEWS.Rd")
compactPDF("NEWS.pdf", gs_quality="printer", verbose = TRUE)
