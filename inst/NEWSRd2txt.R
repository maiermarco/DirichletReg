library(tools)

package_grep <- function(regex){ grep("revdep", dir("..", pattern = regex, recursive = TRUE, full.names = TRUE), invert = TRUE, value = TRUE) }
.aspell_control <- function(wordlist) c("--master=en_US", if(missing(wordlist)) NULL else sprintf("--add-wordlists=%s", wordlist))

aspell_R   <- aspell(package_grep("\\.R$"), filter = "R", control = .aspell_control("WORDLIST_R.txt"))
aspell_Rd  <- aspell(package_grep("\\.Rd$"), filter = list("Rd", drop = "\\references"), control = .aspell_control("WORDLIST_Rd.txt"))
#aspell_pot <- aspell(package_grep("\\.pot$"), filter = "pot", control = .aspell_control())
aspell_dcf <- aspell("../DESCRIPTION", filter = "dcf", control = .aspell_control())
#aspell_md  <- aspell(package_grep("\\.R?md$"), filter = "md", control = .aspell_control())

print(eval(str2expression(sprintf("rbind(%s)", paste(grep("^aspell\\_(R|Rd|pot|dcf|md)$", ls(), value = TRUE), collapse = ", ")))), verbose = TRUE)



# NEWS.Rd --> ../NEWS.md
parsed_rd <- rd2markdown::get_rd(file = "NEWS.Rd")
parsed_rd <- parsed_rd[!sapply(parsed_rd, function(x){ x[[1]] == "UTF-8" })]
#cat(rd2markdown::rd2markdown(parsed_rd, fragments = c()), file = "../NEWS.md")



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
