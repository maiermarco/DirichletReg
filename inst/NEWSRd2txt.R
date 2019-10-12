library("tools")
Rd2txt_options(
  width            = 80L,
  minIndent        = 2L,
  extraIndent      = 2L,
  sectionIndent    = 2L,
  sectionExtra     = 2L,
  underline_titles = FALSE
)

Rd2txt("NEWS.Rd", out="../NEWS")

system("R CMD Rd2pdf --no-preview --encoding=UTF-8 --force NEWS.Rd")

tools:::compactPDF("NEWS.pdf", gs_quality="screen")

#system("gswin64c -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dEmbedAllFonts=true -dPDFSETTINGS=/screen -sColorConversionStrategy=/RGB -sOutputICCProfile=srgb.icc -dCompatibilityLevel=1.5 -dAutoRotatePages=/None -sOutputFile=NEWSx.pdf NEWS.pdf")

#system("qpdf --stream-data=compress --object-streams=generate NEWSx.pdf NEWS.pdf")

#file.remove("NEWSx.pdf")
