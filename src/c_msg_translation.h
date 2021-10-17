// enable message-translation (just in case)
// see 1.8.1 C-level messages
// https://cran.r-project.org/doc/manuals/r-release/R-exts.html#C_002dlevel-messages

// for simple messages and singular/plural forms:
#ifdef ENABLE_NLS
#include<libintl.h>
#define _(String) dgettext("DirichletReg", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif
