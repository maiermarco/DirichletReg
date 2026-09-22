# BayesSpace ()

* GitHub: <https://github.com/maiermarco/DirichletReg>
* Email: <mailto:marco_maier@posteo.de>

Run `revdepcheck::revdep_details(, "BayesSpace")` for more info

## Error before installation

### Devel

```

  There are binary versions available but the source versions are later:
                     binary  source needs_compilation
AnnotationDbi        1.74.0  1.75.2             FALSE
assorthead            1.6.3  1.7.10             FALSE
BiocBaseUtils        1.14.2  1.15.1             FALSE
BiocFileCache         3.2.0   3.3.0             FALSE
BiocGenerics         0.58.1 0.59.12             FALSE
BiocIO               1.22.0  1.23.3             FALSE
biocmake              1.4.0   1.5.1             FALSE
...
Error: package or namespace load failed for 'AnnotationHub' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
 there is no package called 'AnnotationDbi'
In addition: Warning messages:
1: multiple methods tables found for 'transform' 
2: replacing previous import 'BiocGenerics::transform' by 'S4Vectors::transform' when loading 'AnnotationHub' 
Execution halted
Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) : 
  namespace 'S4Vectors' 0.50.3 is being loaded, but >= 0.51.8 is required
Calls: <Anonymous> ... withCallingHandlers -> loadNamespace -> namespaceImport -> loadNamespace
Execution halted


installing the source packages 'AnnotationDbi', 'assorthead', 'BiocBaseUtils', 'BiocFileCache', 'BiocGenerics', 'BiocIO', 'biocmake', 'BiocVersion', 'ComplexHeatmap', 'DelayedArray', 'dir.expiry', 'ExperimentHub', 'GenomicRanges', 'KEGGREST', 'MatrixGenerics', 'ScaledMatrix', 'scater', 'Seqinfo', 'SingleCellExperiment', 'SpatialExperiment', 'spatialLIBD', 'SummarizedExperiment', 'tinytex'

* installing *source* package 'assorthead' ...
** this is package 'assorthead' version '1.7.10'
** package 'assorthead' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
...
ERROR: dependencies 'SpatialExperiment', 'scater', 'ExperimentHub', 'SummarizedExperiment', 'SingleCellExperiment', 'GenomicRanges' are not available for package 'spatialLIBD'
Perhaps try a variation of:
install.packages(c('SpatialExperiment', 'scater', 'ExperimentHub', 'SummarizedExperiment', 'SingleCellExperiment', 'GenomicRanges'))
* removing 'C:/********/DirichletReg/revdep/library/BayesSpace/spatialLIBD'

The downloaded source packages are in
	'C:\********/Temp\RtmpALYwkb\downloaded_packages'
Error in (function (libdir, packages, quiet, repos)  : 
  all(packages %in% rownames(installed.packages(libdir[1]))) is not TRUE
In addition: There were 24 warnings (use warnings() to see them)


```
### CRAN

```

  There are binary versions available but the source versions are later:
                     binary  source needs_compilation
AnnotationDbi        1.74.0  1.75.2             FALSE
assorthead            1.6.3  1.7.10             FALSE
BiocBaseUtils        1.14.2  1.15.1             FALSE
BiocFileCache         3.2.0   3.3.0             FALSE
BiocGenerics         0.58.1 0.59.12             FALSE
BiocIO               1.22.0  1.23.3             FALSE
biocmake              1.4.0   1.5.1             FALSE
...
Error: package or namespace load failed for 'AnnotationHub' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
 there is no package called 'AnnotationDbi'
In addition: Warning messages:
1: multiple methods tables found for 'transform' 
2: replacing previous import 'BiocGenerics::transform' by 'S4Vectors::transform' when loading 'AnnotationHub' 
Execution halted
Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) : 
  namespace 'S4Vectors' 0.50.3 is being loaded, but >= 0.51.8 is required
Calls: <Anonymous> ... withCallingHandlers -> loadNamespace -> namespaceImport -> loadNamespace
Execution halted


installing the source packages 'AnnotationDbi', 'assorthead', 'BiocBaseUtils', 'BiocFileCache', 'BiocGenerics', 'BiocIO', 'biocmake', 'BiocVersion', 'ComplexHeatmap', 'DelayedArray', 'dir.expiry', 'ExperimentHub', 'GenomicRanges', 'KEGGREST', 'MatrixGenerics', 'ScaledMatrix', 'scater', 'Seqinfo', 'SingleCellExperiment', 'SpatialExperiment', 'spatialLIBD', 'SummarizedExperiment', 'tinytex'

* installing *source* package 'assorthead' ...
** this is package 'assorthead' version '1.7.10'
** package 'assorthead' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
...
ERROR: dependencies 'SpatialExperiment', 'scater', 'ExperimentHub', 'SummarizedExperiment', 'SingleCellExperiment', 'GenomicRanges' are not available for package 'spatialLIBD'
Perhaps try a variation of:
install.packages(c('SpatialExperiment', 'scater', 'ExperimentHub', 'SummarizedExperiment', 'SingleCellExperiment', 'GenomicRanges'))
* removing 'C:/********/DirichletReg/revdep/library/BayesSpace/spatialLIBD'

The downloaded source packages are in
	'C:\********/Temp\RtmpALYwkb\downloaded_packages'
Error in (function (libdir, packages, quiet, repos)  : 
  all(packages %in% rownames(installed.packages(libdir[1]))) is not TRUE
In addition: There were 24 warnings (use warnings() to see them)


```
