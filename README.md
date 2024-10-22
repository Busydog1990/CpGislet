remote::install_github("Busydog1990/CpGislet")

require(CpGislet)

You can see workflow of CpGislet in Busydog1990/CpGislet_workflow

If you have any question, contact zzhe@mail.hzau.edu.cn


Package: CpGislet

Type: Package

Title: Detective of CpG islets in the genome and predict their structures

Version: 0.1.7

Author: person(given = 'Zhe', family = 'Zhang', email = 'zzhe@mail.hzau.edu.cn', 
               role = 'aut', comment = c(ORCID = '0000-0001-5265-6069'))
               
Maintainer: Zhang, Zhe <mail@mail.hzau.edu.cn>

Description: A tool for predicting CpG islet structure in the genome. Also including annotation, classification, and statistics of CpG islet.

Depends: R (>= 4.0.0), Rcpp (>= 1.0.4), Matrix (>= 1.3.3), Biostrings,GenomicRanges, stringr, ChIPseeker, GenomicFeatures, mlr3verse, data.table, txdbmaker, clusterProfiler, data.table

License: GPL-3

Encoding: UTF-8

RoxygenNote: 7.3.2

Suggests: rmarkdown, knitr, ggpubr, ggupset, ggplot2, patchwork

VignetteBuilder: knitr

NeedsCompilation: no

Packaged: 2024-10-22 02:15:51 UTC; admin
