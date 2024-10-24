% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CGislet_annot.R
\name{get_methy_level}
\alias{get_methy_level}
\title{Converts genome coordinates from result file of LiftOver}
\usage{
get_methy_level(CGislet, methy_dat, methy_genome, methy_levels)
}
\arguments{
\item{CGislet}{GRange object}

\item{liftover_file}{character;bed filepath obtained from LiftOver software}
}
\value{
position box ("chr4:100,001-100,001") of converted genome coordinates
}
\description{
Converts genome coordinates from result file of LiftOver
}
