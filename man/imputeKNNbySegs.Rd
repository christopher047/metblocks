% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block_methods.R
\name{imputeKNNbySegs}
\alias{imputeKNNbySegs}
\title{getimputeKNNbySegs}
\usage{
imputeKNNbySegs(object, ncores = 16, nb = 5)
}
\arguments{
\item{object}{an object of meth class}

\item{ncores}{the number of cores to use}

\item{nb}{the number of neighhours used in KNN}
}
\value{
an object of meth class
}
\description{
This function imputes missing values by KNN per segment
Stores the imputed relative methylation values (raw_meth/raw_count) in object@imp_mat
}
