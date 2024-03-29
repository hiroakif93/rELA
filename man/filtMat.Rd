% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtMat.R
\name{filtMat}
\alias{filtMat}
\title{Filtering species}
\usage{
filtMat(x = data, th = 10, nsample = nrow(x), minth = 0.05, maxth = 0.95)
}
\arguments{
\item{th}{: th is an abundance by given an integer.}

\item{nsample}{; nsample is total sample.}

\item{minth}{: max.itr is minimum rate of occurence in samples.}

\item{data}{: data is a matrix. Row is sample, column is species.}
}
\description{
Filtering low abundance/occurence species
}
