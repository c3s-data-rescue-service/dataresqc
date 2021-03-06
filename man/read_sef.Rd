% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{read_sef}
\alias{read_sef}
\title{Read data files in Station Exchange Format version 1.0.0}
\usage{
read_sef(file = file.choose(), all = FALSE)
}
\arguments{
\item{file}{Character string giving the path of the SEF file.}

\item{all}{If FALSE (the default), omit the columns 'Period' and 'Meta'
(also 'Hour' and 'Minute' for non-instantaneous data)}
}
\value{
A data frame with up to 9 variables, depending on whether
\code{all} is set to TRUE.
The variables are: variable code, year, month, day, hour, minute,
value, period, metadata.
}
\description{
Read data files in Station Exchange Format version 1.0.0
}
\author{
Yuri Brugnara
}
