% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clarity2sims.R
\name{update_model}
\alias{update_model}
\title{update_model}
\usage{
update_model(mod, moddat, ...)
}
\arguments{
\item{mod}{A compiled cmdstanr stanmodel}

\item{moddat}{Model data for use with mod}

\item{...}{Additional arguments to sample}
}
\value{
Posterior draws
}
\description{
Update the model used in clarity2sims
}
