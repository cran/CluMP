% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump_view.R
\name{CluMP_view}
\alias{CluMP_view}
\title{Cluster profiles' (CluMP results) visualisation}
\usage{
CluMP_view(
  CluMPoutput,
  type = "all",
  nb_intervals = NULL,
  return_table = FALSE,
  title = NULL,
  x_title = NULL,
  y_title = NULL,
  plot_NA = FALSE
)
}
\arguments{
\item{CluMPoutput}{An object (output) from the \code{\link{CluMP}} function.}

\item{type}{String. Indicates which type of graph is required. Possible values for this argument are: \emph{"all"} (plots all data with non-parametric mean trajectories), \emph{"cont"} (only non-parametric mean trajectories) or \emph{"breaks"} (mean trajectories with error bars).}

\item{nb_intervals}{An integer, positive number (scalar) specifying the number of regular timepoints into which should be follow-up period split. This argument works only with graph type = \emph{"breaks"}. In case of other graph types the argument is ignored. The number of error bars is equal to the number of timepoints specified by this argument.}

\item{return_table}{Logical scalar indicating if the summary table of plotted values in the graph of type = \emph{"breaks"} should be returned. Default is \emph{FALSE}.}

\item{title}{String. Optional title for a plot. If undefined, no title will used.}

\item{x_title}{String. An optional title for \emph{x} axis. If undefined, the variable name after ~ in \code{formula} will used.}

\item{y_title}{String. An optional title for \emph{y} axis. If undefined, the variable name before ~ in \code{formula} will used.}

\item{plot_NA}{Plot \emph{NA} cluster if exists. Default is \emph{FALSE}. \emph{NA} cluster gathers improper individuals (< 3 observations) for longitudinal clustering.}
}
\value{
Returns graph for type \emph{"all"} and \emph{"cont"} or (list with) graph and table of mean trajectories (if specified) for type = \emph{"breaks"}.
}
\description{
This graphical function enables to visualise cluster profiles (mean representatives of each cluster). Available are three types of plots: non-parametric (LOESS method for small/medium or GAM method for complex data of large size. Both methods are applied from ggplot2 representatives (mean within-cluster trajectories) with/without all individual (original) trajectories, and nonparametric mean trajectories with error bars.
}
\examples{
set.seed(123)
dataMale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
dataMale$Gender <- "M"
dataFemale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
dataFemale$ID <- dataFemale$ID + 50
dataFemale$Gender <- "F"
data <- rbind(dataMale, dataFemale)

CluMPoutput <- CluMP(formula = Y ~ Time, group = "ID", data = data, cl_numb = 3)
title <- "Plotting clusters' representatives with error bars"
CluMP_view(CluMPoutput, type = "all" , return_table = TRUE)
CluMP_view(CluMPoutput, type = "cont")
CluMP_view(CluMPoutput, type = "breaks", nb_intervals = 5, return_table=TRUE, title = title)

}
\keyword{CluMP}
