#' Cluster profiles' (CluMP results) visualisation
#'
#' This graphical function enables to visualise cluster profiles (mean representatives of each cluster). Available are three types of plots: non-parametric (LOESS method for small/medium or GAM method for complex data of large size. Both methods are applied from ggplot2 representatives (mean within-cluster trajectories) with/without all individual (original) trajectories, and nonparametric mean trajectories with error bars.
#' @param CluMPoutput An object (output) from the \code{\link{CluMP}} function.
#' @param type String. Indicates which type of graph is required. Possible values for this argument are: \emph{"all"} (plots all data with non-parametric mean trajectories), \emph{"cont"} (only non-parametric mean trajectories) or \emph{"breaks"} (mean trajectories with error bars).
#' @param nb_intervals An integer, positive number (scalar) specifying the number of regular timepoints into which should be follow-up period split. This argument works only with graph type = \emph{"breaks"}. In case of other graph types the argument is ignored. The number of error bars is equal to the number of timepoints specified by this argument.
#' @param return_table Logical scalar indicating if the summary table of plotted values in the graph of type = \emph{"breaks"} should be returned. Default is \emph{FALSE}.
#' @param title String. Optional title for a plot. If undefined, no title will used.
#' @param x_title String. An optional title for \emph{x} axis. If undefined, the variable name after ~ in \code{formula} will used.
#' @param y_title String. An optional title for \emph{y} axis. If undefined, the variable name before ~ in \code{formula} will used.
#' @param plot_NA Plot \emph{NA} cluster if exists. Default is \emph{FALSE}. \emph{NA} cluster gathers improper individuals (< 3 observations) for longitudinal clustering.
#' @keywords CluMP
#' @return Returns graph for type \emph{"all"} and \emph{"cont"} or (list with) graph and table of mean trajectories (if specified) for type = \emph{"breaks"}.
#' @export
#' @import ggplot2
#' @import data.table
#' @rawNamespace import(dplyr, except = c(last, first, between))
#' @examples
#' set.seed(123)
#' dataMale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
#' dataMale$Gender <- "M"
#' dataFemale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
#' dataFemale$ID <- dataFemale$ID + 50
#' dataFemale$Gender <- "F"
#' data <- rbind(dataMale, dataFemale)
#'
#' CluMPoutput <- CluMP(formula = Y ~ Time, group = "ID", data = data, cl_numb = 3)
#' title <- "Plotting clusters' representatives with error bars"
#' CluMP_view(CluMPoutput, type = "all" , return_table = TRUE)
#' CluMP_view(CluMPoutput, type = "cont")
#' CluMP_view(CluMPoutput, type = "breaks", nb_intervals = 5, return_table=TRUE, title = title)
#'
CluMP_view <- function(CluMPoutput, type = "all", nb_intervals = NULL, return_table = FALSE,
                       title = NULL, x_title = NULL, y_title = NULL, plot_NA = FALSE) {

  # define global variables
  CluMP_ID = CluMP_X1 = CluMP_Y = ID = Visit = X1 = X1_ann = Y = Y.x = Y.y = Y_ci =
    abs_angle_radian = abs_change = abs_change_ann = angle_radian = best = bestVal =
    cluster = cos_denom = cos_nom  = cosinus = f_up = mean_Time = mean_Y =
    memb_CluMP = nVisit = number = obsah_trojuh = sd_Y = slope =
    slope_first_last = timepoint = value = . = .. = ..colour =
    ..cols = ..cont_vars = ..group = ..scale_cols = Time = NULL


  if(!(type %in% c("all", "cont", "breaks"))) {
    stop('Type should be one of "all", "cont" or "breaks".')
  }

  if(type != "breaks" & !is.null(nb_intervals)) {
    warning('Argument nb_intervals is not possible to combine with types "all" and "cont". nb_intervals ignored.')
  }

  cols <- c(CluMPoutput$variables, CluMPoutput$group, "memb_CluMP")
  PlotData <- setDT(CluMPoutput$data)[, ..cols]
  PlotData$memb_CluMP <- as.factor(PlotData$memb_CluMP)

  if (!plot_NA) {
    PlotData <- setDT(PlotData)[!is.na(memb_CluMP),]
  }

  colnames(PlotData)[which(colnames(PlotData) %in% CluMPoutput$variables)] <- c("Y", "X1")
  colnames(PlotData)[which(colnames(PlotData) %in% CluMPoutput$group)] <- "ID"


  if (type == "all") {

    graf <- ggplot(data = PlotData, aes(x = X1, y = Y, group = ID, colour = memb_CluMP)) +
      geom_line(size = .8, alpha = .5) +
      geom_smooth(aes(x = X1, y = Y, group = memb_CluMP, colour = memb_CluMP), size = 2, method = "auto") +
      scale_colour_discrete(name = "cluster") +
      theme_bw() +
      labs(y = "mean, CI", x = CluMPoutput$variables[2])

  }else if (type == "cont") {

    graf <- ggplot(data = PlotData, aes(x = X1, y = Y, group = memb_CluMP, colour = memb_CluMP)) +
      geom_smooth(method = "auto") +
      scale_colour_discrete(name = "cluster") +
      theme_bw() +
      labs(y = "mean, CI", x = CluMPoutput$variables[2])

  }else if (type == "breaks") {

    if (is.null(nb_intervals)) {
      stop("numb of intervals should be defined.")
    }else if (!is.numeric(nb_intervals) | nb_intervals <= 0) {
      stop("numb of intervals should be a positive number.")
    }

    time_point_intervals <- seq(min(PlotData$X1), max(PlotData$X1), length.out = nb_intervals)
    # create data frame
    PlotData <- setDT(PlotData)[, timepoint := findInterval(X1, time_point_intervals),]
    PlotData <- setDT(PlotData)[, ':='(
      mean_Y = mean(Y, na.rm = TRUE),
      sd_Y = stats::sd(Y, na.rm = TRUE),
      n = .N,
      mean_Time = mean(X1, na.rm = TRUE)
    ), by = .(memb_CluMP, timepoint)]
    PlotData <- setDT(PlotData)[, .(memb_CluMP, timepoint, mean_Y, sd_Y, n, mean_Time)]
    PlotData <- setDT(PlotData)[!duplicated(setDT(PlotData))]
    PlotData <- setDT(PlotData)[, Y_ci := 1.96*sd_Y/sqrt(n), ]

    # create graph
    # * na.rm = T -- suppress warnings
    graf <- ggplot(data = PlotData, aes(x = mean_Time, y = mean_Y, group = memb_CluMP, colour = memb_CluMP)) +
      geom_errorbar(aes(ymin = mean_Y - Y_ci, ymax = mean_Y + Y_ci), width=.1, size = .9, alpha = 0.5, na.rm = TRUE) +
      geom_line(size = .9) +
      geom_point(size = 3, shape = 21, fill = "white") +
      scale_colour_discrete(name = "cluster") +
      theme_bw() +
      scale_x_continuous(name = CluMPoutput$variables[2]) +
      labs(y = "mean, CI", x = CluMPoutput$variables[2])

  }

  # add other optional ggplot objects
  graf <- graf +
    if(!is.null(title)) {ggtitle(title)}
  graf <- graf +
    if(!is.null(x_title)) {xlab(x_title)}
  graf <- graf +
    if(!is.null(y_title)) {ylab(y_title)}
  # Adjust font
  graf <- graf +
    theme(text = element_text(size=16))


  if (type == "breaks" & return_table == TRUE) {
    output <- list(graph = graf, table = as.data.frame(PlotData))
  }else if (type == "cont" | type == "all" | (type == "breaks" & return_table == FALSE)) {
    output <- graf
  }

  return(output)

}
