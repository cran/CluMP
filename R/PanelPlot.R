#' Plot Micro-Panel (longitudinal) Data
#'
#' This function plots micro-panel (longitudinal) data from stored \code{\link[base]{data.frame}} or randomly generated panel data from \code{\link{GeneratePanel}} function.
#' @param data A data frame containing the variables named in \code{formula} and \code{group} arguments.
#' @param formula A two-sided \code{\link[stats]{formula}} object, with a numeric, clustering variable (Y) on the left of a ~ separator and the time (numeric) variable on the right. Time is measured from the start of the follow-up period (baseline).
#' @param group A grouping factor variable (vector), i.e. single identifier for each (trajectory).
#' @param colour Character, which is a variable's name in data. The trajectories are distinguished by colour according to this variable.
#' @param show_legend Logical scalar. It indicates whether to show cluster legend. Default is \emph{TRUE}.
#' @param title String. Is an optional title for a plot. Otherwise no title will used.
#' @param x_title String. Is an optional title for x axis. Otherwise variable name after ~ in \code{formula} will used.
#' @param y_title String. Is an optional title for y axis. Otherwise variable name before ~ in \code{formula} will used.
#' @keywords CLUMP
#' @return Returns plot using package ggplot2.
#' @export
#' @import dplyr ggplot2
#' @examples
#' set.seed(123)
#' dataMale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
#' dataMale$Gender <- "M"
#' dataFemale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
#' dataFemale$ID <- dataFemale$ID + 50
#' dataFemale$Gender <- "F"
#' data <- rbind(dataMale, dataFemale)
#'
#' PanelPlot(data = data, formula = Y ~ Time, group = "ID", colour = "Gender")
#'

PanelPlot <- function(data, formula = Y ~ Time, group = "ID", colour = NA,
                      show_legend = TRUE,
                      title = NULL, x_title = NULL, y_title = NULL){

  # library(dplyr)
  # library(ggplot2)
  # define global variables
  CluMP_ID = CluMP_X1 = CluMP_Y = ID = Visit = X1 = X1_ann = Y = Y.x = Y.y = Y_ci =
    abs_angle_radian = abs_change = abs_change_ann = angle_radian = best = bestVal =
    cluster = cos_denom = cos_nom  = cosinus = f_up = mean_Time = mean_Y =
    memb_CluMP = nVisit = number = obsah_trojuh = sd_Y = slope =
    slope_first_last = timepoint = value = . = NULL

  # Modelframe from formula
  mf <- stats::model.frame(formula = formula, data = data, na.action = NULL)
  variables <- names(mf)
  names(mf) <- c("Y", "X1")
  mf <- cbind(mf, "ID" = data[, group])
  if (!is.numeric(mf$Y)){
    stop('Clustering variable should be numeric')
  }
  # mf_complete <- mf[complete.cases(mf),]
  # mf_complete <- mf_complete %>% group_by(ID) %>% mutate(count = n()) %>% filter(count > 2)
  # mf_non_complete <- mf %>% filter(!(ID %in% mf_complete$ID))
  # mf <- mf_complete
  #
  # if(nrow(mf_non_complete) > 1) {
  #   warning('Some of objects contains NA values.', call. = FALSE)
  # }
  if(!is.na(colour) & !(colour %in% colnames(data))){
    warning(paste0('Trere is no colour variable in the data set named ', colour, '. colour by ID variable.'), call. = FALSE)
    colour <- NA

  }
  if(!is.na(colour)){
    if(length(unique(data[,colour])) > length(unique(data[,group]))) {
    warning(paste0('You are going to assign colour by column ', colour, '. This column has more unique objects than number of observations. colour by ID variable.'), call. = FALSE)
    colour <- "ID"
    }
  }


  if (!is.na(colour)) {
    Cluster = data %>% dplyr::filter((!!as.name(group)) %in% unique(mf$ID)) %>% dplyr::select(colour)
    mf <- cbind(mf, Cluster = Cluster[, colour])
    plotMP <- ggplot(data = mf, aes(x = X1, y = Y, group = ID)) +
      geom_line(aes(colour = factor(Cluster))) +
      labs(x = variables[2], y = variables[1]) +
      scale_colour_discrete(name = colour, labels = unique(mf$Cluster)) +
      theme(legend.position = "right") +
      theme_classic()
    plotMP <- plotMP +
      if(!show_legend) {guides(colour=FALSE)}
  }else {
    plotMP <- ggplot(data = mf, aes(x = X1, y = Y, group = ID)) +
      geom_line(aes(colour = factor(ID))) +
      labs(x = variables[2], y = variables[1]) +
      scale_colour_discrete(labels = unique(mf$ID)) +
      theme_classic() +
      guides(colour=FALSE)
  }


  # add other optional ggplot objects
  plotMP <- plotMP +
    if(!is.null(title)) {ggtitle(title)}
  plotMP <- plotMP +
    if(!is.null(x_title)) {xlab(x_title)}
  plotMP <- plotMP +
    if(!is.null(y_title)) {ylab(y_title)}
  # Adjust font
  plotMP <- plotMP +
    theme(text = element_text(size=16))

  return(plotMP)

}

