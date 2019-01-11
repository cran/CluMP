#' Finding an optimal number of clusters
#'
#' This function finds optimal number of clusters based on evaluation criteria (indices) available from the NbClust package.
#' @param formula A two-sided \code{\link[stats]{formula}} object, with a numeric, clustering variable (Y) on the left of a ~ separator and the time (numeric) variable on the right. Time is measured from the start of the follow-up period (baseline).
#' @param group A grouping factor variable (vector), i.e. single identifier for each individual (trajectory).
#' @param data A data frame containing the variables named in \code{formula} and \code{group} arguments.
#' @param index String vector of indices to be computed. Default is c("silhouette", "ch", "db"). See NbClust package for available indices and their description.
#' @param max_clust An integer, positive number (scalar) defining the maximum number of clusters to check. Default value of this argument is 10 or maximum number of individuals.
#' @param base_val Indicates whether include a value at zero time point as an additional clustering variable. Default is \emph{FALSE} and the standard number (7) of clustering parameters is used.
#' @keywords CluMP
#' @return Determine the optimal number of clusters, returns graphical output (red dot in plot indicates the recommended number of clusters according to that index) and table with indices.
#' @source Malika Charrad, Nadia Ghazzali, Veronique Boiteau, Azam Niknafs (2014). NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. Journal of Statistical Software, 61(6), 1-36. URL http://www.jstatsoft.org/v61/i06/.
#' @export
#' @import dplyr ggplot2 NbClust
#' @examples
#' data <- GeneratePanel(n = 100, Param = ParamLinear, NbVisit = 10)
#' OptiNum(data = data, formula = Y ~ Time, group = "ID")
#'

OptiNum <- function(formula, group, data, index = c("silhouette", "ch", "db"), max_clust = 10, base_val = FALSE) {

  # define global variables
  CluMP_ID = CluMP_X1 = CluMP_Y = ID = Visit = X1 = X1_ann = Y = Y.x = Y.y = Y_ci =
    abs_angle_radian = abs_change = abs_change_ann = angle_radian = best = bestVal =
    cluster = cos_denom = cos_nom  = cosinus = f_up = mean_Time = mean_Y =
    memb_CluMP = nVisit = number = obsah_trojuh = sd_Y = slope =
    slope_first_last = timepoint = value = . =NULL

  # Modelframe from formula
  mf <- stats:: model.frame(formula = formula, data = data, na.action = NULL)
  variables <- names(mf)
  names(mf) <- c("CluMP_Y", "CluMP_X1")
  mf <- cbind(mf, "CluMP_ID" = data[, group])
  mf_complete <- mf[stats:: complete.cases(mf),]
  mf_complete <- mf_complete %>% group_by(CluMP_ID) %>% mutate(count = n()) %>% filter(count > 2)
  mf_non_complete <- mf %>% filter(!(CluMP_ID %in% mf_complete$CluMP_ID))
  mf <- mf_complete

  if (!is.numeric(mf$CluMP_Y)){
    stop('Clustering variable should be numeric')
  }
  if (!is.numeric(mf$CluMP_X1)){
    stop('Time variable should be numeric')
  }
  if(nrow(mf_non_complete) > 1) {
    warning('Some of objects has less then 3 observations or contains NA values. Result data set is reduced.', call. = FALSE)
  }


  if (is.null(max_clust)){
    max_clust <- min(length(unique(mf$CluMP_ID)), 10)
  }else if(max_clust > length(unique(mf$CluMP_ID))){
    warning(paste0('Number of clusters cannot exceed the number of groups in data. The number of cluster  set as ', max(data$ID)), call. = FALSE)
    cl_numb <- length(unique(mf$CluMP_ID))
  }


  # CluMP prepare
  cluster_tmp = mf %>%
    group_by(CluMP_ID) %>%
    arrange(CluMP_X1) %>%
    mutate(Visit = row_number()) %>%
    filter(max(Visit) >= 2) %>%
    mutate(base_val = first(CluMP_Y),
           obsah_trojuh = c(NA, diff(CluMP_Y)/diff(CluMP_X1)/2),
           slope = (CluMP_Y - first(CluMP_Y))/(CluMP_X1 - first(CluMP_X1)),
           slope_first_last = (last(CluMP_Y) - first(CluMP_Y))/(last(CluMP_X1) - first(CluMP_X1)),
           cos_nom = (last(CluMP_X1) - first(CluMP_X1))*(CluMP_X1 - first(CluMP_X1)) + (last(CluMP_Y) - first(CluMP_Y))*(CluMP_Y - first(CluMP_Y)),
           cos_denom = sqrt((last(CluMP_X1) - first(CluMP_X1))^2 + (last(CluMP_Y) - first(CluMP_Y))^2)*sqrt((CluMP_X1 - first(CluMP_X1))^2 + (CluMP_Y - first(CluMP_Y))^2),
           cosinus = ifelse(cos_nom/cos_denom > 1, 1, cos_nom/cos_denom),
           abs_angle_radian = acos(cosinus),
           angle_radian = ifelse(slope_first_last < slope, 1, -1)*abs_angle_radian) %>%
    summarise(base_val = first(base_val),
              mean_triangle = mean(obsah_trojuh, na.rm = T),
              sd_triangle = stats:: sd(obsah_trojuh, na.rm = T),
              mean_abs_triangle = mean(abs(obsah_trojuh), na.rm = T),
              sd_abs_triangle = stats:: sd(abs(obsah_trojuh), na.rm = T),
              k = stats:: cov(CluMP_Y, CluMP_X1)/stats:: var(CluMP_X1),
              rel_pos = sum(diff(CluMP_Y) > 0, na.rm = T)/max(sum(diff(CluMP_Y) < 0, na.rm = T), 0.1),
              angle_radian = ifelse(test = all(is.na(angle_radian)),
                                    yes  = NA,
                                    no   = ifelse(test = max(angle_radian, na.rm = T) == max(abs(angle_radian), na.rm = T),
                                                  yes  = max(abs(angle_radian), na.rm = T),
                                                  no   = -max(abs(angle_radian), na.rm = T))))
  if (base_val == F) cluster_tmp <- select(cluster_tmp, -base_val)
  # Scaling values for clustering
  cluster_norm_tmp <- cbind.data.frame("CluMP_ID" = cluster_tmp$CluMP_ID,
                                       scale(cluster_tmp[, -1], center = TRUE, scale = TRUE))




  #Prepare indices for plotting
  plotData <- data.frame()
  for (ind in index) {

    tmpClust <- NbClust:: NbClust(data     = cluster_norm_tmp[, -1],
                        distance = "euclidean",
                        min.nc   = 2,
                        max.nc   = max_clust,
                        method   = "complete",
                        index    = ind)

    tmpClust <- data.frame(number = (2:max_clust),
                           value  = tmpClust$All.index,
                           best   = as.numeric(tmpClust$Best.nc[1]),
                           index  = ind)

    plotData <- rbind(plotData, tmpClust)

  }


  # jeden vystup - doporuceny pocet clust na zaklade indexu
  # druhy vystup - graf: plotIndex
  # treti vystup- tabulka s indexy: plotData
  bestClust <- plotData %>%
    group_by(index) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(best) %>%
    pull()

  #bestClust <- as.vector(bestClust$best)
  ux <- unique(bestClust)
  modus <- ux[tabulate(match(bestClust, ux)) == max(tabulate(match(bestClust, ux)))]
  indx <- paste(index, collapse = ', ')
  bests <- paste(modus, collapse = ' or ')

  # if there is no mode
  if(length(modus) == length(unique(bestClust)) & length(modus) > 1) {
    bN <- paste("Unable to find best number of clusters (mode does not exists)")
  } else {
    bN <- paste("The best number of clusters based on", indx, ifelse(length(index) > 1, "indices", "index"), "is", bests)
  }
  if(length(modus) == length(unique(bestClust)) & length(modus) > 1) {
    bestNumber <- NA
  } else {
    bestNumber <- modus
  }


  # Plot everything
  bestIndex <- plotData %>%
    filter(number == best) %>%
    rename(bestVal = value) %>%
    dplyr::select(bestVal, index)
  plotData <- left_join(plotData, bestIndex, by = "index")

  plotIndex <- ggplot(data = plotData, aes(x = number, y = value, group = index)) +
    geom_line() + geom_point() +
    geom_vline(aes(xintercept = best, group = index)) +
    geom_point(aes(x = best, y = bestVal),
               size = 2.5, color = "red") +
    labs(x = "number of clusters", y = "coefficient value",
         caption = bN) +
    facet_wrap( ~ index, ncol = max(1,length(index)%/%2), scales="free_y")


  # return(list(plot = plotIndex, indexTable = as.data.frame(plotData), bestNumber = bN))
  # return(list(plot = plotIndex, indexTable = as.data.frame(plotData), bestNumber = as.list(bestNumber)))
  return(list("plot" = plotIndex, "indexTable" = as.data.frame(plotData), "bestNumber" = bestNumber,
              "formula" = formula, "group" = group, "data" = data, "CLUMPdata" = cluster_norm_tmp, "variables" = variables))

}
