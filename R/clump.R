#' Cluster Micro-Panel (longitudinal) Data employing the CluMP algorithm
#'
#' This function clusters Micro-Panel (longitudinal) Data (or trajectories) to a pre-defined number of clusters by employing Feature-Based Clustering of Micro-Panel (longitudinal) Data algorithm called CluMP (see Reference). Currently, only univariate clustering analysis is available.
#' @param formula A two-sided \code{\link[stats]{formula}} object with a numeric clustering variable (Y) on the left of a ~ separator and the time (numeric) variable on the right. Time is measured from the start of the follow-up period (baseline). Any time units are possible.
#' @param group A grouping factor variable (vector), i.e. single identifier for each individual (trajectory).
#' @param data A data frame containing the variables named in the \code{formula} and \code{group} arguments.
#' @param cl_numb An integer, positive number (scalar) specifying the number of clusters. The \code{\link{OptiNum}} function can be used to determine the optimal number of clusters according to common evaluation criteria (indices).
#' @param base_val Indicates whether include a value at zero time point as an additional clustering variable. Default is \emph{FALSE} and the standard number (7) of clustering parameters is used.
#' @keywords CluMP
#' @return Cluster Micro-Panel data. The output is the \code{\link[base]{list}} of 5 components which contain results from clustering.
#' @source Sobisek, L., Stachova, M., Fojtik, J. (2018) Novel Feature-Based Clustering of Micro-Panel Data (CluMP). Working paper version online: www.arxiv.org
#' @export
#' @import dplyr rlang
#' @examples
#' data <- GeneratePanel(n = 100, Param = ParamLinear, NbVisit = 10)
#' CluMP(formula = Y ~ Time, group = "ID", data = data, cl_numb = 3, base_val = FALSE)
#'
#' CluMP(formula = Y ~ Time, group = "ID", data = data, cl_numb = 3, base_val = TRUE)
#'
#'
CluMP <- function(formula, group, data, cl_numb = NA, base_val = FALSE){


  #library(dplyr)
  #library(amap)
  # define global variables
  CluMP_ID = CluMP_X1 = CluMP_Y = ID = Visit = X1 = X1_ann = Y = Y.x = Y.y = Y_ci =
  abs_angle_radian = abs_change = abs_change_ann = angle_radian = best = bestVal =
  cluster = cos_denom = cos_nom  = cosinus = f_up = mean_Time = mean_Y =
  memb_CluMP = nVisit = number = obsah_trojuh = sd_Y = slope =
  slope_first_last = timepoint = value = . = NULL

  if(!any(all.vars(formula) %in% names(data))) {
    varErr <- all.vars(formula)[!(all.vars(formula) %in% names(data))]
    if(length(varErr) == 1) {
      stop(paste0('Variable "', varErr, '" does not present in the data.'))
    }else {
      stop(paste0('Variables "', paste(varErr, collapse = ', '), '" do not present in the data.'))
    }
  }

  # Modelframe from formula
  mf <- stats:: model.frame(formula = formula, data = data, na.action = NULL)
  variables <- names(mf)
  if(length(variables) > 2) {
    stop('It is not possible to cluster multiple variables. Consider using one of them.')
  }
  if(!any(variables %in% names(data))) {
    varErr <- variables[!(variables %in% names(data))]
    stop(paste0('Variable "', varErr, '" does not present in the data.'))
  }
  names(mf) <- c("CluMP_Y", "CluMP_X1")
  mf <- cbind(mf, "CluMP_ID" = data[, group])
  mf_complete <- mf[stats:: complete.cases(mf),]
  mf_complete <- mf_complete %>% group_by(CluMP_ID) %>% mutate(count = n()) %>% filter(count > 2)
  mf_non_complete <- mf %>% filter(!(CluMP_ID %in% mf_complete$CluMP_ID))
  mf <- mf_complete



  if(nrow(mf_non_complete) > 1) {
    warning('Some of objects has less then 3 observations or contains NA values. Result data set is reduced.', call. = FALSE)
  }

  #Not defined entry
  if (!is.numeric(mf$CluMP_Y)){
    stop('Clustering variable should be numeric')
  }
  if (!is.numeric(mf$CluMP_X1)){
    stop('Time variable should be numeric')
  }
  if (is.na(cl_numb)) {
    warning('Number of clusters defined as 3. If you want different number of clusters define parameter ClNumb', call. = FALSE)
    cl_numb <- 3
  }
  if (!is.numeric(cl_numb) || cl_numb <= 1 || (cl_numb %% 1) != 0){
    stop('Number of clusters argument cl_numb should be an integer greater than 1')
  }
  if (cl_numb > length(unique(mf$CluMP_ID))){
    warning(paste0('Number of clusters cannot exceed the number of groups in data. The number of cluster  set as ', max(data$ID)), call. = FALSE)
    cl_numb <- length(unique(mf$CluMP_ID))
  }


  # CluMP prepare
  cluster_tmp = mf %>%
    dplyr::group_by(CluMP_ID) %>%
    dplyr::arrange(CluMP_X1) %>%
    dplyr::mutate(Visit = row_number()) %>%
    dplyr::filter(max(Visit) >= 2) %>%
    dplyr::mutate(base_val = first(CluMP_Y),
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
              sd_triangle = stats::sd(obsah_trojuh, na.rm = T),
              mean_abs_triangle = mean(abs(obsah_trojuh), na.rm = T),
              sd_abs_triangle = stats::sd(abs(obsah_trojuh), na.rm = T),
              k = stats::cov(CluMP_Y, CluMP_X1)/stats::var(CluMP_X1),
              rel_pos = sum(diff(CluMP_Y) > 0, na.rm = T)/max(sum(diff(CluMP_Y) < 0, na.rm = T), 0.1),
              angle_radian = ifelse(test = all(is.na(angle_radian)),
                                    yes  = NA,
                                    no   = ifelse(test = max(angle_radian, na.rm = T) == max(abs(angle_radian), na.rm = T),
                                                  yes  = max(abs(angle_radian), na.rm = T),
                                                  no   = -max(abs(angle_radian), na.rm = T))))
  if (base_val == F) cluster_tmp <- dplyr::select(cluster_tmp, -base_val)
  # Scaling values for clustering
  cluster_norm_tmp <- cbind.data.frame("CluMP_ID" = cluster_tmp$CluMP_ID,
                                       scale(cluster_tmp[, -1], center = TRUE, scale = TRUE))

  # Clustering
  tmp_comb <- amap:: Dist(cluster_norm_tmp[, -1], method = "euclid")
  tmp_comb <- stats:: hclust(tmp_comb, method = "ward.D")
  cluster_norm_tmp$memb_CluMP <- stats:: cutree(tmp_comb, k = cl_numb)

  cluster_tmp <- dplyr::left_join(cluster_tmp, cluster_norm_tmp[, c("CluMP_ID", "memb_CluMP")], by = "CluMP_ID")
  cluster_tmp <- cluster_tmp %>%
    dplyr::select(CluMP_ID, everything())

  names(cluster_tmp)[1] <- group # put back old names

  data <- left_join(data, cluster_tmp[, c(group, "memb_CluMP")], by = group)

  return(list("CluMP" = cluster_tmp, "data" = data, "formula" = formula, "variables" = variables, "group" = group))
}
