
#' Summary characteristics of identified clusters via CluMP
#'
#' The function CluMP_profiles provides a description (profile) for each cluster. The description is in the form of a summary list containing descriptive statistics of a cluster variable, time variable, cluster parameters and other variables (covariates), both continuous and categorical.
#' @param CluMPoutput An object (output) from the \code{\link{CluMP}} function.
#' @param cont_vars An optional single character or a character vector of continuous variables' names (from the original dataset).
#' @param cat_vars An optional single character or a character vector of categorical variables' names (from the original dataset).
#' @param show_NA Logical scalar. Should be calculated and shown descriptive statistics for \emph{NA} cluster if exists? Default is \emph{FALSE}. \emph{NA} cluster gathers improper individuals (trajectories with < 3 not missing observations) for longitudinal clustering.
#' @keywords CluMP
#' @return Returns a \code{\link[base]{list}} with cluster variable (Y) summary, both baseline and changes; time and a summary of the number of observations (visits); clustering parameters summary and optional continuous variables summary (baseline and changes) and categorical variables summary (baseline and end).
#' @export
#' @import dplyr
#' @importFrom stats median sd
#' @examples
#' set.seed(123)
#' dataMale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
#' dataMale$Gender <- "M"
#' dataFemale <- GeneratePanel(n = 50, Param = ParamLinear, NbVisit = 10)
#' dataFemale$ID <- dataFemale$ID + 50
#' dataFemale$Gender <- "F"
#' data <- rbind(dataMale, dataFemale)
#'
#' CLUMP3 <- CluMP(formula = Y ~ Time, group = "ID", data = data, cl_numb = 3)
#' CluMP_profiles(CLUMP3, cat_vars = "Gender")
#'
CluMP_profiles <- function(CluMPoutput, cont_vars = NULL, cat_vars = NULL, show_NA = FALSE) {

  # define global variables
  CluMP_ID = CluMP_X1 = CluMP_Y = ID = Visit = X1 = X1_ann = Y = Y.x = Y.y = Y_ci =
    abs_angle_radian = abs_change = abs_change_ann = angle_radian = best = bestVal =
    cluster = cos_denom = cos_nom  = cosinus = f_up = mean_Time = mean_Y =
    memb_CluMP = nVisit = number = obsah_trojuh = sd_Y = slope =
    slope_first_last = timepoint = value = . = NULL
  
  nb <- CluMPoutput$CluMP %>%
    group_by(memb_CluMP) %>%
    summarise(n = n())

  data <- CluMPoutput$data
  colnames(data)[which(colnames(data) %in% CluMPoutput$variables)] <- c("Y", "X1")
  colnames(data)[which(colnames(data) %in% CluMPoutput$group)] <- "ID"

  if (!show_NA) {
    data <- data %>%
      filter(!is.na(memb_CluMP))
  }

  #Y - baseline, change
  #baseline
  Bl_data <- data %>%
    group_by(ID) %>%
    arrange(X1) %>%
    mutate(Visit = row_number()) %>%
    filter(Visit == 1) %>%
    ungroup() %>%
    mutate(memb_CluMP = as.factor(as.character(memb_CluMP))) %>%
    group_by(memb_CluMP) %>%
    summarise(n = n_distinct(ID), meanY = mean(Y, na.rm = TRUE), sdY = stats:: sd(Y, na.rm = TRUE),
              medianY = median(Y, na.rm = TRUE), Q25Y = stats:: quantile(Y,0.25, na.rm = TRUE), Q75Y = stats:: quantile(Y,0.75, na.rm = TRUE),
              minY = min(Y, na.rm = TRUE), maxY = max(Y, na.rm = TRUE))
  #change
  change_data <- data %>%
    group_by(ID) %>%
    arrange(X1) %>%
    mutate(Visit = row_number()) %>%
    slice(c(1,n())) %>%
    dplyr::select(ID, Visit, X1, Y, memb_CluMP) %>%
    mutate(Visit = row_number()) %>%
    ungroup()
  visit_f <- change_data %>% filter(Visit == 1) %>% dplyr::select(-c(Visit, X1, memb_CluMP))
  visit_l <- change_data %>% filter(Visit == 2) %>% dplyr::select(-Visit)
  change_data <- left_join(visit_f, visit_l, by = "ID")
  change_data <- change_data %>%
    mutate(abs_change = Y.y - Y.x, abs_change_ann = (Y.y - Y.x)/X1) %>%
    group_by(memb_CluMP) %>%
    summarise(n = n_distinct(ID),
              mean_abs_change = mean(abs_change, na.rm = TRUE),
              sd_abs_change = stats:: sd(abs_change, na.rm = TRUE),
              q25_abs_change = stats:: quantile(abs_change, 0.25, na.rm = TRUE),
              q50_abs_change = stats:: quantile(abs_change, 0.5, na.rm = TRUE),
              q75_abs_change = stats:: quantile(abs_change, 0.75, na.rm = TRUE),
              mean_abs_change_ann = mean(abs_change_ann, na.rm = TRUE),
              sd_abs_change_ann = stats:: sd(abs_change_ann, na.rm = TRUE),
              q25_abs_change_ann = stats:: quantile(abs_change_ann, 0.25, na.rm = TRUE),
              q50_abs_change_ann = stats:: quantile(abs_change_ann, 0.5, na.rm = TRUE),
              q75_abs_change_ann = stats:: quantile(abs_change_ann, 0.75, na.rm = TRUE))
  Y <- list("baseline" = as.data.frame(Bl_data), "change" = as.data.frame(change_data))

  #Time
  f_up_data <- data %>%
    group_by(ID) %>%
    arrange(X1) %>%
    summarise(nVisit = n(), f_up = last(X1), memb_CluMP = first(memb_CluMP)) %>%
    ungroup() %>%
    group_by(memb_CluMP) %>%
    summarise(n = n_distinct(ID),
              q25_visit = stats:: quantile(nVisit, 0.25, na.rm = TRUE),
              q50_visit = stats:: quantile(nVisit, 0.5, na.rm = TRUE),
              q75_visit = stats:: quantile(nVisit, 0.75, na.rm = TRUE),
              mean_fup = mean(f_up, na.rm = TRUE),
              sd_fup = stats:: sd(f_up, na.rm = TRUE),
              q25_fup = stats:: quantile(f_up, 0.25, na.rm = TRUE),
              q50_fup = stats:: quantile(f_up, 0.5, na.rm = TRUE),
              q75_fup = stats:: quantile(f_up, 0.75, na.rm = TRUE)) # min max
  time <- as.data.frame(f_up_data)

  if(!is.null(cont_vars)) {
    #cont vars - bl, change
    cont_bl <- data %>%
      group_by(ID) %>%
      arrange(X1) %>%
      mutate(Visit = row_number()) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(cont_vars, memb_CluMP) %>%
      group_by(memb_CluMP) %>%
      summarise_all(list(~mean(., na.rm = TRUE),~stats:: sd(., na.rm = TRUE)))
    cont_bl <- left_join(nb, cont_bl, by = "memb_CluMP")

    cont_change <- data %>%
      group_by(ID) %>%
      arrange(X1) %>%
      mutate(Visit = row_number()) %>%
      slice(c(1,n())) %>%
      mutate(Visit = row_number()) %>%
      ungroup() %>%
      dplyr::select(cont_vars, Visit, ID, X1, memb_CluMP)
    visit_f <- cont_change %>% filter(Visit == 1) %>% dplyr::select(memb_CluMP, cont_vars)
    visit_l <- cont_change %>% filter(Visit == 2) %>% dplyr::select(X1, cont_vars)
    cont_change <- cbind(visit_l[,-1] - visit_f[,-1], visit_f[,1], visit_l[,1])
    cont_change <- cont_change %>%
      group_by(memb_CluMP) %>%
      mutate_all((list(ann = ./X1))) %>%
      dplyr::select(-c(X1, X1_ann)) %>%
      summarise_all(list(~mean(., na.rm = TRUE),~stats:: sd(., na.rm = TRUE)))
    cont_change <- left_join(nb, cont_change, by = "memb_CluMP")
    cont <- list(baseline = as.data.frame(cont_bl), change = as.data.frame(cont_change))
  } else if (is.null(cont_vars)) {
    cont <- NULL
  }

  if(!is.null(cat_vars)) {
    # #baseline
    data_bl <- data %>% group_by(ID) %>% arrange(X1) %>% slice(1)
    cat_bl <- tableone:: CreateCatTable(data = data_bl, strata = "memb_CluMP", vars = cat_vars)
    # #end
    data_end <- data %>% group_by(ID) %>% arrange(X1) %>% slice(n())
    cat_end <- tableone:: CreateCatTable(data = data_end, strata = "memb_CluMP", vars = cat_vars)
    categ <- list(baseline = as.data.frame(print(cat_bl, exact = "ascites", quote = F)),
                  end = as.data.frame(print(cat_end, exact = "ascites", quote = F)))
  } else if (is.null(cat_vars)) {
    categ <- NULL
  }

  #CluMP_params
  param_table <- CluMPoutput$CluMP[,-1] %>%
    group_by(memb_CluMP) %>%
    summarise_all(list(~mean, ~sd))
  param_table <- left_join(nb, param_table, by = "memb_CluMP")
  param_table <- as.data.frame(param_table)

  if(is.null(cont) & is.null(categ)) {
  return(list(Y = Y, time = time, CluMP_params = param_table))
  } else if (is.null(cont) & !is.null(categ)) {
    return(list(Y = Y, time = time, CluMP_params = param_table, cat = categ))
  } else if (!is.null(cont) & is.null(categ)) {
    return(list(Y = Y, time = time, CluMP_params = param_table, cont = cont))
  } else if (!is.null(cont) & !is.null(categ)) {
    return(list(Y = Y, time = time, CluMP_params = param_table, cont = cont, cat = categ))
  }
}
