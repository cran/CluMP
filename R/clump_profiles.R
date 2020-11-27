
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
#' @import data.table
#' @rawNamespace import(rlang, except = ':=')
#' @rawNamespace import(dplyr, except = c(last, first, between))
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
#' CluMPoutput <- CluMP(formula = Y ~ Time, group = "ID", data = data, cl_numb = 3)
#' CluMP_profiles(CluMPoutput, cat_vars = "Gender")
#'
CluMP_profiles <- function(CluMPoutput, cont_vars = NULL, cat_vars = NULL, show_NA = FALSE) {

  # define global variables
  CluMP_ID = CluMP_X1 = CluMP_Y = ID = Visit = X1 = X1_ann = Y = Y.x = Y.y = Y_ci =
    abs_angle_radian = abs_change = abs_change_ann = angle_radian = best = bestVal =
    cluster = cos_denom = cos_nom  = cosinus = f_up = mean_Time = mean_Y =
    memb_CluMP = nVisit = number = obsah_trojuh = sd_Y = slope =
    slope_first_last = timepoint = value = . = .. = ..colour =
    ..cols = ..cont_vars = ..group = ..scale_cols = Time = NULL


  nb <- setDT(CluMPoutput$CluMP)[, n := .N, by = memb_CluMP]
  cols <- c("memb_CluMP", "n")
  nb <- setDT(nb)[, ..cols]
  nb <- nb[!duplicated(nb),]

  data <- CluMPoutput$data
  colnames(data)[which(colnames(data) %in% as.character(CluMPoutput$formula[2]))] <- "Y"
  colnames(data)[which(colnames(data) %in% as.character(CluMPoutput$formula[3]))] <- "X1"
  colnames(data)[which(colnames(data) %in% CluMPoutput$group)] <- "ID"

  if (!show_NA) {
    data <- setDT(data)[!is.na(memb_CluMP), ]
  }

  #Y - baseline, change
  #baseline
  Bl_data <- setDT(data)[,Visit := 1:.N, by = ID]
  Bl_data <- setDT(Bl_data)[Visit == 1,,]
  Bl_data$memb_CluMP = as.factor(as.character(Bl_data$memb_CluMP))
  Bl_data <- setDT(Bl_data)[, ':='(
    n = n_distinct(ID),
    meanY = mean(Y, na.rm = TRUE),
    sdY = stats::sd(Y, na.rm = TRUE),
    medianY = median(Y, na.rm = TRUE),
    Q25Y = stats::quantile(Y,0.25, na.rm = TRUE),
    Q75Y = stats::quantile(Y,0.75, na.rm = TRUE),
    minY = min(Y, na.rm = TRUE),
    maxY = max(Y, na.rm = TRUE)
  ), by = memb_CluMP]
  cols <- c("memb_CluMP", "n", "meanY", "sdY",
            "medianY","Q25Y", "Q75Y", "minY", "maxY")
  Bl_data <- setDT(Bl_data)[, ..cols]
  Bl_data <- Bl_data[!duplicated(Bl_data),]

  #change
  change_data <- setDT(data)[,Visit := 1:.N, by = ID]
  change_data <- setDT(change_data)[, .SD[c(1,.N)], by = ID]
  change_data <- setDT(change_data)[,.(ID, Visit, X1, Y, memb_CluMP),]
  change_data <- setDT(change_data)[,Visit := 1:.N, by = ID]

  visit_f <- setDT(change_data)[Visit == 1,,]
  visit_f$Visit <- NULL
  visit_f$X1 <- NULL
  visit_f$memb_CluMP <- NULL

  visit_l <- setDT(change_data)[Visit == 2,,]
  visit_l$Visit <- NULL
  change_data <- merge(visit_f, visit_l, by = "ID")

  change_data <- setDT(change_data)[, ':='(
    abs_change = Y.y - Y.x,
    abs_change_ann = (Y.y - Y.x)/X1),]
  change_data <- setDT(change_data)[, ':='(
    n = n_distinct(ID),
    mean_abs_change = mean(abs_change, na.rm = TRUE),
    sd_abs_change = stats:: sd(abs_change, na.rm = TRUE),
    q25_abs_change = stats:: quantile(abs_change, 0.25, na.rm = TRUE),
    q50_abs_change = stats:: quantile(abs_change, 0.5, na.rm = TRUE),
    q75_abs_change = stats:: quantile(abs_change, 0.75, na.rm = TRUE),
    mean_abs_change_ann = mean(abs_change_ann, na.rm = TRUE),
    sd_abs_change_ann = stats:: sd(abs_change_ann, na.rm = TRUE),
    q25_abs_change_ann = stats:: quantile(abs_change_ann, 0.25, na.rm = TRUE),
    q50_abs_change_ann = stats:: quantile(abs_change_ann, 0.5, na.rm = TRUE),
    q75_abs_change_ann = stats:: quantile(abs_change_ann, 0.75, na.rm = TRUE)
    ), by = memb_CluMP]
  cols <- c("memb_CluMP", "n", "mean_abs_change", "sd_abs_change",
            "q25_abs_change", "q50_abs_change", "q75_abs_change",
            "mean_abs_change_ann", "sd_abs_change_ann", "q25_abs_change_ann",
            "q50_abs_change_ann", "q75_abs_change_ann")
  change_data <- setDT(change_data)[, ..cols]
  change_data <- change_data[!duplicated(change_data),]

  Y <- list("baseline" = as.data.frame(Bl_data), "change" = as.data.frame(change_data))

  #Time
  f_up_data <- setDT(data)[, ':='(nVisit = .N,
                                  f_up = data.table::last(X1),
                                  memb_CluMP = data.table::first(memb_CluMP))
                           , by = ID][order(X1)]
  f_up_data <- setDT(data)[, ':='(n = n_distinct(ID),
                                  q25_visit = stats:: quantile(nVisit, 0.25, na.rm = TRUE),
                                  q50_visit = stats:: quantile(nVisit, 0.5, na.rm = TRUE),
                                  q75_visit = stats:: quantile(nVisit, 0.75, na.rm = TRUE),
                                  mean_fup = mean(f_up, na.rm = TRUE),
                                  sd_fup = stats:: sd(f_up, na.rm = TRUE),
                                  q25_fup = stats:: quantile(f_up, 0.25, na.rm = TRUE),
                                  q50_fup = stats:: quantile(f_up, 0.5, na.rm = TRUE),
                                  q75_fup = stats:: quantile(f_up, 0.75, na.rm = TRUE))
                           , by = memb_CluMP]
  cols <- c("memb_CluMP", "n", "q25_visit", "q50_visit", "q75_visit",
            "mean_fup", "sd_fup", "q25_fup", "q50_fup", "q75_fup")
  f_up_data <- setDT(f_up_data)[, ..cols]
  f_up_data <- f_up_data[!duplicated(f_up_data),]
  time <- as.data.frame(f_up_data)

  if(!is.null(cont_vars)) {
    #cont vars - bl, change
    cont_bl <- setDT(data)[, Visit := 1:.N, by = ID][order(X1)]
    cont_bl <- setDT(cont_bl)[Visit == 1,,]
    cols <- c(cont_vars, "memb_CluMP")
    cont_bl <- setDT(cont_bl)[, ..cols]
    cont_bl <- cont_bl[!duplicated(cont_bl),]

    cont_bl <- cont_bl[, c(lapply(.SD, mean),
                           lapply(.SD, sd)), by = memb_CluMP]
    names(cont_bl)[-1] <- c(paste0(cont_vars, "_mean"), paste0(cont_vars, "_sd"))

    cont_bl <- merge(nb, cont_bl, by = "memb_CluMP")

    cols <- c("X1", "memb_CluMP", cont_vars)
    visit_f <- setDT(data)[, .SD[1], by = ID][order(X1)][, ..cols]
    visit_l <- setDT(data)[, .SD[.N], by = ID][order(X1)][, ..cols]
    cont_change <- cbind(visit_f[, "memb_CluMP",], visit_l[, !"memb_CluMP",] - visit_f[, !"memb_CluMP",])

    cont_change <- cbind(cont_change[, "memb_CluMP"], cont_change[, ..cont_vars]/cont_change[, "X1"])
    cont_change <- cont_change[, c(lapply(.SD, mean),
                                   lapply(.SD, sd)), by = memb_CluMP]
    names(cont_change)[-1] <- c(paste0(cont_vars, "_ann_mean"), paste0(cont_vars, "_ann_sd"))

    cont_change <- merge(nb, cont_change, by = "memb_CluMP")

    cont <- list(baseline = as.data.frame(cont_bl), change = as.data.frame(cont_change))
  } else if (is.null(cont_vars)) {
    cont <- NULL
  }

  if(!is.null(cat_vars)) {
    # #baseline
    data_bl <- setDT(data)[, .SD[1], by = ID][order(X1)]
    cat_bl <- tableone::CreateCatTable(data = data_bl, strata = "memb_CluMP", vars = cat_vars)
    # #end
    data_end <- setDT(data)[, .SD[.N], by = ID][order(X1)]
    cat_end <- tableone::CreateCatTable(data = data_end, strata = "memb_CluMP", vars = cat_vars)
    categ <- list(baseline = as.data.frame(print(cat_bl, exact = "ascites", quote = F)),
                  end = as.data.frame(print(cat_end, exact = "ascites", quote = F)))
  } else if (is.null(cat_vars)) {
    categ <- NULL
  }

  #CluMP_params
  param_table <- CluMPoutput$CluMP[,-1]
  param_names <- names(param_table)[-length(names(param_table))]
  param_table <- param_table[, c(lapply(.SD, mean),
                                 lapply(.SD, sd)), by = memb_CluMP]
  names(param_table)[-1] <- c(paste0(param_names, "_mean"), paste0(param_names, "_sd"))

  param_table <- merge(nb, param_table, by = "memb_CluMP")
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
