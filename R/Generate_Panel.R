#' Generate an artificial Micro-Panel (longitudinal) Data
#'
#' This function creates artificial linear or non-linear micro-panel (longitudinal) data coming from
#' generating process with a certain function (linear, quadratic, cubic, exponencial) set of parameters
#' (fixed and random (intercept, slope) effects of time).
#' @param n An integer specifying the number of individuals (trajectories) being observed.
#' @param Param Object of \code{\link[base]{data.frame}} containing regression parameters for each cluster.
#' The dimensions are the various number of generating clusters and the fixed number of parameters.
#' The second dimension (the fixed number of parameters) is given by the type of regression model specified by the argument "RegModel".
#' For more information about the parameters, see documentation of: \code{\link{ParamLinear}} for linear model, \code{\link{ParamQuadrat}} for quadratic,
#' \code{\link{ParamCubic}} for cubic model and \code{\link{ParamExpon}} for exponential model.
#' @param NbVisit A positive integer numeric input defining expected number of visits. Option is \emph{Fixed} or \emph{Random}.
#' Number of visits given by the argument \code{VisitFreq}. If \code{VisitFreq} is \emph{Fixed}, the \code{NbVisits} defines exact number of visits for all individuals.
#' If \code{VisitFreq} is \emph{Random} then each individual has different number of visits.
#' The number of visits is then generated from the poisson distribution with the mean (lambda) equal to \code{NbVisits}.
#' @param VisitFreq String that defines the frequency of visits for each individual. Option is \emph{Random} or \emph{Fixed}. If set to \emph{Fixed} or not defined, each individual has the same number of visits given by \code{NbVisits}.
#' If set as \emph{Random} the number of visits is generated from poisson distribution for each individual with the mean equal to the argument \code{NbVisits}.
#' For example if this parameter is set as 5 then the random integer from interval of -5 to 5 is drawned and added to the time variable.
#' Make sure that \code{TimeVar} must be lower then the number of days in parameter \code{units}.
#' @param TimeVar A positive integer representing daily, time variability of the occurrence of repeated measurement (timepoint) from the regular,
#' fixed occurrence (visit) given by the argument units.
#' For example, if this argument is set to 5 then the random integer from interval of -5 to 5 is drawn and added to the time variable.
#' TimeVar must be lower than the regular frequency of repeat measurement given by the argument units.
#' @param RegModel String specifying the mathematical function for generating trajectory for each of n individuals. Options are \emph{linear},
#' \emph{quadratic}, \emph{cubic} or \emph{exponential}.
#' If set to \emph{linear} or not defined, then each trajectory has a linear trend.
#' If set to \emph{quadratic}, then each trajectory has a quadratic development in time.
#' If set to \emph{cubic} then each trajectory has cubic development.
#' If set to \emph{exponential}, then each trajectory has exponential development.
#' @param ClusterProb Numeric scalar (for 2 clusters) or a vector of numbers (for >2 clusters) defining the probability of each cluster.
#' If not defined, then each cluster has the same occurrence probability.
#' @param Rho A numeric scalar specifying autocorrelation parameter with the values from range 0 to 1.
#' If set as 0 or not define then there is no autocorrelation between the within-individual repeated observations.
#' @param units String defining the units of time series. Options are \emph{day}, \emph{week}, \emph{month} or \emph{year}.
#' @keywords CluMP
#' @return Generates artificial panel data.
#' @importFrom MASS mvrnorm
#' @importFrom stats rpois
#' @export
#' @examples
#' set.seed(123)
#' #Simple Linear model where each individual has 10 observations.
#' data <- GeneratePanel(n = 100, Param = ParamLinear, NbVisit = 10)
#'
#' #Exponential model where each individual has 10 observations.
#' data <- GeneratePanel(100, ParamExpon, NbVisit = 10, VisitFreq = "Fixed", RegModel = "exponential")
#' PanelPlot(data)
#'
#' #Cubic model where each individual has random number of observations on daily basis.
#' #Average number of observation is given by parameter NbVisit.
#' data <- GeneratePanel(n = 100, Param = ParamCubic, NbVisit = 100, RegModel = "cubic", units = "day")
#' PanelPlot(data)
#'
#' #Quadratic model where each individual has random number of observations.
#' #Each object is observede weekly with variability 2 days.
#' data <- GeneratePanel(5,ParamQuadrat,NbVisit=50,RegModel="quadratic",units="week",TimeVar=2)
#' PanelPlot(data)
#'
#' #Generate panel data with linear trend with 75% objects in first cluster and 25% in the second.
#' data <- GeneratePanel(n = 100, Param = ParamLinear, NbVisit = 10, ClusterProb = c(0.75, 0.25))
#' PanelPlot(data, colour = "Cluster")

GeneratePanel <- function(n, Param, NbVisit, VisitFreq = NULL, TimeVar = NULL, RegModel = NULL, ClusterProb = NULL, Rho = NULL, units = NULL) {

  # define global variables
  CluMP_ID = CluMP_X1 = CluMP_Y = ID = Visit = X1 = X1_ann = Y = Y.x = Y.y = Y_ci =
    abs_angle_radian = abs_change = abs_change_ann = angle_radian = best = bestVal =
    cluster = cos_denom = cos_nom  = cosinus = f_up = mean_Time = mean_Y =
    memb_CluMP = nVisit = number = obsah_trojuh = sd_Y = slope =
    slope_first_last = timepoint = value = . = .. = ..colour =
    ..cols = ..cont_vars = ..group = ..scale_cols = Time = NULL

  k = nrow(Param) #defines the number of clusters

  #Setting not defined parameters
  if (is.null(ClusterProb))
    ClusterProb <- rep(1/k, k)

  if (is.null(ClusterProb) & sum(ClusterProb) != 1)
    stop("The sum of occurance probabilities in parameter lusterProb has to be equal to 1.")

  if (is.null(Rho))
    Rho <- 0

  if (!is.numeric(NbVisit))
    stop("Parameter NbVisit is not numeric")

  if (is.null(VisitFreq))
    VisitFreq <- "Fixed"

  if (is.null(TimeVar) || !is.numeric(TimeVar))
    TimeVar <- 0

  if (is.null(RegModel)){
    RegModel <- "linear"
  } else if (!RegModel %in% c("linear", "quadratic", "cubic", "exponential")) {
    stop("Incorrect RegModel sellection. Choose between 'linear', 'quadratic' or 'cubic' only.")
  }

  if (is.null(units))
    units <- "year"

  if (units == "day") {
    time.unit <- 1
  } else if (units == "week") {
    time.unit <- 7
  } else if (units == "month") {
    time.unit <- 30
  } else if (units == "year") {
    time.unit <- 365
  }

  if (TimeVar >= time.unit) {
    stop(paste0("Parameter TimeVar must be lower then ", time.unit, ", because the parameter units is '",units,"' (", time.unit," days)"))
  }

  if (VisitFreq == "Fixed") {
    t <- rep(NbVisit, n)
  } else if (VisitFreq == "Random") {
    t <- rpois(n, lambda = NbVisit)
  }

  Visit <- sequence(t)

  Time <- ifelse(Visit == 1, 0, (Visit-1)*time.unit + sample(x = -TimeVar:TimeVar, size = length(Visit), replace = T))/time.unit

  Cluster <- rep(sort(sample(x = 1:k, size = n, replace = T, prob = ClusterProb)), t)
  ID <- rep(1:n, t)

  data <- cbind.data.frame(ID, Cluster, Visit, Time)

  #Fixed and random effects
  for(i in 1:n){
    data.i <- subset(data, ID == i)

    cluster <- unique(data.i$Cluster)

    #random effect
    cov <- Param[cluster, "corr"]*sqrt(Param[cluster, "varU0"])*sqrt(Param[cluster, "varU1"])
    G   <- matrix(c(Param[cluster, "varU0"], cov, cov, Param[cluster, "varU1"]), 2, 2) #new
    U.i <- t(as.matrix(mvrnorm(n = 1, mu = rep(0, nrow(G)), Sigma = G)))

    #residuals
    R <- diag(nrow(data.i))
    R <- Rho^abs(row(R)-col(R))
    R <- R*Param[cluster, "varE"]
    E.i <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(R)), Sigma = R)

    #Regressor X
    Time = data.i[, "Time"]

    #Regression function
    if (RegModel == "linear"){
      data[data$ID == i, "Y"] <- (Param[cluster, "b0"] + U.i[, 1]) + (Param[cluster, "b1"] + U.i[, 2])*Time + E.i
    }else if (RegModel == "quadratic"){
      data[data$ID == i, "Y"] <- (Param[cluster, "b0"] + U.i[, 1]) + (Param[cluster, "b1"] + U.i[, 2])*Time + Param[cluster, "b2"]*Time^2 + E.i
    }else if (RegModel == "cubic"){
      data[data$ID == i, "Y"] <- (Param[cluster, "b0"] + U.i[, 1]) + (Param[cluster, "b1"] + U.i[, 2])*Time + Param[cluster, "b2"]*Time^2 + Param[cluster, "b3"]*Time^3 + E.i
    } else if (RegModel == "exponential") {
      data[data$ID == i, "Y"] <- (Param[cluster, "b0"] + U.i[, 1]) + (Param[cluster, "b1"] + U.i[, 2])*exp(-Time*Param[cluster, "b2"]) + E.i
    }

  }

  data <- as.data.frame(data)
  data$Y <- as.numeric(data$Y)
  return(data)
}
