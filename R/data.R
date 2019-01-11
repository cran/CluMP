#'Parameters of linear model
#'
#' Default parameters to generate micro-panel (longitudinal) data with linear trend. The parameters may differ per each cluster.
#' The parameters of each cluster are in rows. Number of rows denotes the number of clusters.
#' Fixed and random effects are taken from Uher et al. (2017).
#' @format It is adviced to keep parameters in \code{\link[base]{data.frame}}. The Parameters structure is as follows:
#' \describe{
#'   \item{b0}{fixed parameter of intercept}
#'   \item{b1}{fixed parameter of slope}
#'   \item{varU0}{variance of random factor U0 given to fixed parameter b0}
#'   \item{varU1}{variance of random factor U1 given to fixed parameter b1}
#'   \item{corr}{correlation between random factors U0 and U1}
#'   \item{varE}{the variability of the residuals}
#' }
#' @source Uher T, Vaneckova M, Krasensky J, Sobisek L, Tyblova M, Volna J, Seidl Z, Bergsland N, Dwyer MG, Zivadinov R, De Stefano N, Sormani MP, Havrdova EK, Horakova D. Pathological cut-offs of global and regional brain volume loss in multiple sclerosis. Mult Scler. 2017 Nov 1:1352458517742739. doi: 10.1177/1352458517742739.
#'
"ParamLinear"


#'Parameters of quadratic model
#'
#' Parameters to generate panel data with quadratic trend. The parameters may differ per each cluster.
#' The parameters of each cluster are in rows. Number of rows denotes the number of clusters.
#' Fixed effects are taken from Allen et al. (2005), and the source for random effects is Uher et al. (2017).
#' @format It is adviced to keep parameters in \code{\link[base]{data.frame}}. The Parameters structure is as follows:
#' \describe{
#'   \item{b0}{fixed parameter of intercept}
#'   \item{b1}{fixed parameter of slope}
#'   \item{b2}{fixed parameter of defining the quadraticity}
#'   \item{varU0}{variance of random factor U0 given to fixed parameter b0}
#'   \item{varU1}{variance of random factor U1 given to fixed parameter b1}
#'   \item{corr}{correlation between random factors U0 and U1}
#'   \item{varE}{the variability of the residuals}
#' }
#' @source Allen, JS, Bruss, J, Brown, CK, Damasio, H. Normal neuroanatomical variation due to age: the major lobes and a parcellation of the temporal region. Neurobiol Aging. 2005 Oct;26(9):1245-60; discussion 1279-82.
#'
#' Uher T, Vaneckova M, Krasensky J, Sobisek L, Tyblova M, Volna J, Seidl Z, Bergsland N, Dwyer MG, Zivadinov R, De Stefano N, Sormani MP, Havrdova EK, Horakova D. Pathological cut-offs of global and regional brain volume loss in multiple sclerosis. Mult Scler. 2017 Nov 1:1352458517742739. doi: 10.1177/1352458517742739.
#'
"ParamQuadrat"


#'Parameters of cubic model
#'
#' Default parameters to generate micro-panel (longitudinal) data with quadratic trend. The parameters may differ per each cluster.
#' The parameters of each cluster are in rows. Number of rows denotes the number of clusters.
#' Fixed effects are taken from Allen et al. (2005), and the source for random effects is Uher et al. (2017).
#' @format Its adviced to keep parameters in \code{\link[base]{data.frame}}. The Parameters structure is as follows:
#' \describe{
#'   \item{b0}{fixed parameter of intercept}
#'   \item{b1}{fixed parameter of slope}
#'   \item{b2}{fixed parameter of defining the quadraticity}
#'   \item{b3}{fixed parameter of defining the cubicity}
#'   \item{varU0}{variance of random factor U0 given to fixed parameter b0}
#'   \item{varU1}{variance of random factor U1 given to fixed parameter b1}
#'   \item{corr}{correlation between random factors U0 and U1}
#'   \item{varE}{the variability of the residuals}
#' }
#' @source Allen, JS, Bruss, J, Brown, CK, Damasio, H. Normal neuroanatomical variation due to age: the major lobes and a parcellation of the temporal region. Neurobiol Aging. 2005 Oct;26(9):1245-60; discussion 1279-82.
#'
#' Uher T, Vaneckova M, Krasensky J, Sobisek L, Tyblova M, Volna J, Seidl Z, Bergsland N, Dwyer MG, Zivadinov R, De Stefano N, Sormani MP, Havrdova EK, Horakova D. Pathological cut-offs of global and regional brain volume loss in multiple sclerosis. Mult Scler. 2017 Nov 1:1352458517742739. doi: 10.1177/1352458517742739.
#'
"ParamCubic"

#'Parameters of exponential model
#'
#' Default parameters to generate micro-panel (longitudinal) data with exponencial trend. The parameters may differ per each cluster.
#' The parameters of each cluster are in rows. Number of rows denotes the number of clusters.
#' Fixed effects are taken from Jones et al. (2013).
#' @format It is adviced to keep parameters in \code{\link[base]{data.frame}}. The Parameters structure is as follows:
#' \describe{
#'   \item{b0}{fixed parameter of intercept}
#'   \item{b1}{fixed parameter of slope}
#'   \item{b2}{fixed parameter of defining the decay}
#'   \item{varU0}{variance of random factor U0 given to fixed parameter b0}
#'   \item{varU1}{variance of random factor U1 given to fixed parameter b1}
#'   \item{corr}{correlation between random factors U0 and U1}
#'   \item{varE}{the variability of the residuals}
#' }
#' @source Jones BC, Nair G, Shea CD, Crainiceanu CM, Cortese IC, Reich DS. Quantification of multiple-sclerosis-related brain atrophy in two heterogeneous MRI datasets using mixed-effects modeling. Neuroimage Clin. 2013 Aug 13;3:171-9. doi: 10.1016/j.nicl.2013.08.001.
#'
"ParamExpon"
