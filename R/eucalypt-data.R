#' Blue Mountains eucalypt and environmental data
#'
#' This is the `BlueMountains` dataset from the package `ppmlasso`, formatted for use with `scampr` models. Contains the observed presence locations of a Sydney eucalypt, the values of four environmental variables and two variables related to site accessibility throughout the region at a spatial resolution of 500m
#'
#' @format A data frame with 39131 rows and 10 variables:
#' \describe{
#'   \item{x}{UTM Easting coordinates (km)}
#'   \item{y}{UTM Northing coordinates (km)}
#'   \item{FC}{Number of fires since 1943}
#'   \item{D_MAIN_RDS}{Distance from the nearest main road (m)}
#'   \item{D_URBAN}{Distance from the nearest urban area (m)}
#'   \item{RAIN_ANN}{Average annual rainfall (mm)}
#'   \item{TMP_MAX}{Average maximum temperature (degrees Celsius)}
#'   \item{TMP_MIN}{Average minimum temperature (degrees Celsius)}
#'   \item{quad.size}{Size of the quadrat in Northing/Easting units squared.}
#'   \item{pres}{Binary indicating whether the row is a presence record (1) or a quadrat (0).}
#' }
#' @source \url{Library `ppmlasso`}
"eucalypt"
