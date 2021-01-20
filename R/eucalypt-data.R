#' Blue Mountains eucalypt and environmental data
#'
#' This is an extended version of the `BlueMountains` dataset from the package `ppmlasso`, formatted for use with `scampr` models. The list contains both a presence-only and presence/absence data set for a Eucalypt tree in the Greater Blue Mountains World Heritage Area. Data includes the values of four environmental variables and two variables related to site accessibility throughout the region at a spatial resolution of 500m.
#'
#' @format A list containing two data frames. Presence-only (PO) with 10230 rows and 10 variables. Presence/Absence (PA) with 250 rows and 9 variables:
#' \describe{
#'   \item{x}{UTM Easting coordinates (km)}
#'   \item{y}{UTM Northing coordinates (km)}
#'   \item{FC}{Number of fires since 1943}
#'   \item{D_MAIN_RDS}{Distance from the nearest main road (m)}
#'   \item{D_URBAN}{Distance from the nearest urban area (m)}
#'   \item{RAIN_ANN}{Average annual rainfall (mm)}
#'   \item{TMP_MAX}{Average maximum temperature (degrees Celsius)}
#'   \item{TMP_MIN}{Average minimum temperature (degrees Celsius)}
#'   \item{quad.size}{Size of the quadrat in Northing/Easting units squared. (PO data)}
#'   \item{pres}{Binary indicating whether the row is a presence record (1) or a quadrat (0). (PO data)}
#'   \item{Y}{Binary indicating whether species was present (1) or absent (0) at a survey site. (PA data)}
#' }
#' @source \url{Library `ppmlasso`}
"eucalypt"
