#' Blue Mountains flora, environmental and site accessibility data
#'
#' This is an extended version of the `BlueMountains` dataset from the package `ppmlasso`, formatted for use with `scampr` models. The list contains both a presence-only and presence/absence data sets for four species of flora in the Greater Blue Mountains World Heritage Area (GBMWHA). Data includes the values of three environmental variables and two variables related to site accessibility throughout the region at a spatial resolution of 1km. All predictor variables are centered and scaled for ease of use in modelling (center and scale can be recovered from attributes attached to each relevant column of the data frames).
#'
#' @format A list of three elements. First element ('po') is a list of four data frames of presence locations of four species of flora ('sp1', 'sp2', 'sp3', 'sp4') with 242, 38, 11 and 194 rows resp. and each with 9 columns. Second element ('pa') is a data frame of presence/absence (survey) data with 8,223 rows and 11 columns. Third element ('quad') a data frame of quadrature points that describe the GBMWHA at a 1km resolution, 86,227 rows and 9 columns (the region's boundary is attached as an attribute).
#' \describe{
#'   \item{x}{UTM Easting coordinates (scaled km)}
#'   \item{y}{UTM Northing coordinates (scaled km)}
#'   \item{D.Main}{Distance from the nearest main road (scaled km)}
#'   \item{D.Urb}{Distance from the nearest urban area (scaled km)}
#'   \item{Rain}{Average annual rainfall (scaled mm)}
#'   \item{MXT}{Average maximum temperature (scaled degrees Celsius)}
#'   \item{MNT}{Average minimum temperature (scaled degrees Celsius)}
#'   \item{quad.size}{Size of the quadrat in Northing/Easting units squared. ('po' and 'quad')}
#'   \item{pres}{Binary indicating whether the row is a presence record (1) or a quadrat (0). ('po' and 'quad')}
#'   \item{sp1}{Binary indicating whether species 1 was present (1) or absent (0) at a survey site. ('pa')}
#'   \item{sp2}{Binary indicating whether species 2 was present (1) or absent (0) at a survey site. ('pa')}
#'   \item{sp3}{Binary indicating whether species 3 was present (1) or absent (0) at a survey site. ('pa')}
#'   \item{sp4}{Binary indicating whether species 4 was present (1) or absent (0) at a survey site. ('pa')}
#' }
#' @source \url{Library `ppmlasso`}
"flora"
