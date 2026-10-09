#' Example Bathymetric data
#'
#' An example of bathymetric data as a `matrix`
#'
#' @format a `matrix` containing `numerical` values of bathymetry with deminsion of 1:603 and 1:1151.
#'
"bathy_matrix"

#' Event log data
#'
#' An example of psat transmitter event log as a `data.frame`.
#'
#' @format `data.frame` containing 3 rows and 8 variables
#'  \describe{
#'    \item{fish_id}{The id of a given tagged fish}
#'    \item{ptt}{The id of a psat tag}
#'    \item{mrpat}{The id of a mrpat tag}
#'    \item{length}{The length of the invidiual}
#'    \item{date}{The date of an event}
#'    \item{lon}{The longitude of an event}
#'    \item{lat}{The latitude of an event}
#'    \item{event}{The event type}
#' }
"event_log"

#' PSAT log data
#'
#' An example of psat transmitter log file as a `data.frame`.
#'
#' @format `data.frame` containing 19,462 rows and 17 variables
#'  \describe{
#'    \item{fish_id}{The id of a given tagged fish}
#'    \item{ptt}{The id of a psat tag}
#'    \item{depth_sensor}{The type of sensor e.g., (`0.5`)}
#'    \item{source}{Transmission source}
#'    \item{instr}{The intstrument type}
#'    \item{date}{The date recored}
#'    \item{time}{The tiem recorded}
#'    \item{date_time}{The date timestamp}
#'    \item{location_quality}{The quality of the location}
#'    \item{latitude}{The latitude of the transmision}
#'    \item{depth}{The depth of the transmitter}
#'    \item{d_range}{The depth range}
#'    \item{temperature}{The temperature of the the transmitter}
#'    \item{t_range}{The temeprature range}
#'    \item{activity}{An activity measurement}
#'    \item{a_range}{An activity range}
#' }

"pdt"
