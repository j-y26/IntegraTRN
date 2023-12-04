# Purpose: Running IntegraTRN as a Shiny app
# Author: Jielin Yang
# Date: 2023-11-28
# Version: 1.0
# Bugs and Issues: None

#' Launch Shiny App for package IntegraTRN
#'
#' @description This is a function that launches the Shiny app for the package
#'              IntegraTRN. The Shiny app permits the users to perform
#'              differential analysis of the user-input omics data as well as
#'              the visualization of the results. It also further allows the
#'              users to perform the network construction and visualization of
#'              the network based on input and inferred regulatory interactions.
#'
#' @return None, but opens the Shiny app
#'
#' @references
#' \insertRef{csardi2006igraph}{IntegraTRN}
#'
#' \insertRef{huynh2010inferring}{IntegraTRN}
#'
#' \insertRef{chang2020mirnet}{IntegraTRN}
#'
#' @export
#'
#' @importFrom shiny runApp
#' @import shinyBS
#' @import DT
#'
#' @examples
#' \dontrun{
#' runIntegraTRN()
#' }
#'
runIntegraTRN <- function() {
  shiny::runApp(system.file("shiny-scripts", package = "IntegraTRN"),
    display.mode = "normal"
  )
}
