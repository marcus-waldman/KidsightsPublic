#' Score Kidsights
#'
#' @param input A data.frame with columns years and items in lex_kidsights (prefix = AA,BB,CC, or DD)..
#' @param version  Defaults to the most current version ("version 2.0"). Other options include "version 1.0".
#' @param ... Other arguments passed to mirt::fscores
#' @return A data.frame (or list) of Kidsights scores
#' @export

fscores<-function(input, version = "version 2.0",...){

  library(tidyverse)
  # Parameters
  pars_mirt = internals[[version]]$pars

  combined = calibdat %>% dplyr::mutate(wgt = 0) %>% dplyr::bind_rows(input %>% dplyr::mutate(wgt = 1))

  # Fit the empirical histogram
  fit_kidsight = mirt::mirt(
    data = combined %>% dplyr::select(dplyr::any_of(internals[[version]]$codebook$lex_kidsight %>% na.omit())),
    model = 1,
    quadpts = 61*4,
    technical = list(theta_lim = c(-20,15), NCYCLES = 2000),
    optimizer = "NR",
    pars = internals$`version 2.0`$pars,
    dentype = "EH",
    TOL = 1E-5,
    survey.weights = combined$wgt
  )

  scores<-do.call(mirt::fscores,
                  args = list(
                    object = fit_kidsight,
                    use_denstype_estimate = T,
                    response.pattern =  combined %>% dplyr::filter(wgt==1) %>% dplyr::select(dplyr::any_of(internals[[version]]$codebook$lex_kidsight %>% na.omit())),
                    theta_lim = c(-20,15),
                    quadpts = 61*4,
                    ...)
                  )

  return(scores)

}




