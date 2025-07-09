#' Score Kidsights
#'
#' @param input A data.frame with columns years and items in lex_kidsights (prefix = AA,BB,CC, or DD)..
#' @param version  Defaults to the most current version ("version 2.0"). Other options include "version 1.0".
#' @param min.responses Minimum number of responses to items for scoring an individual. Defaults to 5.
#' @param ... Other arguments passed to mirt::fscores
#' @return A data.frame (or list) of Kidsights scores
#' @export

kidsights_scores<-function(input,  version = "version 2.0", min.responses = 1, ...){

  library(tidyverse)
  # Parameters

  if(! ("id" %in% names(input)) ){stop("id must be a variable in input")}

  # Clean up the input
    kidsight_items = internals[[version]]$codebook$lex_kidsight %>% na.omit()
    # Get rid of observations with less than the minimum number of responses
    ians = input %>% dplyr::select(dplyr::any_of(kidsight_items)) %>% apply(1,function(x){sum(!is.na(x))})
    if(sum(ians<min.responses)>0){stop("Some observations have less than the minimum responses.")}

  combined = calibdat %>% dplyr::mutate(wgt = 0) %>% dplyr::bind_rows(input %>% dplyr::mutate(wgt = 1))

  #if(is.null(formula)){
  #  dentype = "EH"
  #  covdata_combined = NULL
  #} else {
    dentype = "Gaussian"
    covdata_combined = combined %>% dplyr::select(years,wgt)

    vals = mirt::mirt(
      data = combined %>% dplyr::select(dplyr::any_of(kidsight_items)),
      model = 1,
      formula = ~ 1 + log(years + .1) +  wgt + wgt*years ,
      covdata = covdata_combined,
      pars = "values",
      dentype = dentype
    )

  vals = vals %>%
      dplyr::left_join(internals[[version]]$pars %>% dplyr::mutate(hat = value) %>% dplyr::select(item,name,hat))

   pars = vals %>% dplyr::mutate(
      value = ifelse(!is.na(hat), hat, value),
      est = !(class %in% c("dich","graded"))
    ) %>% dplyr::select(-hat)

  #}

   # Fit the empirical histogram


  fit_kidsight = mirt::mirt(
    data = combined %>% dplyr::select(any_of(kidsight_items)),
    model = 1,
    quadpts = 61*4,
    technical = list(theta_lim = c(-20,15), NCYCLES = 2000),
    #optimizer = "nloptr",
    #nloptr_args = list(opts = list("algorithm" = "NLOPT_LD_TNEWTON_PRECOND_RESTART")),
    pars = pars,
    dentype = dentype,
    TOL = 1E-3,
    formula = ~ 1 + log(years + .1) +  wgt + wgt*years ,
    covdata = covdata_combined,
    large = F
  )

  print(tail(mirt::mod2values(fit_kidsight)))


  scores<-do.call(mirt::fscores,
                  args = list(
                    object = fit_kidsight,
                    use_denstype_estimate = T,
                    theta_lim = c(-20,15),
                    quadpts = 61*4,
                    ...)
                  ) %>%
    dplyr::bind_cols() %>%
    dplyr::mutate(id = combined$id, wgt = combined$wgt) %>%
    dplyr::filter(wgt==1) %>%
    dplyr::select(-wgt) %>%
    dplyr::relocate(id)

  return(scores)

}




