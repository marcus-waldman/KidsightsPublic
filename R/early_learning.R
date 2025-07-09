#' obtained early learning scores
#'
#' @param input A data.frame with columns years and items in lex_kidsights (prefix = AA,BB,CC, or DD)..
#' @param min.responses Minimum number of responses to items for scoring an individual. Defaults to 5.
#' @param ... Other arguments passed to mirt::fscores
#' @return A data.frame (or list) of Kidsights scores
#' @export

early_learning<-function(input,  min.responses = 1, ...){

  library(tidyverse)
  library(mirt)
  # Parameters


  if(! ("id" %in% names(input)) ){stop("id must be a variable in input")}
  if(min(input$years)<3 | max(input$years>6)){stop("Ages must be between 3 and 6 years old")}

  # Clean up the input
  ee_items = internals[["version 2.0"]]$codebook %>% dplyr::filter(domain_hrtl == "early_learning") %>% purrr::pluck("lex_kidsight")

  # Get rid of observations with less than the minimum number of responses
    ians = input %>% dplyr::select(dplyr::any_of(ee_items)) %>% apply(1,function(x){sum(!is.na(x))})
    if(sum(ians<min.responses)>0){stop("Some observations have less than minimum responses.")}
    dat = input[ians>=min.responses, ]

  combined = hrtl_early_learning_calibdat %>% dplyr::mutate(wgt = 0) %>% dplyr::bind_rows(dat %>% dplyr::mutate(wgt = 1))

  #if(is.null(formula)){
  #  dentype = "EH"
  #  covdata_combined = NULL
  #} else {
    dentype = "Gaussian"
    covdata_combined = combined %>% dplyr::select(years,wgt)

    vals = mirt::mirt(
      data = combined %>% dplyr::select(dplyr::any_of(ee_items)),
      model = 1,
      formula = ~ wgt*years ,
      covdata = covdata_combined,
      pars = "values",
      dentype = dentype
    )

  vals = vals %>%
      dplyr::left_join(hrtl_internals[["2022"]][["Early Learning"]]$pars %>% dplyr::mutate(hat = value) %>% dplyr::select(item,name,hat))

   pars = vals %>% dplyr::mutate(
      value = ifelse(!is.na(hat), hat, value),
      est = !(class %in% c("dich","graded"))
    ) %>% dplyr::select(-hat)

  #}

   # Fit the empirical histogram


  fit_early_learning = mirt::mirt(
    data = combined %>% dplyr::select(any_of(ee_items)),
    model = 1,
    quadpts = 61*4,
    technical = list(theta_lim = c(-10,10), NCYCLES = 2000),
    #optimizer = "nloptr",
    #nloptr_args = list(opts = list("algorithm" = "NLOPT_LD_TNEWTON_PRECOND_RESTART")),
    pars = pars,
    dentype = dentype,
    TOL = 1E-3,
    formula = ~  wgt*years,
    covdata = covdata_combined,
    large = T
  )

  print(tail(mirt::mod2values(fit_early_learning)))


  scores<-do.call(mirt::fscores,
                  args = list(
                    object = fit_early_learning,
                    use_denstype_estimate = T,
                    theta_lim = c(-10,10),
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




