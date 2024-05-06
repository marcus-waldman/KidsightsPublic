#' Score Kidsights
#'
#' @param input A data.frame with columns years and items in lex_kidsights (prefix = AA,BB,CC, or DD)..
#' @param est_hyperpriors Defaults to True. Estimates the hyperpriors meand and variances based on available data.
#' @param plausible_values Defaults to 0, in which case the EAP estimate is returned.
#' @return A data.frame (or list) of EAP scores (or plausible values).
#' @export

fscores<-function(input, est_hyperpriors = T, plausible_values = 0){


  pars_mirt = internals$pars_mirt
  if(!est_hyperpriors){pars_mirt$est = F}

  fit_kidsights <-
    mirt::mirt(
      data = input |> dplyr::select(starts_with("AA"),starts_with("BB"),starts_with("CC"),starts_with("DD")),
      model = 1,
      covdata = input |> dplyr::select(years),
      formula = ~ poly(years,4),
      quadpts = 61*2,
      TOL = 1E-5,
      technical = list(theta_lim = c(-15,10), NCYCLES = 2000),
      pars = internals$pars
    )

  scores<-mirt::fscores(object = fit_kidsights, plausible.draws = plausible_values)

  if(plausible_values==0){
    output = data.frame(id = input$id, scores = as.vector(scores))
  } else{
    output = scores
  }

  return(output)


}




