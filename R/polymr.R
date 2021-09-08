#' Perform polynomial Mendelian randomization
#'
#' @param x A vector containing the individual exposure values
#' @param y A vector containing the individual outcome values
#' @param G The genetic matrix
#' @param others Character vector containing other values to return.
#'     'observational' returns the polynomial approximation of the observed association,
#'     'offset' returns the scaling parameters of x and y,
#'     'means' returns the mean value of y for each bin while 'medians' returns their median.
#' @param bins Number of bins for which to return median or mean values
#' @param reverse_t Threshold to use for reverse causality filtering (T statistic), NULL for no filtering
#' @param ... Other parameters to pass to the non-linear methods
#'
#' @return List of results for PolyMR and the other selected values.
#'     'polymr' is a list containing 'pval0' (the p-value of the full model vs. no causal effect),
#'     'pval1' (the p-value for non-linearity),
#'     'r.squared' (the difference in explained variance full model vs. no causal effect),
#'     'coefficients' (the alpha parameter estimates and standard errors), and
#'     'vcov' (the variance-covariance matrix of the parameter estimates).
#'
#' @importFrom lmtest lrtest
#' @importFrom meta metagen
#' @importFrom stringr str_detect str_match str_replace str_extract
#'
#' @export
#'
polymr  =  function( x,
                     y,
                     G,
                     others  = c( 'observational' ),
                     bins = 100,
                     reverse_t = NULL,
                     ... ) {
  others  =  match.arg( others )
  if ('observational' %in% others) {
    results  =  list( offset = data.frame( mean = c( mean(x), mean(y) ),
                                           sd   = c(   sd(x),   sd(y) ),
                                           row.names = c( 'x', 'y' )) )
  } else {
    results = list()
  }

  x = scale(x)
  y = scale(y)

  # Removing invariant IVs
  G  =  G[ , apply( G, 2, var ) != 0 ]

  bx  =  summary(lm( x~G ))$coefficients[ 1:ncol(G)+1, 1:2 ]
  betax  =  bx[ , 1 ]

  if (!is.null( reverse_t )) {
    sx  =  bx[ , 2 ]
    betay  =  lm( y~G )
    if (any(is.na( betay$coefficients ))) {
      warning( 'Columns of G are not independent, i.e. some IVs have an LD of 1 in the provided sample, discarding duplicates.' )
      G  =  G[ , !is.na( betay$coefficients[-1] ) ]
    }
    betay  =  summary( betay )$coefficients[ -1, 1:2 ]

    to_keep  =  (abs( betax ) - abs( betay[ , 1 ] )) / sqrt( sx^2 + betay[ , 2 ]^2 ) > reverse_t
    G  =  G[ , to_keep, drop = FALSE ]
    betax  =  betax[ to_keep ]
  }

  epsX  =  x - G %*% betax

  results$polymr  =  .polymr_wrapper( x = x,
                                      y = y,
                                      epsX = epsX,
                                      ... )

  if ('medians' %in% others) {
    results$medians  =
      sapply( 1:bins, function( k ){
        bin  =  x > quantile( x, (k-1)/bins ) & x < quantile( x, k/bins )
        c( median(x[ bin ]),
           sd( y[bin] ),
           median(y[ bin ]) )
      } )
    rownames( results$medians )  =  c( 'x', 'y_sd', 'y_median' )
  }

  if ('means' %in% others | 'observational' %in% others) {
    results$means  =
      sapply( 1:bins, function( k ){
        bin  =  x > quantile( x, (k-1)/bins ) & x <= quantile( x, k/bins )
        c( median(x[ bin ]),
           sd( y[bin] ),
           mean( y[ bin ] ) )
      } )
    rownames( results$means )  =  c( 'x', 'y_sd', 'y_mean' )
  }

  if ('observational' %in% others) {
    results$observational  =  .model_wrapper( x = x,
                                              y = y,
                                              ... )
  }

  results
}

.polymr_wrapper  =  function( x,
                              y,
                              epsX,
                              powers = 1:10,
                              confounder_powers = 1:min(max(powers), max_confounder_power),
                              max_confounder_power = NULL,
                              ... ) {
  dat  =  data.frame( y = y,
                      x1 = x,
                      epsX1 = epsX )
  x_terms      =  paste0( 'x', powers )
  dat  =  .create_dat( dat,
                       c( x_terms, paste0( 'epsX', confounder_powers ) ) )

  .loop( dat = dat,
         powers = powers,
         method_function = .polymr_regression,
         coef_names = 'alpha',
         intercept_name = NULL,
         max_confounder_power = max_confounder_power,
         ... )
}

.create_dat  =  function( dat,
                          all_terms ){
  new_terms  =  all_terms[ !(all_terms %in% colnames(dat)) ]
  for (variable in unique( str_extract( new_terms, '^[a-zA-Z]+' ) )) {
    powers  =  str_extract( new_terms[ str_detect( new_terms, paste0('^', variable, '[0-9]+$') ) ], '[0-9]+$' )
    new_dat  =  data.frame( outer( dat[ , paste0(variable, 1) ],
                                   as.numeric( powers ),
                                   '^') )
    colnames( new_dat )  =  paste0( variable, powers )
    dat  =  cbind( dat, new_dat )
  }
  dat
}

.get_formula  =  function( all_terms,
                           remove = NULL,
                           add = NULL,
                           independent_variable = 'y' ){
  if (!is.null( remove )) {
    all_terms  =  all_terms[ !str_detect( all_terms, remove ) ]
  }

  f = paste0( independent_variable, ' ~ ',
              paste( c( add, all_terms ), collapse = ' + ' ) )
  as.formula( f )
}

.polymr_regression  =  function( dat,
                                 powers = 1:10,
                                 x_terms = paste0( 'x', powers ),
                                 other_terms = paste0( 'epsX', confounder_powers ),
                                 max_confounder_power = NULL,
                                 confounder_powers = 1:min(max(powers), max_confounder_power),
                                 ... ){
  if (any( !(c(x_terms, other_terms) %in% colnames(dat)) )) {
    dat  =  .create_dat( dat, c(x_terms, other_terms) )
  }
  model = lm( .get_formula( c(x_terms, other_terms) ), dat )
  to_return  =  list( model = model )

  pval1  =  .get_lrt_pval( lmtest::lrtest( model,
                                           update( model,
                                                   .get_formula( c( 'x1', other_terms ) ),
                                                   data = model$model ) ) )

  model0  =  update( model,
                     .get_formula( other_terms ),
                     data = model$model )

  list( pval0 = .get_lrt_pval( lmtest::lrtest( model, model0 ) ),
        pval1 = pval1,
        r.squared = summary( model )$r.squared - summary( model0 )$r.squared,
        dat = dat,
        model = model )
}

.get_lrt_pval  =  function( lrt ) {
  pchisq( 2 * ( lrt[ 1, 'LogLik' ] - lrt[ 2, 'LogLik' ] ),
          abs( lrt[[ '#Df' ]][ 1 ] - lrt[[ '#Df' ]][ 2 ] ),
          lower.tail = FALSE )
}


.loop  =  function( dat,
                    powers = 1:min(10, max_power),
                    ...,
                    method_function = .polymr_regression,
                    power_step = 2,
                    max_power = 10,
                    drop_lower = TRUE,
                    p_thr_add = .05/power_step,
                    p_thr_drop = .05/length( powers ),
                    return_vcov = TRUE,
                    coef_names = 'alpha',
                    intercept_name = NULL ){
  temp_results  =  final_results  =  method_function( dat = dat,
                                                      powers = powers,
                                                      ... )
  selection_terms  =  paste0( 'x', powers )
  while (min(summary( temp_results$model )$coefficients[ selection_terms, 'Pr(>|t|)' ]) < p_thr_add) {
    final_results  =  temp_results
    if (max(powers) < max_power) {
      new_powers  =  c(max(powers)+1:power_step)[ c(max(powers)+1:power_step) <= max_power ]
      selection_terms  =  paste0( 'x', new_powers )
      powers  =  c( powers, new_powers )
      temp_results  =  method_function( dat = dat,
                                        powers = powers,
                                        ... )
    } else {
      print( 'Max power reached' )
      break
    }
  }

  powers  =  stringr::str_match( names( final_results$model$coefficients ), '^x([0-9]+)$' )
  powers  =  powers[ !is.na(powers[,2]), 2 ]
  powers  =  as.numeric( powers )

  pvals  =  summary( final_results$model )$coefficients
  pvals  =  pvals[ paste0( 'x', powers ), 'Pr(>|t|)' ]

  while (drop_lower & any(pvals > p_thr_drop)) {
    if (length(powers) == 1) {
      return( c( list( coefficients = data.frame( estimate = numeric(), se = numeric() ) ),
                 final_results[ c( 'pval0', 'pval1' )[ c( 'pval0', 'pval1' ) %in% names( final_results ) ] ] ))
    }

    powers  =  powers[ -which.max( pvals ) ]
    final_results  =  method_function( dat = dat,
                                       powers = powers,
                                       ... )
    pvals  =  summary( final_results$model )$coefficients
    pvals  =  pvals[ paste0( 'x', powers ), 'Pr(>|t|)' ]
  }

  final_results$coefficients  =  summary( final_results$model )$coefficients

  final_results$coefficients  =  final_results$coefficients[ paste0( 'x', powers ),
                                                             c('Estimate', 'Std. Error'),
                                                             drop = FALSE ]
  rownames( final_results$coefficients )  =  stringr::str_replace( rownames( final_results$coefficients ),
                                                                   '^x',
                                                                   coef_names )

  colnames( final_results$coefficients )[ 1:2 ]  =  c( 'estimate', 'se' )

  if (return_vcov) {
    final_results$vcov  =  vcov( final_results$model )
    rownames( final_results$vcov )  =
      colnames( final_results$vcov )  =  stringr::str_replace( rownames( final_results$vcov ),
                                                               '^x',
                                                               coef_names )
  }
  final_results$model = final_results$dat = NULL

  final_results
}

.model_wrapper  =  function( y,
                             x,
                             powers = 0:10,
                             ... ){
  dat  =  data.frame( y = y,
                      x1 = x )
  dat  =  .create_dat( dat, paste0( 'x', powers ) )
  .loop( dat = dat,
         powers = unique( c( 0, powers ) ),
         method_function = .model_regression,
         coef_names = 'x',
         ... )
}

.model_regression  =  function( dat,
                                powers = 0:2,
                                null_model_powers = 0:1,
                                ... ) {
  model = lm( y ~ .-1, data = dat )
  to_return  =  list( model = model )

  if (is.null( null_model_powers )) {
    pval1 = NULL
  } else {
    pval1  =  lmtest::lrtest( model,
                              lm( .get_formula(paste0('x', null_model_powers)),
                                  data = model$model ) )
    pval1  =  .get_lrt_pval( pval1 )
  }

  list( model = model,
        pval1 = pval1,
        r.squared = summary( model )$r.squared )
}
















