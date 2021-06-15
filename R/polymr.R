#' Perform polynomial Mendelian randomization
#'
#' @param x A vector containing the individual exposure values
#' @param y A vector containing the individual outcome values
#' @param G The genetic matrix
#' @param others Character vector containing other values to return.
#'     'model' returns the polynomial approximation of the observed association,
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
#' @importFrom stringr str_detect str_match str_replace
#'
#' @export
#'
polymr  =  function( x,
                     y,
                     G,
                     others  = c( 'model', 'means', 'offset' ),
                     bins = 100,
                     reverse_t = NULL,
                     ... ) {

  if ('offset' %in% others) {
    results  =  list( offset = data.frame( mean = c( mean(x), mean(y) ),
                                           sd   = c(   sd(x),   sd(y) ),
                                           row.names = c( 'x', 'y' )) )
  } else {
    results = list()
  }

  x = scale(x)
  y = scale(y)

  G  =  G[ , apply( G, 2, var ) != 0 ]
  bx  =  summary(lm( x~G ))$coefficients[ 1:ncol(G)+1, 1:2 ]
  betax  =  bx[ , 1 ]
  sx  =  bx[ , 2 ]

  if (!is.null( reverse_t )) {
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

  results$polymr  =  .polymr_wrapper( list( x = x, y = y, epsX = epsX ),
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

  if ('means' %in% others) {
    results$means  =
      sapply( 1:bins, function( k ){
        bin  =  x > quantile( x, (k-1)/bins ) & x <= quantile( x, k/bins )
        c( median(x[ bin ]),
           sd( y[bin] ),
           mean( y[ bin ] ) )
      } )
    rownames( results$means )  =  c( 'x', 'y_sd', 'y_mean' )
  }

  if ('model' %in% others) {
    results$model  =  .model_wrapper( list( x = x, y = y ),
                                      ... )
  }

  results
}

.polymr_wrapper  =  function( data_list,
                              ... ) {
  .loop( data_list,
         method_function = .polymr_regression,
         intercept_name = NULL,
         coef_names = 'alpha',
         ... )
}

.polymr_regression  =  function( data_list,
                                 powers = 1:2,
                                 max_confounder_power = max(powers),
                                 confounder_powers = 1:min(max(powers), max_confounder_power),
                                 ... ){
  dat  =  cbind( data_list$y,
                 outer( c(data_list$x), powers, '^' ),
                 outer( c(data_list$epsX), confounder_powers, '^' ) )
  dat  =  data.frame( dat )
  colnames( dat )  =  c( 'y',
                         paste0( 'x', powers ),
                         paste0( 'epsX', confounder_powers ) )
  model = lm( y ~ ., data = dat )
  to_return  =  list( model = model )

  if (!(1 %in% powers)) {
    lrt1  =  lmtest::lrtest( model, lm( y ~ ., data = cbind( dat[ , -c( 1:length(powers)+1 ) ],
                                                             data.frame( x1 = data_list$x ) ) ) )
    pval1  =  .get_lrt_pval( lrt1 )

    if (length( powers ) == 1) {
      pval1  =  max( summary( model )$coefficients[ paste0( 'x', powers ), 'Pr(>|t|)' ],
                     pval1 )
    }
  } else {
    formula1  =  paste( '.~.', paste( '-x', setdiff( powers, 1 ), sep = '', collapse = '' ) )
    pval1  =  .get_lrt_pval( lmtest::lrtest( model, update( model, as.formula( formula1 ) ) ) )
  }

  formula0  =  as.formula( paste( '.~.', paste( '-x', powers, sep = '', collapse = '' ) ) )
  model0  =  update( model, formula0 )

  c( list( pval0 = .get_lrt_pval( lmtest::lrtest( model, model0 ) ),
           pval1 = pval1,
           r.squared = summary( model )$r.squared - summary( model0 )$r.squared ),
     to_return )
}

.get_lrt_pval  =  function( lrt,
                            infer_dfs = TRUE ) {
  if (infer_dfs) {
    df1  =  stringr::str_match( attributes( lrt )$heading[2], 'x([0-9]+)[^x]*Model' )
    df2  =  stringr::str_match( attributes( lrt )$heading[2], 'x([0-9]+)[^x]*$' )
    dfs  =  c(df1[ , 2 ], df2[ , 2 ])
    dfs  =  as.numeric( replace( dfs, is.na( dfs ), 0 ) )
    dfs  =  dfs[1] - dfs[2]
  } else {
    dfs  =  lrt[ 2, 'Df' ]
  }
  pchisq( 2 * ( lrt[ 1, 'LogLik' ] - lrt[ 2, 'LogLik' ] ),
          abs( dfs ),
          lower.tail = FALSE )
}

.loop  =  function( data_list,
                    method_function = .polymr_regression,
                    powers = 1:min(10, max_power),
                    power_step = 2,
                    max_power = 10,
                    drop_lower = TRUE,
                    p_thr_add = .05/power_step,
                    p_thr_drop = .05/length( powers ),
                    return_vcov = TRUE,
                    coef_names = 'alpha',
                    intercept_name = NULL,
                    ... ){
  temp_results  =  final_results  =  method_function( data_list,
                                                      powers = powers,
                                                      ... )

  while (min(summary( temp_results$model )$coefficients[ paste0( 'x', powers ), 'Pr(>|t|)' ]) < p_thr_add) {
    final_results  =  temp_results
    if (max(powers) < max_power) {
      new_powers  =  c(max(powers)+1:power_step)[ c(max(powers)+1:power_step) <= max_power ]
      powers  =  c( powers, new_powers )
      temp_results  =  method_function( data_list,
                                        powers = powers,
                                        ... )
    } else {
      print( 'Max power reached' )
      break
    }
  }

  if (!is.null( final_results$powers )) {
    powers  =  final_results$powers
  } else {
    powers  =  stringr::str_match( names( final_results$model$coefficients ), '^x([0-9]+)$' )
    powers  =  powers[ !is.na(powers[,2]), 2 ]
    powers  =  as.numeric( powers )
  }
  temp_results  =  final_results
  pval  =  summary( final_results$model )$coefficients
  pval  =  pval[ stringr::str_detect( rownames(pval), '^x[0-9]+' ), 'Pr(>|t|)' ]

  while (drop_lower & any(pval > p_thr_drop)) {
    if (length(powers) == 1) {
      return( c( list( coefficients = data.frame( estimate = numeric(), se = numeric() ) ),
                 final_results[ c( 'pval0', 'pval1' )[ c( 'pval0', 'pval1' ) %in% names( final_results ) ] ] ))
    }

    powers  =  powers[ -which.max( pval ) ]
    final_results  =  method_function( data_list,
                                       powers = powers,
                                       ... )
    pval  =  summary( final_results$model )$coefficients
    pval  =  pval[ stringr::str_detect( rownames(pval), '^x[0-9]+' ), 'Pr(>|t|)' ]
  }

  final_results$coefficients  =  summary( final_results$model )$coefficients
  rownames( final_results$coefficients )  =  stringr::str_replace( rownames( final_results$coefficients ),
                                                                   'x',
                                                                   coef_names )
  powers  =  stringr::str_match( names( final_results$model$coefficients ),
                                 '^x([0-9]+)$' )[ , 2 ]
  powers  =  as.numeric( powers[ !is.na(powers) ] )

  final_results$coefficients  =  final_results$coefficients[ paste0( coef_names, powers ),
                                                             c('Estimate', 'Std. Error'),
                                                             drop = FALSE ]

  colnames( final_results$coefficients )[ 1:2 ]  =  c( 'estimate', 'se' )

  if (return_vcov) {
    final_results$vcov  =  vcov( final_results$model )
    rownames( final_results$vcov )  =
      colnames( final_results$vcov )  =  stringr::str_replace( rownames( final_results$vcov ),
                                                               'x',
                                                               coef_names )
  }
  final_results$model = final_results$pval = final_results$powers = NULL

  final_results
}

.model_wrapper  =  function( data_list,
                             powers = 0:2,
                             ... ){
  .loop( data_list,
         method_function = .model_regression,
         powers = unique( c( 0, powers ) ),
         coef_names = 'x',
         ... )
}

.model_regression  =  function( data_list,
                                powers = 0:2,
                                null_model_powers = 0:1,
                                ... ) {
  dat  =  cbind( data_list$y,
                 outer( c(data_list$x), powers, '^' ) )
  dat  =  data.frame( dat )
  colnames( dat )  =  c( 'y',
                         paste0( 'x', powers ) )
  model = lm( y ~ .-1, data = dat )
  to_return  =  list( model = model )

  if (is.null( null_model_powers )) {
    return( c( to_return, list( pval = NULL ) ) )
  }

  if (!all(null_model_powers %in% powers)) {
    pval1  =  lmtest::lrtest( model,
                      lm( y ~ .-1, data = cbind( dat[ , 'y', drop = FALSE ],
                                                 outer( c(data_list$x), null_model_powers, '^' ) ) ) )
    pval1  =  .get_lrt_pval( pval1 )
    if (length( powers ) == length( null_model_powers ) & length( powers ) == 1) {
      pval1  =  max( summary( model )$coefficients[ paste0( 'x', powers ), 'Pr(>|t|)' ],
                     pval1 )
    }
  } else {
    formula1  =  as.formula( paste( '.~.', paste( '-x',
                                                  setdiff( powers, null_model_powers ),
                                                  sep = '',
                                                  collapse = '' ) ) )
    pval1  =  .get_lrt_pval( lmtest::lrtest( model, update( model, formula1 ) ) )
  }

  c( to_return,
     list( pval1 = pval1,
           r.squared = summary( model )$r.squared ) )
}
















