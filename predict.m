NAME
	predict, pointwise, pred_free_mem, pw_free_mem

SYNOPSIS
	#include "loess.h"

	double  *eval, coverage;
	long    m, se;
	struct  loess_struct    *lo;
	struct  predict_struct  *pre;
	struct  ci_struct 	*ci;

	void	predict(eval, m, lo, pre, se)
	
	void	pointwise(pre, m, coverage, ci)

	void	pred_free_mem(pre)

	void	pw_free_mem(ci)

PARAMETERS

	eval	a vector of length m * p specifying the values of the
		predictors at which the evaluation is to be carried out.
		The j-th coordinate of the i-th point is in eval[i+m*j],
		where 0<=j<p, 0<=i<m.

	m	number of evaluations.

	lo	k-d tree and coefficients.

	pre	predicted values; output by predict(), input to pointwise().

	se	logical flag for computing standard errors at eval.

	ci	pointwise confidence limits.

	coverage  (input) confidence level of the limits as a fraction.

DESCRIPTION

	predict() takes all the loess structures from earlier calls to
	loess_setup() and loess(), does an evaluation based on 
	eval and m, and stores the results in the pre structure.
	if se is TRUE, then pre.se_fit are computed along with the 
	surface values, pre.fit.  These returned vectors 
	are vectors of the same structural arrangement as that of eval.

	pointwise() computes the pointwise confidence limits from the
	result of predict().
	
	pred_free_mem() and pw_free_mem() frees up the allocated memory 
	used by the pre and ci structures respectively.

	loess_struct, pred_struct, and ci_struct are defined in loess.h 
	and documented in struct.m. 

NOTES

	The computations of predict() that produce the component se_fit
	are much more costly than those that producing the fit,
	so the number of points at which standard errors are
	computed should be modest compared to those at which we do
	evaluations.  Often this means calling predict() twice,
	once at a large number of points, with se = FALSE,
	to get a thorough description of the surface; and once 
	at a small number of points, with se = TRUE,
	to get standard-error information.

	Suppose the computation method for loess surfaces is
	interpolate, the default for the argument surface. Then the
	evaluation values of a numeric predictor must lie within
	the range of the values of the predictor used in the fit.

SEE ALSO

	loess_setup, loess, loess_summary, loess_free_mem
