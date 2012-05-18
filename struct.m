(All default values mentioned here are set by loess_setup().)

struct  loess_struct    *lo;

in	
	n:		number of observations.

	p:		number of numeric predictors.
	
	y:		vector of response (length n).

	x:		vector of predictors, of length (n * p).
			The j-th coordinate of the i-th point is in x[i+n*j],
			where 0<=j<p, 0<=i<n.

	weights:	vector of weights to be given to individual 
			observations in the sum of squared residuals that 
			forms the local fitting criterion.
			By default, an unweighted fit is carried out.
			If supplied, weights should be a non-negative 
			numeric vector.  If the different observations 
			have non-equal variances, weights should be 
			inversely proportional to the variances.

model
	span:		smoothing parameter. Default is 0.75.

	degree:		overall degree of locally-fitted polynomial. 1 is 
			locally-linear fitting and 2 is locally-quadratic 
			fitting.  Default is 2.

	normalize:	logical that determines if numeric predictors should 
			be normalized.  If TRUE (1) - the default - the 
			standard normalization is used. If FALSE (0), no
			normalization is carried out.

	parametric:	for two or more numeric predictors, this argument
			specifies those variables that should be 
			conditionally-parametric. The argument should be a 
			logical vector of length p, specified in the order 
			of the predictor group ordered in x.
			Default is a vector of 0's of length p.

	drop_square:	for cases with degree = 2, and with two or more 
			numeric predictors, this argument specifies those 
			numeric predictors whose squares should be dropped 
			from the set of fitting variables. The method of 
			specification is the same as for parametric.
			Default is a vector of 0's of length p.

	family:		the assumed distribution of the errors. The values 
			are "gaussian" or "symmetric". The first value is 
			the default.  If the second value is specified, 
			a robust fitting procedure is used.

control
	surface:	determines whether the fitted surface is computed 
			directly at all points ("direct") or whether an 
			interpolation method is used ("interpolate"). 
			The latter, the default, is what most users should 
			use unless special circumstances warrant.

        statistics:	determines whether the statistical quantities are 
			computed exactly ("exact") or approximately 
			("approximate"). The latter is the default. The former
			should only be used for testing the approximation in 
			statistical development and is not meant for routine 
			usage because computation time can be horrendous.

        cell:		if interpolation is used to compute the surface, this
			argument specifies the maximum cell size of the k-d 
			tree.  Suppose k = floor(n*cell*span) where n is the 
			number of observations.  Then a cell is further 
			divided if the number of observations within it
			is greater than or equal to k.

	trace_hat:	when surface is "approximate", determines the 
			computational method used to compute the trace of 
			the hat matrix, which is used in the computation of 
			the statistical quantities.  If "exact", an exact 
			computation is done; normally this goes quite fast 
			on the fastest machines until n, the number of 
			observations is 1000 or more, but for very slow 
			machines, things can slow down at n = 300.  
			If "wait.to.decide" is selected, then a default 
			is chosen in loess();  the default is "exact" for 
			n < 500 and "approximate" otherwise.  If surface 
			is "exact", an exact computation is always done 
			for the trace. Set trace_hat to "approximate" for 
			large dataset will substantially reduce the 
			computation time.

	iterations:	if family is "symmetric", the number of iterations 
			of the robust fitting method.  Default is 0 for
			family being "gaussian" by default.

kd_tree:	k-d tree, an output of loess().

out	
	fitted_values:	fitted values of the local regression model

	fitted_residuals:	residuals of the local regression fit

        enp:		equivalent number of parameters.

        s:		estimate of the scale of the residuals.

        one_delta:	a statistical parameter used in the computation of 
			standard errors.

        two_delta:	a statistical parameter used in the computation of 
			standard errors.

        pseudovalues:	adjusted values of the response when robust 
			estimation is used.

	trace_hat:	trace of the operator hat matrix.

        diagonal:	diagonal of the operator hat matrix.

        robust:		robustness weights for robust fitting.

        divisor:	normalization divisor for numeric predictors.


struct  pred_struct	*pre;

	fit: 		the evaluated loess surface at eval.

	se_fit:		estimates of the standard errors of the surface values.

	residual_scale: estimate of the scale of the residuals.
	
	df:    		the degrees of freedom of the t-distribution used to
		        compute pointwise confidence intervals for the 
			evaluated surface. 


struct  anova_struct	*aov;
	
	dfn:		degrees of freedom of the numerator.

	dfd:		degrees of freedom of the denominator.

	F_values:	F statistic.

	Pr_F:		probability F_value is exceeded if null hypothesis
			is true.


struct	ci_struct	*ci;

	fit:		the evaluated loess surface at eval (see pred_struct).

	upper:		upper limits of pointwise confidence intervals.

	lower:		lower limits of pointwise confidence intervals.






