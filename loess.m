NAME

	loess_setup, loess, loess_summary, loess_free_mem, anova

SYNOPSIS

	#include "loess.h"
	double  *x, *y;
	long    n, p;
        struct  loess_struct    *lo, *lo2;
	struct  anova_struct    *aov;
	
	void 	loess_setup(x, y, n, p, lo)

	void	loess(lo)

	void	loess_summary(lo)

	void	loess_free_mem(lo)

        void    anova(lo, lo2, aov);

PARAMETERS

	x	predictors vector (of length n * p)
		The j-th coordinate of the i-th point is in x[i+n*j],
		where 0<=j<p, 0<=i<n.
 
	y	response vector (of length n).

	n	number of observations.

	p	number of variables/predictors.

	lo	copy of data;  controls;  k-d tree and coefficients.

	aov	results of the F-test in the analysis of variance.

DESCRIPTION

	loess_setup() sets up all default values in loess_struct's in, 
	model, and control structures; it also allocates memory for the 
	kd_tree and out structures based on n and p.  Caller can then
	override any of these parameters by explicitly redefining them
	before the call to loess() (see sample.c).  loess_setup()
	has the side-effect of copying x, y, n, and p into the in
	structure for ease of arguments-passing in subsequent calls to
	other loess and predict routines.

	loess() takes this structure, and does the actual loess
	computation.  It stored the results in the out structure.

	loess_summary() is a simple utility routine that summarizes the
	results of the loess computation.  Since it takes in the whole 
	loess structure as its argument, it has the potential of printing 
	out any parameter of interest with only a slight modification to 
	the code.
	
	loess_free_mem() frees up all dynamically allocated memory 
	used by the loess structure.

	anova() performs an analysis of variance on two loess models, and
	stores the results of the F-test in the anova_struct structure.

	loess_struct and anova_struct are defined in loess.h and documented 
	in struct.m.  Although the internal arrays are allocated by
	loess_setup(), the struct arguments (lo, lo2, aov) should be
	allocated by the caller.  Thus a typical call would be
		struct loess_struct lo;
		loess_setup(x,y,n,p.&lo);

SEE ALSO

	predict, pointwise, pred_free_mem, pw_free_mem

