#include <stdio.h>
#include "loess.h"

struct  loess_struct    ethanol, ethanol_cp;
struct  pred_struct     ethanol_pred, ethanol_grid;
struct  ci_struct       ethanol_ci;
double	NOx[] = {3.741, 2.295, 1.498, 2.881, 0.76, 3.12, 0.638, 1.17, 2.358, 
		 0.606, 3.669, 1, 0.981, 1.192, 0.926, 1.59, 1.806, 1.962, 
		 4.028, 3.148, 1.836, 2.845, 1.013, 0.414, 0.812, 0.374, 3.623,
		 1.869, 2.836, 3.567, 0.866, 1.369, 0.542, 2.739, 1.2, 1.719, 
		 3.423, 1.634, 1.021, 2.157, 3.361, 1.39, 1.947, 0.962, 0.571,
		 2.219, 1.419, 3.519, 1.732, 3.206, 2.471, 1.777, 2.571, 3.952,
		 3.931, 1.587, 1.397, 3.536, 2.202, 0.756, 1.62, 3.656, 2.964,
		 3.76, 0.672, 3.677, 3.517, 3.29, 1.139, 0.727, 2.581, 0.923, 
		 1.527, 3.388, 2.085, 0.966, 3.488, 0.754, 0.797, 2.064, 3.732,
		 0.586, 0.561, 0.563, 0.678, 0.37, 0.53, 1.9};
double	C_E[] = {12, 12, 12, 12, 12, 9, 9, 9, 12, 12, 12, 12, 15, 18, 7.5, 12, 
	       12, 15, 15, 9, 9, 7.5, 7.5, 18, 18, 15, 15, 7.5, 7.5, 9, 15, 15,
	       15, 15, 15, 9, 9, 7.5, 7.5, 7.5, 18, 18, 18, 18, 9, 9, 9, 9, 
	       7.5, 7.5, 7.5, 15, 18, 18, 15, 15, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5,
	       7.5, 18, 18, 18, 12, 12, 9, 9, 9, 15, 15, 15, 15, 15, 7.5, 7.5,
	       9, 7.5, 18, 18, 7.5, 9, 12, 15, 18, 18,
	       0.907, 0.761, 1.108, 1.016, 1.189, 1.001, 1.231, 1.123, 1.042,
               1.215, 0.93, 1.152, 1.138, 0.601, 0.696, 0.686, 1.072, 1.074,
               0.934, 0.808, 1.071, 1.009, 1.142, 1.229, 1.175, 0.568, 0.977,
               0.767, 1.006, 0.893, 1.152, 0.693, 1.232, 1.036, 1.125, 1.081,
               0.868, 0.762, 1.144, 1.045, 0.797, 1.115, 1.07, 1.219, 0.637,
               0.733, 0.715, 0.872, 0.765, 0.878, 0.811, 0.676, 1.045, 0.968,
               0.846, 0.684, 0.729, 0.911, 0.808, 1.168, 0.749, 0.892, 1.002,
               0.812, 1.23, 0.804, 0.813, 1.002, 0.696, 1.199, 1.03, 0.602,
               0.694, 0.816, 1.037, 1.181, 0.899, 1.227, 1.18, 0.795, 0.99,
               1.201, 0.629, 0.608, 0.584, 0.562, 0.535, 0.655};
double  newdata[] = {7.5, 9.0, 12.0, 15.0, 18.0, 0.6, 0.8, 1.0, 0.8, 0.6};
double	Cmin = 7.5, Cmax = 18.0, Emin = 0.535, Emax = 1.232;
double  Cm[7], Em[16], grid[224];
double  tmp, coverage = .99;
long    n = 88, p = 2, m = 5, se_fit = FALSE;
int     i, j, k;

main() {
        printf("\nloess(&ethanol): (span = 0.5)\n");
        loess_setup(C_E, NOx, n, p, &ethanol);
        ethanol.model.span = 0.5;
        loess(&ethanol);
        loess_summary(&ethanol);        
	
	printf("\nloess(&ethanol): (span = 0.25)\n");
        ethanol.model.span = 0.25;
        loess(&ethanol);
        loess_summary(&ethanol);

	printf("\nloess(&ethanol_cp): (span = 0.25)\n");
        loess_setup(C_E, NOx, n, p, &ethanol_cp);
        ethanol_cp.model.span = 0.25;
	ethanol_cp.model.parametric[0] = TRUE;
	ethanol_cp.model.drop_square[0] = TRUE;
        loess(&ethanol_cp);
        loess_summary(&ethanol_cp);

	printf("\nloess(&ethanol_cp): (span = 0.5)\n");
        ethanol_cp.model.span = 0.5;
        loess(&ethanol_cp);
        loess_summary(&ethanol_cp);

        printf("\npredict(newdata, m, &ethanol, &ethanol_pred, %d):\n", se_fit);
	predict(newdata, m, &ethanol_cp, &ethanol_pred, se_fit);
	for(i = 0; i < m; i++)
	        printf("%g ", ethanol_pred.fit[i]);
	printf("\n");

	m = 112;
	se_fit = TRUE;
	tmp = (Cmax - Cmin) / 6;
	for(i = 0; i < 7; i++) 
		Cm[i] = Cmin + tmp * i;
	tmp = (Emax - Emin) / 15;
	for(i = 0; i < 16; i++)
		Em[i] = Emin + tmp * i;
	for(i = 0; i < 16; i++) {
		k = i * 7;
		for(j = 0; j < 7; j++) {
			grid[k + j] = Cm[j];
			grid[m + k + j] = Em[i];
		}
	}
	predict(grid, m, &ethanol_cp, &ethanol_grid, se_fit);
	pointwise(&ethanol_grid, m, coverage, &ethanol_ci);

        loess_free_mem(&ethanol);
	loess_free_mem(&ethanol_cp);
        pred_free_mem(&ethanol_pred);     
        pred_free_mem(&ethanol_grid);     
}
