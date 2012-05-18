/* sample program for the gas data using loess */

#include <stdio.h>
#include "loess.h"

struct  loess_struct    gas, gas_null;
struct  pred_struct 	gas_pred;
struct  ci_struct       gas_ci;
struct  anova_struct    gas_anova;
double  NOx[] = {4.818, 2.849, 3.275, 4.691, 4.255, 5.064, 2.118, 4.602,
		 2.286, 0.97, 3.965, 5.344, 3.834, 1.99, 5.199, 5.283,
		 3.752, 0.537, 1.64, 5.055, 4.937, 1.561};
double  E[] = {0.831, 1.045, 1.021, 0.97, 0.825, 0.891, 0.71, 0.801, 
	       1.074, 1.148, 1, 0.928, 0.767, 0.701, 0.807, 0.902, 
	       0.997, 1.224, 1.089, 0.973, 0.98, 0.665};
double	gas_fit_E[] = {0.665, 0.949, 1.224};
double  newdata[] = {0.6650000, 0.7581667, 0.8513333, 0.9445000,
                     1.0376667, 1.1308333, 1.2240000};
double  coverage = .99;
long    i, n = 22, p = 1, m = 3, se_fit = FALSE;

main() {
	printf("\nloess(&gas):\n");
        loess_setup(E, NOx, n, p, &gas);
        gas.model.span = 2.0 / 3.0;
        loess(&gas);
	loess_summary(&gas);

	printf("\nloess(&gas_null):\n");
	loess_setup(E, NOx, n, p, &gas_null);
        gas_null.model.span = 1.0;
        loess(&gas_null);
	loess_summary(&gas_null);

	printf("\npredict(gas_fit_E, m, &gas, &gas_pred, %d):\n", se_fit);
	predict(gas_fit_E, m, &gas, &gas_pred, se_fit);
	for(i = 0; i < m; i++)
              printf("%g ", gas_pred.fit[i]);
	printf("\n");

        m = 7;
        se_fit = TRUE;
        predict(newdata, m, &gas, &gas_pred, se_fit);
	printf("\npointwise(&gas_pred, m, coverage, &gas_ci):\n");
        pointwise(&gas_pred, m, coverage, &gas_ci);
        for(i = 0; i < m; i++)
              printf("%g ", gas_ci.upper[i]);
        printf("\n");
        for(i = 0; i < m; i++)
              printf("%g ", gas_ci.fit[i]);
        printf("\n");
        for(i = 0; i < m; i++)
              printf("%g ", gas_ci.lower[i]);
        printf("\n");

	printf("\nanova(&gas_null, &gas, &gas_anova):\n");
	anova(&gas_null, &gas, &gas_anova);
	printf("%g %g %g %g\n", gas_anova.dfn, gas_anova.dfd,
	        gas_anova.F_value, gas_anova.Pr_F);

        loess_free_mem(&gas);
	loess_free_mem(&gas_null);
	pred_free_mem(&gas_pred);	
	pw_free_mem(&gas_ci);
}




