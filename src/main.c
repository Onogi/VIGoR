#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "header1.h"
#include "header2.h"
#include "Rdynload.h"

/*
VIGoR function for R package.

Ver.1.1.0

=====Modifications from Ver. 1.0=====
==Major modifications
(1) GBLUP was added.
(2) Multiple regression methods can be included in the same regression models. That is,
   Y = f(Xf) + g(Xg) + h(Xh) + e
   when the number of methods is 3.

==Minor things
(1) Response variables (Y) should be standardized in the wrapper R function
(2) N (number of records or Y) and Nn (number of non-missing records) were changed to Ny.
(3) Y and X/K are linked with YtoX and XtoY. Individuals in Y and X/K should be one to one correspondence, but the order is arbitrary.

Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata

Released under the MIT license
http://opensource.org/licenses/mit-license.php

*/
//__declspec(dllexport)
static void vigor_c(int *Priortype, int *Methodcode, int *Nm, int *P, int *Ny, int *Nx, int *YtoX, int *XtoY,
                    int *Maxiteration, int *RandomIni, int *Division_G, int *Division_H, int *Division_E, int *Division_V,
                    double *Thresholdvalue, double *Observations, double *Genotypes,  double *Hyperparameters,
                    double *Tau0, double *LB, double *Rmonitor, double *Expectation, double *Uncertainty, double *Variance, double *Gamma,
                    double *Fittedvalue)
{
	/* CAUTION: In XtoY, enter Missingvalue for elements which are not linked with Y.*/

	/* Observations, genotypes, covariates, and hyperparameters */
	Ystruct *Y = NULL;
	Xstruct *X = NULL;
	Hstruct *H = NULL;

	/* Repeat statement */
	int	method, record;

	/* Calculation time */
	double	FittingTime = 0.0;

	/* Used seed */
	unsigned long usedseed;

	/* Allocate memory */
	Y = (Ystruct*)calloc(1, sizeof(Ystruct));
	Y[0].stobs = (double*)calloc(Ny[0], sizeof(double));
	memcpy(Y[0].stobs, Observations, sizeof(double) * Ny[0]);
	Y[0].expErrors = (double*)calloc(Ny[0], sizeof(double));
	memcpy(Y[0].expErrors, Y[0].stobs, sizeof(double) * Ny[0]);

	X = (Xstruct*)calloc(Nm[0], sizeof(Xstruct));
	H = (Hstruct*)calloc(Nm[0], sizeof(Hstruct));
	for (method = 0; method < Nm[0]; method++)
	{
		AllocateMemory(Priortype[method], Methodcode[method], Nx[0], P[method], X + method, Genotypes + Division_G[method]);
		AssignHyperparameters(Methodcode[method], H + method, Hyperparameters + Division_H[method]);
	}

	/* Initialize the seed */
	usedseed = (unsigned long)time(NULL);
	init_genrand(usedseed);

	/* Calculation */
	GenomeWideRegression(Methodcode, Nm[0], P, Ny[0], Nx[0], YtoX, XtoY, Thresholdvalue[0], Maxiteration[0],
		Y, X, H, Tau0, LB, Rmonitor, &FittingTime, RandomIni[0]);

	/* Copy results */
	for (method = 0; method < Nm[0]; method++)
		CopyResults(Methodcode[method], P[method], Nx[0], X + method, Expectation + Division_E[method], Uncertainty + Division_E[method],
			Variance + Division_V[method], Gamma + Division_E[method]);

	/* Fitted values */
	for (record = 0; record < Ny[0]; record++)
		Fittedvalue[record] = Y[0].stobs[record] - Y[0].expErrors[record];

	/* free */
	free(Y[0].stobs);
	free(Y[0].expErrors);
	free(Y);
	for (method = 0; method < Nm[0]; method++)
		FreeMemory(Priortype[method], Methodcode[method], X + method);
	free(X);
	free(H);
}

static const
R_CMethodDef cMethods[] = {
    { "vigor_c", (DL_FUNC)&vigor_c, 26 },
    {NULL,0}
};

void R_init_VIGoR(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

