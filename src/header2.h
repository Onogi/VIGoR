#pragma once
/*

Ver. 1.1

Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata

Released under the MIT license
http://opensource.org/licenses/mit-license.php

*/

/*----define global objects------------------------------------------------------------------------------------------------------*/
#define	Missingvalue		99999999			/* Missing value used in this program */
#define log2pi				1.83787706640935

/*----define structure-----------------------------------------------------------------------------------------------------------*/
typedef struct {

	double	*stobs;			/* standardized observations */
	double	*expErrors;		/* expected errors */
	char	name[101];		/* name of trait */

} Ystruct;

typedef struct {

	double	*covariates;	/* covariates */
	double	*x2;			/* sum of squared covariates */
	double	*expEffect;		/* expectation of marker effect. When BayesC, this indicates E[B|gamma=1] */
	double	*varEffect;		/* variance of marker effect. When BayesC, this indicates V[B|gamma=1] */
	double	*exp2Effect;	/* 2nd moment (for wBSR, BayesC, SSVS, MIX, and GBLUP). When BayesC, this indicates E[B^2|gamma=1] */
	double	*expGamma;		/* expectation of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	*exp2Gamma;		/* 2nd moment (for wBSR) */
	double	*expTau2;		/* expectation of tau (for BL and EBL) */
	double	*expInTau2;		/* expectation of the inverse of tau (for BL and EBL) */
	double  expDelta2;      /* expectation of delta (for BL and EBL). This is also used to store the determinant of K for GBLUP*/
	double	*expEta2;		/* expectation of eta (for EBL) */
	double	*expSigma2;		/* expectation of sigma (for wBSR, BayesC, SSVS, MIX, and GBLUP).*/
	double	*S2;			/* parameter of the posterior distribution of marker effect variance (for wBSR, BayesC, SSVS, and MIX) */
	double	*a3;			/* shape parameter of the posterior gamma distribution of eta2 (for EBL) */
	double	*b3;			/* rate parameter of the posterior gamma distribution of eta2 (for EBL) */
	double	*evalue;		/* eigen values of iK (inverse of the relationship matrix) (for GBLUP) */
	double	*evector;		/* eigen vectors of iK (inverse of the relationship matrix) (for GBLUP) */
	double	*tevector;		/* transposed eigen vectors of iK (for GBLUP) */

} Xstruct;

typedef struct {

	double	deltaShape;		/* for BL and EBL */
	double	deltaRate;		/* for BL and EBL */
	double	etaShape;		/* for EBL */
	double	etaRate;		/* for EBL */
	double	v;				/* for wBSR, BayesC, SSVS, MIX, and GBLUP */
	double	S2;				/* for wBSR, BayesC, SSVS, MIX, and GBLUP */
	double	Pi;				/* for wBSR, SSVS, MIX, and BayesC */
	double	c;				/* for SSVS and MIX */

} Hstruct;

/*----generate a random integer -------------------------------------------------------------------------------------------------*/
int RandomInteger(int low, int high)
{
	int k;
	double value;

	value = genrand_real3();
	k = (int)(value * (double)(high - low + 1));
	return (low + k);
}

/*----generate a real number from Uniform (0,1)------------------------------------------------------------------------------------*/
double rnd(void)
{
	return(genrand_real3());
}

/*----generate real numbers from N (0,1) ( from STRUCTURE )------------------------------------------------------------------------*/
double snorm(void)
{
	static double a[32] = {
		0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
		0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
		0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
		1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
		1.862732,2.153875
	};
	static double d[31] = {
		0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
		0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
		0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
		0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
	};
	static double t[31] = {
		7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
		1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
		2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
		4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
		9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
	};
	static double h[31] = {
		3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
		4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
		4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
		5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
		8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
	};
	static long inr;
	static double snorm, u, ss, ustar, anr, w, y, tt;
	u = rnd();   /* was ranf--JKP*/
	ss = 0.0;
	if (u > 0.5) ss = 1.0;
	u += (u - ss);
	u = 32.0*u;
	inr = (long)(u);
	if (inr == 32) inr = 31;
	if (inr == 0) goto S100;
	/*
	START CENTER
	*/
	ustar = u - (double)inr;
	anr = *(a + inr - 1);
S40:
	if (ustar <= *(t + inr - 1)) goto S60;
	w = (ustar - *(t + inr - 1))**(h + inr - 1);
S50:
	/*
	EXIT   (BOTH CASES)
	*/
	y = anr + w;
	snorm = y;
	if (ss == 1.0) snorm = -y;
	return snorm;
S60:
	/*
	CENTER CONTINUED
	*/
	u = rnd();                /*was ranf--JKP*/
	w = u * (*(a + inr) - anr);
	tt = (0.5*w + anr)*w;
	goto S80;
S70:
	tt = u;
	ustar = rnd();                /*was ranf--JKP*/
S80:
	if (ustar > tt) goto S50;
	u = rnd();               /*was ranf--JKP*/
	if (ustar >= u) goto S70;
	ustar = rnd();               /*was ranf--JKP*/
	goto S40;
S100:
	/*
	START TAIL
	*/
	inr = 6;
	anr = *(a + 31);
	goto S120;
S110:
	anr += *(d + inr - 1);
	inr += 1;
S120:
	u += u;
	if (u < 1.0) goto S110;
	u -= 1.0;
S140:
	w = u * *(d + inr - 1);
	tt = (0.5*w + anr)*w;
	goto S160;
S150:
	tt = u;
S160:
	ustar = rnd();               /*was ranf--JKP*/
	if (ustar > tt) goto S50;
	u = rnd();               /*was ranf--JKP*/
	if (ustar >= u) goto S150;
	u = rnd();               /*was ranf--JKP*/
	goto S140;
}

/*----generate real number from N(mu, sd) (from STRUCTURE)-------------------------------------------------------------------------*/
double RNormal(double mu, double sd)
{
	return (mu + sd * snorm());
}

/*----calculate log of the gamma function (from STRUCTURE)-------------------------------------------------------------------------*/
double mylgamma(double zpo)
{
	/* LGAMMA function

	double_value = lgamma(<double_value > 0.>)

	returns the natural log of the gamma function

	Uses Lanczos-type approximation to ln(gamma) for zpo > 0.
	Reference:
	Lanczos, C. 'A precision approximation of the gamma
	function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.

	Original was in FORTRAN
	Accuracy: About 14 significant digits except for small regions
	in the vicinity of 1 and 2.
	Programmer: Alan Miller
	CSIRO Division of Mathematics & Statistics
	Latest revision - 17 April 1988

	Translated and modified into C by Peter Beerli 1997
	Tested against Mathematica's Log[Gamma[x]]
	*/

	double alg[9] = { 0.9999999999995183, 676.5203681218835,
		-1259.139216722289, 771.3234287757674, -176.6150291498386,
		12.50734324009056, -0.1385710331296526, 9.934937113930748e-6,
		1.659470187408462e-7 };
	double lnsqrt2pi = 0.9189385332046727;
	double result;
	long jlg;
	double tmp;

	result = 0.;
	tmp = zpo + 7.;
	for (jlg = 9; jlg >= 2; --jlg)
	{
		result += alg[jlg - 1] / tmp;
		tmp -= 1.;
	}
	result += alg[0];
	result = log(result) + lnsqrt2pi - (zpo + 6.5) + (zpo - 0.5) * log(zpo + 6.5);
	return result;
}

/*----Shuffle an array------------------------------------------------------------------------------------------------------------*/
void Shuffle(int *asl, int arraylength, int timesofshuffle)
{
	int ts, content, target1, target2;

	for (ts = 0; ts<timesofshuffle; ts++)
	{
		target1 = RandomInteger(0, arraylength - 1);
		content = asl[target1];
		target2 = RandomInteger(0, arraylength - 1);
		asl[target1] = asl[target2];
		asl[target2] = content;
	}
}

/*----digamma function-----------------------------------------------------------------------------------------------------------------------*/
/* taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.C. But now cannot be found */
double Digamma(const double x)
{
	double c = 8.5;
	double d1 = -0.5772156649;
	double r;
	double s = 0.00001;
	double s3 = 0.08333333333;
	double s4 = 0.0083333333333;
	double s5 = 0.003968253968;
	double value;
	double y;

	//  if ( x <= 0.0 ){printf ( "Error: negative value in Digamma\n" );}
	y = x;
	value = 0.0;
	if (y <= s)
	{
		value = d1 - 1.0 / y;
		return (value);
	}
	while (y < c)
	{
		value = value - 1.0 / y;
		y = y + 1.0;
	}
	r = 1.0 / y;
	value = value + log(y) - 0.5 * r;
	r = r * r;
	value = value - r * (s3 - r * (s4 - r * s5));

	return (value);
}

/*----Inner product---------------------------------------------------------------------------------------------------------------------------*/
/* C = t(A)%*%B
A: m x n
B: m x k
C: n x k
The order of elements in each matrix is [row1column1, row2column1, ..., row1column2, row2column2, ...]
*/
void Innerproduct_tAB(double *A, double *B, int n, int m, int k, double *C)
{
	int nn, mm, kk, pos1;

	for (nn = 0; nn < n; nn++)
	{
		for (kk = 0; kk < k; kk++)
		{
			pos1 = kk * n + nn;
			C[pos1] = 0.0;
			for (mm = 0; mm < m; mm++)
			{
				C[pos1] += A[nn*m + mm] * B[kk*m + mm];
			}
		}
	}
}

/* C = t(A)%*%D%*%A where D is diagonal
A: m x n
D: m x m
C: n x n
The order of elements in each matrix is [row1column1, row2column1, ..., row1column2, row2column2, ...]
*/
void Innerproduct_tADA(double *A, double *D, int n, int m, double *C)
{
	int nn1, nn2, mm, pos1, pos2;

	for (nn1 = 0; nn1 < n; nn1++)
	{
		for (nn2 = nn1; nn2 < n; nn2++)
		{
			pos1 = nn2 * n + nn1;
			pos2 = nn1 * n + nn2;
			C[pos1] = 0.0;
			for (mm = 0; mm < m; mm++)
			{
				C[pos1] += A[nn1 * m + mm] * A[nn2 * m + mm] * D[mm];
			}
			C[pos2] = C[pos1];
		}
	}
}

//test
/*__declspec(dllexport) void Innerproduct_tADA_dll(double *A, double *D, int *n, int *m, double *C)
{
	int nn1, nn2, mm, pos1, pos2;

	for (nn1 = 0; nn1 < n[0]; nn1++)
	{
		for (nn2 = nn1; nn2 < n[0]; nn2++)
		{
			pos1 = nn2 * n[0] + nn1;
			pos2 = nn1 * n[0] + nn2;
			C[pos1] = 0.0;
			for (mm = 0; mm < m[0]; mm++)
			{
				C[pos1] += A[nn1 * m[0] + mm] * A[nn2 * m[0] + mm] * D[mm];
			}
			C[pos2] = C[pos1];
		}
	}
}*/

/*---Allocate and free memory--------------------------------------------------------------------------------------------------------------------------------------*/
void AllocateMemory(int Priortype, int Methodcode, int Nx, int P, Xstruct *X, double *Genotypes)
{
	switch (Priortype) {
	case 1: /* BL and EBL */
		X[0].covariates = (double*)calloc(Nx * P, sizeof(double));
		memcpy(X[0].covariates, Genotypes, sizeof(double) * Nx * P);
		X[0].x2 = (double*)calloc(P, sizeof(double));
		X[0].expEffect = (double*)calloc(P, sizeof(double));
		X[0].exp2Effect = (double*)calloc(P, sizeof(double));
		X[0].varEffect = (double*)calloc(P, sizeof(double));
		X[0].expTau2 = (double*)calloc(P, sizeof(double));
		X[0].expInTau2 = (double*)calloc(P, sizeof(double));
		X[0].expEta2 = (double*)calloc(P, sizeof(double));
		if (Methodcode == 2)
		{
			X[0].a3 = (double*)calloc(P, sizeof(double));
			X[0].b3 = (double*)calloc(P, sizeof(double));
		}
		break;
	case 2: /* wBSR, BayesC, SSVS, MIX, and BayesB */
		X[0].covariates = (double*)calloc(Nx * P, sizeof(double));
		memcpy(X[0].covariates, Genotypes, sizeof(double) * Nx * P);
		X[0].x2 = (double*)calloc(P, sizeof(double));
		X[0].expEffect = (double*)calloc(P, sizeof(double));
		X[0].exp2Effect = (double*)calloc(P, sizeof(double));
		X[0].varEffect = (double*)calloc(P, sizeof(double));
		X[0].expGamma = (double*)calloc(P, sizeof(double));
		X[0].exp2Gamma = (double*)calloc(P, sizeof(double));
		if (Methodcode == 3 || Methodcode == 7)/* wBSR and BayesB */
		{
			X[0].expSigma2 = (double*)calloc(P, sizeof(double));
			X[0].S2 = (double*)calloc(P, sizeof(double));
		}
		else
		{
			X[0].expSigma2 = (double*)calloc(2, sizeof(double));
			X[0].S2 = (double*)calloc(2, sizeof(double));
		}
		break;
	case 3: /* GBLUP */
		X[0].evalue = (double*)calloc(Nx, sizeof(double));
		X[0].evector = (double*)calloc(Nx * Nx, sizeof(double));
		X[0].tevector = (double*)calloc(Nx * Nx, sizeof(double));
		memcpy(X[0].evalue, Genotypes, sizeof(double) * Nx);
		memcpy(X[0].evector, Genotypes + Nx, sizeof(double) * Nx * Nx);
		memcpy(X[0].tevector, Genotypes + Nx + Nx * Nx, sizeof(double) * Nx * Nx);
		X[0].expEffect = (double*)calloc(Nx, sizeof(double));
		X[0].exp2Effect = (double*)calloc(1, sizeof(double));
		X[0].varEffect = (double*)calloc(Nx * Nx, sizeof(double));
		X[0].expSigma2 = (double*)calloc(1, sizeof(double));
		X[0].S2 = (double*)calloc(1, sizeof(double));
		break;
	case 4: /* FIXED */
		X[0].covariates = (double*)calloc(Nx * P, sizeof(double));
		memcpy(X[0].covariates, Genotypes, sizeof(double) * Nx * P);
		X[0].x2 = (double*)calloc(P, sizeof(double));
		X[0].expEffect = (double*)calloc(P, sizeof(double));
		X[0].exp2Effect = (double*)calloc(P, sizeof(double));
		X[0].varEffect = (double*)calloc(P, sizeof(double));
	}
}

void FreeMemory(int Priortype, int Methodcode, Xstruct *X)
{
	switch (Priortype) {
	case 1: /* BL and EBL */
		free(X[0].covariates);
		free(X[0].x2);
		free(X[0].expEffect);
		free(X[0].exp2Effect);
		free(X[0].varEffect);
		free(X[0].expTau2);
		free(X[0].expInTau2);
		free(X[0].expEta2);
		if (Methodcode == 2)
		{
			free(X[0].a3);
			free(X[0].b3);
		}
		break;
	case 2: /* wBSR, BayesC, SSVS, MIX, and BayesB */
		free(X[0].covariates);
		free(X[0].x2);
		free(X[0].expEffect);
		free(X[0].exp2Effect);
		free(X[0].varEffect);
		free(X[0].expGamma);
		free(X[0].exp2Gamma);
		free(X[0].expSigma2);
		free(X[0].S2);
		break;
	case 3: /* GBLUP */
		free(X[0].evalue);
		free(X[0].evector);
		free(X[0].tevector);
		free(X[0].expEffect);
		free(X[0].exp2Effect);
		free(X[0].varEffect);
		free(X[0].expSigma2);
		free(X[0].S2);
		break;
	case 4: /* FIXED */
		free(X[0].covariates);
		free(X[0].x2);
		free(X[0].expEffect);
		free(X[0].exp2Effect);
		free(X[0].varEffect);
		break;
	}
}


/*---Assign hyperparameter values------------------------------------------------------------------------------------------------------------------------------*/
void AssignHyperparameters(int Methodcode, Hstruct *H, double *Hyperparameters)
{
	switch (Methodcode) {
	case 1:  /* BL */
		H[0].deltaShape = Hyperparameters[0];
		H[0].deltaRate = Hyperparameters[1];
		break;
	case 2: /* EBL */
		H[0].deltaShape = Hyperparameters[0];
		H[0].deltaRate = Hyperparameters[1];
		H[0].etaShape = Hyperparameters[2];
		H[0].etaRate = Hyperparameters[3];
		break;
	case 3: /* wBSR */
		H[0].v = Hyperparameters[0];
		H[0].S2 = Hyperparameters[1];
		H[0].Pi = Hyperparameters[2];
		break;
	case 4: /* BayesC */
		H[0].v = Hyperparameters[0];
		H[0].S2 = Hyperparameters[1];
		H[0].Pi = Hyperparameters[2];
		break;
	case 5: /* SSVS */
		H[0].c = Hyperparameters[0];
		H[0].v = Hyperparameters[1];
		H[0].S2 = Hyperparameters[2];
		H[0].Pi = Hyperparameters[3];
		break;
	case 6: /* MIX */
		H[0].c = Hyperparameters[0];
		H[0].v = Hyperparameters[1];
		H[0].S2 = Hyperparameters[2];
		H[0].Pi = Hyperparameters[3];
		break;
	case 7: /* BayesB */
		H[0].v = Hyperparameters[0];
		H[0].S2 = Hyperparameters[1];
		H[0].Pi = Hyperparameters[2];
		break;
	case 8: /* GBLUP */
		H[0].v = Hyperparameters[0];
		H[0].S2 = Hyperparameters[1];
		break;
	}
}


/*--- Copy the results of GenomeWideRegression to the arguments of vigor -------------------------------------------------------------------------*/
void CopyResults(int Methodcode, int P, int Nx, Xstruct *X, double *Expectation, double *Uncertainty, double *Variance, double *Gamma)
{
	int locus, level;
	double temp;

	switch (Methodcode) {
	case 1: /* BL */
		memcpy(Expectation, X[0].expEffect, sizeof(double) * P);
		for (locus = 0; locus<P; locus++)
			Uncertainty[locus] = sqrt(X[0].varEffect[locus]);
		break;
	case 2: /* EBL */
		memcpy(Expectation, X[0].expEffect, sizeof(double) * P);
		for (locus = 0; locus<P; locus++)
			Uncertainty[locus] = sqrt(X[0].varEffect[locus]);
		break;
	case 3: /* wBSR */
		memcpy(Gamma, X[0].expGamma, sizeof(double) * P);
		memcpy(Variance, X[0].expSigma2, sizeof(double) * P);
		for (locus = 0; locus<P; locus++)
		{
			Expectation[locus] = X[0].expEffect[locus] * X[0].expGamma[locus];
			temp = X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]) * X[0].varEffect[locus]
				+ X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]) * pow(X[0].expEffect[locus], 2.0)
				+ pow(X[0].expGamma[locus], 2.0) * X[0].varEffect[locus];
			Uncertainty[locus] = sqrt(temp);
		}
		break;
	case 4: /* BayesC */
		memcpy(Gamma, X[0].expGamma, sizeof(double) * P);
		Variance[0] = X[0].expSigma2[0];
		for (locus = 0; locus<P; locus++)
		{
			Expectation[locus] = X[0].expEffect[locus] * X[0].expGamma[locus];
			temp = pow(X[0].expEffect[locus], 2.0) * X[0].expGamma[locus]  * (1.0 - X[0].expGamma[locus]) + X[0].varEffect[locus] * X[0].expGamma[locus];
			Uncertainty[locus] = sqrt(temp);
		}
		break;
	case 5: /* SSVS */
		memcpy(Expectation, X[0].expEffect, sizeof(double) * P);
		memcpy(Gamma, X[0].expGamma, sizeof(double) * P);
		Variance[0] = X[0].expSigma2[0];
		for (locus = 0; locus<P; locus++)
			Uncertainty[locus] = sqrt(X[0].varEffect[locus]);
		break;
	case 6: /* MIX */
		memcpy(Expectation, X[0].expEffect, sizeof(double) * P);
		memcpy(Gamma, X[0].expGamma, sizeof(double) * P);
		Variance[0] = X[0].expSigma2[0];
		Variance[1] = X[0].expSigma2[1];
		for (locus = 0; locus<P; locus++)
			Uncertainty[locus] = sqrt(X[0].varEffect[locus]);
		break;
	case 7: /* BayesB */
		memcpy(Gamma, X[0].expGamma, sizeof(double) * P);
		memcpy(Variance, X[0].expSigma2, sizeof(double) * P);
		for (locus = 0; locus<P; locus++)
		{
			Expectation[locus] = X[0].expEffect[locus] * X[0].expGamma[locus];
			temp = pow(X[0].expEffect[locus], 2.0) * X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]) + X[0].varEffect[locus] * X[0].expGamma[locus];
			Uncertainty[locus] = sqrt(temp);
		}
		break;
	case 8: /* GBLUP */
		memcpy(Expectation, X[0].expEffect, sizeof(double) * Nx);
		for (level = 0; level<Nx; level++)
			Uncertainty[level] = sqrt(X[0].varEffect[level * Nx + level]);
		Variance[0] = X[0].expSigma2[0];
		break;
	case 9: /* FIXED */
		memcpy(Expectation, X[0].expEffect, sizeof(double) * P);
		for (locus = 0; locus<P; locus++)
			Uncertainty[locus] = sqrt(X[0].varEffect[locus]);
		break;
	}
}


/*---Lower bound of marginal log probability of data-----------------------------------------------------------------------------------------------------------*/
double LowerBound(int Methodcode, int N, int P, Xstruct *X, double a2, double b2, Hstruct *H, double Tau0)
{
	int		locus, level;
	double	LB = 0.0, temp, temp2, vj, vS2, Pi2, v2, vN;
	double	Gammalowerbound = 1e-100, Gammaupperbound;

	Gammaupperbound = 1.0 - Gammalowerbound;
	switch (Methodcode)
	{
	case 1: /* BL */
		for (locus = 0; locus<P; locus++)
		{
			LB -= 0.5 * log(X[0].expDelta2);
			LB += 0.5 * log(X[0].varEffect[locus]);
			LB -= 0.5 * (X[0].expDelta2 * X[0].expInTau2[locus]);
		}
		LB -= (a2 * (log(b2) - 1.0) - mylgamma(a2));
		LB -= H[0].deltaRate * X[0].expDelta2;
		LB += ((double)(P) * 0.5 * log2pi);
		if (H[0].deltaShape>0.0) LB -= mylgamma(H[0].deltaShape);			/* When informative prior is used */
		if (H[0].deltaRate>0.0) LB += H[0].deltaShape * log(H[0].deltaRate);
		LB += ((double)(2 * P) * 0.5 - (double)P * log(2.0));
		break;
	case 2: /* EBL */
		for (locus = 0; locus<P; locus++)
		{
			LB -= 0.5*log(X[0].expDelta2 * X[0].expEta2[locus]);
			LB += 0.5*log(X[0].varEffect[locus]);
			LB -= (X[0].a3[locus] * (log(X[0].b3[locus]) - 1.0) - mylgamma(X[0].a3[locus]));
			LB -= 0.5*(X[0].expDelta2 * X[0].expEta2[locus] * X[0].expInTau2[locus]);
			LB -= H[0].etaRate * X[0].expEta2[locus];
		}
		LB -= (a2 * (log(b2) - 1.0) - mylgamma(a2));
		LB -= H[0].deltaRate * X[0].expDelta2;
		LB += ((double)(P) * 0.5 * log2pi);
		if (H[0].deltaShape>0.0)	LB -= mylgamma(H[0].deltaShape);				/* When informative prior is used */
		if (H[0].deltaRate>0.0)		LB += H[0].deltaShape * log(H[0].deltaRate);
		if (H[0].etaShape>0.0)		LB -= (double)P * mylgamma(H[0].etaShape);		/* When informative prior is used */
		if (H[0].etaRate>0.0)		LB += (double)P * (H[0].etaShape * log(H[0].etaRate));
		LB += ((double)(2 * P) * 0.5 - (double)P * log(2.0));
		break;
	case 3:	/* wBSR */
		vj = 1.0 + H[0].v;
		vS2 = H[0].v * H[0].S2;
		Pi2 = (1.0 - H[0].Pi);
		if (H[0].Pi<1.0) /* BayesB */
		{
			for (locus = 0; locus<P; locus++)
			{
				LB -= 0.5 * log(vj * X[0].S2[locus] * 0.5) * vj;
				LB += 0.5 * log(X[0].varEffect[locus]);
				if (X[0].expGamma[locus] > Gammalowerbound && X[0].expGamma[locus] < Gammaupperbound)
				{
					LB += X[0].expGamma[locus] * log(H[0].Pi / X[0].expGamma[locus]);
					LB += (1.0 - X[0].expGamma[locus])*log(Pi2 / (1.0 - X[0].expGamma[locus]));
				}
				else
				{
					if (X[0].expGamma[locus] > Gammalowerbound)
					{
						LB += log(H[0].Pi); /* expGamma==1.0 */
					}
					else
					{
						LB += log(Pi2);/* expGamma==0.0 */
					}
				}
			}
		}
		else
		{	/* BayesA */
			for (locus = 0; locus<P; locus++)
			{
				LB -= 0.5 * log(vj * X[0].S2[locus] * 0.5) * vj;
				LB += 0.5 * log(X[0].varEffect[locus]);
			}
		}
		LB += (double)P * (0.5 * H[0].v * log(0.5 * vS2) - mylgamma(0.5 * H[0].v) + mylgamma(0.5 * vj));
		LB += (double)P * 0.5;
		break;
	case 4:	/* BayesC */
		if (H[0].Pi<1.0) /* BayesC */
		{
			vS2 = H[0].v * H[0].S2;
			Pi2 = (1.0 - H[0].Pi);
			for (locus = 0, temp = 0.0; locus<P; locus++)
			{
				LB += 0.5 * X[0].expGamma[locus] * log(X[0].varEffect[locus]);
				if (X[0].expGamma[locus] > Gammalowerbound && X[0].expGamma[locus] < Gammaupperbound)
				{
					LB += X[0].expGamma[locus] * log(H[0].Pi / X[0].expGamma[locus]);
					LB += (1.0 - X[0].expGamma[locus])*log(Pi2 / (1.0 - X[0].expGamma[locus]));
				}
				else
				{
					if (X[0].expGamma[locus] > Gammalowerbound)
					{
						LB += log(H[0].Pi); /* expGamma==1.0 */
					}
					else
					{
						LB += log(Pi2);/* expGamma==0.0 */
					}
				}
				temp += X[0].expGamma[locus];
			}
			vj = H[0].v + temp;
			LB -= (0.5 * log(vj * X[0].S2[0] * 0.5) * vj - mylgamma(0.5 * vj));
			LB += (0.5 * H[0].v * log(0.5 * vS2) - mylgamma(0.5 * H[0].v));
			LB += temp * 0.5;
		}
		else
		{	/* BRR */
			vj = (double)P + H[0].v;
			vS2 = H[0].v * H[0].S2;
			for (locus = 0; locus<P; locus++) { LB += 0.5 * log(X[0].varEffect[locus]); }
			LB -= (0.5 * log(vj * X[0].S2[0] * 0.5) * vj - mylgamma(0.5 * vj));
			LB += (0.5 * H[0].v * log(0.5 * vS2) - mylgamma(0.5 * H[0].v));
			LB += (double)P * 0.5;
		}
		break;
	case 5:	/* SSVS */
		vj = (double)P + H[0].v;
		vS2 = H[0].v * H[0].S2;
		Pi2 = (1.0 - H[0].Pi);
		for (locus = 0, temp = 0.0; locus<P; locus++)
		{
			LB += 0.5 * log(X[0].varEffect[locus]);
			if (X[0].expGamma[locus] > Gammalowerbound && X[0].expGamma[locus] < Gammaupperbound)
			{
				LB += X[0].expGamma[locus] * log(H[0].Pi / X[0].expGamma[locus]);
				LB += (1.0 - X[0].expGamma[locus])*log(Pi2 / (1.0 - X[0].expGamma[locus]));
			}
			else
			{
				if (X[0].expGamma[locus] > Gammalowerbound)
				{
					LB += log(H[0].Pi);
				}
				else
				{
					LB += log(Pi2);
				}
			}
			temp += X[0].expGamma[locus];
		}
		LB -= 0.5 * ((double)P - temp) * log(H[0].c);
		LB -= (0.5 * log(vj * X[0].S2[0] * 0.5) * vj - mylgamma(0.5 * vj));
		LB += (0.5 * H[0].v * log(0.5 * vS2) - mylgamma(0.5 * H[0].v));
		LB += (double)P * 0.5;
		break;
	case 6: /* MIX */
		vS2 = H[0].v * H[0].S2;
		Pi2 = (1.0 - H[0].Pi);
		for (locus = 0, vj = H[0].v, v2 = H[0].v + (double)P; locus<P; locus++)
		{
			LB += 0.5 * log(X[0].varEffect[locus]);
			if (X[0].expGamma[locus] > Gammalowerbound && X[0].expGamma[locus] < Gammaupperbound)
			{
				LB += X[0].expGamma[locus] * log(H[0].Pi / X[0].expGamma[locus]);
				LB += (1.0 - X[0].expGamma[locus])*log(Pi2 / (1.0 - X[0].expGamma[locus]));
			}
			else
			{
				if (X[0].expGamma[locus] > Gammalowerbound)
				{
					LB += log(H[0].Pi); /* expGamma==1.0 */
				}
				else
				{
					LB += log(Pi2);/* expGamma==0.0 */
				}
			}
			vj += X[0].expGamma[locus];
			v2 -= X[0].expGamma[locus];
		}
		LB -= (0.5 * log(vj * X[0].S2[0] * 0.5) * vj - mylgamma(0.5 * vj));
		LB -= (0.5 * log(v2 * X[0].S2[1] * 0.5) * v2 - mylgamma(0.5 * v2));
		LB += (H[0].v * log(0.5 * vS2) - 2.0 * mylgamma(0.5 * H[0].v) + 0.5 * H[0].v * log(H[0].c));
		LB += (double)P * 0.5;
		break;
	case 7:	/* BayesB */
		if (H[0].Pi<1.0) /* BayesB */
		{
			vS2 = H[0].v * H[0].S2;
			Pi2 = (1.0 - H[0].Pi);
			for (locus = 0, temp = 0.0; locus<P; locus++)
			{
				LB += 0.5 * X[0].expGamma[locus] * log(X[0].varEffect[locus]);
				if (X[0].expGamma[locus] > Gammalowerbound && X[0].expGamma[locus] < Gammaupperbound)
				{
					LB += X[0].expGamma[locus] * log(H[0].Pi / X[0].expGamma[locus]);
					LB += (1.0 - X[0].expGamma[locus])*log(Pi2 / (1.0 - X[0].expGamma[locus]));
				}
				else
				{
					if (X[0].expGamma[locus] > Gammalowerbound)
					{
						LB += log(H[0].Pi); /* expGamma==1.0 */
					}
					else
					{
						LB += log(Pi2);/* expGamma==0.0 */
					}
				}
				vj = H[0].v + X[0].expGamma[locus];
				LB -= (0.5 * log(vj * X[0].S2[locus] * 0.5) * vj - mylgamma(0.5 * vj));
				temp += X[0].expGamma[locus];
			}
			LB += (double)P * (0.5 * H[0].v * log(0.5 * vS2) - mylgamma(0.5 * H[0].v));
			LB += (double)temp * 0.5;
		}
		else
		{	/* BayesA */
			vS2 = H[0].v * H[0].S2;
			vj = 1.0 + H[0].v;
			for (locus = 0; locus<P; locus++) {
				LB += 0.5 * log(X[0].varEffect[locus]);
				LB -= (0.5 * log(vj * X[0].S2[locus] * 0.5) * vj);
			}
			LB += (double)P * (0.5 * H[0].v * log(0.5 * vS2) - mylgamma(0.5 * H[0].v) + mylgamma(0.5 * vj));
			LB += (double)P * 0.5;
		}
		break;
	case 8: /* GBLUP */
		vN = N + H[0].v;
		temp2 = log(X[0].S2[0]);
		LB += -0.5 * X[0].expDelta2 - 0.5 * X[0].exp2Effect[0] / X[0].S2[0];
		LB += 0.5 * H[0].v * log(0.5 * H[0].v * H[0].S2) - mylgamma(0.5 * H[0].v) - 0.5 * H[0].v * H[0].S2 / X[0].S2[0];
		for (level = 0, temp = 0.0; level < N; level++)
		{
			temp += temp2 - log (X[0].evalue[level] + X[0].S2[0] * Tau0);
		}
		LB += 0.5 * (temp + N);
		LB += (-0.5 * vN * log(0.5 * vN * X[0].S2[0])) + mylgamma(0.5 * vN) + 0.5 * vN;
		break;
	case 9: /* FIXED */
		LB += ((double) P * 0.5 * log2pi);
		for (locus = 0; locus<P; locus++) { LB += 0.5 * log(X[0].varEffect[locus]); }
		LB += (double)P * 0.5;
		break;
	}
	return(LB);
}

void	Initialize_BL(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, int RandomIni)
{
	/* For repeat statement */
	int		record, locus;

	/* Temporary objects */
	double	temp;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	temp = 1.0;
	X[0].expDelta2 = 1.0;
	if (RandomIni)
	{
		for (locus = 0; locus<P; locus++)
		{
			X[0].expTau2[locus] = (double)P / temp;
			X[0].expInTau2[locus] = temp / (double)P;
			X[0].expEta2[locus] = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[0].expEffect[locus] = RNormal(0.0, sqrt(1.0 / X[0].expTau2[locus]));
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
		}
	}
	else
	{
		for (locus = 0; locus<P; locus++)
		{
			X[0].expTau2[locus] = (double)P / temp;
			X[0].expInTau2[locus] = temp / (double)P;
			X[0].expEta2[locus] = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[0].expEffect[locus] = 0.0;
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
		}
	}

	for (locus = 0; locus<P; locus++)
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus];
}


void	UpdateB_BL(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *sumVarB, double *sumTau2B2, double *a2, double *b2,
	int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2, temp3;

	/* New values */
	double	prop, prop2;

	/* for update of Delta2 */
	double	sumEta2InTau2;

	temp3 = Tau0[0];
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];

		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target]);

		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0] + X[0].expTau2[target] * temp3);

		prop = temp * temp2;
		prop2 = prop * prop + temp2;

		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * (X[0].expEffect[target] - prop));
		sumVarB[0] += (X[0].x2[target] * temp2);

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);

		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}

	/* update of Tau2 or Sigma2*/
	temp = Tau0[0];
	for (locus = 0, sumEta2InTau2 = 0.0; locus<P; locus++)
	{
		target = Order[locus];

		prop = sqrt(X[0].expDelta2 * X[0].expEta2[target] / (X[0].exp2Effect[target] * temp));

		Check1[0] += pow((prop - X[0].expTau2[target]), 2.0);
		Check2[0] += pow(prop, 2.0);

		sumTau2B2[0] += X[0].exp2Effect[target] * prop;
		X[0].expTau2[target] = prop;

		prop = 1.0 / X[0].expTau2[target] + 1.0 / (X[0].expDelta2 * X[0].expEta2[target]);
		sumEta2InTau2 += prop * X[0].expEta2[target];

		X[0].expInTau2[target] = prop;
	}

	/* update of Delta2 */
	a2[0] = (double)P + H[0].deltaShape;
	b2[0] = 0.5 * sumEta2InTau2 + H[0].deltaRate;
	prop = a2[0] / b2[0];
	Check1[0] += pow((prop - X[0].expDelta2), 2.0);
	Check2[0] += pow(prop, 2.0);
	X[0].expDelta2 = prop;
}


void	Initialize_EBL(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, int RandomIni)
{
	/* For repeat statement */
	int		record, locus;

	/* Temporary objects */
	double	temp;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	temp = 1.0;
	X[0].expDelta2 = 1.0;
	if (RandomIni)
	{
		for (locus = 0; locus<P; locus++)
		{
			X[0].expTau2[locus] = (double)P / temp;
			X[0].expInTau2[locus] = temp / (double)P;
			X[0].expEta2[locus] = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[0].expEffect[locus] = RNormal(0.0, sqrt(1.0 / X[0].expTau2[locus]));
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			X[0].a3[locus] = 0.0;
			X[0].b3[locus] = 0.0;
		}
	}
	else
	{
		for (locus = 0; locus<P; locus++)
		{
			X[0].expTau2[locus] = (double)P / temp;
			X[0].expInTau2[locus] = temp / (double)P;
			X[0].expEta2[locus] = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[0].expEffect[locus] = 0.0;
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			X[0].a3[locus] = 0.0;
			X[0].b3[locus] = 0.0;
		}
	}

	for (locus = 0; locus<P; locus++)
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus];
}


void	UpdateB_EBL(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *sumVarB, double *sumTau2B2, double *a2, double *b2,
	int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2, temp3;

	/* New values */
	double	prop, prop2;

	/* for update of Delta2 */
	double	sumEta2InTau2;

	temp3 = Tau0[0];
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];

		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target]);

		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0] + X[0].expTau2[target] * temp3);

		prop = temp * temp2;
		prop2 = prop * prop + temp2;

		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * (X[0].expEffect[target] - prop));
		sumVarB[0] += (X[0].x2[target] * temp2);

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);

		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}

	/* update of Tau2 or Sigma2*/
	temp = Tau0[0];
	for (locus = 0, sumEta2InTau2 = 0.0; locus<P; locus++)
	{
		target = Order[locus];
		prop = sqrt(X[0].expDelta2 * X[0].expEta2[target] / (X[0].exp2Effect[target] * temp));

		Check1[0] += pow((prop - X[0].expTau2[target]), 2.0);
		Check2[0] += pow(prop, 2.0);

		sumTau2B2[0] += X[0].exp2Effect[target] * prop;

		X[0].expTau2[target] = prop;

		prop = 1.0 / X[0].expTau2[target] + 1.0 / (X[0].expDelta2 * X[0].expEta2[target]);
		sumEta2InTau2 += prop * X[0].expEta2[target];

		X[0].expInTau2[target] = prop;
	}

	/* update of Delta2 */
	a2[0] = (double)P + H[0].deltaShape;
	b2[0] = 0.5 * sumEta2InTau2 + H[0].deltaRate;
	prop = a2[0] / b2[0];

	Check1[0] += pow((prop - X[0].expDelta2), 2.0);
	Check2[0] += pow(prop, 2.0);
	X[0].expDelta2 = prop;

	/* update of Eta2 */
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];
		X[0].a3[target] = 1.0 + H[0].etaShape;
		X[0].b3[target] = 0.5 * X[0].expDelta2 * X[0].expInTau2[target] + H[0].etaRate;
		prop = X[0].a3[target] / X[0].b3[target];

		Check1[0] += pow((prop - X[0].expEta2[target]), 2.0);
		Check2[0] += pow(prop, 2.0);
		X[0].expEta2[target] = prop;
	}
}


void	Initialize_wBSR(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, Hstruct *H, int RandomIni)
{
	/* For repeat statement */
	int		record, locus;

	/* Temporary objects */
	double	temp;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	temp = 1.0;
	if (H[0].Pi<1.0)
	{/* BayesB */
		if (RandomIni)
		{
			for (locus = 0; locus<P; locus++)
			{
				X[0].expSigma2[locus] = temp / ((double)P * H[0].Pi);
				X[0].S2[locus] = (double)P * H[0].Pi / temp;
				X[0].expGamma[locus] = 0.5;
				X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
				X[0].expEffect[locus] = RNormal(0.0, sqrt(X[0].expSigma2[locus]));
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			}
		}
		else {
			for (locus = 0; locus<P; locus++)
			{
				X[0].expEffect[locus] = 0.0;
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
				X[0].expSigma2[locus] = temp / ((double)P * H[0].Pi);
				X[0].S2[locus] = (double)P * H[0].Pi / temp;
				X[0].expGamma[locus] = 0.5;
				X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
			}
		}
	}
	else
	{/* BayesA */
		if (RandomIni)
		{
			for (locus = 0; locus<P; locus++)
			{
				X[0].expSigma2[locus] = temp / (double)P;
				X[0].S2[locus] = (double)P / temp;
				X[0].expGamma[locus] = 1.0;
				X[0].exp2Gamma[locus] = 1.0;
				X[0].expEffect[locus] = RNormal(0.0, sqrt(X[0].expSigma2[locus]));
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			}
		}
		else {
			for (locus = 0; locus<P; locus++)
			{
				X[0].expEffect[locus] = 0.0;
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
				X[0].expSigma2[locus] = temp / (double)P;
				X[0].S2[locus] = (double)P / temp;
				X[0].expGamma[locus] = 1.0;
				X[0].exp2Gamma[locus] = 1.0;
			}
		}
	}

	for (locus = 0; locus<P; locus++)
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus] * X[0].expGamma[locus];
}


void	UpdateB_wBSR(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *sumVarB, int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2;

	/* New values */
	double	prop, prop2;

	/* Product of hyperparameters */
	double	vS2;

	/* For variable selection */
	double	logPi, log1minusPi;

	/* Probability of gamma */
	double	ProbInclude, ProbExclude;

	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi);
	log1minusPi = log(1.0 - H[0].Pi);

	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];
		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target] * X[0].expGamma[target]);
		temp *= X[0].expGamma[target];
		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0] * X[0].exp2Gamma[target] + 1.0 / X[0].S2[target]);

		/* prop: E[B], prop2: E[B^2], temp2: V[B] */
		prop = temp * temp2;
		prop2 = pow(prop, 2.0) + temp2;

		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * X[0].expGamma[target] * (X[0].expEffect[target] - prop));
		if ((int)H[0].Pi == 1) /* when BayesA */
			sumVarB[0] += (X[0].x2[target] * temp2);

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);
		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}

	/* update of Tau2 or Sigma2*/
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];
		temp = X[0].exp2Effect[target] + vS2;

		X[0].expSigma2[target] = temp / (H[0].v - 1.0);
		prop = temp / (H[0].v + 1.0);

		Check1[0] += pow((prop - X[0].S2[target]), 2.0);
		Check2[0] += pow(prop, 2.0);

		X[0].S2[target] = prop;
	}

	if (H[0].Pi<1.0)
	{ /*when BayesB */
		for (locus = 0; locus<P; locus++)
		{
			target = Order[locus];

			for (record = 0, ProbInclude = 0.0, ProbExclude = 0.0; record<Ny; record++)
			{
				ProbInclude += pow((Y[0].expErrors[record] + (X[0].expGamma[target] - 1.0) * X[0].expEffect[target] * X[0].covariates[target * Nx + YtoX[record]]), 2.0);
				ProbExclude += pow((Y[0].expErrors[record] + X[0].expGamma[target] * X[0].expEffect[target] * X[0].covariates[target * Nx + YtoX[record]]), 2.0);
			}

			ProbInclude += X[0].x2[target] * X[0].varEffect[target];
			ProbInclude *= -0.5 * Tau0[0];
			ProbExclude *= -0.5 * Tau0[0];
			ProbInclude += logPi;
			ProbExclude += log1minusPi;

			temp = ProbInclude;
			if (temp < ProbExclude) temp = ProbExclude;
			ProbInclude -= temp;
			ProbExclude -= temp;

			ProbInclude = exp(ProbInclude);
			ProbExclude = exp(ProbExclude);
			prop = ProbInclude / (ProbInclude + ProbExclude);
			prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

			// update sumVarB and expErrors
			sumVarB[0] += (X[0].x2[target] * (prop2 * X[0].exp2Effect[target] - pow(prop * X[0].expEffect[target], 2.0)));
			for (record = 0; record<Ny; record++)
				Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target] * (X[0].expGamma[target] - prop));

			// for convergence check
			Check1[0] += pow((prop - X[0].expGamma[target]), 2.0);
			Check2[0] += pow(prop, 2.0);

			X[0].expGamma[target] = prop;
			X[0].exp2Gamma[target] = prop2;
		}
	}
}


void	Initialize_BayesC(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, Hstruct *H, int RandomIni)
{

	/* For repeat statement */
	int		record, locus;

	/* Temporary objects */
	double	temp;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	temp = 1.0;
	X[0].expSigma2[0] = temp;
	X[0].S2[0] = 1.0 / temp;

	if (H[0].Pi<1.0)
	{/* BayesC */
		if (RandomIni)
		{
			for (locus = 0; locus<P; locus++)
			{
				X[0].expGamma[locus] = 0.5;
				X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
				X[0].expEffect[locus] = RNormal(0.0, sqrt(temp / ((double)P * H[0].Pi)));
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			}
		}
		else {
			for (locus = 0; locus<P; locus++)
			{
				X[0].expEffect[locus] = 0.0;
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
				X[0].expGamma[locus] = 0.5;
				X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
			}
		}
	}
	else
	{/* BRR */
		if (RandomIni)
		{
			for (locus = 0; locus<P; locus++)
			{
				X[0].expGamma[locus] = 1.0;
				X[0].exp2Gamma[locus] = 1.0;
				X[0].expEffect[locus] = RNormal(0.0, sqrt(temp / (double)P));
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			}
		}
		else {
			for (locus = 0; locus<P; locus++)
			{
				X[0].expEffect[locus] = 0.0;
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
				X[0].expGamma[locus] = 1.0;
				X[0].exp2Gamma[locus] = 1.0;
			}
		}
	}

	for (locus = 0; locus<P; locus++)
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus] * X[0].expGamma[locus];
}


void	UpdateB_BayesC(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *sumVarB, double *sumGamma, double *sumGammaB2,
	int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2, temp3;

	/* New values */
	double	prop, prop2;

	/* Product of hyperparameters */
	double	vS2;

	/* For variable selection */
	double	logPi;

	/* Probability of gamma */
	double	ProbInclude, ProbIncludeConstant;

	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi);

	/* update of B */
	if (H[0].Pi<1.0)
	{	/* BayesC */
		sumGammaB2[1] = 0.0;
		sumGamma[1] = 0.0;
		ProbIncludeConstant = 0.5 * Digamma(0.5*(H[0].v + sumGamma[0])) - 0.5*log(0.5*(sumGammaB2[0] + vS2)) + logPi;
	}
	else
	{	/* BRR */
		sumGammaB2[0] = 0.0;
	}

	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];

		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target] * X[0].expGamma[target]);

		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0] + 1.0 / X[0].S2[0]);

		/* prop: E[B], prop2: E[B^2], temp2: V[B] */
		prop = temp * temp2;
		prop2 = pow(prop, 2.0) + temp2;

		if (H[0].Pi<1.0)
		{	/* BayesC */
			/* update Gamma */
			ProbInclude = 0.5 * temp2 * temp * temp + 0.5 * log(temp2);
			ProbInclude += ProbIncludeConstant;
			if (ProbInclude>20.0) ProbInclude = 20.0; /* to avoid overflow */
			ProbInclude = exp(ProbInclude);
			temp3 = ProbInclude / (ProbInclude + 1.0 - H[0].Pi);

			/* update residuals */
			for (record = 0; record<Ny; record++)
			{
				Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * X[0].expGamma[target] * X[0].expEffect[target]);
				Y[0].expErrors[record] -= (X[0].covariates[target * Nx + YtoX[record]] * temp3 * prop);
			}
			X[0].expGamma[target] = temp3;
			X[0].exp2Gamma[target] = pow(X[0].expGamma[target], 2.0) + X[0].expGamma[target] * (1.0 - X[0].expGamma[target]);
			sumVarB[0] += (X[0].x2[target] * X[0].expGamma[target] * (prop2 - X[0].expGamma[target] * prop * prop));
			sumGammaB2[1] += prop2 * X[0].expGamma[target];
			sumGamma[1] += X[0].expGamma[target];
			/*---- Note -----------------------------------------------------------------------------------------------------
			sumGamma[0] and sumGammaB2[0] are not updated here, because these are the posterior parameters of X[0].expSigma.
			These values are replaced by sumGamma[1] and sumGammaB2[1] at the update of X[0].expSigma.
			---------------------------------------------------------------------------------------------------------------*/
		}
		else
		{	/* BRR */
			for (record = 0; record<Ny; record++)
				Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * (X[0].expEffect[target] - prop));
			sumVarB[0] += (X[0].x2[target] * temp2);
			sumGammaB2[0] += prop2;
		}

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);
		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}

	/* update of Tau2 or Sigma2*/
	if (H[0].Pi<1.0)
	{	/* BayesC */
		sumGammaB2[0] = sumGammaB2[1];
		sumGamma[0] = sumGamma[1];
		temp = sumGammaB2[0] + vS2;
		X[0].expSigma2[0] = temp / (H[0].v + sumGamma[0] - 2.0);
		prop = temp / (H[0].v + sumGamma[0]);
	}
	else
	{	/* BRR */
		temp = sumGammaB2[0] + vS2;
		X[0].expSigma2[0] = temp / (H[0].v + (double)P - 2.0);
		prop = temp / (H[0].v + (double)P);
	}
	Check1[0] += pow((prop - X[0].S2[0]), 2.0);
	Check2[0] += pow(prop, 2.0);
	X[0].S2[0] = prop;
}


void	Initialize_SSVS(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, Hstruct *H, int RandomIni)
{
	/* For repeat statement */
	int		record, locus;

	/* Temporary objects */
	double	temp;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	temp = 1.0;
	X[0].expSigma2[0] = temp / 2.0;
	X[0].S2[0] = 2.0 / temp;
	if (RandomIni)
	{
		for (locus = 0; locus<2; locus++)
		{
			X[0].expGamma[locus] = 0.5;
			X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
			X[0].expEffect[locus] = RNormal(0.0, sqrt(X[0].expSigma2[0] / (double)P));
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
		}
	}
	else
	{
		for (locus = 0; locus<2; locus++)
		{
			X[0].expEffect[locus] = 0.0;
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			X[0].expGamma[locus] = 0.5;
			X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
		}
	}

	for (locus = 0; locus<P; locus++)
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus];
}


void	UpdateB_SSVS(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *sumVarB, double *sumGammaB2,
	int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4;

	/* New values */
	double	prop, prop2;

	/* Product of hyperparameters */
	double	vS2;

	/* For variable selection */
	double	logPi, log1minusPi;

	/* Probability of gamma */
	double	ProbInclude, ProbExclude;

	/* Used when SSVS */
	double	invC, OneMinusInvC, logC;

	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi);
	log1minusPi = log(1.0 - H[0].Pi);
	invC = 1.0 / H[0].c;
	OneMinusInvC = 1.0 - invC;
	logC = log(H[0].c);

	sumGammaB2[0] = 0.0;
	sumGammaB2[1] = 0.0;
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];

		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target]);

		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0] + (X[0].expGamma[target] * OneMinusInvC + invC) / X[0].S2[0]);

		/* prop: E[B], prop2: E[B^2], temp2: V[B] */
		prop = temp * temp2;
		prop2 = pow(prop, 2.0) + temp2;

		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * (X[0].expEffect[target] - prop));

		sumVarB[0] += (X[0].x2[target] * temp2);
		temp3 = prop2;
		sumGammaB2[0] += temp3 * X[0].expGamma[target];
		sumGammaB2[1] += temp3 * (1.0 - X[0].expGamma[target]);

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);
		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}

	/* update of Tau2 or Sigma2*/
	temp = sumGammaB2[0] + sumGammaB2[1] * invC + vS2;
	X[0].expSigma2[0] = temp / (H[0].v + (double)P - 2.0);
	prop = temp / (H[0].v + (double)P);

	Check1[0] += pow((prop - X[0].S2[0]), 2.0);
	Check2[0] += pow(prop, 2.0);
	X[0].S2[0] = prop;

	/* Update of Gamma (for wBSR, SSVS, and MIX) */
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];
		temp4 = X[0].exp2Effect[target];

		ProbInclude = -0.5 * temp4 / X[0].S2[0] + logPi;
		ProbExclude = -0.5 * temp4 / X[0].S2[0] * invC + log1minusPi - 0.5 * logC;

		temp = ProbInclude;
		if (temp < ProbExclude) temp = ProbExclude;
		ProbInclude -= temp;
		ProbExclude -= temp;

		ProbInclude = exp(ProbInclude);
		ProbExclude = exp(ProbExclude);
		prop = ProbInclude / (ProbInclude + ProbExclude);
		prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

		/* for convergence check */
		Check1[0] += pow((prop - X[0].expGamma[target]), 2.0);
		Check2[0] += pow(prop, 2.0);

		X[0].expGamma[target] = prop;
		X[0].exp2Gamma[target] = prop2;
	}
}


void	Initialize_MIX(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, Hstruct *H, int RandomIni)
{

	/* For repeat statement */
	int		record, locus;

	/* Temporary objects */
	double	temp;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	temp = 1.0;
	X[0].expSigma2[0] = temp / 2.0;
	X[0].S2[0] = 2.0 / temp;
	X[0].expSigma2[1] = temp / 2.0;
	X[0].S2[1] = 2.0 / temp;
	if (RandomIni)
	{
		for (locus = 0; locus<2; locus++)
		{
			X[0].expGamma[locus] = 0.5;
			X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
			X[0].expEffect[locus] = RNormal(0.0, sqrt(X[0].expSigma2[locus] / (double)P));
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
		}
	}
	else
	{
		for (locus = 0; locus<2; locus++)
		{
			X[0].expEffect[locus] = 0.0;
			X[0].varEffect[locus] = 0.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			X[0].expGamma[locus] = 0.5;
			X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
		}
	}

	for (locus = 0; locus<P; locus++)
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus];
}


void	UpdateB_MIX(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *sumVarB, double *sumGamma, double *sumGammaB2,
	int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4;

	/* New values */
	double	prop, prop2;

	/* Product of hyperparameters */
	double	vS2, vcS2;

	/* For variable selection */
	double	logPi, log1minusPi;

	/* Probability of gamma */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	vS2 = H[0].v * H[0].S2;
	vcS2 = H[0].c * H[0].v * H[0].S2;
	logPi = log(H[0].Pi);
	log1minusPi = log(1.0 - H[0].Pi);

	sumGammaB2[0] = 0.0;
	sumGammaB2[1] = 0.0;
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];

		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target]);
		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0] + X[0].expGamma[target] / X[0].S2[0] + (1.0 - X[0].expGamma[target]) / X[0].S2[1]);

		/* prop: E[B], prop2: E[B^2], temp2: V[B] */
		prop = temp * temp2;
		prop2 = pow(prop, 2.0) + temp2;

		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * (X[0].expEffect[target] - prop));

		sumVarB[0] += (X[0].x2[target] * temp2);
		temp3 = prop2;

		sumGammaB2[0] += temp3 * X[0].expGamma[target];
		sumGammaB2[1] += temp3 * (1.0 - X[0].expGamma[target]);

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);
		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}

	/* update of Tau2 or Sigma2*/
	temp = sumGammaB2[0] + vS2;
	X[0].expSigma2[0] = temp / (H[0].v + sumGamma[0] - 2.0);
	prop = temp / (H[0].v + sumGamma[0]);

	Check1[0] += pow((prop - X[0].S2[0]), 2.0);
	Check2[0] += pow(prop, 2.0);
	X[0].S2[0] = prop;

	temp = sumGammaB2[1] + vcS2;
	X[0].expSigma2[1] = temp / (H[0].v + (double)P - sumGamma[0] - 2.0);
	prop = temp / ((double)P - sumGamma[0]);

	Check1[0] += pow((prop - X[0].S2[1]), 2.0);
	Check2[0] += pow(prop, 2.0);
	X[0].S2[1] = prop;


	/* Update of Gamma */
	ProbIncludeConstant = Digamma(0.5*(H[0].v + sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[0] + vS2)) + logPi;
	ProbExcludeConstant = Digamma(0.5*(H[0].v + (double)P - sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[1] + vS2)) + log1minusPi;
	/*---- Note -----------------------------------------------------------------------------------------------------
	sumGamma[0] and sumGammaB2 are not updated here, because these are the posterior parameters of X[0].expSigma.
	sumGamma[0] is replaced by sumGamma[1] after the update of all gamma.
	sumGammaB2 are updated when marker effects are updated.
	---------------------------------------------------------------------------------------------------------------*/
	for (locus = 0, sumGamma[1] = 0.0; locus<P; locus++)
	{
		target = Order[locus];
		temp4 = X[0].exp2Effect[target];

		ProbInclude = ProbIncludeConstant - 0.5 * temp4 / X[0].S2[0];
		ProbExclude = ProbExcludeConstant - 0.5 * temp4 / X[0].S2[1];

		temp = ProbInclude;
		if (temp < ProbExclude) temp = ProbExclude;
		ProbInclude -= temp;
		ProbExclude -= temp;

		ProbInclude = exp(ProbInclude);
		ProbExclude = exp(ProbExclude);
		prop = ProbInclude / (ProbInclude + ProbExclude);
		prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

		/* for convergence check */
		Check1[0] += pow((prop - X[0].expGamma[target]), 2.0);
		Check2[0] += pow(prop, 2.0);

		X[0].expGamma[target] = prop;
		X[0].exp2Gamma[target] = prop2;
		sumGamma[1] += X[0].expGamma[target];
	}
	sumGamma[0] = sumGamma[1];
}


void	Initialize_BayesB(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, Hstruct *H, int RandomIni)
{
	/* For repeat statement */
	int		record, locus;

	/* Temporary objects */
	double	temp;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	temp = 1.0;
	if (H[0].Pi<1.0)
	{/* BayesB */
		if (RandomIni)
		{
			for (locus = 0; locus<P; locus++)
			{
				X[0].expSigma2[locus] = temp / ((double)P * H[0].Pi);
				X[0].S2[locus] = (double)P * H[0].Pi / temp;
				X[0].expGamma[locus] = 0.5;
				X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
				X[0].expEffect[locus] = RNormal(0.0, sqrt(X[0].expSigma2[locus]));
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			}
		}
		else {
			for (locus = 0; locus<P; locus++)
			{
				X[0].expEffect[locus] = 0.0;
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
				X[0].expSigma2[locus] = temp / ((double)P * H[0].Pi);
				X[0].S2[locus] = (double)P * H[0].Pi / temp;
				X[0].expGamma[locus] = 0.5;
				X[0].exp2Gamma[locus] = pow(X[0].expGamma[locus], 2.0) + X[0].expGamma[locus] * (1.0 - X[0].expGamma[locus]);
			}
		}
	}
	else
	{/* BayesA */
		if (RandomIni)
		{
			for (locus = 0; locus<P; locus++)
			{
				X[0].expSigma2[locus] = temp / (double)P;
				X[0].S2[locus] = (double)P / temp;
				X[0].expGamma[locus] = 1.0;
				X[0].exp2Gamma[locus] = 1.0;
				X[0].expEffect[locus] = RNormal(0.0, sqrt(X[0].expSigma2[locus]));
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			}
		}
		else {
			for (locus = 0; locus<P; locus++)
			{
				X[0].expEffect[locus] = 0.0;
				X[0].varEffect[locus] = 0.0;
				X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
				X[0].expSigma2[locus] = temp / (double)P;
				X[0].S2[locus] = (double)P / temp;
				X[0].expGamma[locus] = 1.0;
				X[0].exp2Gamma[locus] = 1.0;
			}
		}
	}

	for (locus = 0; locus<P; locus++)
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus] * X[0].expGamma[locus];
}


void	UpdateB_BayesB(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *sumVarB, int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2, temp3;

	/* New values */
	double	prop, prop2;

	/* Product of hyperparameters */
	double	vS2;

	/* For variable selection */
	double	logPi;

	/* Probability of gamma */
	double	ProbInclude, ProbIncludeConstant;

	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi);

	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];

		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target] * X[0].expGamma[target]);

		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0] + 1.0 / X[0].S2[target]);

		/* prop: E[B], prop2: E[B^2], temp2: V[B] */
		prop = temp * temp2;
		prop2 = pow(prop, 2.0) + temp2;

		if (H[0].Pi<1.0)
		{	/* BayesB */
			/* update Gamma */
			ProbIncludeConstant = 0.5 * Digamma(0.5*(H[0].v + X[0].expGamma[target])) - 0.5*log(0.5*(X[0].expGamma[target] * X[0].exp2Effect[target] + vS2)) + logPi;
			ProbInclude = 0.5 * temp2 * temp * temp + 0.5 * log(temp2);
			ProbInclude += ProbIncludeConstant;
			if (ProbInclude>20.0) ProbInclude = 20.0; /* to avoid overflow */
			ProbInclude = exp(ProbInclude);
			temp3 = ProbInclude / (ProbInclude + 1.0 - H[0].Pi);

			/* update residuals */
			for (record = 0; record<Ny; record++)
			{
				Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * X[0].expGamma[target] * X[0].expEffect[target]);
				Y[0].expErrors[record] -= (X[0].covariates[target * Nx + YtoX[record]] * temp3 * prop);
			}
			X[0].expGamma[target] = temp3;
			X[0].exp2Gamma[target] = pow(X[0].expGamma[target], 2.0) + X[0].expGamma[target] * (1.0 - X[0].expGamma[target]);
			sumVarB[0] += (X[0].x2[target] * X[0].expGamma[target] * (prop2 - X[0].expGamma[target] * prop * prop));
		}
		else
		{	/* BayesA */
			for (record = 0; record<Ny; record++)
				Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * (X[0].expEffect[target] - prop));
			sumVarB[0] += (X[0].x2[target] * temp2);
		}

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);
		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}

	/* update of Tau2 or Sigma2*/
	if (H[0].Pi<1.0)
	{	/* BayesB */
		for (locus = 0; locus<P; locus++)
		{
			target = Order[locus];
			temp = X[0].expGamma[target] * X[0].exp2Effect[target] + vS2;
			X[0].expSigma2[target] = temp / (H[0].v + X[0].expGamma[target] - 2.0);
			prop = temp / (H[0].v + X[0].expGamma[target]);
			Check1[0] += pow((prop - X[0].S2[target]), 2.0);
			Check2[0] += pow(prop, 2.0);
			X[0].S2[target] = prop;
		}
	}
	else
	{	/* BayesA */
		for (locus = 0; locus<P; locus++)
		{
			target = Order[locus];
			temp = X[0].exp2Effect[target] + vS2;
			X[0].expSigma2[target] = temp / (H[0].v - 1.0);
			prop = temp / (H[0].v + 1.0);
			Check1[0] += pow((prop - X[0].S2[target]), 2.0);
			Check2[0] += pow(prop, 2.0);
			X[0].S2[target] = prop;
		}
	}
}


void	Initialize_GBLUP(int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, double Tau0, int RandomIni)
{
	/* For repeat statement */
	int		record, level;

	/* Temporary objects */
	double	temp;

	/* Diagonal of posterior variance */
	double *D;

	D = (double*)calloc(Nx, sizeof(double));

	X[0].expSigma2[0] = 1.0;
	X[0].S2[0] = 1.0;
	X[0].exp2Effect[0] = 0.0;
	X[0].expDelta2 = 0.0;//used for determinant

	for (level = 0; level < Nx; level++)
		D[level] = X[0].S2[0] / (X[0].evalue[level] + X[0].S2[0] * Tau0);

	Innerproduct_tADA(X[0].tevector, D, Nx, Nx, X[0].varEffect);

	temp = sqrt(X[0].expSigma2[0]);
	for (level = 0; level < Nx; level++)
	{
		if (RandomIni) X[0].expEffect[level] = RNormal(0.0, temp); else X[0].expEffect[level] = 0.0;
		X[0].expDelta2 += log(X[0].evalue[level]);
	}

	for (record = 0; record<Ny; record++)
		Y[0].expErrors[record] -= X[0].expEffect[YtoX[record]];

	free(D);
}



void	UpdateB_GBLUP(int Nx, int Ny, int *XtoY, Ystruct *Y, Xstruct *X, Hstruct *H, double Tau0, double *sumVarB,  double *Check1, double *Check2)
{
	/* For repeat statement */
	int		level, level2;

	/* Temporary objects */
	double	temp, temp2, a;

	/* New values */
	double	*prop;

	/* Summation of hyperparameters */
	double	vN;

	/* Other variables */
	double	*Yr, *D;

	prop = (double*)calloc(Nx, sizeof(double));
	Yr = (double*)calloc(Nx, sizeof(double));
	D = (double*)calloc(Nx, sizeof(double));
	vN = H[0].v + (double)Nx;
	a = X[0].S2[0] * Tau0;

	/* Random effect */
	for (level = 0; level < Nx; level++)
	{
		D[level] = X[0].S2[0] / (X[0].evalue[level] + a);
		Yr[level] = Y[0].expErrors[XtoY[level]] + X[0].expEffect[level];
	}

	Innerproduct_tADA(X[0].tevector, D, Nx, Nx, X[0].varEffect);
	Innerproduct_tAB(X[0].varEffect, Yr, Nx, Nx, 1, prop);
	for (level = 0; level < Nx; level++)
	{
		prop[level] *= Tau0;
		Y[0].expErrors[XtoY[level]] += (X[0].expEffect[level] - prop[level]);
		Check1[0] += pow((prop[level] - X[0].expEffect[level]), 2.0);
		Check2[0] += pow(prop[level], 2.0);
		X[0].expEffect[level] = prop[level];
		sumVarB[0] += X[0].varEffect[level * Nx + level];
	}

	/* Sigma2 */
	for (level = 0, X[0].exp2Effect[0] = 0.0; level < Nx; level++)
	{
		for (level2 = 0, temp = 0.0; level2 < Nx; level2++)
		{
			temp += X[0].expEffect[level2] * X[0].evector[level * Nx + level2];
		}
		X[0].exp2Effect[0] += X[0].evalue[level] * temp * temp;
		X[0].exp2Effect[0] += (X[0].evalue[level] * X[0].S2[0]) / (X[0].S2[0] * Tau0 + X[0].evalue[level]);
	}//X[0].exp2Effect[0] is the expectation of U'iKU
	temp2 = X[0].exp2Effect[0] + H[0].v * H[0].S2;
	X[0].S2[0] = temp2 / vN;
	prop[0] = vN * X[0].S2[0] / (vN - 2.0);
	Check1[0] += pow((prop[0] - X[0].expSigma2[0]), 2.0);
	Check2[0] += pow(prop[0], 2.0);
	X[0].expSigma2[0] = prop[0];

	free(prop);  free(Yr); free(D);
}


void	Initialize_FIXED(int P, int Ny, int Nx, int *YtoX, Xstruct *X, Ystruct *Y, int RandomIni)
{
	/* For repeat statement */
	int		record, locus;

	for (locus = 0; locus<P; locus++)
	{
		X[0].x2[locus] = 0.0;
		for (record = 0; record<Ny; record++)
			X[0].x2[locus] += pow(X[0].covariates[locus * Nx + YtoX[record]], 2.0);
	}

	if (RandomIni)
	{
		for (locus = 0; locus<P; locus++)
		{
			X[0].expEffect[locus] = RNormal(0.0, 1.0);
			X[0].varEffect[locus] = 1.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			for (record = 0; record<Ny; record++)
				Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus];
		}
	}
	else
	{
		for (locus = 0; locus<P; locus++)
		{
			X[0].expEffect[locus] = 0.0;
			X[0].varEffect[locus] = 1.0;
			X[0].exp2Effect[locus] = pow(X[0].expEffect[locus], 2.0) + X[0].varEffect[locus];
			for (record = 0; record<Ny; record++)
				Y[0].expErrors[record] -= X[0].covariates[locus * Nx + YtoX[record]] * X[0].expEffect[locus];
		}
	}

}


void	UpdateB_FIXED(int P, int Ny, int Nx, int *YtoX, Ystruct *Y, Xstruct *X, double *Tau0, double *sumVarB, int *Order, double *Check1, double *Check2)
{
	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp, temp2;

	/* New values */
	double	prop, prop2;

	// update of fixed effects
	for (locus = 0; locus<P; locus++)
	{
		target = Order[locus];

		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += X[0].covariates[target * Nx + YtoX[record]] * (Y[0].expErrors[record] + X[0].covariates[target * Nx + YtoX[record]] * X[0].expEffect[target]);

		temp *= Tau0[0];
		temp2 = 1.0 / (X[0].x2[target] * Tau0[0]);

		prop = temp * temp2;
		prop2 = prop * prop + temp2;
		for (record = 0; record<Ny; record++)
			Y[0].expErrors[record] += (X[0].covariates[target * Nx + YtoX[record]] * (X[0].expEffect[target] - prop));

		sumVarB[0] += (X[0].x2[target] * temp2);

		Check1[0] += pow((prop - X[0].expEffect[target]), 2.0);
		Check2[0] += pow(prop, 2.0);
		X[0].expEffect[target] = prop;
		X[0].exp2Effect[target] = prop2;
		X[0].varEffect[target] = temp2;
	}
}


/*---Perform genome-wide regression----------------------------------------------------------------------------------------------------------------------------*/
void GenomeWideRegression(int *Methodcode, int Nm, int *P, int Ny, int Nx, int *YtoX, int *XtoY, double Thresholdvalue, int Maxiteration,
	Ystruct *Y, Xstruct *X, Hstruct *H, double *Tau0, double *LB, double *Rmonitor, double *FittingTime, int RandomIni)
{
	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		record, locus, method, method_s;

	/* Total number of covariates */
	int		*TotalP = NULL;

	/* Update orders, timse of shuffling orders, locus to be updated*/
	int		*Order = NULL, *Times = NULL;
	int		*Morder = NULL;

	/* Temporary objects */
	double	temp;

	/* New values */
	double	prop;

	/* For update of Re2 */
	double	sumVarB = 0.0, sumTau2B2 = 0.0;

	/* For update of Gamma */
	double sumGamma[2] = { 0.0, 0.0 }, sumGammaB2[2] = { 0.0, 0.0 };

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1 = 0.0, b1 = 0.0, *a2 = NULL, *b2 = NULL;
	int		*UseBL = NULL;

	/* check convergence */
	double	Check1 = 0.0, Check2 = 0.0;

	/* calculate computation time */
	clock_t start, end;

	/* initialization and calculatio of residual errors */
	UseBL = (int*)calloc(Nm, sizeof(int));
	a2 = (double*)calloc(Nm, sizeof(double));
	b2 = (double*)calloc(Nm, sizeof(double));
	temp = 1.0;
	Tau0[0] = 100.0 / temp;
	Tau0[1] = temp;
	for (method = 0; method < Nm; method++)
	{
		switch (Methodcode[method])
		{
		case 1:
			Initialize_BL(P[method], Ny, Nx, YtoX, X + method, Y, RandomIni);
			break;
		case 2:
			Initialize_EBL(P[method], Ny, Nx, YtoX, X + method, Y, RandomIni);
			break;
		case 3:
			Initialize_wBSR(P[method], Ny, Nx, YtoX, X + method, Y, H + method, RandomIni);
			break;
		case 4:
			Initialize_BayesC(P[method], Ny, Nx, YtoX, X + method, Y, H + method, RandomIni);
			break;
		case 5:
			Initialize_SSVS(P[method], Ny, Nx, YtoX, X + method, Y, H + method, RandomIni);
			break;
		case 6:
			Initialize_MIX(P[method], Ny, Nx, YtoX, X + method, Y, H + method, RandomIni);
			break;
		case 7:
			Initialize_BayesB(P[method], Ny, Nx, YtoX, X + method, Y, H + method, RandomIni);
			break;
		case 8:
			Initialize_GBLUP(Ny, Nx, YtoX, X + method, Y, Tau0[0], RandomIni);
			break;
		case 9:
			Initialize_FIXED(P[method], Ny, Nx, YtoX, X + method, Y, RandomIni);
		}

		//check whether BL or EBL is used
		if (Methodcode[method] == 1 || Methodcode[method] == 2) UseBL[method] = 1;
	}

	/* update order of marker effects */
	TotalP = (int*)calloc(Nm + 1, sizeof(int));
	Times = (int*)calloc(Nm, sizeof(int));
	for (method = 1; method <= Nm; method++) TotalP[method] += TotalP[method-1] + P[method-1];
	Order = (int*)calloc(TotalP[Nm], sizeof(int));
	for (method = 0; method < Nm; method++)
	{
		for (locus = 0; locus<P[method]; locus++) Order[TotalP[method] + locus] = locus;
		Shuffle(Order + TotalP[method], P[method], P[method]);
		Times[method] = P[method] / 10;
	}

	/* update order of regression methods */
	Morder = (int*)calloc(Nm, sizeof(int));
	for (method = 0; method < Nm; method++) Morder[method] = method;
	Shuffle(Morder, Nm, Nm);
	start = clock();

	/* start optimization */
	for (ite = 1; ite <= Maxiteration; ite++)
	{
		// To check convergence
		Check1 = 0.0; Check2 = 0.0;

		// For update of residual variance
		sumVarB = 0.0;
		sumTau2B2 = 0.0;

		// Update of B
		Shuffle(Morder, Nm, Nm);
		for (method = 0; method < Nm; method++)
		{
			method_s = Morder[method];
			Shuffle(Order + TotalP[method_s], P[method_s], Times[method_s]);

			switch (Methodcode[method_s]) {
			case 1:
				UpdateB_BL(P[method_s], Ny, Nx, YtoX, Y, X + method_s, H + method_s, Tau0, &sumVarB, &sumTau2B2, a2 + method_s, b2 + method_s, Order + TotalP[method_s], &Check1, &Check2);
				break;
			case 2:
				UpdateB_EBL(P[method_s], Ny, Nx, YtoX, Y, X + method_s, H + method_s, Tau0, &sumVarB, &sumTau2B2, a2 + method_s, b2 + method_s, Order + TotalP[method_s], &Check1, &Check2);
				break;
			case 3:
				UpdateB_wBSR(P[method_s], Ny, Nx, YtoX, Y, X + method_s, H + method_s, Tau0, &sumVarB, Order + TotalP[method_s], &Check1, &Check2);
				break;
			case 4:
				UpdateB_BayesC(P[method_s], Ny, Nx, YtoX, Y, X + method_s, H + method_s, Tau0, &sumVarB, sumGamma, sumGammaB2, Order + TotalP[method_s], &Check1, &Check2);
				break;
			case 5:
				UpdateB_SSVS(P[method_s], Ny, Nx, YtoX, Y, X + method_s, H + method_s, Tau0, &sumVarB, sumGammaB2, Order + TotalP[method_s], &Check1, &Check2);
				break;
			case 6:
				UpdateB_MIX(P[method_s], Ny, Nx, YtoX, Y, X + method_s, H + method_s, Tau0, &sumVarB, sumGamma, sumGammaB2, Order + TotalP[method_s], &Check1, &Check2);
				break;
			case 7:
				UpdateB_BayesB(P[method_s], Ny, Nx, YtoX, Y, X + method_s, H + method_s, Tau0, &sumVarB, Order + TotalP[method_s], &Check1, &Check2);
				break;
			case 8:
				UpdateB_GBLUP(Nx, Ny, XtoY, Y, X + method_s, H + method_s, Tau0[0], &sumVarB, &Check1, &Check2);
				break;
			case 9:
				UpdateB_FIXED(P[method_s], Ny, Nx, YtoX, Y, X + method_s, Tau0, &sumVarB, Order + TotalP[method_s], &Check1, &Check2);
				break;
			}
		}

		// Update of residual precision
		for (record = 0, temp = 0.0; record<Ny; record++)
			temp += pow(Y[0].expErrors[record], 2.0);

		a1 = (double)Ny * 0.5;
		for (method = 0; method < Nm; method++)
			if (UseBL[method]) a1 += 0.5 * (double)P[method];

		b1 = 0.5*(temp + sumVarB + sumTau2B2);

		prop = a1 / b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		// Lower bound
		Rmonitor[ite - 1] = 1.0 / Tau0[0];
		temp = Check1 / Check2;
		/*
		LB[ite - 1] = -0.5 * (double)Ny * log2pi - a1 * log(b1) + mylgamma(a1);
		for (method = 0; method < Nm; method++)
		{
			LB[ite - 1] += LowerBound(Methodcode[method], Ny, P[method], X + method, a2[method], b2[method], H + method, Tau0[0]);
		}
		*/
		if (ite == Maxiteration || temp<Thresholdvalue)
		{
			LB[0] = - 0.5 * (double) Ny * log2pi - a1 * log(b1) + mylgamma(a1);
			for (method = 0; method < Nm; method++)
			{
				LB[0] += LowerBound(Methodcode[method], Ny, P[method], X + method, a2[method], b2[method], H + method, Tau0[0]);
			}
			break;
		}
	}

	end = clock();
	FittingTime[0] = (double)(end - start) / CLOCKS_PER_SEC;

	free(Order); free(Morder); free(TotalP); free(Times); free(UseBL); free(a2); free(b2);
}
