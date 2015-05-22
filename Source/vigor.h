/*
	Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata

	Released under the MIT license
	http://opensource.org/licenses/mit-license.php
*/

/*----define global objects------------------------------------------------------------------------------------------------------*/
#define Version				"1.0 (May 2015)"
#define	Missingvalue		99999999		/* Missing value used in this program */
#define	MaxBufSize			1048576			/* Maximum buffer size (32^4) */

/*----define structure-----------------------------------------------------------------------------------------------------------*/
typedef struct {

	double	*observations;	/* observations */
	double	*stobs;			/* standardized observations */
	double	*expErrors;		/* expected errors */
	int		nmiss;			/* number of missing records */
	int		*use;			/* used records (non-missing record) */
	char	name[101];		/* name of trait */

} Ystruct;

typedef struct {

	double	*covariates;	/* covariates */
	double	x2;				/* sum of squared covariates */
	double	expEffect;		/* expectation of marker effect. When BayesC and BayesB, this indicates E[B|gamma=1] */
	double	varEffect;		/* variance of marker effect. When BayesC and BayesB, this indicates V[B|gamma=1] */
	double	exp2Effect;		/* 2nd moment (for wBSR, BayesB, BayesC, SSVS, and MIX). When BayesC and BayesB, this indicates E[B^2|gamma=1] */
	double	expGamma;		/* expectation of gamma (for wBSR, SSVS, MIX, BayesC, and BayesB) */
	double	exp2Gamma;		/* 2nd moment (for wBSR) */
	double	expTau2;		/* expectation of tau (for BL and EBL) */
	double	expInTau2;		/* expectation of the inverse of tau (for BL and EBL) */
//	double	varTau2;		/* variance of tau (for BL and EBL) */
	double	expEta2;		/* expectation of eta (for EBL) */
	double	expSigma2;		/* expectation of sigma (for wBSR, BayesB, BayesC, SSVS, and MIX).*/
	double	S2;				/* parameter of the posterior distribution of marker effect variance (for wBSR, BayesB, BayesC, SSVS, and MIX) */
	double	a3;				/* shape parameter of the posterior gamma distribution of eta2 (for EBL) */
	double	b3;				/* rate parameter of the posterior gamma distribution of eta2 (for EBL) */

} Xstruct;

typedef struct {

	double	deltaShape;		/* for BL and EBL */
	double	deltaRate;		/* for BL and EBL */
	double	etaShape;		/* for EBL */
	double	etaRate;		/* for EBL */
	double	v;				/* for wBSR, BayesB, BayesC, SSVS, and MIX */
	double	S2;				/* for wBSR, BayesB, BayesC, SSVS, and MIX */
	double	Pi;				/* for wBSR, BayesB, BayesC, SSVS, and MIX */
	double	c;				/* for SSVS and MIX */

} Hstruct;

typedef struct {

	int ntest;				/* Number of observations in the test set */
	int *obstest;			/* Observation in the test set */

} Astruct;

/*----generate a random integer -------------------------------------------------------------------------------------------------*/
int RandomInteger ( int low, int high )
{
	int k;
	double value;

	value = genrand_real3();
	k = (int) ( value * (double)(high - low + 1) );
	return (low + k);
}

/*----generate a real number from Uniform (0,1)------------------------------------------------------------------------------------*/
double rnd()
{
	return(genrand_real3());
}

/*----generate a real number from Gamma (n, lambda) (from STRUCTURE), E[x] = n/lambda, V[x] = n/lambda^2 --------------------------*/
/* lambda = E[x]/V[x], and n = E[x]*lambda = E[x]^2/V[x] */
double RGamma(double n,double lambda)
{
  double ar;
  double w;

	double x = 0.0;
	if ( n < 1 )
	{
		const double E = 2.71828182;
		const double b = ( n + E ) / E;
		double p = 0.0;
		one: 
		p = b * rnd ();
		if ( p > 1 ) goto two;
		x = exp ( log ( p ) / n );
		if ( x > -log ( rnd() ) ) goto one;
		goto three;
		two: 
		x = -log ( ( b - p ) / n );
		if ( ( ( n - 1 ) * log ( x ) ) < log ( rnd() ) ) goto one;
		three:;	
	}
	else 
	{
	    if ( n == 1.0 )
		{                         /* exponential random variable, from Ripley, 1987, P230  */
		    double a = 0.0;
		    double u, u0, ustar;
            ten:
		    u = rnd();
	     	u0 = u;
         	twenty:
		    ustar = rnd();
		    if ( u < ustar) goto thirty;
	    	u = rnd();
		    if ( u < ustar) goto twenty;
	    	a++;
		    goto ten;
        	thirty:
		    return ( a + u0 ) / lambda;
		}
	    else
		{
	    	double static nprev = 0.0;
	     	double static c1 = 0.0;
	    	double static c2 = 0.0;
	    	double static c3 = 0.0;
	    	double static c4 = 0.0;
	    	double static c5 = 0.0;
	    	double u1;
	    	double u2;
	    	if(n != nprev)
			{
		    	c1 = n - 1.0;
		    	ar = 1.0 / c1;
		    	c2 = ar * ( n - 1 / ( 6 * n ) );
		    	c3 = 2 * ar;
		    	c4 = c3 + 2;
		    	if( n > 2.5 ) c5 = 1 / sqrt(n);
			}
		    four:
	       	u1 = rnd();
	       	u2 = rnd();
	       	if ( n <= 2.5 ) goto five;
	       	u1 = u2 + c5 * ( 1 - 1.86 * u1);
	       	if ( ( u1 <= 0 ) || ( u1 >= 1 ) ) goto four;
	       	five:
	       	w = c2 * u2 / u1;
	        if( c3 * u1 + w + 1.0 /w < c4) goto six;
		    if( c3 * log ( u1 ) - log ( w ) + w >= 1) goto four;
	       	six:
	       	x = c1 * w;		
	       	nprev = n;
		}
	}

	return x/lambda;
}

/*----generate real numbers from N (0,1) ( from STRUCTURE )------------------------------------------------------------------------*/
double snorm()    /*was snorm(void) -- JKP*/
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
static double snorm,u,ss,ustar,anr,w,y,tt;
    u = rnd();   /* was ranf--JKP*/
    ss = 0.0;
    if(u > 0.5) ss = 1.0;
    u += (u-ss);
    u = 32.0*u;
    inr = (long) (u);
    if(inr == 32) inr = 31;
    if(inr == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(double)inr;
    anr = *(a+inr-1);
S40:
    if(ustar <= *(t+inr-1)) goto S60;
    w = (ustar-*(t+inr-1))**(h+inr-1);
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = anr+w;
    snorm = y;
    if(ss == 1.0) snorm = -y;
    return snorm;
S60:
/*
                                CENTER CONTINUED
*/
    u = rnd();                /*was ranf--JKP*/
    w = u*(*(a+inr)-anr);
    tt = (0.5*w+anr)*w;
    goto S80;
S70:
    tt = u;
    ustar = rnd();                /*was ranf--JKP*/
S80:
    if(ustar > tt) goto S50;
    u = rnd();               /*was ranf--JKP*/
    if(ustar >= u) goto S70;
    ustar = rnd();               /*was ranf--JKP*/
    goto S40;
S100:
/*
                                START TAIL
*/
    inr = 6;
    anr = *(a+31);
    goto S120;
S110:
    anr += *(d+inr-1);
    inr += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+inr-1);
    tt = (0.5*w+anr)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = rnd();               /*was ranf--JKP*/
    if(ustar > tt) goto S50;
    u = rnd();               /*was ranf--JKP*/
    if(ustar >= u) goto S150;
    u = rnd();               /*was ranf--JKP*/
    goto S140;
}

/*----generate real number from N(mu, sd) (from STRUCTURE)-------------------------------------------------------------------------*/
double RNormal(double mu, double sd) 
     /* Returns Normal rv with mean mu, variance sigsq.   
        Uses snorm function of Brown and Lovato.  By JKP*/
{
  return (mu + sd*snorm());
}

/*----perform Bernoulli ----------------------------------------------------------------------------------------------------------*/
double RBernoulli (double prop)
{
	double flag=0.0;
	if(rnd()<prop) flag=1.0;
	return(flag);
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
  result = log (result) + lnsqrt2pi - (zpo + 6.5) + (zpo - 0.5) * log (zpo + 6.5);
  return result;
}

/*----Shuffle an array------------------------------------------------------------------------------------------------------------*/
void Shuffle ( int *asl, int arraylength, int timesofshuffle )
{
	int ts, content, target1, target2;

	for ( ts=0; ts<timesofshuffle; ts++ )
	{
		target1 = RandomInteger ( 0, arraylength-1);
		content = asl [target1];
		target2 = RandomInteger ( 0, arraylength-1 );
		asl [target1] = asl [target2];
		asl [target2] = content;
	}
}
void ShuffleDouble ( double *asl, int arraylength, int timesofshuffle )
{
	int ts, target1, target2;
	double content;

	for ( ts=0; ts<timesofshuffle; ts++ )
	{
		target1 = RandomInteger ( 0, arraylength-1);
		content = asl [target1];
		target2 = RandomInteger ( 0, arraylength-1 );
		asl [target1] = asl [target2];
		asl [target2] = content;
	}
}

/*----Sort an array----------------------------------------------------------------------------------------------------------------*/
/* decreasing */
void Sortarray (double *a, int n)
{
	int i,j,lastsort,count;
	double	temp;

	lastsort=n-2;
	do{
		for(i=0,count=0;i<=lastsort;i++)
			if(a[i]<a[i+1]){temp=a[i+1]; a[i+1]=a[i]; a[i]=temp; j=i; count++;}
		lastsort=j;
	}while(count>0);

}

/*----return MSE between a and b. Missingvalues are ignored-------------------------------------------------------------------------*/
double MSE (double *a, double *b, int n)
{
	int ii, count=0;
	double mse;

	for(ii=0, mse=0.0; ii<n; ii++)
		if(a[ii] != Missingvalue && b[ii] != Missingvalue)
		{
			mse += pow (a[ii] - b[ii], 2.0);
			count ++;
		}

	mse /= (double)count;
	return(mse);
}

/*----calculate Pearson correlation coefficient ----------------------------------------------------------------------------------------*/
double Pearson (double *a, double *b, int n)
{
	double	meana=0.0, meanb=0.0, va=0.0, vb=0.0, cov=0.0;
	int		ii, count=0;

	for(ii=0; ii<n; ii++)
		if(a[ii] != Missingvalue && b[ii] != Missingvalue)
		{
			meana+=a[ii]; meanb+=b[ii];
			count ++;
		}
	meana/=(double)count; meanb/=(double)count;

	for(ii=0; ii<n; ii++)
		if(a[ii] != Missingvalue && b[ii] != Missingvalue)
		{
			va+= pow(a[ii]-meana,2.0);	vb+= pow(b[ii]-meanb,2.0);
			cov+=(a[ii]-meana)*(b[ii]-meanb);
		}
	return(cov/sqrt(va)/sqrt(vb));
}

/*----calculate Kendall correlation coefficient ----------------------------------------------------------------------------------------*/
double Kendall (double *a, double *b, int n)
{
	int all=n, count=0, ii, jj;
	
	for(ii=0;ii<n;ii++)
		if(a[ii]==Missingvalue||b[ii]==Missingvalue) all--;

	for(ii=0;ii<(n-1);ii++)
		for(jj=ii+1;jj<n;jj++)
			if(a[ii] != Missingvalue && a[jj] != Missingvalue && b[ii] != Missingvalue && b[jj] != Missingvalue)
				if((a[ii]-a[jj])*(b[ii]-b[jj])>0) {count++;} else {count--;}

	return(2.0*(double)count/(double)(all*(all-1)));
}

/*----calculate NDCG------------------------------------------------------------------------------------------------------------*/
double NDCG (int k, double *y, double *prediction, int n){

	double *discount, loge2=log(2.0), *orderedy1, *orderedpre, *orderedy2, max=-1e100, sum1=0.0, sum2=0.0;
	int	i, j, m;

	discount=(double*)calloc(k,sizeof(double));
	orderedy1=(double*)calloc(k,sizeof(double));
	orderedy2=(double*)calloc(k,sizeof(double));
	orderedpre=(double*)calloc(k,sizeof(double));
	for(j=0;j<k;j++)
	{
		discount[j] = loge2/log((double)(j+2));
		orderedy1[j] = max; orderedpre[j] = max;
	}

	for(i=0;i<n;i++)
		if(y[i]!=Missingvalue)
			for(j=0;j<k;j++)
				if(y[i]>orderedy1[j]){if(j==k-1){orderedy1[j]=y[i];}else{for(m=k-2;m>=j;m--){orderedy1[m+1]=orderedy1[m];} orderedy1[j]=y[i];}break;}

	for(i=0;i<n;i++)
		if(y[i]!=Missingvalue)
			for(j=0;j<k;j++)
				if(prediction[i]>orderedpre[j])
				{
					if(j==k-1)
					{
						orderedpre[j]=prediction[i];orderedy2[j]=y[i];
					}
					else
					{
						for(m=k-2;m>=j;m--){orderedpre[m+1]=orderedpre[m];orderedy2[m+1]=orderedy2[m];} 
						orderedpre[j]=prediction[i];orderedy2[j]=y[i];
					}
					break;
				}

	for(j=0;j<k;j++){
		sum1+=orderedy1[j]*discount[j]; sum2+=orderedy2[j]*discount[j];
	}

	free(discount);	free(orderedy1); free(orderedpre); free(orderedy2);
	return(sum2/sum1);
}

/*----calculate log likelihood for record i-------------------------------------------------------------------------------------------*/
/* Not used */
double LogLikelihood (Ystruct *Y, Xstruct *Q, int F, Xstruct *X, int P, int i, int Methodcode, double re)
{
	/* i indicates the position in Q and X. re indicates the precision of errors */
	int locus;
	double	temp, loglike, log2pi=1.83787706640935;

	temp=Y[0].observations[i];
	for (locus=0; locus<F; locus++)
		temp -= Q[locus].covariates[i] * Q[locus].expEffect;

	if (Methodcode==3||Methodcode==4||Methodcode==7)
	{
		for (locus=0; locus<P; locus++)
			temp -= X[locus].covariates[i] * X[locus].expEffect * X[locus].expGamma;
	}
	else
	{
		for (locus=0; locus<P; locus++)
			temp -= X[locus].covariates[i] * X[locus].expEffect;
	}
	loglike = (-0.5*re*temp*temp);
	loglike += 0.5*(log(re)-log2pi);
	return(loglike);
}

/*----calculate LPDS ----------------------------------------------------------------------------------------------------------------*/
double LPDS (Ystruct *Y, Xstruct *Q, int F, Xstruct *X, int P, Astruct *A2, int fold2, int Methodcode, double *Tau0, int N, int Nsamples)
{
	double density=0.0, *simQ, *simX, *simGamma, simTau, temp, log2pi=1.83787706640935;
	int	sample, locus, record, target;
	
	simX = (double*)calloc(P, sizeof(double));
	simQ = (double*)calloc(F, sizeof(double));
	if(Methodcode==3||Methodcode==4||Methodcode==7) simGamma = (double*)calloc(P, sizeof(double));

	for(sample=0;sample<Nsamples;sample++)
	{
		for(locus=0;locus<F;locus++) {
			simQ[locus] = RNormal(Q[locus].expEffect, sqrt(Q[locus].varEffect));
		}
		for(locus=0;locus<P;locus++) {
			simX[locus] = RNormal(X[locus].expEffect, sqrt(X[locus].varEffect));
		}
		if(Methodcode==3||Methodcode==4||Methodcode==7)
			for(locus=0;locus<P;locus++) simGamma[locus] = RBernoulli(X[locus].expGamma);
		simTau = RGamma(pow(Tau0[0],2.0)/Tau0[1], Tau0[0]/Tau0[1]);

		for(record=0;record<N; record++)
		{
			target=A2[fold2].obstest[record];
			if(Y[0].observations[target]!=Missingvalue)
			{
				temp=Y[0].observations[target];
				for (locus=0; locus<F; locus++)
					temp -= Q[locus].covariates[target] * simQ[locus];

				if (Methodcode==3||Methodcode==4||Methodcode==7)
				{
					for (locus=0; locus<P; locus++)
						temp -= X[locus].covariates[target] * simX[locus] * simGamma[locus];
				}
				else
				{
					for (locus=0; locus<P; locus++)
						temp -= X[locus].covariates[target] * simX[locus];
				}
				density += (-0.5*simTau*temp*temp);
				density += 0.5*(log(simTau)-log2pi);
			}
		}
	}
	free(simX); free(simQ); if(Methodcode==3||Methodcode==4||Methodcode==7) free(simGamma);
	return(density/(double)Nsamples);
}

/*----digamma function-----------------------------------------------------------------------------------------------------------------------*/
/* from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.C. But now cannot be found */
double Digamma ( const double x )
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

  if ( x <= 0.0 ){fprintf (stderr,"Negative value in Digamma\n" );}
  y = x;
  value = 0.0;
  if ( y <= s )
  {
    value = d1 - 1.0 / y;
    return (value);
  }
   while ( y < c )
  {
    value = value - 1.0 / y;
    y = y + 1.0;
  }
  r = 1.0 / y;
  value = value + log ( y ) - 0.5 * r;
  r = r * r;
  value = value - r * ( s3 - r * ( s4 - r * s5 ) );

  return (value);
}

/*---Lower bound of marginal log probability of data-----------------------------------------------------------------------------------------------------------*/
double LowerBound (int Methodcode, int N, int P, Xstruct *X, int F, Xstruct *Q, double a1, double b1, 
	double a2, double b2, double expDelta2, Hstruct *H, int CondResidual)
{
	int		locus;
	double	LB=0.0, log2pi=1.83787706640935, temp, vj, vS2, Pi2, v2;
	double	Gammalowerbound=1e-100, Gammaupperbound;

	Gammaupperbound = 1.0-Gammalowerbound;
	switch(Methodcode)
	{
		case 1: /* BL */
			for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
			for(locus=0;locus<P;locus++)
			{
//				LB -= 0.5*(3.0*log(X[locus].expTau2)-log(X[locus].varTau2+pow(X[locus].expTau2,2.0)));
				LB -= 0.5*log(expDelta2);
				LB += 0.5*log(X[locus].varEffect);
				LB -= 0.5*(expDelta2*X[locus].expInTau2);
				if(CondResidual==0)	LB -= 0.5*X[locus].expTau2*X[locus].exp2Effect;
			}
			LB -= (a1*log(b1) - mylgamma(a1));
			LB -= (a2*(log(b2)-1.0) - mylgamma(a2));
			LB -= H[0].deltaRate*expDelta2;
			LB += ((double)(P-N+F)*0.5*log2pi);
			if(H[0].deltaShape>0.0) LB-=mylgamma(H[0].deltaShape);			/* When informative prior is used */
			if(H[0].deltaRate>0.0) LB+=H[0].deltaShape*log(H[0].deltaRate);
			LB += ((double)(2*P+F)*0.5 - (double)P*log(2.0));
			break;
		case 2: /* EBL */
			for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
			for(locus=0;locus<P;locus++)
			{
//				LB -= 0.5*(3.0*log(X[locus].expTau2)-log(X[locus].varTau2+pow(X[locus].expTau2,2.0)));
				LB -= 0.5*log(expDelta2*X[locus].expEta2);
				LB += 0.5*log(X[locus].varEffect);
				LB -= (X[locus].a3*(log(X[locus].b3)-1.0) - mylgamma(X[locus].a3));
				LB -= 0.5*(expDelta2*X[locus].expEta2*X[locus].expInTau2);
				LB -= H[0].etaRate*X[locus].expEta2;
				if(CondResidual==0)	LB -= 0.5*X[locus].expTau2*X[locus].exp2Effect;
			}
			LB -= (a1*log(b1) - mylgamma(a1));
			LB -= (a2*(log(b2)-1.0) - mylgamma(a2));
			LB -= H[0].deltaRate*expDelta2;
			LB += ((double)(P-N+F)*0.5*log2pi);
			if(H[0].deltaShape>0.0)		LB-=mylgamma(H[0].deltaShape);				/* When informative prior is used */
			if(H[0].deltaRate>0.0)		LB+=H[0].deltaShape*log(H[0].deltaRate);
			if(H[0].etaShape>0.0)		LB-=(double)P*mylgamma(H[0].etaShape);		/* When informative prior is used */
			if(H[0].etaRate>0.0)		LB+=(double)P*(H[0].etaShape*log(H[0].etaRate));
			LB += ((double)(2*P+F)*0.5 - (double)P*log(2.0));
			break;
		case 3:	/* wBSR */
			vj=1.0+H[0].v;
			vS2=H[0].v*H[0].S2;
			Pi2=(1.0-H[0].Pi);
			for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
			if(H[0].Pi<1.0) /* BayesB */
			{
				for(locus=0;locus<P;locus++)
				{
					LB -= 0.5*log(vj*X[locus].S2*0.5)*vj;
					LB += 0.5*log(X[locus].varEffect);
					if(X[locus].expGamma>Gammalowerbound&&X[locus].expGamma<Gammaupperbound)
					{
						LB += X[locus].expGamma*log(H[0].Pi/X[locus].expGamma);
						LB += (1.0-X[locus].expGamma)*log(Pi2/(1.0-X[locus].expGamma));
					}
					else
					{
						if(X[locus].expGamma>Gammalowerbound) {LB += log(H[0].Pi); /* expGamma==1.0 */} else {LB += log(Pi2);/* expGamma==0.0 */}
					}
				}
			}
			else
			{	/* BayesA */
				for(locus=0;locus<P;locus++)
				{
					LB -= 0.5*log(vj*X[locus].S2*0.5)*vj;
					LB += 0.5*log(X[locus].varEffect);
				}
			}
			LB -= (a1*log(b1) - mylgamma(a1));
			LB += (double)P * (0.5*H[0].v*log(0.5*vS2)-mylgamma(0.5*H[0].v)+mylgamma(0.5*vj));
			LB -= ((double)(N-F)*0.5*log2pi - (double)(F+P)*0.5);		
			break;
		case 4:	/* BayesC */
			if(H[0].Pi<1.0) /* BayesC */
			{
				vS2=H[0].v*H[0].S2;
				Pi2=(1.0-H[0].Pi);
				for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
				for(locus=0,temp=0.0;locus<P;locus++)
				{
					LB += 0.5*X[locus].expGamma*log(X[locus].varEffect);
					if(X[locus].expGamma>Gammalowerbound&&X[locus].expGamma<Gammaupperbound)
					{
						LB += X[locus].expGamma*log(H[0].Pi/X[locus].expGamma);
						LB += (1.0-X[locus].expGamma)*log(Pi2/(1.0-X[locus].expGamma));
					}
					else
					{	
						if(X[locus].expGamma>Gammalowerbound) {LB += log(H[0].Pi); /* expGamma==1.0 */} else {LB += log(Pi2);/* expGamma==0.0 */}
					}
					temp+=X[locus].expGamma;
				}
				vj=H[0].v+temp;
				LB -= (0.5*log(vj*X[0].S2*0.5)*vj - mylgamma(0.5*vj));
				LB += (0.5*H[0].v*log(0.5*vS2) - mylgamma(0.5*H[0].v));
				LB -= (a1*log(b1) - mylgamma(a1));
				LB -= ((double)(N-F)*0.5*log2pi - ((double)F+temp)*0.5);
			}
			else
			{	/* BRR */
				vj=(double)P+H[0].v;
				vS2=H[0].v*H[0].S2;
				for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
				for(locus=0;locus<P;locus++){LB += 0.5*log(X[locus].varEffect);}
				LB -= (0.5*log(vj*X[0].S2*0.5)*vj - mylgamma(0.5*vj));
				LB += (0.5*H[0].v*log(0.5*vS2) - mylgamma(0.5*H[0].v));
				LB -= (a1*log(b1) - mylgamma(a1));			
				LB -= ((double)(N-F)*0.5*log2pi - (double)(F+P)*0.5);
			}
			break;
		case 5:	/* SSVS */
			vj=(double)P+H[0].v;
			vS2=H[0].v*H[0].S2;
			Pi2=(1.0-H[0].Pi);
			for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
			for(locus=0,temp=0.0;locus<P;locus++)
			{
				LB += 0.5*log(X[locus].varEffect);
				if(X[locus].expGamma>Gammalowerbound&&X[locus].expGamma<Gammaupperbound)
				{
					LB += X[locus].expGamma*log(H[0].Pi/X[locus].expGamma);
					LB += (1.0-X[locus].expGamma)*log(Pi2/(1.0-X[locus].expGamma));
				}
				else
				{
					if(X[locus].expGamma>Gammalowerbound) {LB += log(H[0].Pi);} else {LB += log(Pi2);}
				}
				temp+=X[locus].expGamma;
			}
			LB -= 0.5*((double)P-temp)*log(H[0].c);
			LB -= (0.5*log(vj*X[0].S2*0.5)*vj - mylgamma(0.5*vj));
			LB += (0.5*H[0].v*log(0.5*vS2) - mylgamma(0.5*H[0].v));
			LB -= (a1*log(b1) - mylgamma(a1));			
			LB -= ((double)(N-F)*0.5*log2pi - (double)(F+P)*0.5);
			break;
		case 6: /* MIX */
			vS2=H[0].v*H[0].S2;
			Pi2=(1.0-H[0].Pi);
			for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
			for(locus=0, vj=H[0].v, v2=H[0].v+(double)P;locus<P;locus++)
			{
				LB += 0.5*log(X[locus].varEffect);
				if(X[locus].expGamma>Gammalowerbound&&X[locus].expGamma<Gammaupperbound)
				{
					LB += X[locus].expGamma*log(H[0].Pi/X[locus].expGamma);
					LB += (1.0-X[locus].expGamma)*log(Pi2/(1.0-X[locus].expGamma));
				}
				else
				{
					if(X[locus].expGamma>Gammalowerbound) {LB += log(H[0].Pi); /* expGamma==1.0 */} else {LB += log(Pi2);/* expGamma==0.0 */}
				}
				vj += X[locus].expGamma;
				v2 -= X[locus].expGamma;
			}
			LB -= (0.5*log(vj*X[0].S2*0.5)*vj - mylgamma(0.5*vj));
			LB -= (0.5*log(v2*X[1].S2*0.5)*v2 - mylgamma(0.5*v2));
			LB += (H[0].v*log(0.5*vS2) - 2.0*mylgamma(0.5*H[0].v) + 0.5*H[0].v*log(H[0].c));
			LB -= (a1*log(b1) - mylgamma(a1));			
			LB -= ((double)(N-F)*0.5*log2pi - (double)(F+P)*0.5);
			break;
		case 7:	/* BayesB */
			if(H[0].Pi<1.0) /* BayesB */
			{
				vS2=H[0].v*H[0].S2;
				Pi2=(1.0-H[0].Pi);
				for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
				for(locus=0,temp=0.0;locus<P;locus++)
				{
					LB += 0.5*X[locus].expGamma*log(X[locus].varEffect);
					if(X[locus].expGamma>Gammalowerbound&&X[locus].expGamma<Gammaupperbound)
					{
						LB += X[locus].expGamma*log(H[0].Pi/X[locus].expGamma);
						LB += (1.0-X[locus].expGamma)*log(Pi2/(1.0-X[locus].expGamma));
					}
					else
					{	
						if(X[locus].expGamma>Gammalowerbound) {LB += log(H[0].Pi); /* expGamma==1.0 */} else {LB += log(Pi2);/* expGamma==0.0 */}
					}
					vj=H[0].v+X[locus].expGamma;
					LB -= (0.5*log(vj*X[locus].S2*0.5)*vj - mylgamma(0.5*vj));
					temp+=X[locus].expGamma;
				}			
				LB += (double)P*(0.5*H[0].v*log(0.5*vS2) - mylgamma(0.5*H[0].v));
				LB -= (a1*log(b1) - mylgamma(a1));
				LB -= ((double)(N-F)*0.5*log2pi - ((double)F+temp)*0.5);
			}
			else
			{	/* BayesA */
				vS2=H[0].v*H[0].S2;
				vj=1.0+H[0].v;
				for(locus=0;locus<F;locus++){LB += 0.5*log(Q[locus].varEffect);}
				for(locus=0;locus<P;locus++){
					LB += 0.5*log(X[locus].varEffect);
					LB -= (0.5*log(vj*X[locus].S2*0.5)*vj);
				}
				LB += (double)P*(0.5*H[0].v*log(0.5*vS2)-mylgamma(0.5*H[0].v)+mylgamma(0.5*vj));
				LB -= (a1*log(b1) - mylgamma(a1));			
				LB -= ((double)(N-F)*0.5*log2pi - (double)(F+P)*0.5);
			}
			break;
	}
	return(LB);
}

/*---Judge significance----------------------------------------------------------------------------------------------------------------------------------------*/
int	JudgeSig (double abseffect, double sd)
{
	int	Sig;
	if(abseffect>2.575829*sd){Sig = 1;}else{ if(abseffect>1.959964*sd){ Sig=5;} else {Sig=0;}}
	return(Sig);
}

int	JudgeSigPermutation (double effect, double *PermutationThreshold)
{
	int	Sig;

	/* These scripts should be active when PermutationTest is used for permutation*/
	//if(effect>PermutationThreshold[0]){Sig = 1;}else{ if(effect>PermutationThreshold[1]){ Sig=5;}else{ Sig=0;}}
	//if(effect<PermutationThreshold[3]){Sig = 1;}else{ if(effect<PermutationThreshold[2]){ Sig=5;}else{ Sig=0;}}

	if(fabs(effect)>PermutationThreshold[0]){Sig = 1;}else{ if(fabs(effect)>PermutationThreshold[1]){ Sig=5;}else{ Sig=0;}}

	return(Sig);
}

int	BL_VB (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* for update of Re2 */
	double	sumVarB;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* for update of Delta2 */
	double	sumEta2InTau2;

	/* Used when CondResidual == 1 */
	double	sumTau2B2;

	/* check convergence */
	double	Check1, Check2;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	expDelta2[0] = 1.0;
	if(RandomIni)
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = RNormal(0.0,sqrt(1.0/X[locus].expTau2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			sumVarB += (Q[target].x2*temp2);/* used when VB */

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		if(CondResidual) { temp3 = Tau0[0]; } else { temp3 = 1.0;}
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];

			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + X[target].expTau2 * temp3);

			prop = temp * temp2;
			prop2 = prop * prop + temp2;

			for(record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));
			sumVarB += (X[target].x2*temp2); /* used when VB */

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expEffect  = prop;
			X[target].exp2Effect = prop2;/* used when VB */
			X[target].varEffect  = temp2;/* used when VB */
		}

		/* update of Tau2 or Sigma2*/
		if (CondResidual) {temp = Tau0[0];} else {temp = 1.0;}	
		for(locus=0, sumEta2InTau2=0.0, sumTau2B2=0.0; locus<P; locus++)
		{
			target = Order[locus];

			prop = sqrt( expDelta2[0] * X[target].expEta2 / (X[target].exp2Effect * temp) ); 
					
			Check1 += pow((prop - X[target].expTau2), 2.0);
			Check2 += pow(prop, 2.0);

			if(CondResidual)
				sumTau2B2 += X[target].exp2Effect * prop;
					
			X[target].expTau2 = prop;

			prop = 1.0/X[target].expTau2 + 1.0/(expDelta2[0] * X[target].expEta2);
			sumEta2InTau2 += prop * X[target].expEta2;

			X[target].expInTau2 = prop;
		}

		/* update of Delta2 (for BL and EBL)*/
		a2=(double)P + H[0].deltaShape;
		b2=0.5 * sumEta2InTau2 + H[0].deltaRate;
		prop = a2/b2;
		Check1 += pow((prop-expDelta2[0]), 2.0);
		Check2 += pow(prop, 2.0);
		expDelta2[0] = prop;

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);
	
		if (CondResidual==1) 
		{
			a1=(double) (P + Nn)*0.5;
			b1=0.5*(temp + sumVarB + sumTau2B2);
		}
		else
		{
			a1=(double)Nn*0.5;
			b1=0.5*(temp + sumVarB);
		}
		prop = a1/b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			LBmonitor[ite-1] = LowerBound(Methodcode,Nn,P,X,F,Q,a1,b1,a2,b2,expDelta2[0],H,CondResidual);
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f LB: %f ", ite, Rmonitor[ite-1], log10(temp), LBmonitor[ite-1]); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

int	BL_EM (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* for update of Delta2 */
	double	sumEta2InTau2;

	/* Used when CondResidual == 1 */
	double	sumTau2B2;

	/* check convergence */
	double	Check1, Check2;

	/* Limits of estimates */
	double	Lowesteffect=1.0e-150;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	expDelta2[0] = 1.0;
	if(RandomIni)
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = RNormal(0.0,sqrt(1.0/X[locus].expTau2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
		}

		/* update of B */
		Shuffle (Order, P, Times);
		if(CondResidual) { temp3 = Tau0[0]; } else { temp3 = 1.0;}
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];

			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + X[target].expTau2 * temp3);

			prop = temp * temp2; if(fabs(prop)<Lowesteffect) {prop=Lowesteffect;}
			prop2 = prop * prop + temp2;

			for(record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expEffect  = prop;
		}

		/* update of Tau2 or Sigma2*/
		if (CondResidual) {temp = Tau0[0];} else {temp = 1.0;}
		for(locus=0, sumEta2InTau2=0.0, sumTau2B2=0.0; locus<P; locus++)
		{
			target = Order[locus];

			prop = sqrt( expDelta2[0] * X[target].expEta2 / temp)/ fabs(X[target].expEffect);

			Check1 += pow((prop - X[target].expTau2), 2.0);
			Check2 += pow(prop, 2.0);

			if(CondResidual) 
				sumTau2B2 += pow(X[target].expEffect, 2.0) * prop;
						
			X[target].expTau2 = prop;

			prop = 1.0/X[target].expTau2 + 1.0/(expDelta2[0] * X[target].expEta2);
			sumEta2InTau2 += prop * X[target].expEta2;

			X[target].expInTau2 = prop;
		}

		/* update of Delta2 (for BL and EBL)*/
		a2=(double)P + H[0].deltaShape;
		b2=0.5 * sumEta2InTau2 + H[0].deltaRate;
		prop = a2/b2;

		Check1 += pow((prop-expDelta2[0]), 2.0);
		Check2 += pow(prop, 2.0);
		expDelta2[0] = prop;

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		if (CondResidual==1) { prop = (double) (P + Nn - 2)/ (temp + sumTau2B2); }else{ prop = (double)(Nn - 2) / temp;}

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f ", ite, Rmonitor[ite-1], log10(temp)); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}


int	EBL_VB (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* for update of Re2 */
	double	sumVarB;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* for update of Delta2 */
	double	sumEta2InTau2;

	/* Used when CondResidual == 1 */
	double	sumTau2B2;

	/* check convergence */
	double	Check1, Check2;

	/* Limits of estimates */
	double	Lowesteffect=1.0e-150;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	expDelta2[0] = 1.0;
	if(RandomIni)
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = RNormal(0.0,sqrt(1.0/X[locus].expTau2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			sumVarB += (Q[target].x2*temp2);/* used when VB */

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		if(CondResidual) { temp3 = Tau0[0]; } else { temp3 = 1.0;}
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];

			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + X[target].expTau2 * temp3);

			prop = temp * temp2;
			prop2 = prop * prop + temp2;

			for(record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));
			sumVarB += (X[target].x2*temp2); /* used when VB */

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expEffect  = prop;
			X[target].exp2Effect = prop2;/* used when VB */
			X[target].varEffect  = temp2;/* used when VB */
		}


		/* update of Tau2 or Sigma2*/
		if (CondResidual) {temp = Tau0[0];} else {temp = 1.0;}
		for(locus=0, sumEta2InTau2=0.0, sumTau2B2=0.0; locus<P; locus++)
		{
			target = Order[locus];
			prop = sqrt( expDelta2[0] * X[target].expEta2 / (X[target].exp2Effect * temp) ); 						

			Check1 += pow((prop - X[target].expTau2), 2.0);
			Check2 += pow(prop, 2.0);

			if(CondResidual) 
				sumTau2B2 += X[target].exp2Effect * prop;
						
			X[target].expTau2 = prop;

			prop = 1.0/X[target].expTau2 + 1.0/(expDelta2[0] * X[target].expEta2);
			sumEta2InTau2 += prop * X[target].expEta2;

			X[target].expInTau2 = prop;
		}

		/* update of Delta2 (for BL and EBL)*/
		a2=(double)P + H[0].deltaShape;
		b2=0.5 * sumEta2InTau2 + H[0].deltaRate;
		prop = a2/b2;

		Check1 += pow((prop-expDelta2[0]), 2.0);
		Check2 += pow(prop, 2.0);
		expDelta2[0] = prop;

		/* update of Eta2 (for EBL)*/
		Shuffle (Order, P, Times);
		for(locus=0; locus<P; locus++)
		{
			target = Order[locus];
			X[target].a3=1.0 + H[0].etaShape;
			X[target].b3=0.5 * expDelta2[0] * X[target].expInTau2 + H[0].etaRate;
			prop = X[target].a3/X[target].b3;

			Check1 += pow((prop - X[target].expEta2), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEta2 = prop;
		}

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		if (CondResidual==1)
		{
			a1=(double) (P + Nn)*0.5;
			b1=0.5*(temp + sumVarB + sumTau2B2);
		}
		else
		{
			a1=(double)Nn*0.5;
			b1=0.5*(temp + sumVarB);
		}
		prop = a1/b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			LBmonitor[ite-1] = LowerBound(Methodcode,Nn,P,X,F,Q,a1,b1,a2,b2,expDelta2[0],H,CondResidual);
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f LB: %f ", ite, Rmonitor[ite-1], log10(temp), LBmonitor[ite-1]); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

int	EBL_EM (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* for update of Re2 */
	double	sumVarB;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* for update of Delta2 */
	double	sumEta2InTau2;

	/* Used when CondResidual == 1 */
	double	sumTau2B2;

	/* check convergence */
	double	Check1, Check2;

	/* Limits of estimates */
	double	Lowesteffect=1.0e-150;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}

	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;

	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	expDelta2[0] = 1.0;
	if(RandomIni)
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = RNormal(0.0,sqrt(1.0/X[locus].expTau2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expTau2   = (double)P/VarTemp;
			X[locus].expInTau2 = VarTemp/(double)P;
			X[locus].expEta2   = 1.0; /* when expEta2 for all markers are always fixed to 1.0, EBL becomes BL */
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
		}

		/* update of B */
		Shuffle (Order, P, Times);
		if(CondResidual) { temp3 = Tau0[0]; } else { temp3 = 1.0;}
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];

			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + X[target].expTau2 * temp3);

			prop = temp * temp2;	if(fabs(prop)<Lowesteffect) {prop=Lowesteffect;}
			prop2 = prop * prop + temp2;

			for(record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expEffect  = prop;
		}

		/* update of Tau2 or Sigma2*/
		if (CondResidual) {temp = Tau0[0];} else {temp = 1.0;}
		for(locus=0, sumEta2InTau2=0.0, sumTau2B2=0.0; locus<P; locus++)
		{
			target = Order[locus];
			prop = sqrt( expDelta2[0] * X[target].expEta2 / temp)/ fabs(X[target].expEffect);

			Check1 += pow((prop - X[target].expTau2), 2.0);
			Check2 += pow(prop, 2.0);

			if(CondResidual) 
				sumTau2B2 += pow(X[target].expEffect, 2.0) * prop;
						
			X[target].expTau2 = prop;

			prop = 1.0/X[target].expTau2 + 1.0/(expDelta2[0] * X[target].expEta2);
			sumEta2InTau2 += prop * X[target].expEta2;

			X[target].expInTau2 = prop;
		}

		/* update of Delta2 (for BL and EBL)*/
		a2=(double)P + H[0].deltaShape;
		b2=0.5 * sumEta2InTau2 + H[0].deltaRate;
		prop = a2/b2;

		Check1 += pow((prop-expDelta2[0]), 2.0);
		Check2 += pow(prop, 2.0);
		expDelta2[0] = prop;

		/* update of Eta2 (for EBL)*/
		Shuffle (Order, P, Times);
		for(locus=0; locus<P; locus++)
		{
			target = Order[locus];
			X[target].a3=1.0 + H[0].etaShape;
			X[target].b3=0.5 * expDelta2[0] * X[target].expInTau2 + H[0].etaRate;
			prop = X[target].a3/X[target].b3;

			Check1 += pow((prop - X[target].expEta2), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEta2 = prop;
		}

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		if (CondResidual==1) { prop = (double) (P + Nn - 2)/ (temp + sumTau2B2); }else{ prop = (double)(Nn - 2) / temp;}

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f ", ite, Rmonitor[ite-1], log10(temp));
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

/*---Perform genome-wide regression----------------------------------------------------------------------------------------------------------------------------*/
int	wBSR_VB (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* for update of Re2 */
	double	sumVarB;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi, log1minusPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* check convergence */
	double	Check1, Check2;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	if (H[0].Pi<1.0){logPi = log(H[0].Pi); log1minusPi = log(1.0-H[0].Pi);} /* when BayesB  */

	if(H[0].Pi<1.0)
	{/* BayesB */
		if(RandomIni)
		{
			for(locus=0; locus<P; locus++)
			{
				X[locus].expSigma2 = VarTemp/((double)P*H[0].Pi);	X[locus].S2 = (double)P*H[0].Pi/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
				X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp/((double)P*H[0].Pi);	X[locus].S2 = (double)P*H[0].Pi/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
			}
		}
	}
	else
	{/* BayesA */
		if(RandomIni)
		{
			for (locus=0; locus<P; locus++)
			{	
				X[locus].expSigma2 = VarTemp/(double)P;	X[locus].S2 = (double)P/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
				X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp/(double)P;	X[locus].S2 = (double)P/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
			}
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect * X[locus].expGamma;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			sumVarB += (Q[target].x2*temp2);/* used when VB */

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect * X[target].expGamma);
			temp*=X[target].expGamma;
			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] * X[target].exp2Gamma + 1.0/X[target].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;

			for(record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * X[target].expGamma * (X[target].expEffect - prop));
			if((int)H[0].Pi==1) /* when BayesA */
				sumVarB += (X[target].x2 * temp2); /* used when VB */

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
			X[target].exp2Effect = prop2;	/* used when VB */
			X[target].varEffect = temp2;	/* used when VB */
		}

		/* update of Tau2 or Sigma2*/
		for(locus=0; locus<P; locus++)
		{
			target = Order[locus];
			temp = X[target].exp2Effect + vS2;
						
			X[target].expSigma2 = temp/(H[0].v - 1.0);
			prop = temp/(H[0].v + 1.0);

			Check1 += pow((prop - X[target].S2), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].S2 = prop;
		}

		if(H[0].Pi<1.0)
		{ /*when BayesB */
			Shuffle (Order, P, Times);
			for(locus=0; locus<P; locus++)
			{
				target = Order[locus];

				for (record=0, ProbInclude=0.0, ProbExclude=0.0; record<Nn; record++)
				{
					ProbInclude += pow((Y[0].expErrors [Y[0].use[record]] + (X[target].expGamma - 1.0) * X[target].expEffect * X[target].covariates[Y[0].use[record]]), 2.0);
					ProbExclude += pow((Y[0].expErrors [Y[0].use[record]] + X[target].expGamma * X[target].expEffect * X[target].covariates[Y[0].use[record]]), 2.0);
				}

				ProbInclude += X[target].x2 * X[target].varEffect; /* when VB */
				ProbInclude *= -0.5 * Tau0[0];	ProbExclude *= -0.5 * Tau0[0];
				ProbInclude += logPi;				ProbExclude += log1minusPi;

				temp = ProbInclude;
				if(temp < ProbExclude) temp = ProbExclude;
				ProbInclude -= temp;	ProbExclude -= temp;

				ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
				prop = ProbInclude/(ProbInclude + ProbExclude);
				prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

				/* update sumVarB and expErrors */
				sumVarB += (X[target].x2 * (prop2*X[target].exp2Effect - pow(prop*X[target].expEffect, 2.0))); /*used when VB */
				for(record=0; record<Nn; record++)
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * X[target].expEffect * (X[target].expGamma - prop));

				/* for convergence check */
				Check1 += pow((prop - X[target].expGamma), 2.0);
				Check2 += pow(prop, 2.0);

				X[target].expGamma = prop;
				X[target].exp2Gamma = prop2; /* used when variational Bayes */
			}
		}

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		a1=(double)Nn*0.5;
		b1=0.5*(temp + sumVarB);

		prop = a1/b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			LBmonitor[ite-1] = LowerBound(Methodcode,Nn,P,X,F,Q,a1,b1,a2,b2,expDelta2[0],H,CondResidual);
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f LB: %f ", ite, Rmonitor[ite-1], log10(temp), LBmonitor[ite-1]); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}



/*---Perform genome-wide regression----------------------------------------------------------------------------------------------------------------------------*/
int	wBSR_EM (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi, log1minusPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* check convergence */
	double	Check1, Check2;

	/* Limits of estimates */
	double	Lowesteffect=1.0e-150;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	if (H[0].Pi<1.0){logPi = log(H[0].Pi); log1minusPi = log(1.0-H[0].Pi);} /* when BayesB  */

	if(H[0].Pi<1.0)
	{/* BayesB */
		if(RandomIni)
		{
			for(locus=0; locus<P; locus++)
			{
				X[locus].expSigma2 = VarTemp/((double)P*H[0].Pi);	X[locus].S2 = (double)P*H[0].Pi/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
				X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp/((double)P*H[0].Pi);	X[locus].S2 = (double)P*H[0].Pi/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
			}
		}
	}
	else
	{/* BayesA */
		if(RandomIni)
		{
			for (locus=0; locus<P; locus++)
			{	
				X[locus].expSigma2 = VarTemp/(double)P;	X[locus].S2 = (double)P/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
				X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp/(double)P;	X[locus].S2 = (double)P/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
			}
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect * X[locus].expGamma;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect * X[target].expGamma);
			temp*=X[target].expGamma;
			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] * pow(X[target].expGamma,2.0) + 1.0/X[target].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;

			for(record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * X[target].expGamma * (X[target].expEffect - prop));

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
		}


		/* update of Tau2 or Sigma2*/
		for(locus=0; locus<P; locus++)
		{
			target = Order[locus];
			temp = pow(X[target].expEffect, 2.0) + vS2;
						
			X[target].expSigma2 = temp/(H[0].v - 1.0);
			prop = temp/(H[0].v + 1.0);

			Check1 += pow((prop - X[target].S2), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].S2 = prop;
		}

		/* Update of Gamma (for wBSR, SSVS, and MIX) */
		if(H[0].Pi<1.0)
		{ /*when BayesB */
			Shuffle (Order, P, Times);
			for(locus=0; locus<P; locus++)
			{
				target = Order[locus];

				for (record=0, ProbInclude=0.0, ProbExclude=0.0; record<Nn; record++)
				{
					ProbInclude += pow((Y[0].expErrors [Y[0].use[record]] + (X[target].expGamma - 1.0) * X[target].expEffect * X[target].covariates[Y[0].use[record]]), 2.0);
					ProbExclude += pow((Y[0].expErrors [Y[0].use[record]] + X[target].expGamma * X[target].expEffect * X[target].covariates[Y[0].use[record]]), 2.0);
				}

				ProbInclude *= -0.5 * Tau0[0];	ProbExclude *= -0.5 * Tau0[0];
				ProbInclude += logPi;				ProbExclude += log1minusPi;

				temp = ProbInclude;
				if(temp < ProbExclude) temp = ProbExclude;
				ProbInclude -= temp;	ProbExclude -= temp;

				ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
				prop = ProbInclude/(ProbInclude + ProbExclude);
				prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

				/* update sumVarB and expErrors */
				for(record=0; record<Nn; record++)
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * X[target].expEffect * (X[target].expGamma - prop));

				/* for convergence check */
				Check1 += pow((prop - X[target].expGamma), 2.0);
				Check2 += pow(prop, 2.0);

				X[target].expGamma = prop;
			}
		}

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		prop = (double)(Nn - 2) / temp;
		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f ", ite, Rmonitor[ite-1], log10(temp)); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

int	BayesB_VB (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* for update of Re2 */
	double	sumVarB;

	/* For variable selection */
	double	logPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* check convergence */
	double	Check1, Check2;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	if(H[0].Pi<1.0) /* BayesB */{logPi = log(H[0].Pi);}

	if(H[0].Pi<1.0)
	{/* BayesB */
		if(RandomIni)
		{
			for(locus=0; locus<P; locus++)
			{
				X[locus].expSigma2 = VarTemp/((double)P*H[0].Pi);	X[locus].S2 = (double)P*H[0].Pi/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
				X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp/((double)P*H[0].Pi);	X[locus].S2 = (double)P*H[0].Pi/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
			}
		}
	}
	else
	{/* BayesA */
		if(RandomIni)
		{
			for (locus=0; locus<P; locus++)
			{	
				X[locus].expSigma2 = VarTemp/(double)P;	X[locus].S2 = (double)P/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
				X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp/(double)P;	X[locus].S2 = (double)P/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
			}
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect * X[locus].expGamma;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			sumVarB += (Q[target].x2*temp2);/* used when VB */

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect * X[target].expGamma);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + 1.0/X[target].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;

			if(H[0].Pi<1.0)
			{	/* BayesB */
				/* update Gamma */
				ProbIncludeConstant= 0.5 * Digamma (0.5*(H[0].v+X[target].expGamma)) - 0.5*log(0.5*(X[target].expGamma*X[target].exp2Effect+vS2)) + logPi;
				ProbInclude = 0.5*temp2*temp*temp+0.5*log(temp2);
				ProbInclude += ProbIncludeConstant;
				if(ProbInclude>20.0) ProbInclude=20.0; /* to avoid overflow */
				ProbInclude = exp(ProbInclude);
				temp3 = ProbInclude/(ProbInclude + 1.0 - H[0].Pi);

				/* update residuals */
				for(record=0; record<Nn; record++)
				{
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * X[target].expGamma * X[target].expEffect);
					Y[0].expErrors[Y[0].use[record]] -= (X[target].covariates[Y[0].use[record]] * temp3 * prop);
				}
				X[target].expGamma = temp3;
				X[target].exp2Gamma = pow(X[target].expGamma, 2.0) + X[target].expGamma * (1.0 - X[target].expGamma);
				sumVarB += (X[target].x2*X[target].expGamma*(prop2-X[target].expGamma*prop*prop)); /* used when VB */
			}
			else
			{	/* BayesA */
				for(record=0; record<Nn; record++)
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));
				sumVarB += (X[target].x2 * temp2); /* used when VB */
			}

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
			X[target].exp2Effect = prop2;	/* used when VB */
			X[target].varEffect = temp2;	/* used when VB */
		}

		/* update of Tau2 or Sigma2*/
		if(H[0].Pi<1.0)
		{	/* BayesB */
			for(locus=0; locus<P; locus++)
			{
				target=Order[locus];
				temp = X[target].expGamma*X[target].exp2Effect + vS2;
				X[target].expSigma2 = temp/(H[0].v + X[target].expGamma - 2.0);
				prop = temp/(H[0].v + X[target].expGamma);
				Check1 += pow((prop - X[target].S2), 2.0);
				Check2 += pow(prop, 2.0);
				X[target].S2 = prop;
			}
		}
		else
		{	/* BayesA */
			for(locus=0; locus<P; locus++)
			{
				target=Order[locus];
				temp = X[target].exp2Effect + vS2;
				X[target].expSigma2 = temp/(H[0].v - 1.0);
				prop = temp/(H[0].v + 1.0);
				Check1 += pow((prop - X[target].S2), 2.0);
				Check2 += pow(prop, 2.0);
				X[target].S2 = prop;
			}
		}


		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		a1=(double)Nn*0.5;
		b1=0.5*(temp + sumVarB);

		prop = a1/b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			LBmonitor[ite-1] = LowerBound(Methodcode,Nn,P,X,F,Q,a1,b1,a2,b2,expDelta2[0],H,CondResidual);
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f LB: %f ", ite, Rmonitor[ite-1], log10(temp), LBmonitor[ite-1]); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

int	BayesC_VB (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* for update of Re2 */
	double	sumVarB;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi, log1minusPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* Sum of marker effects. Used when BayesC, SSVS, and MIX */
	double	sumGammaB2[2];

	/* Used when MIX and BayesC (sumGamma, sum of gamma), and when Mix (vsS2, product of hyperparameters). */
	double	sumGamma[2], vcS2;

	/* check convergence */
	double	Check1, Check2;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	if(H[0].Pi<1.0) /* BayesC */
	{
		logPi = log(H[0].Pi); log1minusPi = log(1.0-H[0].Pi);
		sumGamma[0]=0.0; sumGamma[1]=0.0;
		sumGammaB2[0]=0.0; sumGammaB2[1]=0.0;
	}
	else
	{	/* BRR */
		sumGammaB2[0]=0.0; sumGammaB2[1]=0.0;
	}

	if(H[0].Pi<1.0)
	{/* BayesC */
		if(RandomIni)
		{
			for(locus=0; locus<P; locus++)
			{
				X[locus].expSigma2 = VarTemp;	X[locus].S2 = 1.0/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
				X[locus].expEffect = RNormal(0.0,sqrt(VarTemp/((double)P*H[0].Pi)));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				sumGamma[0] += X[locus].expGamma;
				sumGammaB2[0] += X[locus].expGamma * X[locus].exp2Effect;			
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp;	X[locus].S2 = 1.0/VarTemp;
				X[locus].expGamma = 0.5;
				X[locus].exp2Gamma = pow(X[locus].expGamma, 2.0) + X[locus].expGamma * (1.0-X[locus].expGamma);
				sumGamma[0] += X[locus].expGamma;
				sumGammaB2[0] += X[locus].expGamma * X[locus].exp2Effect;		
			}
		}
	}
	else
	{/* BRR */
		if(RandomIni)
		{
			for (locus=0; locus<P; locus++)
			{	
				X[locus].expSigma2 = VarTemp;	X[locus].S2 = 1.0/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
				X[locus].expEffect = RNormal(0.0,sqrt(VarTemp/(double)P));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
				sumGammaB2[0] += X[locus].expGamma * X[locus].exp2Effect;			
			}
		}else{
			for (locus=0; locus<P; locus++)
			{
				X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
				X[locus].expSigma2 = VarTemp;	X[locus].S2 = 1.0/VarTemp;
				X[locus].expGamma = 1.0;
				X[locus].exp2Gamma = 1.0;
				sumGammaB2[0] += X[locus].expGamma * X[locus].exp2Effect;			
			}
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect * X[locus].expGamma;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			sumVarB += (Q[target].x2*temp2);/* used when VB */

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		if(H[0].Pi<1.0)
		{	/* BayesC */
			sumGammaB2[1]=0.0; sumGamma[1]=0.0;
			//ProbIncludeConstant=Digamma (0.5*(H[0].v+sumGamma[0])) - 0.5*log(0.5*(sumGammaB2[0]+vS2)) + logPi - log1minusPi;
			ProbIncludeConstant= 0.5 * Digamma (0.5*(H[0].v+sumGamma[0])) - 0.5*log(0.5*(sumGammaB2[0]+vS2)) + logPi;
		}
		else
		{	/* BRR */
			sumGammaB2[0]=0.0;
		}

		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect * X[target].expGamma);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + 1.0/X[0].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;

			if(H[0].Pi<1.0)
			{	/* BayesC */
				/* update Gamma */
				ProbInclude = 0.5*temp2*temp*temp+0.5*log(temp2);
				ProbInclude += ProbIncludeConstant;
				if(ProbInclude>20.0) ProbInclude=20.0; /* to avoid overflow */
				ProbInclude = exp(ProbInclude);	
				//ProbExclude = 1.0;
				//temp3 = ProbInclude/(ProbInclude + ProbExclude);
				temp3 = ProbInclude/(ProbInclude + 1.0 - H[0].Pi);

				/* update residuals */
				for(record=0; record<Nn; record++)
				{
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * X[target].expGamma * X[target].expEffect);
					Y[0].expErrors[Y[0].use[record]] -= (X[target].covariates[Y[0].use[record]] * temp3 * prop);
				}
				X[target].expGamma = temp3;
				X[target].exp2Gamma = pow(X[target].expGamma, 2.0) + X[target].expGamma * (1.0 - X[target].expGamma);
				sumVarB += (X[target].x2*X[target].expGamma*(prop2-X[target].expGamma*prop*prop)); /* used when VB */
				sumGammaB2[1] += prop2 * X[target].expGamma;
				sumGamma[1] += X[target].expGamma;
				/*---- Note -----------------------------------------------------------------------------------------------------
				sumGamma[0] and sumGammaB2[0] are not updated here, because these are the posterior parameters of X[0].expSigma. 
				These values are replaced by sumGamma[1] and sumGammaB2[1] at the update of X[0].expSigma.
				---------------------------------------------------------------------------------------------------------------*/
			}
			else
			{	/* BRR */
				for(record=0; record<Nn; record++)
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));
				sumVarB += (X[target].x2 * temp2); /* used when VB */
				sumGammaB2[0] += prop2;
			}

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
			X[target].exp2Effect = prop2;	/* used when VB */
			X[target].varEffect = temp2;	/* used when VB */
		}

		/* update of Tau2 or Sigma2*/
		if(H[0].Pi<1.0)
		{	/* BayesC */
			sumGammaB2[0] = sumGammaB2[1];
			sumGamma[0] = sumGamma[1];
			temp = sumGammaB2[0] + vS2;
			X[0].expSigma2 = temp/(H[0].v + sumGamma[0] - 2.0);
			prop = temp/(H[0].v + sumGamma[0]);
		}
		else
		{	/* BRR */
			temp = sumGammaB2[0] + vS2;
			X[0].expSigma2 = temp/(H[0].v + (double) P - 2.0);
			prop = temp/(H[0].v + (double)P);
		}
		Check1 += pow((prop - X[0].S2), 2.0);
		Check2 += pow(prop, 2.0);
		X[0].S2 = prop;


		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		a1=(double)Nn*0.5;
		b1=0.5*(temp + sumVarB);

		prop = a1/b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			LBmonitor[ite-1] = LowerBound(Methodcode,Nn,P,X,F,Q,a1,b1,a2,b2,expDelta2[0],H,CondResidual);
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f LB: %f ", ite, Rmonitor[ite-1], log10(temp), LBmonitor[ite-1]); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

int	SSVS_VB (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* for update of Re2 */
	double	sumVarB;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi, log1minusPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* Sum of marker effects. Used when BayesC, SSVS, and MIX */
	double	sumGammaB2[2];

	/* Used when SSVS */
	double	invC, OneMinusInvC, logC;

	/* check convergence */
	double	Check1, Check2;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi); log1minusPi = log(1.0-H[0].Pi);
	invC = 1.0/H[0].c; OneMinusInvC = 1.0 - invC; logC = log(H[0].c);

	if(RandomIni)/* Although we assign initial values for all loci, only the values for the first and second loci are used*/
	{
		for (locus=0; locus<P; locus++)
		{			
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
			X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2/(double)P));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;

	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			sumVarB += (Q[target].x2*temp2);/* used when VB */

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		sumGammaB2[0]=0.0; sumGammaB2[1]=0.0;
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + (X[target].expGamma * OneMinusInvC + invC)/X[0].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;
			for(record=0; record<Nn; record++)
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));
			sumVarB += (X[target].x2 * temp2); /* used when VB */
			temp3 = prop2;
			sumGammaB2[0] += temp3 * X[target].expGamma;
			sumGammaB2[1] += temp3 * (1.0-X[target].expGamma);

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
			X[target].exp2Effect = prop2;	/* used when VB */
			X[target].varEffect = temp2;	/* used when VB */
		}

		/* update of Tau2 or Sigma2*/
		temp = sumGammaB2[0] + sumGammaB2[1]*invC + vS2;
		X[0].expSigma2 = temp/(H[0].v + (double) P - 2.0);
		prop = temp/(H[0].v + (double)P);

		Check1 += pow((prop - X[0].S2), 2.0);
		Check2 += pow(prop, 2.0);
		X[0].S2 = prop;

		/* Update of Gamma (for wBSR, SSVS, and MIX) */
		Shuffle (Order, P, Times);
		for(locus=0; locus<P; locus++)
		{
			target = Order[locus];
			temp4 = X[target].exp2Effect;

			ProbInclude = - 0.5 * temp4 / X[0].S2 + logPi;
			ProbExclude = - 0.5 * temp4 / X[0].S2 * invC + log1minusPi - 0.5 * logC;

			temp = ProbInclude;
			if(temp < ProbExclude) temp = ProbExclude;
			ProbInclude -= temp;	ProbExclude -= temp;

			ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
			prop = ProbInclude/(ProbInclude + ProbExclude);
			prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

			/* for convergence check */
			Check1 += pow((prop - X[target].expGamma), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expGamma = prop;	
			X[target].exp2Gamma = prop2;
		}

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);
		a1=(double)Nn*0.5;
		b1=0.5*(temp + sumVarB);
		prop = a1/b1;
		Tau0[1] = a1 / pow(b1, 2.0);
		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			LBmonitor[ite-1] = LowerBound(Methodcode,Nn,P,X,F,Q,a1,b1,a2,b2,expDelta2[0],H,CondResidual);
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f LB: %f ", ite, Rmonitor[ite-1], log10(temp), LBmonitor[ite-1]); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

/*---Perform genome-wide regression----------------------------------------------------------------------------------------------------------------------------*/
int	SSVS_EM (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi, log1minusPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* Sum of marker effects. Used when BayesC, SSVS, and MIX */
	double	sumGammaB2[2];

	/* Used when SSVS */
	double	invC, OneMinusInvC, logC;

	/* check convergence */
	double	Check1, Check2;

	/* Limits of estimates */
	double	Lowesteffect=1.0e-150;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}

	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi); log1minusPi = log(1.0-H[0].Pi);
	invC = 1.0/H[0].c; OneMinusInvC = 1.0 - invC; logC = log(H[0].c);

	if(RandomIni)/* Although we assign initial values for all loci, only the values for the first and second loci are used*/
	{
		for (locus=0; locus<P; locus++)
		{			
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
			X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2/(double)P));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
		}

		/* update of B */
		Shuffle (Order, P, Times);
		sumGammaB2[0]=0.0; sumGammaB2[1]=0.0;
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
					temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);
			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + (X[target].expGamma * OneMinusInvC + invC)/X[0].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;

			for(record=0; record<Nn; record++)
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));

			temp3 = pow(prop, 2.0);
			sumGammaB2[0] += temp3 * X[target].expGamma;
			sumGammaB2[1] += temp3 * (1.0-X[target].expGamma);

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
			X[target].exp2Effect = prop2;	/* used when VB */
			X[target].varEffect = temp2;	/* used when VB */
		}

		/* update of Tau2 or Sigma2*/
		temp = sumGammaB2[0] + sumGammaB2[1]*invC + vS2;
		X[0].expSigma2 = temp/(H[0].v + (double) P - 2.0);
		prop = temp/(H[0].v + (double)P);

		Check1 += pow((prop - X[0].S2), 2.0);
		Check2 += pow(prop, 2.0);
		X[0].S2 = prop;

		/* Update of Gamma (for wBSR, SSVS, and MIX) */
		Shuffle (Order, P, Times);
		for(locus=0; locus<P; locus++)
		{
			target = Order[locus];
			temp4 = pow(X[target].expEffect,2.0);

			ProbInclude = - 0.5 * temp4 / X[0].S2 + logPi;
			ProbExclude = - 0.5 * temp4 / X[0].S2 * invC + log1minusPi - 0.5 * logC;

			temp = ProbInclude;
			if(temp < ProbExclude) temp = ProbExclude;
			ProbInclude -= temp;	ProbExclude -= temp;

			ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
			prop = ProbInclude/(ProbInclude + ProbExclude);
			prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

			/* for convergence check */
			Check1 += pow((prop - X[target].expGamma), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expGamma = prop;	
			X[target].exp2Gamma = prop2;
		}

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		prop = (double)(Nn - 2) / temp;
		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f ", ite, Rmonitor[ite-1], log10(temp)); 
			break;

		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}

int	MIX_VB (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* for update of Re2 */
	double	sumVarB;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi, log1minusPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* Sum of marker effects. Used when BayesC, SSVS, and MIX */
	double	sumGammaB2[2];

	/* Used when MIX and BayesC (sumGamma, sum of gamma), and when Mix (vsS2, product of hyperparameters). */
	double	sumGamma[2], vcS2;

	/* check convergence */
	double	Check1, Check2;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi); log1minusPi = log(1.0-H[0].Pi);
	sumGamma[0]=0.0; sumGamma[1]=0.0;
	vcS2 = H[0].c * H[0].v * H[0].S2;

	if(RandomIni)/* Although we assign initial values for all loci, only the values for the first and second loci are used*/
	{
		for (locus=0; locus<P; locus++)
		{			
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
			X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2/(double)P));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
			sumGamma[0] += X[locus].expGamma;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
			sumGamma[0] += X[locus].expGamma;
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;
		
		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			sumVarB += (Q[target].x2*temp2);/* used when VB */

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
			Q[target].exp2Effect = prop2;/* used when VB */
			Q[target].varEffect  = temp2;/* used when VB */
		}

		/* update of B */
		Shuffle (Order, P, Times);
		sumGammaB2[0]=0.0; sumGammaB2[1]=0.0;
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);
			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + X[target].expGamma/X[0].S2 + (1.0 - X[target].expGamma)/X[1].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;

			for(record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));
			sumVarB += (X[target].x2 * temp2); /* used when VB */
			temp3 = prop2;

			sumGammaB2[0] += temp3 * X[target].expGamma;
			sumGammaB2[1] += temp3 * (1.0-X[target].expGamma);

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
			X[target].exp2Effect = prop2;	/* used when VB */
			X[target].varEffect = temp2;	/* used when VB */
		}


		/* update of Tau2 or Sigma2*/
		temp = sumGammaB2[0] + vS2;
		X[0].expSigma2 = temp/(H[0].v + sumGamma[0] - 2.0);
		prop = temp/(H[0].v + sumGamma[0]);

		Check1 += pow((prop - X[0].S2), 2.0);
		Check2 += pow(prop, 2.0);
		X[0].S2 = prop;

		temp = sumGammaB2[1] + vcS2;
		X[1].expSigma2 = temp/(H[0].v + (double)P - sumGamma[0] - 2.0);
		prop = temp/((double)P - sumGamma[0]);

		Check1 += pow((prop - X[1].S2), 2.0);
		Check2 += pow(prop, 2.0);
		X[1].S2 = prop;


		/* Update of Gamma */
		Shuffle (Order, P, Times);
		ProbIncludeConstant = Digamma (0.5*(H[0].v+sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[0]+vS2)) + logPi;
		ProbExcludeConstant = Digamma (0.5*(H[0].v+(double)P-sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[1]+vS2)) + log1minusPi;
		/*---- Note -----------------------------------------------------------------------------------------------------
			sumGamma[0] and sumGammaB2 are not updated here, because these are the posterior parameters of X[0].expSigma. 
			sumGamma[0] is replaced by sumGamma[1] after the update of all gamma.
			sumGammaB2 are updated when marker effects are updated.
		---------------------------------------------------------------------------------------------------------------*/
		for(locus=0, sumGamma[1]=0.0; locus<P; locus++)
		{
			target = Order[locus];
			temp4 = X[target].exp2Effect;

			ProbInclude = ProbIncludeConstant - 0.5 * temp4 / X[0].S2;
			ProbExclude = ProbExcludeConstant - 0.5 * temp4 / X[1].S2;

			temp = ProbInclude;
			if(temp < ProbExclude) temp = ProbExclude;
			ProbInclude -= temp;	ProbExclude -= temp;

			ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
			prop = ProbInclude/(ProbInclude + ProbExclude);
			prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

			/* for convergence check */
			Check1 += pow((prop - X[target].expGamma), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expGamma = prop;	
			X[target].exp2Gamma = prop2;
			sumGamma[1] += X[target].expGamma;
		}
		sumGamma[0] = sumGamma[1];


		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		a1=(double)Nn*0.5;
		b1=0.5*(temp + sumVarB);

		prop = a1/b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			LBmonitor[ite-1] = LowerBound(Methodcode,Nn,P,X,F,Q,a1,b1,a2,b2,expDelta2[0],H,CondResidual);
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f LB: %f ", ite, Rmonitor[ite-1], log10(temp), LBmonitor[ite-1]); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);



	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}


int	MIX_EM (int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo, int RandomIni)
{
	/* Number of non-missing record */
	int		Nn;

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		ii;
	int		record, locus;

	/* Update orders, Timse of shuffling orders, locus to be updated*/
	int		*Order, Times, target;

	/* Temporary objects */
	double	temp, temp2, temp3, temp4; 
		
	/* New values */
	double	prop, prop2;

	/* mean, sd and Var of observations */
	double	Mean, Sd, Var, VarTemp;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi, log1minusPi;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude, ProbExclude, ProbIncludeConstant, ProbExcludeConstant;

	/* Shape (a) ans rate (b) parameters of posterior gamma distributions */
	double	a1=0.0, b1=0.0, a2=0.0, b2=0.0;

	/* Sum of marker effects. Used when BayesC, SSVS, and MIX */
	double	sumGammaB2[2];

	/* Used when MIX and BayesC (sumGamma, sum of gamma), and when Mix (vsS2, product of hyperparameters). */
	double	sumGamma[2], vcS2;

	/* check convergence */
	double	Check1, Check2;

	/* Limits of estimates */
	double	Lowesteffect=1.0e-150;

	/* calculate computation time */
	clock_t start, end;

	/* Check missing records */
	for (record=0, Y[0].nmiss=0; record<N; record++)
		if (Y[0].observations[record] == Missingvalue)
			Y[0].nmiss ++; 
	
	Nn = N - Y[0].nmiss;

	Y[0].use = (int*) calloc (Nn, sizeof(int));
	for(record=0, ii=0; record<N; record++)
		if(Y[0].observations[record] != Missingvalue) { Y[0].use[ii] = record; ii++;}


	/* Standardize observations */
	for (record=0, Mean=0.0; record<Nn; record++)
		Mean += Y[0].observations[Y[0].use[record]];
	Mean /= (double) Nn;

	for (record=0, Var=0.0; record<Nn; record++)
		Var += pow( (Y[0].observations[Y[0].use[record]] - Mean), 2.0);
	Var /= (double) Nn;
	Sd = sqrt(Var);

	Y[0].stobs =(double*) calloc (N, sizeof(double));
	memcpy (Y[0].stobs, Y[0].observations, sizeof(double) * N);
	for (record=0; record<Nn; record++)
		Y[0].stobs[Y[0].use[record]] = (Y[0].observations[Y[0].use[record]] - Mean)/Sd;


	/* Calculate sum of squared covariates for non-missing records */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			Q[locus].x2 += pow(Q[locus].covariates[Y[0].use[record]], 2.0);
	}

	for (locus=0; locus<P; locus++)
	{
		X[locus].x2 = 0.0;
		for (record=0; record<Nn; record++)
			X[locus].x2 += pow(X[locus].covariates[Y[0].use[record]], 2.0);
	}

	/* allocation and initialization */
	for (locus=0; locus<F; locus++)
	{
		Q[locus].expEffect = 0.0;	Q[locus].varEffect = 0.0;	Q[locus].exp2Effect = pow(Q[locus].expEffect,2.0) + Q[locus].varEffect;
	}

	VarTemp=1.0;
	Tau0[0] = 100.0/VarTemp;	Tau0[1]=1.0;	vS2 = H[0].v * H[0].S2;
	logPi = log(H[0].Pi); log1minusPi = log(1.0-H[0].Pi);
	sumGamma[0]=0.0; sumGamma[1]=0.0;
	vcS2 = H[0].c * H[0].v * H[0].S2;

	if(RandomIni)/* Although we assign initial values for all loci, only the values for the first and second loci are used*/
	{
		for (locus=0; locus<P; locus++)
		{			
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
			X[locus].expEffect = RNormal(0.0,sqrt(X[locus].expSigma2/(double)P));	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;
			sumGamma[0] += X[locus].expGamma;
		}
	}
	else
	{
		for (locus=0; locus<P; locus++)
		{
			X[locus].expEffect = 0.0;	X[locus].varEffect = 0.0;	X[locus].exp2Effect = pow(X[locus].expEffect,2.0) + X[locus].varEffect;	
			X[locus].expSigma2 = VarTemp/2.0;	X[locus].S2 = 2.0/VarTemp;
			X[locus].expGamma = 0.5;
			sumGamma[0] += X[locus].expGamma;
		}
	}

	/* Calcualte residual errors */
	Y[0].expErrors = (double*) calloc (N, sizeof(double));
	for(record=0; record<Nn; record++)
	{
		Y[0].expErrors[Y[0].use[record]] = Y[0].stobs[Y[0].use[record]];
		for(locus=0; locus<F; locus++)
			Y[0].expErrors[Y[0].use[record]] -= Q[locus].covariates [Y[0].use[record]] * Q[locus].expEffect;

		for(locus=0; locus<P; locus++)
			Y[0].expErrors[Y[0].use[record]] -= X[locus].covariates [Y[0].use[record]] * X[locus].expEffect;
	}

	/* update order of marker effects */
	Order = (int*) calloc (P, sizeof(int));
	for(locus=0; locus<P; locus++) Order[locus] = locus;
	Shuffle (Order, P, P);
	Times = P / 10;

	if(Printinfo) printf ("Start (observation mean: %f SD %f)\n", Mean, Sd);
	start = clock ();

	/* start optimization */
	for (ite=1; ite<=Maxiteration; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;

		/* update of fixed effects */
		for (target=0; target<F; target++)
		{
			for (record=0, temp=0.0; record<Nn; record++)
					temp += Q[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + Q[target].covariates[Y[0].use[record]] * Q[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(Q[target].x2 * Tau0[0]);

			prop  = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record=0; record<Nn; record++)
				Y[0].expErrors[Y[0].use[record]] += (Q[target].covariates[Y[0].use[record]] * (Q[target].expEffect - prop));

			Check1 += pow((prop - Q[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			Q[target].expEffect  = prop;
		}

		/* update of B */
		Shuffle (Order, P, Times);
		sumGammaB2[0]=0.0; sumGammaB2[1]=0.0;
		for (locus=0; locus<P; locus++)
		{
			target = Order[locus];
			for(record=0, temp=0.0; record<Nn; record++)
				temp += X[target].covariates[Y[0].use[record]] * (Y[0].expErrors[Y[0].use[record]] + X[target].covariates[Y[0].use[record]] * X[target].expEffect);

			temp *= Tau0[0];
			temp2 = 1.0/(X[target].x2 * Tau0[0] + X[target].expGamma/X[0].S2 + (1.0 - X[target].expGamma)/X[1].S2);

			/* prop: E[B], prop2: E[B^2], temp2: V[B] */
			prop  = temp * temp2;
			prop2 = pow(prop, 2.0) + temp2;

			for(record=0; record<Nn; record++)
					Y[0].expErrors[Y[0].use[record]] += (X[target].covariates[Y[0].use[record]] * (X[target].expEffect - prop));

			temp3 = pow(prop, 2.0); 
			sumGammaB2[0] += temp3 * X[target].expGamma;
			sumGammaB2[1] += temp3 * (1.0-X[target].expGamma);

			Check1 += pow((prop - X[target].expEffect), 2.0);
			Check2 += pow(prop, 2.0);
			X[target].expEffect = prop;
		}

		/* update of Tau2 or Sigma2*/
		temp = sumGammaB2[0] + vS2;
		X[0].expSigma2 = temp/(H[0].v + sumGamma[0] - 2.0);
		prop = temp/(H[0].v + sumGamma[0]);

		Check1 += pow((prop - X[0].S2), 2.0);
		Check2 += pow(prop, 2.0);
		X[0].S2 = prop;

		temp = sumGammaB2[1] + vcS2;
		X[1].expSigma2 = temp/(H[0].v + (double)P - sumGamma[0] - 2.0);
		prop = temp/((double)P - sumGamma[0]);

		Check1 += pow((prop - X[1].S2), 2.0);
		Check2 += pow(prop, 2.0);
		X[1].S2 = prop;


		/* Update of Gamma (for wBSR, SSVS, and MIX) */
		Shuffle (Order, P, Times);
		ProbIncludeConstant = Digamma (0.5*(H[0].v+sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[0]+vS2)) + logPi;
		ProbExcludeConstant = Digamma (0.5*(H[0].v+(double)P-sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[1]+vS2)) + log1minusPi;
		/*---- Note -----------------------------------------------------------------------------------------------------
			sumGamma[0] and sumGammaB2 are not updated here, because these are the posterior parameters of X[0].expSigma. 
			sumGamma[0] is replaced by sumGamma[1] after the update of all gamma.
			sumGammaB2 are updated when marker effects are updated.
		---------------------------------------------------------------------------------------------------------------*/
		for(locus=0, sumGamma[1]=0.0; locus<P; locus++)
		{
			target = Order[locus];
			temp4 = pow(X[target].expEffect,2.0);
			ProbInclude = ProbIncludeConstant - 0.5 * temp4 / X[0].S2;
			ProbExclude = ProbExcludeConstant - 0.5 * temp4 / X[1].S2;

			temp = ProbInclude;
			if(temp < ProbExclude) temp = ProbExclude;
			ProbInclude -= temp;	ProbExclude -= temp;

			ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
			prop = ProbInclude/(ProbInclude + ProbExclude);
			prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

			/* for convergence check */
			Check1 += pow((prop - X[target].expGamma), 2.0);
			Check2 += pow(prop, 2.0);

			X[target].expGamma = prop;	
			X[target].exp2Gamma = prop2;
			sumGamma[1] += X[target].expGamma;
		}
		sumGamma[0] = sumGamma[1];

		/* update of Re */
		for (record=0, temp=0.0; record<Nn; record++)
			temp += pow(Y[0].expErrors[Y[0].use[record]], 2.0);

		prop = (double)(Nn - 2) / temp;
		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		/* Lower bound and residual variance */
		Rmonitor  [ite-1] = 1.0/Tau0[0];
		temp = Check1/Check2;

		if (Printinfo&&ite%100==0)	printf ( "%d Re2:%f Conv:%f\n", ite, Rmonitor[ite-1], log10(temp));
		if (ite==Maxiteration||temp<Thresholdvalue)
		{
			if(Printinfo) printf ( "Finish at %d (Re2: %f Conv: %f ", ite, Rmonitor[ite-1], log10(temp)); 
			break;
		}

	}	/* ite */

	end = clock ();
	FittingTime[0] = (double)(end-start)/CLOCKS_PER_SEC;
	if(Printinfo) printf ( "Time: %f sec.)\n\n", FittingTime[0]);

	/* Rescale estimates */
	for (locus=0; locus<F; locus++) { Q[locus].expEffect *= Sd; Q[locus].varEffect *= Var;}
	Q[0].expEffect += Mean;	/* the first effect is recognized as the overall mean */
	for (locus=0; locus<P; locus++) { X[locus].expEffect *= Sd; X[locus].varEffect *= Var;}
	Tau0[0] /= Var;
	Tau0[1] /= (Var*Var);

	free (Y[0].expErrors); free (Order); free (Y[0].use); free(Y[0].stobs);

	return(0);
}


/*---Perform genome-wide regression----------------------------------------------------------------------------------------------------------------------------*/
void GenomeWideRegression (int Algorithm, int Priortype, int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Y, Xstruct *X, Xstruct *Q, Hstruct *H, double *expDelta2, double *Tau0, double *LBmonitor, double *Rmonitor, double *FittingTime, int Printinfo,int RandomIni)
{
	switch (Methodcode){
	case 1:
		switch (Algorithm){
		case 1: BL_VB (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break; 
		case 2: BL_EM (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break;
		} 
		break;
	case 2:
		switch (Algorithm){
		case 1: EBL_VB (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break; 
		case 2: EBL_EM (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break;
		} 
		break;
	case 3:
		switch (Algorithm){
		case 1: wBSR_VB (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break; 
		case 2: wBSR_EM (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break;
		} 
		break;
	case 4:
		BayesC_VB (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);
		break;
	case 5:
		switch (Algorithm){
		case 1: SSVS_VB (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break; 
		case 2: SSVS_EM (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break;
		} 
		break;
	case 6:
		switch (Algorithm){
		case 1: MIX_VB (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break; 
		case 2: MIX_EM (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);break;
		} 
		break;
	case 7:
		BayesB_VB (Methodcode,CondResidual,P,F,N,Thresholdvalue,Maxiteration,Y,X,Q,H,expDelta2,Tau0,LBmonitor,Rmonitor,FittingTime,Printinfo,RandomIni);
		break;
	}
}

/*----make a prediction for record i ------------------------------------------------------------------------------------------------------------------*/
double PredictYhat (Xstruct *Q, int F, Xstruct *X, int P, int i, int Methodcode)
{
	/* i indicates the position in Q and X */

	int locus;
	double	temp=0.0;

	for (locus=0; locus<F; locus++)
		temp += Q[locus].covariates[i] * Q[locus].expEffect;

	if (Methodcode==3||Methodcode==4||Methodcode==7)
	{
		for (locus=0; locus<P; locus++)
			temp += X[locus].covariates[i] * X[locus].expEffect * X[locus].expGamma;
	}
	else
	{
		for (locus=0; locus<P; locus++)
			temp += X[locus].covariates[i] * X[locus].expEffect;
	}

	return(temp);
}

double PredictBV (Xstruct *X, int P, int i, int Methodcode)
{
	/* i indicates the position in Q and X */

	int locus;
	double	temp=0.0;

	if (Methodcode==3||Methodcode==4||Methodcode==7)
	{
		for (locus=0; locus<P; locus++)
			temp += X[locus].covariates[i] * X[locus].expEffect * X[locus].expGamma;
	}
	else
	{
		for (locus=0; locus<P; locus++)
			temp += X[locus].covariates[i] * X[locus].expEffect;
	}

	return(temp);
}

/*----Permutation tests--------------------------------------------------------------------------------------------------------------------------------*/
/* PermutationTest uses the distribution of all marker effects. This generated too much false positives.
	PermutationTestMax uses the distribution of maximum absolute marker effects. */
int PermutationTest (int Nperm, double *PermutationThreshold, int Algorithm, int Priortype, int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Yperm, Xstruct *Xperm, Xstruct *Qperm, Hstruct *H, int *Currentrun, int Nruns, int set, int Printinfo, int RandomIni)
{
	double	expDelta2, Tau0[2], *LBmonitor, *Rmonitor, FittingTime; /* these objects are required by GenomeWideRegression, but have no importance in the permutation tests */
	double	LimitCoef=2.0, Mean, Sd, Effect, Increment;
	int		n, record, locus, *CountPositive, *CountNegative, Nbreaks=10000, Nn, Position, Nsamples, Nsamples1, Nsamples5, Count, Check;
	
	/* allocate memory */
	LBmonitor = (double*) calloc (Maxiteration, sizeof(double));
	Rmonitor = (double*) calloc (Maxiteration, sizeof(double));
	CountPositive = (int*) calloc (Nbreaks+1, sizeof(double));
	CountNegative = (int*) calloc (Nbreaks+1, sizeof(double));

	/* Calculate the SD of observations. The limit of the effects is LimitCoef times Sd of observations */
	for (record=0, Mean=0.0, Nn=0; record<N; record++)
		if(Yperm[0].observations[record]!=Missingvalue){ Mean += Yperm[0].observations[record]; Nn++;}
	Mean /= (double) Nn;
	for (record=0, Sd=0.0; record<N; record++)
		if(Yperm[0].observations[record]!=Missingvalue)
			Sd += pow( (Yperm[0].observations[record] - Mean), 2.0);
	Sd /= (double) Nn;
	Sd = sqrt(Sd);

	/* Increment */
	Increment = LimitCoef * Sd / (double)Nbreaks;
	
	/* Perform permutation tests */
	for(n=0; n<Nperm; n++)
	{
		ShuffleDouble (Yperm[0].observations, N, N);
		if(Printinfo) printf ("****Progress:%d/%d trait:%s set:%d **** permutation %d/%d\n", Currentrun[0], Nruns, Yperm[0].name, set+1, n+1, Nperm);
		GenomeWideRegression (Algorithm, Priortype, Methodcode, CondResidual, P, F, N, Thresholdvalue, Maxiteration, Yperm, Xperm, Qperm, H, &expDelta2, Tau0, LBmonitor, Rmonitor, &FittingTime, Printinfo, RandomIni);
		Currentrun[0] ++;
		/* Count */
		for(locus=0; locus<P; locus++)
		{
			if(Methodcode==3||Methodcode==4||Methodcode==7) {Effect=Xperm[locus].expEffect*Xperm[locus].expGamma;} else {Effect=Xperm[locus].expEffect;}
			Position = (int)(fabs(Effect)/Increment);
			if(Position>Nbreaks) Position = Nbreaks;
			if(Effect>0){CountPositive[Position]++;} else {CountNegative[Position]++;}
		}
	}

	/* Calculate the PermutationThreshold */
	Nsamples = P * Nperm;
	Nsamples1 = (int)((double)Nsamples*0.01);
	Nsamples5 = (int)((double)Nsamples*0.05);
	for(n=Nbreaks, Count=0, Check=0; n>0; n--)
	{
		Count += CountPositive[n];
		if(Count>Nsamples1&&Check==0) {PermutationThreshold[0]=Increment*(double)(n+1); Check++;}
		if(Count>Nsamples5) {PermutationThreshold[1]=Increment*(double)(n+1); break;}
	}
	for(n=Nbreaks, Count=0, Check=0; n>0; n--)
	{
		Count += CountNegative[n];
		if(Count>Nsamples1&&Check==0) {PermutationThreshold[3]=-1.0*Increment*(double)(n+1); Check++;}
		if(Count>Nsamples5) {PermutationThreshold[2]=-1.0*Increment*(double)(n+1); break;}
	}
	free(LBmonitor); free(Rmonitor); free(CountPositive); free(CountNegative);
	return(0);
}

int PermutationTestMax (int Nperm, double *PermutationThreshold, int Algorithm, int Priortype, int Methodcode, int CondResidual, int P, int F, int N, double Thresholdvalue, int Maxiteration, 
	Ystruct *Yperm, Xstruct *Xperm, Xstruct *Qperm, Hstruct *H, int *Currentrun, int Nruns, int set, int Printinfo, int RandomIni)
{
	double	expDelta2, Tau0[2], *LBmonitor, *Rmonitor, FittingTime; /* these objects are required by GenomeWideRegression, but have no importance in the permutation tests */
	double	*MaximumEffects, Effect;
	int		n, locus, Quantile;
	
	/* allocate memory */
	LBmonitor = (double*) calloc (Maxiteration, sizeof(double));
	Rmonitor = (double*) calloc (Maxiteration, sizeof(double));
	MaximumEffects = (double*) calloc (Nperm, sizeof(double));
	
	/* Perform permutation tests */
	for(n=0; n<Nperm; n++)
	{
		ShuffleDouble (Yperm[0].observations, N, N);
		if(Printinfo) printf ("****Progress:%d/%d trait:%s set:%d **** permutation %d/%d\n", Currentrun[0], Nruns, Yperm[0].name, set+1, n+1, Nperm);
		GenomeWideRegression (Algorithm, Priortype, Methodcode, CondResidual, P, F, N, Thresholdvalue, Maxiteration, Yperm, Xperm, Qperm, H, &expDelta2, Tau0, LBmonitor, Rmonitor, &FittingTime, Printinfo, RandomIni);
		/* Search for maximum values */
		for(locus=0; locus<P; locus++)
		{
			if(Methodcode==3||Methodcode==4||Methodcode==7) {Effect=fabs(Xperm[locus].expEffect)*Xperm[locus].expGamma;} else {Effect=fabs(Xperm[locus].expEffect);}
			if(Effect>MaximumEffects[n]) MaximumEffects[n]=Effect;
		}
		Currentrun[0]++;
	}

	/* Calculate the PermutationThreshold */
	Sortarray(MaximumEffects, Nperm);
	Quantile=(int)(Nperm*0.01); /* 99% quantile */
	Quantile-=2; if(Quantile<0) Quantile=0;
	PermutationThreshold[0] = MaximumEffects[Quantile];
	Quantile=(int)(Nperm*0.05); /* 95% quantile */
	Quantile-=2; if(Quantile<0) Quantile=0;
	PermutationThreshold[1] = MaximumEffects[Quantile];

	//for(n=0;n<Nperm;n++)printf("%f ",MaximumEffects[n]); putchar('\n');
	//printf("%f %f\n", PermutationThreshold[0], PermutationThreshold[1]);

	free(LBmonitor); free(Rmonitor); free(MaximumEffects);
	return(0);
}

/*----Output results-----------------------------------------------------------------------------------------------------------------------------------*/
/* used when AnalysisType is 0 or 1. */
void WriteResults (int Algorithm, int Methodcode, char *Phenotypefile, char *Genotypefile, char *Covariatefile, char *Partitionfile,
					int AnalysisType, int Ntraits, int P, int F, int N, double Threshold, int Maxiteration, int usedseed, Ystruct *Y, Hstruct *H, Xstruct *X, Xstruct *Q, double *Rmonitor, double *LBmonitor, double expDelta2, 
					int trait, int set, int CVFoldTuning, int Nset, double *ValueList, double *PredictedValues, double *BreedingValues, double FittingTime, int MPP, int NDCGrank, int Nperm, double *PermutationThreshold, int RandomIni)
{

	FILE	*fp;
	char	Output[301]="", Method[11]="", AlgorithmLabel[11]="", MeasureLabel[21]="";
	int		locus, ii, record, SigSD, SigPerm=-9;
	double	temp;

	switch (Methodcode) {
		case 1: sprintf(Method, "BL");		break;
		case 2: sprintf(Method, "EBL");		break;
		case 3: sprintf(Method, "wBSR");	break;
		case 4: sprintf(Method, "BayesC");	break;
		case 5: sprintf(Method, "SSVS");	break;
		case 6: sprintf(Method, "MIX");		break;
		case 7: sprintf(Method, "BayesB");	break;
	}

	switch (Algorithm) {
		case 1: sprintf(AlgorithmLabel, "VB");	break;
		case 2: sprintf(AlgorithmLabel, "EM");	break;
	}

	switch (MPP) {
		case 1: sprintf(MeasureLabel, "MSE");	break;
		case 2: sprintf(MeasureLabel, "Cor");	break;
		case 3: sprintf(MeasureLabel, "NDCG@%d", NDCGrank);	break;
		case 4: sprintf(MeasureLabel, "LPDS");	break;
		case 5: sprintf(MeasureLabel, "MeanNDCG@%d", NDCGrank);	break;
		case 6: sprintf(MeasureLabel, "Kendall");	break;
	}

//	sprintf ( Output, "%s_%s_%s_set%d.fitting", Method, AlgorithmLabel, Y[trait].name, set+1 );
	sprintf ( Output, "%s_%s_set%d.fitting", Method, Y[trait].name, set+1 );
	fp = fopen ( Output, "w" );

	/* Marker effects and others*/
	fprintf (fp, "====Markers============================================================\n");
	switch (Algorithm) {
		case 1:	/* VB */
			switch (Methodcode) {
				case 1:	/* BL */
					fprintf (fp, "Effect SD Tau2 Lambda2 Sig(SD) Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						SigSD = JudgeSig(fabs(X[locus].expEffect),sqrt(X[locus].varEffect));
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						if(X[locus].expTau2>100000.0)
						{
							fprintf (fp, "%f %f %.6e %f %d %d\n", X[locus].expEffect, sqrt(X[locus].varEffect), X[locus].expTau2, expDelta2, SigSD, SigPerm);
						}
						else
						{
							fprintf (fp, "%f %f %f %f %d %d\n", X[locus].expEffect, sqrt(X[locus].varEffect), X[locus].expTau2, expDelta2, SigSD, SigPerm);
						}
					}
					break;
				case 2:	/* EBL */
					fprintf (fp, "Effect SD Tau2 Delta2 Eta2 Sig(SD) Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						SigSD = JudgeSig(fabs(X[locus].expEffect),sqrt(X[locus].varEffect));
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						if(X[locus].expTau2>100000.0)
						{
							fprintf (fp, "%f %f %.6e %f %f %d %d\n", X[locus].expEffect, sqrt(X[locus].varEffect), X[locus].expTau2, expDelta2, X[locus].expEta2, SigSD, SigPerm);
						}
						else
						{
							fprintf (fp, "%f %f %f %f %f %d %d\n", X[locus].expEffect, sqrt(X[locus].varEffect), X[locus].expTau2, expDelta2, X[locus].expEta2, SigSD, SigPerm);
						}
					}
					break;
				case 3:	/* wBSR */
					fprintf (fp, "EffectxGamma SD(EffectxGamma) Effect SD Sigma2 Gamma Sig(SD) Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						temp = X[locus].expGamma * (1.0-X[locus].expGamma) * X[locus].varEffect 
							+ X[locus].expGamma * (1.0-X[locus].expGamma) * pow(X[locus].expEffect, 2.0) 
							+ pow(X[locus].expGamma,2.0) * X[locus].varEffect;
						SigSD = JudgeSig(fabs(X[locus].expEffect*X[locus].expGamma),sqrt(temp));
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect*X[locus].expGamma,PermutationThreshold);
						fprintf (fp, "%f %f %f %f %f %f %d %d\n", X[locus].expEffect*X[locus].expGamma, sqrt(temp), X[locus].expEffect, sqrt(X[locus].varEffect), X[locus].expSigma2, X[locus].expGamma, SigSD, SigPerm);
					}
					break;
				case 4:	/* BayesC. "Effect" indicates E[B], i.e., E[B|gamma=1]*E[gamma=1].
						Similarly, "SD" indicates sqrt(Var[B]) where Var[B]=E[B|gamma=1]^2*Var[gamma]+Var[B|gamma=1]*E[gamma=1]*/
					fprintf (fp, "Effect SD Rho Sig(SD) Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						temp = pow(X[locus].expEffect,2.0)*X[locus].expGamma*(1.0-X[locus].expGamma)+X[locus].varEffect*X[locus].expGamma;
						SigSD = JudgeSig(fabs(X[locus].expEffect*X[locus].expGamma),sqrt(temp));
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect*X[locus].expGamma,PermutationThreshold);
						fprintf (fp, "%f %f %f %d %d\n", X[locus].expEffect*X[locus].expGamma, sqrt(temp), X[locus].expGamma, SigSD, SigPerm);
					}
					break;
				case 5: /* SSVS */
					fprintf (fp, "Effect SD Rho Sig(SD) Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						SigSD = JudgeSig(fabs(X[locus].expEffect),sqrt(X[locus].varEffect));
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						fprintf (fp, "%f %f %f %d %d\n", X[locus].expEffect, sqrt(X[locus].varEffect), X[locus].expGamma, SigSD, SigPerm);
					}
					break;
				case 6:	/* MIX */
					fprintf (fp, "Effect SD Rho Sig(SD) Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						SigSD = JudgeSig(fabs(X[locus].expEffect),sqrt(X[locus].varEffect));
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						fprintf (fp, "%f %f %f %d %d\n", X[locus].expEffect, sqrt(X[locus].varEffect), X[locus].expGamma, SigSD, SigPerm);
					}
					break;
				case 7:	/* BayesB. "Effect" indicates E[B], i.e., E[B|gamma=1]*E[gamma=1].
						Similarly, "SD" indicates sqrt(Var[B]) where Var[B]=E[B|gamma=1]^2*Var[gamma]+Var[B|gamma=1]*E[gamma=1]*/
					fprintf (fp, "Effect SD Rho Sigma2 Sig(SD) Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						temp = pow(X[locus].expEffect,2.0)*X[locus].expGamma*(1.0-X[locus].expGamma)+X[locus].varEffect*X[locus].expGamma;
						SigSD = JudgeSig(fabs(X[locus].expEffect*X[locus].expGamma),sqrt(temp));
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect*X[locus].expGamma,PermutationThreshold);
						fprintf (fp, "%f %f %f %f %d %d\n", X[locus].expEffect*X[locus].expGamma, sqrt(temp), X[locus].expGamma, X[locus].expSigma2, SigSD, SigPerm);
					}
					break;
			}
			break;
		case 2:	/* EM */
			switch (Methodcode) {
				case 1:	/* BL */
					fprintf (fp, "Effect Tau2 Lambda2 Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						if(X[locus].expTau2>100000.0)
						{
							fprintf (fp, "%f %.6e %f %d\n", X[locus].expEffect, X[locus].expTau2, expDelta2, SigPerm);
						}
						else
						{
							fprintf (fp, "%f %f %f %d\n", X[locus].expEffect, X[locus].expTau2, expDelta2, SigPerm);
						}
					}
					break;
				case 2:	/* EBL */
					fprintf (fp, "Effect Tau2 Delta2 Eta2 Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						if(X[locus].expTau2>100000.0)
						{
							fprintf (fp, "%f %.6e %f %f %d\n", X[locus].expEffect, X[locus].expTau2, expDelta2, X[locus].expEta2, SigPerm);
						}
						else
						{
							fprintf (fp, "%f %f %f %f %d\n", X[locus].expEffect, X[locus].expTau2, expDelta2, X[locus].expEta2, SigPerm);
						}
					}
					break;
				case 3:	/* wBSR */
					fprintf (fp, "EffectxGamma Effect Sigma2 Gamma, Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect*X[locus].expGamma,PermutationThreshold);
						fprintf (fp, "%f %f %f %f %d\n", X[locus].expEffect*X[locus].expGamma, X[locus].expEffect, X[locus].expSigma2, X[locus].expGamma, SigPerm);
					}
					break;
				case 4:	/* BayesC. Not supported. */
					fprintf (fp, "Effect Rho Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect*X[locus].expGamma,PermutationThreshold);
						fprintf (fp, "%f %f %d\n", X[locus].expEffect*X[locus].expGamma, X[locus].expGamma, SigPerm);
					}
					break;
				case 5: /* SSVS */
					fprintf (fp, "Effect Rho Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						fprintf (fp, "%f %f %d\n", X[locus].expEffect, X[locus].expGamma, SigPerm);
					}
					break;
				case 6:	/* MIX */
					fprintf (fp, "Effect Rho Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect,PermutationThreshold);
						fprintf (fp, "%f %f %d\n", X[locus].expEffect, X[locus].expGamma, SigPerm);
					}
					break;
				case 7:	/* BayesB. Not supported. */
					fprintf (fp, "Effect Rho Sigma2 Sig(Perm)\n");
					for (locus=0; locus<P; locus++)
					{
						if(Nperm>0) SigPerm = JudgeSigPermutation(X[locus].expEffect*X[locus].expGamma,PermutationThreshold);
						fprintf (fp, "%f %f %f %d\n", X[locus].expEffect*X[locus].expGamma, X[locus].expGamma, X[locus].expSigma2, SigPerm);
					}
					break;
			}
			break;
	}
	fprintf (fp, "\n");

	/* Fitted values*/
	fprintf (fp, "====Fitted values======================================================\n");
	fprintf (fp, "Yhat BV\n");
	for (record=0; record<N; record++)
		fprintf (fp, "%f %f\n", PredictedValues[record], BreedingValues[record]);
	fprintf (fp, "\n");

	/* Fixed effects */
	switch(Algorithm){
		case 1:/*VB*/
			fprintf (fp, "====Covariates=========================================================\n");
			fprintf (fp, "Effect SD\n");
			for (locus=0; locus<F; locus++)
				fprintf (fp, "%f %f\n", Q[locus].expEffect, sqrt(Q[locus].varEffect));
			fprintf (fp, "\n");
		break;
		case 2:/*EM*/
			fprintf (fp, "====Covariates=========================================================\n");
			fprintf (fp, "Effect\n");
			for (locus=0; locus<F; locus++)
				fprintf (fp, "%f\n", Q[locus].expEffect);
			fprintf (fp, "\n");
		break;
	}

	/* Lower bound of log likelihood of data */
	if(Algorithm==1)
	{
		fprintf (fp, "====logP(Y)============================================================\n");
		for ( ii=(Maxiteration-1); ii>=0; ii-- )if (LBmonitor[ii]<0.0) break;
		fprintf (fp, "%f\n", LBmonitor[ii]);
		fprintf (fp, "\n");
	}

	/* std.Sigma*/
	switch (Methodcode){
		case 4:
			fprintf (fp, "====Marker effect variance==============================================\n");
			fprintf (fp, "Sigma2: %f\n", X[0].expSigma2);
			fprintf (fp, "\n");
			break;
		case 5:
			fprintf (fp, "====Marker effect variance==============================================\n");
			fprintf (fp, "Sigma2: %f\n", X[0].expSigma2);
			fprintf (fp, "\n");
			break;
		case 6:
			fprintf (fp, "====Marker effect variances=============================================\n");
			fprintf (fp, "Sigma2A: %f\n", X[0].expSigma2);
			fprintf (fp, "Sigma2B: %f\n", X[1].expSigma2);
			fprintf (fp, "\n");
			break;
	}

	/* Monitor */
	fprintf (fp, "====Iterations=========================================================\n");
	fprintf (fp, "Iteration ResidualVariance\n");
	for ( ii=0; ii<Maxiteration; ii++ )
	{
		if (Rmonitor[ii]==0.0) break;
		fprintf (fp, "%d %f\n", ii+1, Rmonitor[ii]);
	}
	fprintf (fp, "\n");

	/* MSE etc. */
	if(AnalysisType==1&&Nset>1)
	{
		fprintf (fp, "====%s under each hyperparameter set====================================\n", MeasureLabel);
		switch (Methodcode) {
			case 1:	/* BL */
				fprintf (fp, "Set Phi Omega %s\n", MeasureLabel);
				for (ii=0; ii<Nset; ii++)
					fprintf (fp, "%d %f %f %f\n", ii+1, H[ii].deltaShape, H[ii].deltaRate, ValueList[ii]);
				break;
			case 2:	/* EBL */
				fprintf (fp, "Set Phi Omega Psi Theta %s\n", MeasureLabel);
				for (ii=0; ii<Nset; ii++)
					fprintf (fp, "%d %f %f %f %f %f\n", ii+1, H[ii].deltaShape, H[ii].deltaRate, H[ii].etaShape, H[ii].etaRate, ValueList[ii]);
				break;
			case 3:	/* wBSR */
				fprintf (fp, "Set Nu S2 Kappa %s\n", MeasureLabel);
				for (ii=0; ii<Nset; ii++)
					fprintf (fp, "%d %f %f %f %f\n", ii+1, H[ii].v, H[ii].S2, H[ii].Pi, ValueList[ii]);
				break;
			case 4:	/* BayesC */
				fprintf (fp, "Set Nu S2 Kappa %s\n", MeasureLabel);
				for (ii=0; ii<Nset; ii++)
					fprintf (fp, "%d %f %f %f %f\n", ii+1, H[ii].v, H[ii].S2, H[ii].Pi, ValueList[ii]);
				break;
			case 5:	/* SSVS */
				fprintf (fp, "Set c Nu S2 Kappa %s\n", MeasureLabel);
				for (ii=0; ii<Nset; ii++)
					fprintf (fp, "%d %f %f %f %f %f\n", ii+1, H[ii].c, H[ii].v, H[ii].S2, H[ii].Pi, ValueList[ii]);
				break;
			case 6:	/* MIX */
				fprintf (fp, "Set c Nu S2 Kappa %s\n", MeasureLabel);
				for (ii=0; ii<Nset; ii++)
					fprintf (fp, "%d %f %f %f %f %f\n", ii+1, H[ii].c, H[ii].v, H[ii].S2, H[ii].Pi, ValueList[ii]);
				break;
			case 7:	/* BayesB */
				fprintf (fp, "Set Nu S2 Kappa %s\n", MeasureLabel);
				for (ii=0; ii<Nset; ii++)
					fprintf (fp, "%d %f %f %f %f\n", ii+1, H[ii].v, H[ii].S2, H[ii].Pi, ValueList[ii]);
				break;
		}
		fprintf (fp, "Set%d was chosen\n", set+1);
		fprintf (fp, "\n");
	}

	/* PermutationThreshold */

	if(Nperm>0)
	{
	/* When PermutationTest is used
		fprintf (fp, "====Significance level=================================================\n");
		fprintf (fp, "1%: %f %f\n", PermutationThreshold[0], PermutationThreshold[3]);
		fprintf (fp, "5%: %f %f\n", PermutationThreshold[1], PermutationThreshold[2]);
		fprintf (fp, "\n");
	*/
		fprintf (fp, "====Significance level=================================================\n");
		fprintf (fp, "1%%: %f\n", PermutationThreshold[0]);
		fprintf (fp, "5%%: %f\n", PermutationThreshold[1]);
		fprintf (fp, "\n");
	}


	/* Parameters */
	fprintf ( fp, "====Parameters=========================================================\n");
	fprintf ( fp, "Phenotypefile	: %s\n", Phenotypefile );
	fprintf ( fp, "Genotypefile     : %s\n", Genotypefile);
	fprintf ( fp, "Covariatefile    : %s\n", Covariatefile);
	fprintf ( fp, "Partitionfile    : %s",	 Partitionfile);

		if(AnalysisType!=-9) { fprintf(fp, " (ignored)\n");}else{ fprintf(fp, "\n");}

	fprintf ( fp, "Option           : %d",   AnalysisType);

		if(AnalysisType==0)  fprintf (fp, " (fitting)\n");
		if(AnalysisType==1&&Nset>1)  fprintf (fp, " (fitting with the hyperparameters tuned by %d-fold CV)\n", CVFoldTuning);
		if(AnalysisType==1&&Nset==1)  fprintf (fp, " (fitting)\n");

	fprintf ( fp, "Method           : %s\n",   Method );
	fprintf ( fp, "Nperm            : %d\n",   Nperm );

		switch (Methodcode) {
			case 1: /* BL */
				fprintf (fp, "Phi              : %f\n", H[set].deltaShape);   
				fprintf (fp, "Omega            : %f\n", H[set].deltaRate);			
				break;
			case 2: /* EBL */
				fprintf (fp, "Phi              : %f\n", H[set].deltaShape);
				fprintf (fp, "Omega            : %f\n", H[set].deltaRate);			
				fprintf (fp, "Psi              : %f\n", H[set].etaShape);	
				fprintf (fp, "Theta            : %f\n", H[set].etaRate);
				break;
			case 3: /* wBSR */
				fprintf (fp, "Nu               : %f\n", H[set].v); 
				fprintf (fp, "S2               : %f\n", H[set].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[set].Pi); 
				break;
			case 4: /* BayesC */
				fprintf (fp, "Nu               : %f\n", H[set].v); 
				fprintf (fp, "S2               : %f\n", H[set].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[set].Pi); 			
				break;
			case 5: /* SSVS */
				fprintf (fp, "c                : %f\n", H[set].c); 
				fprintf (fp, "Nu               : %f\n", H[set].v); 
				fprintf (fp, "S2               : %f\n", H[set].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[set].Pi);
				break;
			case 6: /* MIX */
				fprintf (fp, "c                : %f\n", H[set].c); 
				fprintf (fp, "Nu               : %f\n", H[set].v); 
				fprintf (fp, "S2               : %f\n", H[set].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[set].Pi);
				break;
			case 7: /* BayesB */
				fprintf (fp, "Nu               : %f\n", H[set].v); 
				fprintf (fp, "S2               : %f\n", H[set].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[set].Pi); 			
				break;
		}

	fprintf ( fp, "Ntraits          : %d\n", Ntraits);
	fprintf ( fp, "K                : %d\n", trait+1);
	fprintf ( fp, "TraitName        : %s\n", Y[trait].name);	
	fprintf ( fp, "P                : %d\n", P);
	fprintf ( fp, "F                : %d\n", F);
	fprintf ( fp, "N                : %d\n", N);
	fprintf ( fp, "Threshold        : %f\n", Threshold);
	fprintf ( fp, "Maxiteration     : %d\n", Maxiteration);
//	fprintf ( fp, "Algorithm        : ");
//	if(Algorithm==1) {fprintf (fp, "VB\n");} else {fprintf (fp, "EM\n");}
	fprintf ( fp, "CVFoldTuning     : %d\n", CVFoldTuning);
	fprintf ( fp, "Randomize        : ");
	if(RandomIni) {fprintf (fp, "Yes\n");} else {fprintf (fp, "No\n");}
	fprintf ( fp, "seed             : %d\n", usedseed);
	fprintf ( fp, "Comp.time (sec.) : %f\n", FittingTime);
//	fprintf ( fp, "Measure          : %s\n", MeasureLabel);
	fprintf ( fp, "\n");

	/* Version information*/
	fprintf ( fp, "====Version============================================================\n");
	fprintf ( fp, "VIGoR C script version %s\n", Version);

	fclose(fp);
}


void WriteResultsPrediction (int Algorithm, int Methodcode, char *Phenotypefile, char *Genotypefile, char * Covariatefile, char *Partitionfile,
							int AnalysisType, int Ntraits, int P, int F, int N, double Threshold, int Maxiteration, int usedseed, 
							Ystruct *Y, Hstruct *H, Xstruct *X, Xstruct *Q, double *Rmonitor, double *LBmonitor, double expDelta2, double *PredictedValues, double *BreedingValues,
							int trait, int Fold, int *BestSet, double *BestValue, int Nset, double CVTime, int CVFoldTuning, Astruct *A, int MPP, int NDCGrank, int Nperm, int RandomIni)
{

	FILE *fp;
	char	Output[301]="", Method[11]="", AlgorithmLabel[11]="", MeasureLabel[21]="";
	int		record, fold, CumNtest;
//	double	*a; /* to check Pearson function */
//	int		ii; /* to check Pearson function */

	switch (Methodcode) {
		case 1: sprintf(Method, "BL");		break;
		case 2: sprintf(Method, "EBL");		break;
		case 3: sprintf(Method, "wBSR");	break;
		case 4: sprintf(Method, "BayesC");	break;
		case 5: sprintf(Method, "SSVS");	break;
		case 6: sprintf(Method, "MIX");		break;
		case 7: sprintf(Method, "BayesB");	break;
	}

	switch (Algorithm) {
		case 1: sprintf(AlgorithmLabel, "VB");	break;
		case 2: sprintf(AlgorithmLabel, "EM");	break;
	}

	switch (MPP) {
		case 1: sprintf(MeasureLabel, "MSE");	break;
		case 2: sprintf(MeasureLabel, "Cor");	break;
		case 3: sprintf(MeasureLabel, "NDCG@%d", NDCGrank);	break;
		case 4: sprintf(MeasureLabel, "LPDS");	break;
		case 5: sprintf(MeasureLabel, "MeanNDCG@%d", NDCGrank);	break;
		case 6: sprintf(MeasureLabel, "Kendall");	break;
	}

//	sprintf ( Output, "%s_%s_%s.crossvalidation", Method, AlgorithmLabel, Y[trait].name );
	sprintf ( Output, "%s_%s.crossvalidation", Method, Y[trait].name );
	fp = fopen ( Output, "w" );

	/* Predicted values */
	fprintf (fp, "====Predicted values==================================================\n");
	fprintf (fp, "Test Y Yhat BV\n");

	for ( fold=0, CumNtest=0; fold<Fold; fold++)
	{
		for ( record=0; record<A[fold].ntest; record++ )
		{
			if(Y[trait].observations[A[fold].obstest[record]]==Missingvalue)
			{
				fprintf (fp, "%d NA NA NA\n", A[fold].obstest[record]+1);
			} else {
				fprintf (fp, "%d %f %f %f\n", A[fold].obstest[record]+1, Y[trait].observations[A[fold].obstest[record]], PredictedValues[record+CumNtest], BreedingValues[record+CumNtest]);
			}
		}
		CumNtest += A[fold].ntest;
	}
	fprintf (fp, "\n");
/*
	a=(double*)calloc(N,sizeof(double));
	for ( fold=0, ii=0; fold<Fold; fold++)
		for ( record=0; record<A[fold].ntest; record++ )
		{
			a[ii]=Y[trait].observations[A[fold].obstest[record]];
			ii++;
		}
	printf("Cor.%f\n", Pearson(a, PredictedValues, N));
	free(a);
*/

	/* MSE etc. */
	if(Nset>1)
	{
		fprintf (fp, "====%s===============================================================\n",MeasureLabel);
		switch (Methodcode) {
			case 1:	/* BL */
				fprintf (fp, "Fold ChosenSet Phi Omega %s\n", MeasureLabel);
				for (fold=0; fold<Fold; fold++)
					fprintf (fp, "%d %d %f %f %f\n", fold+1, BestSet[fold]+1, H[BestSet[fold]].deltaShape, H[BestSet[fold]].deltaRate, BestValue[fold]);
				break;
			case 2:	/* EBL */
				fprintf (fp, "Fold ChosenSet Phi Omega Psi Theta %s\n", MeasureLabel);
				for (fold=0; fold<Fold; fold++)
					fprintf (fp, "%d %d %f %f %f %f %f\n", fold+1, BestSet[fold]+1, H[BestSet[fold]].deltaShape, H[BestSet[fold]].deltaRate, 
							H[BestSet[fold]].etaShape, H[BestSet[fold]].etaRate, BestValue[fold]);
				break;
			case 3:	/* wBSR */
				fprintf (fp, "Fold ChosenSet Nu S2 Kappa %s\n", MeasureLabel);
				for (fold=0; fold<Fold; fold++)
					fprintf (fp, "%d %d %f %f %f %f\n", fold+1, BestSet[fold]+1, H[BestSet[fold]].v, H[BestSet[fold]].S2, H[BestSet[fold]].Pi, BestValue[fold]);
				break;
			case 4:	/* BayesC */
				fprintf (fp, "Fold ChosenSet Nu S2 Kappa %s\n", MeasureLabel);
				for (fold=0; fold<Fold; fold++)
					fprintf (fp, "%d %d %f %f %f %f\n", fold+1, BestSet[fold]+1, H[BestSet[fold]].v, H[BestSet[fold]].S2, H[BestSet[fold]].Pi, BestValue[fold]);
				break;
			case 5:	/* SSVS */
				fprintf (fp, "Fold ChosenSet c Nu S2 Kappa %s\n", MeasureLabel);
				for (fold=0; fold<Fold; fold++)
					fprintf (fp, "%d %d %f %f %f %f %f\n", fold+1, BestSet[fold]+1, H[BestSet[fold]].c, H[BestSet[fold]].v, H[BestSet[fold]].S2, H[BestSet[fold]].Pi, BestValue[fold]);
				break;
			case 6:	/* MIX */
				fprintf (fp, "Fold ChosenSet c Nu S2 Kappa %s\n", MeasureLabel);
				for (fold=0; fold<Fold; fold++)
					fprintf (fp, "%d %d %f %f %f %f %f\n", fold+1, BestSet[fold]+1, H[BestSet[fold]].c, H[BestSet[fold]].v, H[BestSet[fold]].S2, H[BestSet[fold]].Pi, BestValue[fold]);
				break;
			case 7:	/* BayesB */
				fprintf (fp, "Fold ChosenSet Nu S2 Kappa %s\n", MeasureLabel);
				for (fold=0; fold<Fold; fold++)
					fprintf (fp, "%d %d %f %f %f %f\n", fold+1, BestSet[fold]+1, H[BestSet[fold]].v, H[BestSet[fold]].S2, H[BestSet[fold]].Pi, BestValue[fold]);
				break;
		}
		fprintf (fp, "\n");
	}

	/* Parameters */
	fprintf ( fp, "====Parameters=========================================================\n");
	fprintf ( fp, "Phenotypefile	: %s\n", Phenotypefile );
	fprintf ( fp, "Genotypefile     : %s\n", Genotypefile);
	fprintf ( fp, "Covariatefile    : %s\n", Covariatefile);
	fprintf ( fp, "Partitionfile    : %s",	 Partitionfile);

		if(AnalysisType!=-9) { fprintf(fp, " (ignored)\n");}else{ fprintf(fp, "\n");}

	fprintf ( fp, "Option           : %d",   AnalysisType);

		if(AnalysisType==-9) fprintf (fp, " (%d-fold CV)\n", Fold);
		if(AnalysisType==-1) fprintf (fp, " (Leave-one-out CV)\n");
		if(AnalysisType>1)   fprintf (fp, " (%d-fold CV)\n", Fold);

	fprintf ( fp, "Method           : %s\n",   Method );
	fprintf ( fp, "Nperm            : %d\n",   Nperm );

	if(Nset==1)
		switch (Methodcode) {
			case 1: 
				fprintf (fp, "Phi              : %f\n", H[0].deltaShape);   
				fprintf (fp, "Omega            : %f\n", H[0].deltaRate);			
				break;
			case 2:
				fprintf (fp, "Phi              : %f\n", H[0].deltaShape);
				fprintf (fp, "Omega            : %f\n", H[0].deltaRate);			
				fprintf (fp, "Psi              : %f\n", H[0].etaShape);	
				fprintf (fp, "Theta            : %f\n", H[0].etaRate);
				break;
			case 3: 
				fprintf (fp, "Nu               : %f\n", H[0].v); 
				fprintf (fp, "S2               : %f\n", H[0].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[0].Pi); 
				break;
			case 4: 
				fprintf (fp, "Nu               : %f\n", H[0].v); 
				fprintf (fp, "S2               : %f\n", H[0].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[0].Pi); 			
				break;
			case 5:
				fprintf (fp, "c                : %f\n", H[0].c); 
				fprintf (fp, "Nu               : %f\n", H[0].v); 
				fprintf (fp, "S2               : %f\n", H[0].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[0].Pi);
				break;
			case 6:
				fprintf (fp, "c                : %f\n", H[0].c); 
				fprintf (fp, "Nu               : %f\n", H[0].v); 
				fprintf (fp, "S2               : %f\n", H[0].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[0].Pi);
				break;
			case 7: 
				fprintf (fp, "Nu               : %f\n", H[0].v); 
				fprintf (fp, "S2               : %f\n", H[0].S2); 			
				fprintf (fp, "Kappa            : %f\n", H[0].Pi); 
				break;
		}

	fprintf ( fp, "Ntraits          : %d\n", Ntraits);
	fprintf ( fp, "K                : %d\n", trait+1);
	fprintf ( fp, "TraitName        : %s\n", Y[trait].name);
	fprintf ( fp, "P                : %d\n", P);
	fprintf ( fp, "F                : %d\n", F);
	fprintf ( fp, "N                : %d\n", N);
	fprintf ( fp, "Threshold	: %f\n", Threshold);
	fprintf ( fp, "Maxiteration     : %d\n", Maxiteration);
//	fprintf ( fp, "Algorithm        : ");
//	if(Algorithm==1) {fprintf (fp, "VB\n");} else {fprintf (fp, "EM\n");}
	fprintf ( fp, "CVFoldTuning     : %d\n", CVFoldTuning);
	fprintf ( fp, "Randomize        : ");
	if(RandomIni) {fprintf (fp, "Yes\n");} else {fprintf (fp, "No\n");}
	fprintf ( fp, "seed             : %d\n", usedseed);
	fprintf ( fp, "Comp.time (sec.) : %f\n", CVTime);
//	fprintf ( fp, "Measure          : %s\n", MeasureLabel);
	fprintf ( fp, "\n");

	/* Version information*/
	fprintf ( fp, "====Version============================================================\n");
	fprintf ( fp, "VIGoR C script version %s\n", Version);

	fclose(fp);
}

/* Count the number of rows and columns of input files------------------------------------------------------------------------------------------------*/
int CountRowsColumns (int *N, int *P, char *Filename)
{
	char	v[MaxBufSize]={0};
	char	*temp, *CR;
	int		n, p, back;
	FILE	*fp;

	/* count the number of columns from the first row */
	if((fp=fopen(Filename,"r"))==NULL){fprintf(stderr,"Cannot open %s\n", Filename); return(1);}
	fgets(v,sizeof(v),fp);
	CR=strchr(v, '\n');
	if (CR != NULL) {
		v[strlen(v) - 1] = '\0';
	} else {
		fprintf(stderr,"Buffer overflow\n"); return(1);
	}

	if(strncmp(CR-1,"\r",1)==0){back=2;}else{back=1;}	
	while (strncmp(CR-back," ",1)==0||strncmp(CR-back,"\t",1)==0)
	{		
		CR=CR-back;
		strcpy(CR,"\n");
		back=1;
	}
	temp=strtok(v," \t");
	p=1;
	while((temp=strtok(NULL," \t"))!=NULL){p++;}
	//printf("%d\n",p);
	*P=p;
	fclose(fp);

	fp=fopen(Filename,"r");
	n=0;
    while(fgets(v,sizeof(v),fp)) {n++;}
	*N=n;
	fclose(fp);
	return(0);
}
