/*
	Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata

	Released under the MIT license
	http://opensource.org/licenses/mit-license.php
*/

/*----define global objects------------------------------------------------------------------------------------------------------*/
#define Version				"1.0 (Apr. 2015)"
#define	MaxBufSize			1048576			/* Maximum buffer size (32^4) */

/*----define structure-----------------------------------------------------------------------------------------------------------*/
typedef struct {

	double	*covariates;	/* covariates */
	double	x2;				/* sum of squared covariates */
	double	expEffect;		/* expectation of marker effect. When BayesC, this indicates E[B|gamma=1] */
	double	varEffect;		/* variance of marker effect. When BayesC, this indicates V[B|gamma=1] */
	double	exp2Effect;		/* 2nd moment (for wBSR, BayesC, SSVS, and MIX). When BayesC, this indicates E[B^2|gamma=1] */
	double	expGamma;		/* expectation of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	exp2Gamma;		/* 2nd moment (for wBSR) */
	double	expTau2;		/* expectation of tau (for BL and EBL) */
	double	expInTau2;		/* expectation of the inverse of tau (for BL and EBL) */
//	double	varTau2;		/* variance of tau (for BL and EBL) */
	double	expEta2;		/* expectation of eta (for EBL) */
	double	expSigma2;		/* expectation of sigma (for wBSR, BayesC, SSVS, and MIX).*/
	double	S2;				/* parameter of the posterior distribution of marker effect variance (for wBSR, BayesC, SSVS, and MIX) */
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

	double	*deltaShape;	/* 0 (for BL and EBL) */ 
	double	*deltaRate;		/* 1 (for BL and EBL) */
	double	*etaShape;		/* 2 (for EBL) */
	double	*v;				/* 3 (for wBSR, BayesB, BayesC, SSVS, and MIX) */
	double	*Pi;			/* 4 (for all the methods) */
	double	*A;				/* 5 (for SSVS and MIX) */

} Vstruct;


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

/* Variance ----------------------------------------------------------------------------------------------------------------------------------------------*/
double Variance (int N, double *Array)
{
	int i;
	double mean=0.0, sum=0.0;

	for(i=0; i<N; i++){mean += Array[i];}
	mean /= (double)N;
	for(i=0; i<N; i++){sum += pow((Array[i]-mean),2.0);}
	return(sum/(double)(N-1));
}