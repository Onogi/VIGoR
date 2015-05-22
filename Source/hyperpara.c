#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "hyperpara.h"

/*	Calculate hyperparameter values for VIGoR

	Seven regression models are supported:
		1.BL		: Bayeian Lasso (Park and Casella 2008)
		2.EBL		: Extended Bayesian Lasso (Mutshinda and Sillanpaa 2010)
		3.wBSR		: weighted Bayesian Shrinkage Regression (Hayashi and Iwata 2010)
		4.BayesC	: BayesC (Habier et al. 2011)
		5.SSVS		: Stchastic search variable selection (George and McCulloch 1997)
		6.MIX		: Bayesian mixture regression (Luan et al. 2009)
		7.BayesB    : BayesB (Meuwissen et al. 2001)
	
	Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata

	Released under the MIT license
	http://opensource.org/licenses/mit-license.php

*/


int main (int argc, char *argv[])
{

	/* Values defined from arguments, options, or parameter files */
	char	Genotypefile[401]="";		/* Genotype file name */
	double	Mvar;
	double	Incoef;
	int     Methodcode;					/* 1,BL; 2,EBL; 3,wBSR; 4,BayesC; 5,SSVS; 6,MIX; 7,BayesB */
	int		P;							/* Number of markers */
	int		N;							/* Number of individuals */
	int		UsePNmatrix;				/* Use P x N matrix for Genotypefile (1) or not (0); */
	int		Xtype;
	int		*Nvalue;
	Vstruct	*GivenValues;
	
	/* Values defined by the developer */
	int		Printinfo=1;				/* Whether print out run information (1) or not (0); */
	int		Mnset=1000;					/* Maximum number of hyperparameter sets */
	int		MaxNumHyperpara=4;			/* Maximum number of hyperparameter values */
	int		Fileformat=1;				/* File format. 1: Default, 2: Plink, 3: Beagle */
	
	/* Number of values that are given for calculation and the upper limit of the number of given values*/
	int		Ngivenvalues=6;
	int		Lgivenvalues=100;
	int		Limitover=0;

	/* Whether values are initialized (1) or not (0) */
	int		IniA=0, IniBLPhi=0, IniEBLPhi=0, IniOmega=0, IniPsi=0, IniNu=0;
	
	/* Number of arguments */
	int		Narg;

	/* for reading PED file */
	int		*AlleleA, allele;

	/* For repeat statement */
	int		ii, jj;
	int		record, locus, set;

	/* Generic use */
	double	temp;

	/* Genotypes */
	Xstruct	*X;

	/* For specification of hyperparameters */
	Hstruct	*H;
	int		Nh, Nset;

	/* Hyperparameter calculation */
	double	SumVar=0.0, Af;
	int		Neach, Ntimes, each, time;

	/* To print results */
	double	*Print;

	/* Inputs and outputs */
	char    Output[301]="", Chartemp[101]="", Method[11]="", AlgorithmLabel[11]="";
	char	v[MaxBufSize]=""; /* for Beagle dosage file */
	FILE    *fp;
	
	/* output help */
	if(getopt(argc, argv, "k:a:f:b:p:g:s:n:toh")=='h')
	{
		printf("Usage\n\n"); 
		printf(" hyperpara <Genotype file> <Mvar>\n");
		printf("       <Method BL|EBL|wBSR|BayesB|BayesC|SSVS|MIX>\n");
		printf("       <Kappa> [options]\n\n");
		printf("Genotype file consists of a N x P (number of markers) matrix. P x N matrix can be used by using option -o.\n\n");
		printf("NOTE: the maximum size of one row is 32^4=1MB. Use P x N matrix for genotype file when P is large and buffer overflow occurs\n\n");
		printf("PED files of Plink and dosage files of Beagle can be used as Genotype file\n");
		printf("PED files include the first 6 mandantory columns. Missing genotypes are not allowed\n");
		printf("In PED files, two alleles of a genotyps are stored in two columns. Alleles can be coded by any single character\n\n");
		printf("Mvar is a scalar representing the assumed proportion of phenotypic variance that markers can explain\n");
		printf("Mvar is < 1.0 for BL and EBL, and <= 1.0 for the other methods.\n");
		printf("Methods are\n");
		printf("  BL    : Bayesian lasso\n");
		printf("  EBL   : Extended Bayesian lasso\n");
		printf("  wBSR  : Weighted Bayesian shrinkage regression\n");
		printf("  BayesB: BayesB\n");
		printf("  BayesC: BayesC\n");
		printf("  SSVS  : Stochastic search variable selection\n");
		printf("  MIX   : Bayesian mixture model\n\n");
		printf("Kappa is a scaler representing the assumed proportion of markers that have non-zero effects.\n");
		printf("Kappa is 0<Kappa<=1. For MIX and SSVS, Kappa is 0<Kappa<1.\n");
		printf("See the pdf manual of VIGoR for the details of hyperpara\n\n");
		printf("Options that take arguments\n");
		printf("-k: Additional Kappa values. Multiple Kappa values can be specified by using this option repeatedly\n");
		printf("-a: Determines A value. A is used when SSVS and MIX, and represents the proportion of Mvar that\n");
		printf("    markers assigned to the normal prior distribution with the larger variance can explain.\n");
		printf("    For example, if Mvar is 0.5, A is 0.9, and Kappa is 0.01, 0.45 (= 0.5 x 0.9) of phenotypic variance is assumed to\n");
		printf("    be explained by markers, and 0.05 (= 0.5 x (1-0.9)) of variance is assumed to be explained by of markers\n");
		printf("-f: A scalar representing the inbreeding coefficient. Enter 1 for inbred species.\n");
		printf("-b: A scaler representing Phi values for BL.\n");
		printf("-p: A scaler representing Phi values for EBL.\n");
		printf("-g: A scaler representing Omega values for EBL.\n");
		printf("-s: A scaler representing Psi values for EBL.\n");
		printf("-n: A scaler representing Nu values for wBSR, BayesB, BayesC, SSVS, and MIX.\n");
		printf("****Note****\n");
		printf("   Options -k, -a, -b, -p, -g, -s, and -n can be used repeatedly, which makes multiple combinations (sets) of hyperparameter values\n");
		printf("************\n");
		printf("Options\n");
		printf("-t: Treat variables in Genotype file as general variables\n");
		printf("-o: Use P x N matrix for Genotype file. Recommended when P>>N. \n");
		printf("-h: help\n\n");
		return(0);
	}

	/* Read arguments */
	Narg=1;
	strcpy(Genotypefile, argv[Narg]);Narg++; //Narg==2
	if(strcmp(Genotypefile+(strlen(Genotypefile)-4),".ped")==0)
	{
		Fileformat=2; /* PED file for PLINK */
	} 
	else
	{
		if(strcmp(Genotypefile+(strlen(Genotypefile)-5),".dose")==0) Fileformat=3; /* Beagle dosage file */
	}

	Mvar = atof(argv[Narg]);Narg++; //Narg==3
	strcpy(Method,argv[Narg]);Narg++; //Narg==4

	GivenValues = (Vstruct*)calloc(1,sizeof(Vstruct));
	GivenValues[0].Pi=(double*)calloc(Lgivenvalues,sizeof(double));
	if(strcmp(Method,"BL")==0)	{
		Methodcode=1; Nh = 2; 
		GivenValues[0].deltaShape=(double*)calloc(Lgivenvalues,sizeof(double));
	}
	if(strcmp(Method,"EBL")==0)	{
		Methodcode=2; Nh = 4;		
		GivenValues[0].deltaShape=(double*)calloc(Lgivenvalues,sizeof(double));
		GivenValues[0].deltaRate =(double*)calloc(Lgivenvalues,sizeof(double));
		GivenValues[0].etaShape  =(double*)calloc(Lgivenvalues,sizeof(double));
	}
	if(strcmp(Method,"wBSR")==0)	{
		Methodcode=3; Nh = 3;
		GivenValues[0].v=(double*)calloc(Lgivenvalues,sizeof(double));
	}
	if(strcmp(Method,"BayesB")==0)	{
		Methodcode=7; Nh = 3;
		GivenValues[0].v=(double*)calloc(Lgivenvalues,sizeof(double));
	}
	if(strcmp(Method,"BayesC")==0)	{
		Methodcode=4; Nh = 3;
		GivenValues[0].v=(double*)calloc(Lgivenvalues,sizeof(double));
	}
	if(strcmp(Method,"SSVS")==0)	{
		Methodcode=5; Nh = 4;
		GivenValues[0].v=(double*)calloc(Lgivenvalues,sizeof(double));
		GivenValues[0].A=(double*)calloc(Lgivenvalues,sizeof(double));
	}
	if(strcmp(Method,"MIX")==0)	{
		Methodcode=6; Nh = 4;
		GivenValues[0].v=(double*)calloc(Lgivenvalues,sizeof(double));
		GivenValues[0].A=(double*)calloc(Lgivenvalues,sizeof(double));
	}
	if(strcmp(Method,"BL")!=0&&strcmp(Method,"EBL")!=0&&strcmp(Method,"wBSR")!=0&&strcmp(Method,"BayesB")!=0&&strcmp(Method,"BayesC")!=0&&strcmp(Method,"SSVS")!=0&&strcmp(Method,"MIX")!=0) 
	{fprintf(stderr,"Number of arguments or method specification\n");return(0);}

	Nvalue = (int*)calloc(Ngivenvalues,sizeof(int));
	GivenValues[0].Pi[0] = atof(argv[Narg]); Narg++; //Narg==5.
	if(Methodcode==5||Methodcode==6){
		if(GivenValues[0].Pi[0]<=0.0 && GivenValues[0].Pi[0]>=1.0) {fprintf(stderr, "Kappa should be 0<Kappa<1\n"); return(0);}
	}
	else
	{
		if(GivenValues[0].Pi[0]<=0.0 && GivenValues[0].Pi[0]>1.0) {fprintf(stderr, "Kappa should be 0<Kappa<=1\n"); return(0);}
	}
	Nvalue[4]=1;

	/*Default values*/
	if(Methodcode==1){
		GivenValues[0].deltaShape[0]=1.0;
		Nvalue[0]=1;
	}
	if(Methodcode==2){
		GivenValues[0].deltaShape[0]=0.1;	Nvalue[0]=1;
		GivenValues[0].deltaRate[0]=0.1;	Nvalue[1]=1;
		GivenValues[0].etaShape[0]=1.0;		Nvalue[2]=1;
	}
	if(Methodcode==3||Methodcode==4||Methodcode==7){
		GivenValues[0].v[0]=5.0;	Nvalue[3]=1;

	}
	if(Methodcode==5||Methodcode==6){
		GivenValues[0].v[0]=5.0;	Nvalue[3]=1;
		GivenValues[0].A[0]=0.9;	Nvalue[5]=1;
	}
	Incoef=0.0;
	UsePNmatrix=0;
	Xtype=0;

	/*Read options*/
	optind=Narg;
	while(((ii=getopt(argc, argv, "k:a:f:b:p:g:s:n:toh"))!=-1)||(Limitover)){
		switch(ii){
		case 'k':
			GivenValues[0].Pi[Nvalue[4]] = atof(optarg);
			if(Methodcode==5||Methodcode==6){
				if(GivenValues[0].Pi[Nvalue[4]]<=0.0 && GivenValues[0].Pi[Nvalue[4]]>=1.0) {fprintf(stderr, "Kappa should be 0<Kappa<1\n"); return(0);}
			}
			else
			{
				if(GivenValues[0].Pi[Nvalue[4]]<=0.0 && GivenValues[0].Pi[Nvalue[4]]>1.0) {fprintf(stderr, "Kappa should be 0<Kappa<=1\n"); return(0);}
			}
			Nvalue[4]++; if(Nvalue[4]==Lgivenvalues) Limitover=1;
			break;
		case 'a':
			if(Methodcode==5||Methodcode==6)
			{
				if(IniA==0){Nvalue[5]=0;IniA=1;}
				GivenValues[0].A [Nvalue[5]] = atof(optarg);
				if(GivenValues[0].A [Nvalue[5]]<=0||GivenValues[0].A [Nvalue[5]]>=1.0) {fprintf(stderr,"A should be 0<A<1\n"); return(0);}
				Nvalue[5]++; if(Nvalue[5]==Lgivenvalues) Limitover=1;				
			}
			break;
		case 'f':
			Incoef = atof(optarg);
			if(Incoef>1.0||Incoef<0.0){fprintf(stderr,"Inbreeding coefficient should be 0<=f<=1\n"); return(0);}
			break;
		case 'b':
			if(Methodcode==1){
				if(IniBLPhi==0){Nvalue[0]=0;IniBLPhi=1;}
				GivenValues[0].deltaShape [Nvalue[0]] = atof(optarg);
				if(GivenValues[0].deltaShape [Nvalue[0]]<=0) {fprintf(stderr,"Phi of BL should be positive\n"); return(0);}
				Nvalue[0]++; if(Nvalue[0]==Lgivenvalues) Limitover=1;
			}
			break;				
		case 'p':
			if(Methodcode==2){
				if(IniEBLPhi==0){Nvalue[0]=0;IniEBLPhi=1;}
				GivenValues[0].deltaShape [Nvalue[0]] = atof(optarg); 
				if(GivenValues[0].deltaShape [Nvalue[0]]<=0) {fprintf(stderr,"Phi of EBL should be positive\n"); return(0);}
				Nvalue[0]++; if(Nvalue[0]==Lgivenvalues) Limitover=1;
			}
			break;
		case 'g':
			if(Methodcode==2){
				if(IniOmega==0){Nvalue[1]=0;IniOmega=1;}
				GivenValues[0].deltaRate [Nvalue[1]] = atof(optarg); 
				if(GivenValues[0].deltaRate [Nvalue[1]]<=0) {fprintf(stderr,"Omega of EBL should be positive\n"); return(0);}
				Nvalue[1]++; if(Nvalue[1]==Lgivenvalues) Limitover=1;
			}
			break;
		case 's':
			if(Methodcode==2){
				if(IniPsi==0){Nvalue[2]=0;IniPsi=1;}
				GivenValues[0].etaShape [Nvalue[2]] = atof(optarg); 
				if(GivenValues[0].etaShape [Nvalue[2]]<=0) {fprintf(stderr,"Psi of EBL should be positive\n"); return(0);}
				Nvalue[2]++; if(Nvalue[2]==Lgivenvalues) Limitover=1;
			}
			break;
		case 'n':
			if(!(Methodcode==1||Methodcode==2)){
				if(IniNu==0){Nvalue[3]=0;IniNu=1;}
				GivenValues[0].v [Nvalue[3]] = atof(optarg);
				if(GivenValues[0].v [Nvalue[3]]<=0) {fprintf(stderr,"Nu should be positive\n"); return(0);}
				Nvalue[3]++; if(Nvalue[3]==Lgivenvalues) Limitover=1;
			}
			break;
		case 't':
			Xtype=1;
			break;
		case 'o':
			UsePNmatrix=1;
			break;
		}
	}

	if(Printinfo){
		printf("Genotype file   : %s\n", Genotypefile);
		printf("Mvar            : %f\n", Mvar);
		printf("Method          : %s\n", Method);
		printf("Kappa           : ");
		for(ii=0; ii<Nvalue[4]; ii++) {printf("%f ",GivenValues[0].Pi[ii]);} printf("\n");
		if(Methodcode==1){
			printf("Phi             : ");
			for(ii=0; ii<Nvalue[0]; ii++){printf("%f ", GivenValues[0].deltaShape[ii]);} printf("\n");

		}
		if(Methodcode==2){
			printf("Phi             : ");
			for(ii=0; ii<Nvalue[0]; ii++){printf("%f ", GivenValues[0].deltaShape[ii]);} printf("\n");
			printf("Omega           : ");
			for(ii=0; ii<Nvalue[1]; ii++){printf("%f ", GivenValues[0].deltaRate[ii]);} printf("\n");
			printf("Psi             : ");
			for(ii=0; ii<Nvalue[2]; ii++){printf("%f ", GivenValues[0].etaShape[ii]);} printf("\n");

		}
		if(Methodcode==3||Methodcode==4||Methodcode==7){
			printf("Nu              : ");
			for(ii=0; ii<Nvalue[3]; ii++){printf("%f ", GivenValues[0].v[ii]);} printf("\n");

		}
		if(Methodcode==5||Methodcode==6){
			printf("Nu              : ");
			for(ii=0; ii<Nvalue[3]; ii++){printf("%f ", GivenValues[0].v[ii]);} printf("\n");
			printf("A               : ");
			for(ii=0; ii<Nvalue[5]; ii++){printf("%f ", GivenValues[0].A[ii]);} printf("\n");

		}
		printf("Inbreeding coef.: %f\n", Incoef);
		if(UsePNmatrix) printf("Genotype file includes P x N matrix\n");
		if(Xtype)		printf("Treat predictor variables as general variables\n");
		if(Limitover) printf("The number(s) of given values for parameter(s) exceeds the upper limit (%d). The first %d values are used\n", Lgivenvalues,Lgivenvalues); 
	}

	/* Open the genotype files */
	switch(Fileformat){
		case 1:
			if(UsePNmatrix)
			{
				jj=CountRowsColumns(&P, &N, Genotypefile);if(jj)return(0);
				X = (Xstruct*) calloc ( P, sizeof (Xstruct) );

				fp = fopen( Genotypefile, "r" );
				for (locus=0; locus<P; locus++)
				{
					X[locus].covariates = (double*) calloc (N, sizeof(double));
				}

				for (locus=0; locus<P; locus++)
					for (record=0; record<N; record++)
					{
						if (fscanf (fp, "%lf", X[locus].covariates+record) == -1) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
						if (X[locus].covariates[record]<0.0||X[locus].covariates[record]>2.0) 
						{
							fprintf (stderr,"Genotypes should be coded as 0(AA), 1(AB), and 2(BB). Doubles between 0 and 2 are also allowed\n"); return(0);
						}
					}

				if (fscanf (fp, "%lf", &temp) > 0) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
				fclose (fp);
			}
			else
			{
				jj=CountRowsColumns(&N, &P, Genotypefile);if(jj)return(0);
				X = (Xstruct*) calloc ( P, sizeof (Xstruct) );

				fp = fopen( Genotypefile, "r" );
				for (locus=0; locus<P; locus++)
				{
					X[locus].covariates = (double*) calloc (N, sizeof(double));
				}

				for (record=0; record<N; record++)
					for (locus=0; locus<P; locus++)
					{
						if (fscanf (fp, "%lf", X[locus].covariates+record) == -1) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
						if (X[locus].covariates[record]<0.0||X[locus].covariates[record]>2.0) 
						{
							fprintf (stderr,"Genotypes should be coded as 0(AA), 1(AB), and 2(BB). Doubles between 0 and 2 are also allowed\n"); return(0);
						}
					}

				if (fscanf (fp, "%lf", &temp) > 0) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
				fclose (fp);
			}
			break;
		case 2: /* PED file */
			jj=CountRowsColumns(&N, &P, Genotypefile);if(jj)return(0);
			P=(P-6)/2;

			X = (Xstruct*) calloc (P, sizeof (Xstruct) );
			AlleleA = (int*) calloc (P, sizeof(int));
			for (locus=0; locus<P; locus++){ X[locus].covariates = (double*) calloc (N, sizeof(double));}

			fp = fopen( Genotypefile, "r" );
			for (record=0; record<N; record++)
			{
				fscanf (fp, "%*s%*s%*s%*s%*s%s", Chartemp);
				for (locus=0; locus<P; locus++)
				{	
					for(ii=0; ii<2; ii++)
					{
						if (fscanf (fp, "%s", Chartemp) == -1) { fprintf (stderr,"Incorrect ped file format\n"); return(0);}
						allele = (int)(unsigned char)Chartemp[0];
						if (allele == 48) { fprintf (stderr, "Incorrect ped file format\n"); return(0);}
						if (record==0&&ii==0) {
							AlleleA[locus] = allele;
						}
						if (allele!=AlleleA[locus])
							X[locus].covariates[record] += 1.0;
					}
				}
			}
			if (fscanf (fp, "%s", Chartemp) > 0) { fprintf (stderr,"Incorrect ped file format\n"); return(0);}
			fclose (fp);
			free(AlleleA);
			break;
		case 3: /* Dosage file */
			jj=CountRowsColumns(&P, &N, Genotypefile);if(jj)return(0);
			P--;
			N=N-3;

			fp = fopen( Genotypefile, "r" );
			X = (Xstruct*) calloc ( P, sizeof (Xstruct) );
			for (locus=0; locus<P; locus++)
			{
				X[locus].covariates = (double*) calloc (N, sizeof(double));
			}

			fgets(v,sizeof(v),fp); /*skip the header*/
			for (locus=0; locus<P; locus++)
			{
				fscanf (fp, "%*s%*s%*s"); /*skip the first 3 columns*/
				for (record=0; record<N; record++)
				{
					if (fscanf (fp, "%lf", X[locus].covariates+record) == -1) { fprintf (stderr,"Incorrect dosage file format\n"); return(0);}
				}
			}

			if (fscanf (fp, "%lf", &temp) > 0) { fprintf (stderr,"Incorrect dosage file format\n"); return(0);}
			fclose (fp);
			break;
	}

	if(Printinfo) { 
		printf("# individuals   : %d\n", N);
		printf("# markers       : %d\n", P);
		printf("\n");		
	}

	
	/* Sum of variances of variables */
	if(Xtype){
		for(locus=0; locus<P; locus++)
		{
			SumVar += Variance(N,X[locus].covariates);
		}
	}else{
		for(locus=0; locus<P; locus++)
		{
			Af=0.0;
			for(record=0; record<N; record++){Af+=X[locus].covariates[record];}
			Af/=((double)N*2.0);
			SumVar += 2.0*(1.0-Af)*Af;
		}
		SumVar*=(1.0+Incoef);
	}

	/* Calculate hyperparameter values */
	Nset = 1;
	if(Methodcode==1){Nset=Nvalue[4]*Nvalue[0];}
	if(Methodcode==2){Nset=Nvalue[4]*Nvalue[0]*Nvalue[1]*Nvalue[2];}
	if(Methodcode==3||Methodcode==4||Methodcode==7){Nset=Nvalue[4]*Nvalue[3];}
	if(Methodcode==5||Methodcode==6){Nset=Nvalue[4]*Nvalue[3]*Nvalue[5];}

	H = (Hstruct*) calloc (Nset, sizeof(Hstruct));
	Print = (double*) calloc(Nset, sizeof(double));
	if(Methodcode==1){

		Neach=Nvalue[4];
		Ntimes=1;
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[0]; ii++)
				for(each=0; each<Neach; each++){
					H[set].deltaShape = GivenValues[0].deltaShape[ii];
					set++;
				}

		Neach=1;
		Ntimes=Nvalue[0];
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[4]; ii++)
				for(each=0; each<Neach; each++){
					H[set].deltaRate = H[set].deltaShape/(2.0*GivenValues[0].Pi[ii]*SumVar*(1.0/Mvar-1.0));
					Print[set] = GivenValues[0].Pi[ii];
					set++;
				}
	}
	if(Methodcode==2){
		Neach=Nvalue[1]*Nvalue[2]*Nvalue[4];
		Ntimes=1;
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[0]; ii++)
				for(each=0; each<Neach; each++){
					H[set].deltaShape = GivenValues[0].deltaShape[ii];
					set++;
				}
		Neach=Nvalue[2]*Nvalue[4];
		Ntimes=Nvalue[0];
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[1]; ii++)
				for(each=0; each<Neach; each++){
					H[set].deltaRate = GivenValues[0].deltaRate[ii];
					set++;
				}
		Neach=Nvalue[4];
		Ntimes=Nvalue[0]*Nvalue[1];
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[2]; ii++)
				for(each=0; each<Neach; each++){
					H[set].etaShape = GivenValues[0].etaShape[ii];
					set++;
				}
		Neach=1;
		Ntimes=Nvalue[0]*Nvalue[1]*Nvalue[2];
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[4]; ii++)
				for(each=0; each<Neach; each++){
					H[set].etaRate = H[set].deltaShape/H[set].deltaRate*H[set].etaShape/(2.0*GivenValues[0].Pi[ii]*SumVar*(1.0/Mvar-1.0));
					Print[set] = GivenValues[0].Pi[ii];
					set++;
				}
	}
	if(Methodcode==3||Methodcode==4||Methodcode==7){
		Neach=Nvalue[4];
		Ntimes=1;
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[3]; ii++)
				for(each=0; each<Neach; each++){
					H[set].v = GivenValues[0].v[ii];
					set++;
				}
		Neach=1;
		Ntimes=Nvalue[3];
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[4]; ii++)
				for(each=0; each<Neach; each++){
					H[set].Pi = GivenValues[0].Pi[ii];
					H[set].S2 = (H[set].v-2.0)/H[set].v*Mvar/(H[set].Pi*SumVar);
					set++;
				}
	}
	if(Methodcode==5||Methodcode==6){
		Neach=Nvalue[4]*Nvalue[5];
		Ntimes=1;
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[3]; ii++)
				for(each=0; each<Neach; each++){
					H[set].v = GivenValues[0].v[ii];
					set++;
				}
		Neach=1;
		Ntimes=Nvalue[3]*Nvalue[5];
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[4]; ii++)
				for(each=0; each<Neach; each++){
					H[set].Pi = GivenValues[0].Pi[ii];
					set++;
				}
		Neach=Nvalue[4];
		Ntimes=Nvalue[3];
		for(time=0, set=0; time<Ntimes; time++)
			for(ii=0; ii<Nvalue[5]; ii++)
				for(each=0; each<Neach; each++){
					H[set].c = (1.0-GivenValues[0].A[ii])/GivenValues[0].A[ii]*H[set].Pi/(1.0-H[set].Pi);
					H[set].S2 = (H[set].v-2.0)/H[set].v*Mvar/((H[set].Pi+H[set].c*(1-H[set].Pi))*SumVar);
					Print[set] = GivenValues[0].A[ii];
					set++;
				}
	}

	/* Print results */
	if(Printinfo)
	{
		printf("Hyperparameters\n");
		if(Methodcode==1){
			printf("Phi Omega (Kappa)\n");
			for(set=0; set<Nset; set++){
				printf("%f %f (%f)\n", H[set].deltaShape, H[set].deltaRate, Print[set]);
			}
		}
		if(Methodcode==2){
			printf("Phi Omega Psi Theta (Kappa)\n");
			for(set=0; set<Nset; set++){
				printf("%f %f %f %f (%f)\n", H[set].deltaShape, H[set].deltaRate, H[set].etaShape, H[set].etaRate, Print[set]);
			}
		}
		if(Methodcode==3||Methodcode==4||Methodcode==7){
			printf("Nu S2 Kappa\n");
			for(set=0; set<Nset; set++){
				printf("%f %f %f\n", H[set].v, H[set].S2, H[set].Pi);
			}
		}
		if(Methodcode==5||Methodcode==6){
			printf("c Nu S2 Kappa (A)\n");
			for(set=0; set<Nset; set++){
				printf("%f %f %f %f (%f)\n", H[set].c, H[set].v, H[set].S2, H[set].Pi, Print[set]);
			}
		}
	}

	/* free memory */
	for (locus=0; locus<P; locus++) free(X[locus].covariates);
	free (X);
	free (Nvalue);
	free (GivenValues[0].Pi);
	switch(Methodcode){
		case 1:
			free(GivenValues[0].deltaShape);
			break;
		case 2:
			free(GivenValues[0].deltaShape);free(GivenValues[0].deltaRate);free(GivenValues[0].etaShape);
			break;
		case 3:
			free(GivenValues[0].v);
			break;
		case 4:
			free(GivenValues[0].v);
			break;
		case 5:
			free(GivenValues[0].v);free(GivenValues[0].A);
			break;
		case 6:
			free(GivenValues[0].v);free(GivenValues[0].A);
			break;
		case 7:
			free(GivenValues[0].v);
			break;
	}
	free(GivenValues);
	free(H);
	free(Print);
	return(0);
}
