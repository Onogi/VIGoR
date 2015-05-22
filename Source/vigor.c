#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "MT.h"
#include "vigor.h"

/*	VIGoR (Variational Bayesian Inference for Genome-wide Regression)

	Seven regression models are available:
		1.BL		: Bayeian Lasso (Park and Casella 2008)
		2.EBL		: Extended Bayesian Lasso (Mutshinda and Sillanpaa 2010)
		3.wBSR		: Weighted Bayesian Shrinkage Regression (Hayashi and Iwata 2010)
		4.BayesC	: BayesC (Habier et al. 2011)
		5.SSVS		: Stchastic search variable selection (George and McCulloch 1997)
		6.MIX		: Bayesian mixture regression (Luan et al. 2009)
		7.BayesB    : BayesB (Meuwissen et al. 2001)

	Prior distribution types of these models are
		1.Laplace	: BL and EBL
		2.Student t	: The others

	In cross-validation, the following six metrics are available (fixed to 1 in the current version)
		1. MSE
		2. Pearson
		3. NDCG
		4. Log predictive density score (LPDS)
		5. MeanNDCG
		6. Kedall

	The basic linear model is 

		Y = ZA + XGB + e

	where Y is the vector of phenotypic values, Z is the covariate matrix, A is the vector of fixed effect, X is the genotype data, 
	G is the matrix of indicator variable (only wBSR), B is the vector of marekr effects and e is the vector of residuals.

	Priors
		A is proportional to constant
		e ~ Normal (0, 1/Re2)
		Re2 is proportional to 1/Re2

			BL: 
				bj     ~ Normal (0, 1/tau2j)
				tau2j  ~ Inv-Gamma (1, Lambda2/2)
				Lambda2~ Gamma (Phi, Omega) (Lambda2 is denoted as Delta2 in the script)
			EBL:
				bj     ~ Normal (0, 1/tau2j)
				tau2j  ~ Inv-Gamma (1, (Delta2*Eta2j)/2)
				Delta2 ~ Gamma (Phi, Omega)
				Eta2j  ~ Gamma (Psi, Theta)
			wBSR:
				bj     ~ Normal (0, Sigma2j)
				Sigma2j~ Scaled-inv-chi square (Nu, S2)
				Gammaj ~ Bernoulli (Kappa) (Kappa is denoted as Pi in the script)
			BayesC:
				bj     ~ Normal (0, Sigma2) if Rhoj=1, 0 if Rhoj=0
				Sigma2 ~ inv-chi squared (Nu, S2)
				Rhoj   ~ Bernoulli (Kappa) (Rho and Kappa are denoted as Gamma and Pi in the script)
			SSVS:
				bj     ~ Normal (0, Sigma2) if Rhoj=1, Normal (0, c*Sigma2) if Rhoj=0
				Sigma2 ~ Scaled inv-chi square (Nu, S2)
				Rhoj   ~ Bernoulli (Kappa) (Rho and Kappa are denoted as Gamma and Pi in the script)
			MIX:
				bj     ~ Normal (0, Sigma21) if Rhoj=1, Normal (0, Sigma22) if Rhoj=0
				Sigma21~ Scaled inv-chi square (Nu, S2)
				Sigma22~ Scaled inv-chi square (Nu, cS2)
				Rhoj   ~ Bernoulli (Kappa) (Rho and Kappa are denoted as Gamma and Pi in the script)
			BayesB:
				bj     ~ Normal (0, Sigma2j) if Rhoj=1, 0 if Rhoj=0
				Sigma2j~ inv-chi squared (Nu, S2)
				Rhoj   ~ Bernoulli (Kappa) (Rho and Kappa are denoted as Gamma and Pi in the script)

	The numbers of hyperparameters (Nh) are
			BL			:	Nh=2 (Phi and Omega)
			EBL			:	Nh=4 (Phi,Omega, Psi, and Theta)
			wBSR		:	Nh=3 (Nu, S2, and Kappa)
			BayesC		:	Nh=3 (Nu, S2, and Kappa)
			SSVS		:	Nh=4 (c, Nu, S2, and Kappa)
			MIX			:	Nh=4 (c, Nu, S2, and Kappa)
			BayesB		:	Nh=3 (Nu, S2, and Kappa)

	Input files are
		Phenotypefile	: (N (number of individuals) + 1(header)) x Ntraits (number of traits) matrix. The header contains the trait names. Both doubles and integers are allowed.
						Missing values should be specified as 'NA'. 'NA' is automatically converted to Missingvalue which is set to 99999999 by default.
						The maximum length of the trait name is 100.
		Genotypefile	: N x P (number of markers) matrix. Both doubles and integers are allowed. Missing values are not allowed.
		Covariatefile	: N x F (number of fixed effects) matrix. Both doubles and integers are allowed. 
						 The first column is recognized as the intercept.
		Partitionfile	: Partition of individuals in cross-validation (CV). Used only when AnalysisType = -9.
						 TESTING individuals are specified in the file.
						 The file format is, for example,

						 1 1 3 2
						 2 4 4 4
						 3 5 5 -9

						 The fold number is 4, and the number of testing individuals at each fold is 3.
						 In the first fold, individuals 1, 2 and 3 are removed from training and tested,
						 in the second, individuals 1, 4 and 5 are removed.
						 In the last fold, only individuals 2 and 4 are removed; -9 is used to fill the gap.

	The maximum length of input file names is 400.
	The maximum length of the folder name is 100.

	In the current implementation, convergence is checked as
		
		|E[parameters'] - E[parameters]|^2 / |E[parameters']|^2 < 10^-(Threshold)

	where E[parameters'] are expectaions of newly updated parameters, and || is the Euclidean norm.

	The pdf manual of VIGoR is provided at https://github.com/Onogi/VIGoR
	
	==============================================================================================================
	Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata

	Released under the MIT license
	http://opensource.org/licenses/mit-license.php
	==============================================================================================================
*/


int main (int argc, char *argv[])
{

	/* Values defined from arguments, options, or parameter files */
	char	Phenotypefile[401]="";		/* Phenotype file name */
	char	Genotypefile[401]="";		/* Genotype file name */
	char	Covariatefile[401]="";		/* Covariate file. When '0', only overall mean is added to the model automatically */
	char	Partitionfile[401]="";		/* Partition of individuals for CV. Used when AnalysisType = -9 */
	int		AnalysisType;				/* If 0, CV is not performed (only model fitting). If 1, CV is performed with each set of hyperparameters and model fitting is conducted with the best set. 
											If -1 or n>1, leave-one-out CV or n-fold CV are performed. If -9, the assignments of individuals are read from the 'Partitionfile'
											When -1, n or -9, nested CV is performed at each CV procedure to select the best set from the hyperparameter sets */
	int     Methodcode;					/* 1,BL; 2,EBL; 3,wBSR; 4,BayesC; 5,SSVS; 6,MIX; 7,BayesB */
	int		Nperm;						/* Number of permutations. 0 if not performed. Effective when Analysis type is 0.*/
	int		Ntraits;					/* Total number of traits included in the Phenotypefile */
	int		K;							/* Trait position to be analysed in the Phenotypefile. If 0, all traits are analyzed */
	int		P;							/* Number of markers */
	int		F;							/* Number of fixed effects. Ignored when Covariatefile = '0' */
	int		N;							/* Number of individuals */
	double	Threshold;					/* Threshold to quit iterations. The actural threshold value is 10^-Threshold */
	int		Maxiteration;				/* Max number of iterations */
	int		CVFoldTuning;				/* Number of folds of the CV for hyperparameter tuning. Used only when AnalysisType != 0 and multiple hyperparameter value sets are given */
	int		UsePNmatrix;				/* Use P x N matrix for Genotypefile (1) or not (0); */
	int		Printinfo;					/* Whether print out run information (1) or not (0); */
	int		RandomIni;					/* Whether randomize initial values (1) or not (0); */
	
	/* Values defined by the developer */
	int		MPP=1;						/* MPP determines the metrics for predictive performance. 1: MSE, 2: Pearson, 3: NDCG@NDCGrank, 4: LPDS, 5: MeanNDCG@NDCGrank, 6: Kendall */
	int		NDCGrank=1;					/* Rank of NDCG */
	int		NsampleLPDS=1000;			/* Number of samples when LPDS is calculated */
	int		Mnset=1000;					/* Maximum number of hyperparameter sets */
	int		MaxNumHyperpara=4;			/* Maximum number of hyperparameter values */
	int		Algorithm=1;				/* Variational Bayes: 1, MAP estimation with EM: 2 */

	/* Seed used */
	int		usedseed;

	/* Number of arguments */
	int		Narg;

	/* File format. 1: Default, 2: Plink, 3: Beagle */
	int		Fileformat=1;

	/* Inputs and outputs */
	char    Output[301]="", Chartemp[101]="", Method[11]="", AlgorithmLabel[11]="";
	char	v[MaxBufSize]=""; /* for Beagle dosage file */
	FILE    *fp;

	/* for reading PED file */
	int		*AlleleA, allele;

	/* For repeat statement */
	int		ii;
	int		record, locus, set, fold, fold2, trait;

	/* Generic use */
	double	temp; 
	int	jj;

	/* specify prior distribution types (Laplace or Student t)*/
	int		Priortype, CondResidual;

	/* Phenotype */
	Ystruct	*Y, *Ycv, *Ycv2, *Yperm;

	/* Genotypes and fixed effects (or covariates) */
	Xstruct	*X, *Q, *Xperm, *Qperm;

	/* For cross validation */
	Astruct	*A, *A2;
	int		Fold, Maxntest, *BestSet, *Order;
	int		Nmissing, Nobserved;
	double	*BestValue, *ValueList;
	double	ValueAtFold;

	/* Delta2 and  residual variance*/
	double	expDelta2=0.0;
	double	Tau0[2];

	/* Monitor iterations */
	double	*Rmonitor, *LBmonitor;

	/* Predicted and fitted values */
	double	*PredictedValues, *PredictedValues2, *BreedingValues;

	/* Threshold */
	double	Thresholdvalue;
	int	ThresholdDetermined=0;

	/* Calculation Time */
	double	FittingTime;
	clock_t	CVstart, CVend;

	/* For specification of hyperparameters */
	Hstruct	*H;
	int		Nh, Nset;

	/* For the progress bar */
	int		Nruns, Currentrun;

	/* Total number of tested samples in the assignment file */
	int		TotalNtest;
	int		CumNtest;

	/* Number of traits to be analyzed*/
	int		NtraitsAnalyzed;
	int		ntrait;

	/* Permutation test */
	double	PermutationThreshold[4];
	
	/* output help */
	if(getopt(argc, argv, "a:c:p:u:n:k:s:i:v:toqrh")=='h')
	{
		printf("Usage\n\n"); 
		printf(" vigor <Phenotype (response variable) file> <Genotype file>\n");
		printf("       <Method BL|EBL|wBSR|BayesB|BayesC|SSVS|MIX>\n");
		printf("       <hyperparameter values> [options]\n\n");
		printf("Phenotype file consists of a header for trait names, and a subsquent N (number of individuals) x Ntraits (number of traits) matrix. ");
		printf("Genotype file consists of a N x P (number of markers) matrix. P x N matrix can be used by using the -o option.\n\n");
		printf("NOTE: the maximum size of one row is 32^4=1MB. Use P x N matrix for genotype file when P is large and buffer overflow occurs\n\n");
		printf("When a PED file of Plink is input, the usage is\n");
		printf(" vigor <PED file>\n");
		printf("       <Method BL|EBL|wBSR|BayesB|BayesC|SSVS|MIX>\n");
		printf("       <hyperparameter values> [options]\n\n");
		printf("The first 6 columns of the PED files are mandantory. Missing genotypes are not allowed\n");
		printf("In the PED files, two alleles of a genotype are stored in two columns. Alleles can be encoded by any single character\n\n");
		printf("The regression methods are\n");
		printf("  BL    : Bayesian lasso\n");
		printf("  EBL   : Extended Bayesian lasso\n");
		printf("  wBSR  : Weighted Bayesian shrinkage regression\n");
		printf("  BayesB: BayesB\n");
		printf("  BayesC: BayesC\n");
		printf("  SSVS  : Stochastic search variable selection\n");
		printf("  MIX   : Bayesian mixture model\n\n");
		printf("All the hyperparameter values required should be given according to the following orders\n");
		printf("  BL    : Phi Omega\n");
		printf("  EBL   : Phi Omega Psi Theta\n");
		printf("  wBSR  : Nu S2 Kappa\n");
		printf("  BayesB: Nu S2 Kappa\n");
		printf("  BayesC: Nu S2 Kappa\n");
		printf("  SSVS  : c, Nu S2 Kappa\n");
		printf("  MIX   : c, Nu S2 Kappa\n\n");
		printf("See the pdf manual of VIGoR for the explanation of these hyperparameters\n\n");
		printf("Options that take arguments\n");
		printf("-a: Covariate file name (including predictor variables other than the markers)\n");
		printf("    The intercept is automatically added to the model when this option is not specified\n");
		printf("-c: Conduct the cross-validation (CV) option\n");
		printf("   -1      : leave one out CV\n");
		printf("   -9      : CV based on Partition file\n");
		printf("    n (>=2): n-fold CV\n");
		printf("-p: Partition file name which defines the partition of individuals in CV.\n");
		printf("    Used when the argument of the -c option is -9\n");
		printf("-u: Fold numbers of CV for hyperaparameters tuning.\n");
		printf("    Used when the -c and -t options are specified, and multiple hyperparameter sets are given\n");
		printf("    The default value is 5\n");
		printf("-n: Numer of permutations. Used only when the fitting function is selected. The default is 0\n");
		printf("-k: Column number (trait) of Phenotype file to be analyzed.\n");
		printf("    The default value is 1. When 0 is input, all the traits are analyzed sequentially.\n");
		printf("-s: Convergence threshold. Larger values indicate stricter thresholds. The default is 2+log10(P) (P is the number of markers)\n");
		printf("-i: Maximum number of iterations. The default value is 1000\n");
		printf("-v: Additional hyperparameter value sets\n");
		printf("    For example, when the regression method is BL, \n");
		printf("        -v 1 0.005\n");
		printf("    indicates that Phi = 1 and Omega = 0.005. \n");
		printf("    Multiple hyperparameter sets can be specified by using this option repeatedly.\n\n");
		printf("Options\n");
		printf("-t: Execute model fitting after hyperparameter tuning\n");
		printf("-o: Use P x N genotype file matrix. Recommended when P >> N. \n");
		printf("-q: Quiet. Information is not output\n");
		printf("-r: Randomize initial values\n");
		printf("-h: Help\n\n");
		return(0);
	}

	Narg=1;
	strcpy(Phenotypefile, argv[Narg]);Narg++; //Narg==2
	if(strcmp(Phenotypefile+(strlen(Phenotypefile)-4),".ped")==0)
	{
		Fileformat=2; /* PED file for PLINK */
		strcpy(Genotypefile, "None");
	} 
	else
	{
		strcpy(Genotypefile, argv[Narg]);Narg++; //Narg==3
		if(strcmp(Genotypefile+(strlen(Genotypefile)-5),".dose")==0) Fileformat=3; /* Beagle dosage file */
	}

	strcpy(Method,argv[Narg]);Narg++; //Narg==3 or 4

	if(strcmp(Method,"BL")==0)		{Methodcode=1; Nh = 2; Priortype = 1; CondResidual=1;}
	if(strcmp(Method,"EBL")==0)		{Methodcode=2; Nh = 4; Priortype = 1; CondResidual=1;}
	if(strcmp(Method,"wBSR")==0)	{Methodcode=3; Nh = 3; Priortype = 2; CondResidual=0;}
	if(strcmp(Method,"BayesB")==0)	{Methodcode=7; Nh = 3; Priortype = 2; CondResidual=0;}
	if(strcmp(Method,"BayesC")==0)	{Methodcode=4; Nh = 3; Priortype = 2; CondResidual=0;}
	if(strcmp(Method,"SSVS")==0)	{Methodcode=5; Nh = 4; Priortype = 2; CondResidual=0;}
	if(strcmp(Method,"MIX")==0)		{Methodcode=6; Nh = 4; Priortype = 2; CondResidual=0;}
	if(strcmp(Method,"BL")!=0&&strcmp(Method,"EBL")!=0&&strcmp(Method,"wBSR")!=0&&strcmp(Method,"BayesB")!=0&&strcmp(Method,"BayesC")!=0&&strcmp(Method,"SSVS")!=0&&strcmp(Method,"MIX")!=0) 
	{fprintf(stderr,"Number of arguments or method specification\n");return(0);}

	H = (Hstruct*) calloc (Mnset, sizeof(Hstruct));
	switch(Methodcode){
		case 1:	/*BL*/
			H[0].deltaShape = atof(argv[Narg]); H[0].deltaRate = atof(argv[Narg+1]);
			if(H[0].deltaShape<=0.0||H[0].deltaRate<=0.0) {fprintf(stderr,"Hyperparameters should be positive\n"); return(0);}
		break;
		case 2: /*EBL*/
			H[0].deltaShape = atof(argv[Narg]); H[0].deltaRate = atof(argv[Narg+1]); H[0].etaShape = atof(argv[Narg+2]); H[0].etaRate = atof(argv[Narg+3]);
			if(H[0].deltaShape<=0.0||H[0].deltaRate<=0.0||H[0].etaShape<=0.0||H[0].etaRate<=0.0) {fprintf(stderr,"Hyperparameters should be positive\n"); return(0);}
		break;
		case 3: /*wBSR*/
			H[0].v = atof(argv[Narg]); H[0].S2 = atof(argv[Narg+1]); H[0].Pi = atof(argv[Narg+2]);
			if(H[0].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
			if(H[0].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
			if(H[0].Pi<=0.0||H[0].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
		break;
		case 4: /*BayesC*/
			H[0].v = atof(argv[Narg]); H[0].S2 = atof(argv[Narg+1]); H[0].Pi = atof(argv[Narg+2]);
			if(H[0].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
			if(H[0].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
			if(H[0].Pi<=0.0||H[0].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
		break;
		case 5: /*SSVS*/
			H[0].c = atof(argv[Narg]); H[0].v = atof(argv[Narg+1]); H[0].S2 = atof(argv[Narg+2]); H[0].Pi = atof(argv[Narg+3]);
			if(H[0].c<=0.0) {fprintf(stderr,"c should be > 0\n"); return(0);}
			if(H[0].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
			if(H[0].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
			if(H[0].Pi<=0.0||H[0].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
		break;
		case 6: /*MIX*/
			H[0].c = atof(argv[Narg]); H[0].v = atof(argv[Narg+1]); H[0].S2 = atof(argv[Narg+2]); H[0].Pi = atof(argv[Narg+3]);
			if(H[0].c<=0.0) {fprintf(stderr,"c should be > 0\n"); return(0);}
			if(H[0].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
			if(H[0].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
			if(H[0].Pi<=0.0||H[0].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
		break;
		case 7: /*BayesB*/
			H[0].v = atof(argv[Narg]); H[0].S2 = atof(argv[Narg+1]); H[0].Pi = atof(argv[Narg+2]);
			if(H[0].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
			if(H[0].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
			if(H[0].Pi<=0.0||H[0].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
		break;
	}
	Nset = 1;
	Narg += Nh;

	/*Default values*/
	strcpy(Covariatefile, "0");
	strcpy(Partitionfile, "0");
	AnalysisType=0; /* fitting */
	CVFoldTuning=5;
	Nperm=0;
	K=1;
	Maxiteration=1000;
	strcpy(AlgorithmLabel,"VB");
	UsePNmatrix=0;
	Printinfo=1;
	RandomIni=0;

	/*Read options*/
	optind=Narg;
	while((ii=getopt(argc, argv, "a:c:p:u:n:k:s:i:v:toqrh"))!=-1){
		switch(ii){
		case 'a':
			strcpy(Covariatefile, optarg);
			break;
		case 'c':
			AnalysisType=atoi(optarg);
			if(AnalysisType<-9||(AnalysisType>-8&&AnalysisType<2)) {fprintf(stderr,"Misspecification in -c\n"); return(0);}
			break;
		case 'p':
			strcpy(Partitionfile, optarg);
			break;
		case 'u':
			CVFoldTuning=atoi(optarg);
			if(CVFoldTuning<=1) {fprintf(stderr,"Misspecification in -u\n"); return(0);}
			break;				
		case 'n':
			Nperm=atoi(optarg);
			if(Nperm<0) {fprintf(stderr,"Misspecification in -n\n"); return(0);}
			break;
		case 'k':
			K=atoi(optarg);
			if(K<0) {fprintf (stderr,"Misspecification in -k\n"); return(0);}
			break;
		case 's':
			Threshold=atof(optarg);
			ThresholdDetermined=1;
			break;
		case 'i':
			Maxiteration=atoi(optarg);
			if(Maxiteration<=0) {fprintf(stderr,"Misspecification in -i\n"); return(0);}
			break;
		case 'v':
			switch(Methodcode){
			case 1:
				H[Nset].deltaShape=atof(argv[optind-1]); H[Nset].deltaRate=atof(argv[optind]);
				if(H[Nset].deltaShape<=0.0||H[Nset].deltaRate<=0.0) {fprintf(stderr,"Hyperparameters should be positive\n"); return(0);}
				optind=optind+1;
				break;
			case 2:
				H[Nset].deltaShape=atof(argv[optind-1]); H[Nset].deltaRate=atof(argv[optind]); 
				H[Nset].etaShape=atof(argv[optind+1]);   H[Nset].etaRate=atof(argv[optind+2]);
				if(H[Nset].deltaShape<=0.0||H[Nset].deltaRate<=0.0||H[Nset].etaShape<=0.0||H[Nset].etaRate<=0.0) {fprintf(stderr,"Hyperparameters should be positive\n"); return(0);}
				optind=optind+3;
				break;
			case 3:
				H[Nset].v=atof(argv[optind-1]); H[Nset].S2=atof(argv[optind]); H[Nset].Pi=atof(argv[optind+1]);
				if(H[Nset].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
				if(H[Nset].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
				if(H[Nset].Pi<=0.0||H[Nset].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
				optind=optind+2;
				break;
			case 4:
				H[Nset].v=atof(argv[optind-1]); H[Nset].S2=atof(argv[optind]); H[Nset].Pi=atof(argv[optind+1]);
				if(H[Nset].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
				if(H[Nset].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
				if(H[Nset].Pi<=0.0||H[Nset].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
				optind=optind+2;
				break;
			case 5:
				H[Nset].c=atof(argv[optind-1]); H[Nset].v=atof(argv[optind]); H[Nset].S2=atof(argv[optind+1]); H[Nset].Pi=atof(argv[optind+2]);
				if(H[Nset].c<=0.0) {fprintf(stderr,"c should be > 0\n"); return(0);}
				if(H[Nset].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
				if(H[Nset].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
				if(H[Nset].Pi<=0.0||H[Nset].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
				optind=optind+3;
				break;
			case 6:
				H[Nset].c=atof(argv[optind-1]); H[Nset].v=atof(argv[optind]); H[Nset].S2=atof(argv[optind+1]); H[Nset].Pi=atof(argv[optind+2]);
				if(H[Nset].c<=0.0) {fprintf(stderr,"c should be > 0\n"); return(0);}
				if(H[Nset].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
				if(H[Nset].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
				if(H[Nset].Pi<=0.0||H[Nset].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
				optind=optind+3;
				break;
			case 7:
				H[Nset].v=atof(argv[optind-1]); H[Nset].S2=atof(argv[optind]); H[Nset].Pi=atof(argv[optind+1]);
				if(H[Nset].v<=2.0) {fprintf(stderr,"Nu should be > 2\n"); return(0);}
				if(H[Nset].S2<0.0) {fprintf(stderr,"S2 should be >= 0\n"); return(0);}
				if(H[Nset].Pi<=0.0||H[Nset].Pi>1.0) {fprintf(stderr,"Kappa should be 0<Kappa<=1\n"); return(0);}
				optind=optind+2;
				break;
			}
			Nset++;
			break;
		case 't':
			AnalysisType=1;
			break;
		case 'o':
			UsePNmatrix=1;
			break;
		case 'q':
			Printinfo=0;
			break;
		case 'r':
			RandomIni=1;
			break;
		}
	}

	if(AnalysisType!=0) Nperm=0; if(Nperm<0) Nperm=0;
	if(Methodcode==7&&Algorithm==2) { fprintf(stderr,"BayesB can be conducted only with VB\n"); return(0);}
	if(Methodcode==4&&Algorithm==2) { fprintf(stderr,"BayesC can be conducted only with VB\n"); return(0);}
	if(Nset==1&&AnalysisType==1) AnalysisType=0;
	if(AnalysisType==-9&&(strcmp(Partitionfile,"0")==0)) { fprintf(stderr, "Specify Partition file\n"); return(0);}

	if(Printinfo){
		printf("Phenotype file	: %s\n", Phenotypefile);
		printf("Genotype file   : %s\n", Genotypefile);
		printf("Method          : %s\n", Method);
		printf("Hyperparameter values:\n");
		switch (Methodcode){
		case 1:
			printf(" Phi Omega\n");
			for(set=0; set<Nset; set++){printf(" %f %f\n", H[set].deltaShape, H[set].deltaRate);}
			break;
		case 2:
			printf(" Phi Omega Psi Theta\n");
			for(set=0; set<Nset; set++){printf(" %f %f %f %f\n", H[set].deltaShape, H[set].deltaRate, H[set].etaShape, H[set].etaRate);}
			break;
		case 3:
			printf(" Nu S2 Kappa\n");
			for(set=0; set<Nset; set++){printf(" %f %f %f\n", H[set].v, H[set].S2, H[set].Pi);}
			break;
		case 4:
			printf(" Nu S2 Kappa\n");
			for(set=0; set<Nset; set++){printf(" %f %f %f\n", H[set].v, H[set].S2, H[set].Pi);}
			break;
		case 5:
			printf(" c Nu S2 Kappa\n");
			for(set=0; set<Nset; set++){printf(" %f %f %f %f\n", H[set].c, H[set].v, H[set].S2, H[set].Pi);}
			break;
		case 6:
			printf(" c Nu S2 Kappa\n");
			for(set=0; set<Nset; set++){printf(" %f %f %f %f\n", H[set].c, H[set].v, H[set].S2, H[set].Pi);}
			break;
		case 7:
			printf(" Nu S2 Kappa\n");
			for(set=0; set<Nset; set++){printf(" %f %f %f\n", H[set].v, H[set].S2, H[set].Pi);}
			break;
		}
		printf("\n");
		if(strcmp(Covariatefile, "0")!=0) printf("Covariatefile   : %s\n", Covariatefile);
		printf("Option        : ");
		if(AnalysisType==0) printf("Fitting\n");
		if(AnalysisType==1) printf("Tuning&fitting\n");
		if(AnalysisType==-1) printf("Leave-one-out CV\n");
		if(AnalysisType==-9) printf("CV defined by the file\n");
		if(AnalysisType>1) printf("%d-fold CV\n", AnalysisType);
		if(AnalysisType==-9) printf("Partitionfile   : %s\n", Partitionfile);
		if(AnalysisType!=0&&Nset>1) printf("CVFoldTuning   : %d\n", CVFoldTuning);
		if(AnalysisType==0) printf("Nperm           : %d\n", Nperm);
		printf("K               : %d\n", K);
		printf("Maxiteration    : %d\n", Maxiteration);
	//	printf("Algorithm       : %s\n", AlgorithmLabel);
		printf("UsePNmatrix     : "); if(UsePNmatrix) printf("Yes\n"); else printf("No\n");
		printf("\n");
	}


	/* Open the phenotype and genotype files */
	if(Fileformat==2) /* read Ped file */
	{
		jj=CountRowsColumns(&N, &P, Phenotypefile);if(jj)return(0);
		P=(P-6)/2;
		K=1; Ntraits=1; NtraitsAnalyzed=1; 

		if(Printinfo) { 
			printf("# individuals: %d\n", N);
			printf("# traits     : %d\n", Ntraits);
			printf("# markers    : %d\n\n", P);
		}

		Y = (Ystruct*) calloc ( Ntraits, sizeof(Ystruct));
		Y[0].observations = (double*) calloc (N, sizeof(double));
		strcpy(Y[0].name, "Trait");
		X = (Xstruct*) calloc (P, sizeof (Xstruct) );
		AlleleA = (int*) calloc (P, sizeof(int));
		for (locus=0; locus<P; locus++){ X[locus].covariates = (double*) calloc (N, sizeof(double));}

		fp = fopen( Phenotypefile, "r" );
		for (record=0; record<N; record++)
		{
			fscanf (fp, "%*s%*s%*s%*s%*s%s", Chartemp);
			if (strcmp(Chartemp, "-9") == 0) 
			{ 
				Y[0].observations[record] = Missingvalue;
			} 
			else 
			{ 
				Y[0].observations[record] = atof(Chartemp);
			}
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

		/*--------------------------
		for(record=0; record<N; record++)
		{
			printf("%f ", Y[0].observations[record]);
			for(locus=0; locus<P; locus++)
			{
				printf("%f ", X[locus].covariates[record]);
			}
			printf("\n");
		}
		return(0);
		----------------------------*/
	}
	else
	{
		/* Phenotype file */
		jj=CountRowsColumns(&N, &Ntraits, Phenotypefile);if(jj)return(0);
		N --;
	
		fp = fopen( Phenotypefile, "r" );
		Y = (Ystruct*) calloc ( Ntraits, sizeof(Ystruct));
		for(trait=0; trait<Ntraits; trait++)
		{
			Y[trait].observations = (double*) calloc (N, sizeof(double));
			fscanf (fp, "%s", Y[trait].name);
		}

		for (record=0; record<N; record++)
			for(trait=0; trait<Ntraits; trait++)
			{
				if (fscanf (fp, "%s", Chartemp) == -1) { fprintf (stderr,"Number of elements in the phenotype file\n"); return(0);}
				if (strcmp(Chartemp, "NA") == 0) 
				{ 
					Y[trait].observations[record] = Missingvalue;
				} 
				else 
				{ 
					Y[trait].observations[record] = atof(Chartemp);
				}
			}

		if (fscanf (fp, "%lf", &temp) > 0) { fprintf (stderr,"Number of elements in the phenotype file\n"); return(0);}
		fclose (fp);

		if(Printinfo) { 
			printf("# individuals: %d\n", N);
			printf("# traits     : %d\n", Ntraits);
		}

		if(K==0) {NtraitsAnalyzed=Ntraits;} else {NtraitsAnalyzed=1;}
		if(K>Ntraits) {fprintf(stderr, "K is larger than the number of trais in phenotype file\n"); return(0);}

		/* Genotype file */
		if(Fileformat==1)/* read default format file */
		{
			if(UsePNmatrix)
			{
				jj=CountRowsColumns(&P, &ii, Genotypefile);if(jj)return(0);
				if(ii!=N) {fprintf(stderr,"The number of individuals differs between the genotype and the phenotype files\n"); return(0);}

				fp = fopen( Genotypefile, "r" );
				X = (Xstruct*) calloc ( P, sizeof (Xstruct) );
				for (locus=0; locus<P; locus++)
				{
					X[locus].covariates = (double*) calloc (N, sizeof(double));
				}

				for (locus=0; locus<P; locus++)
					for (record=0; record<N; record++)
					{
						if (fscanf (fp, "%lf", X[locus].covariates+record) == -1) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
					}

				if (fscanf (fp, "%lf", &temp) > 0) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
				fclose (fp);
			}
			else
			{
				jj=CountRowsColumns(&ii, &P, Genotypefile);if(jj)return(0);
				if(ii!=N) {fprintf(stderr,"The number of individuals differs between the genotype and phenotype files\n"); return(0);}

				fp = fopen( Genotypefile, "r" );
				X = (Xstruct*) calloc ( P, sizeof (Xstruct) );
				for (locus=0; locus<P; locus++)
				{
					X[locus].covariates = (double*) calloc (N, sizeof(double));
				}

				for (record=0; record<N; record++)
					for (locus=0; locus<P; locus++)
					{
						if (fscanf (fp, "%lf", X[locus].covariates+record) == -1) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
					}

				if (fscanf (fp, "%lf", &temp) > 0) { fprintf (stderr,"Number of elements in the genotype file\n"); return(0);}
				fclose (fp);
			}
		}

		if(Fileformat==3) /* read Beagle dosage file */
		{
			jj=CountRowsColumns(&P, &ii, Genotypefile);if(jj)return(0);
			P--;
			ii=ii-3;
			if(ii!=N) {fprintf(stderr,"The number of individuals differs between the genotype and phenotype files\n"); return(0);}

			/*---------
			N=ii;
			-----------*/

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

			/*------------------
			printf("P%d N%d\n", P,ii);
			for(locus=401; locus<405; locus++)
			{
				for(record=371; record<381; record++)
				{
					printf("%f ", X[locus].covariates[record]);
				}
				printf("\n");
			}
			return(0);
			--------------------*/
		}

		if(Printinfo) { 
			printf("# markers    : %d\n", P);
		}
	}

	/* Threshold value */
	if(ThresholdDetermined==0){Threshold=2.0+log10((double)P);}
	Thresholdvalue = 1.0/pow(10,Threshold);
	printf("Threshold    : %f\n", Threshold);

	/* Open the covariate file */
	if (strcmp(Covariatefile, "0") == 0) /* only overall mean is included in the model */
	{
		F = 1;
		Q = (Xstruct*) calloc ( F, sizeof(Xstruct));
		Q[0].covariates = (double*) calloc (N, sizeof (double));
		for (record=0; record<N; record++)
			Q[0].covariates[record] = 1.0;
	}
	else	/* Covariates are read from the covariate file */
	{
		jj=CountRowsColumns(&ii, &F, Covariatefile);if(jj)return(0);
		if(ii!=N) {fprintf(stderr,"The number of individuals differs between the covariate and the phenotype/genotype files\n"); return(0);}

		fp = fopen( Covariatefile, "r" );
		Q = (Xstruct*) calloc ( F, sizeof (Xstruct) );
		for (locus=0; locus<F; locus++)
		{
			Q[locus].covariates = (double*) calloc (N, sizeof(double));
		}

		for (record=0; record<N; record++)
			for (locus=0; locus<F; locus++)
			{
				if (fscanf (fp, "%lf", Q[locus].covariates+record) == -1) { fprintf (stderr,"Number of elements in the covariate file\n"); return(0);}
			}

		if (fscanf (fp, "%lf", &temp) > 0) { fprintf (stderr,"Number of elements in the covariate file\n"); return(0);}

		fclose (fp);
		if(Printinfo) { 
			printf("# covariates : %d\n", F);
		}
	}

	/* Objects for permutation tests */
	if(Nperm>0)
	{
		Yperm = (Ystruct*) calloc ( Ntraits, sizeof(Ystruct));
		for(trait=0; trait<Ntraits; trait++)
		{
			Yperm[trait].observations = (double*) calloc (N, sizeof(double));
			memcpy(Yperm[trait].observations,Y[trait].observations,sizeof(double)*N);
			strcpy(Yperm[trait].name, Y[trait].name);
		}

		Xperm = (Xstruct*) calloc (P, sizeof (Xstruct) );
		for(locus=0; locus<P; locus++)
		{
			Xperm[locus].covariates = (double*) calloc (N, sizeof(double));
			memcpy(Xperm[locus].covariates, X[locus].covariates, sizeof(double)*N);
		}

		Qperm = (Xstruct*) calloc ( F, sizeof(Xstruct));
		for (locus=0; locus<F; locus++)
		{
			Qperm[locus].covariates = (double*) calloc (N, sizeof(double));
			memcpy(Qperm[locus].covariates, Q[locus].covariates,sizeof(double)*N);
		}
	}

	/*Cross validation schemes. When AnalysisType>1, random partitions vary among traits*/
	if(AnalysisType==0) 
	{
		if(Nperm>0)
			Nruns = Nset*NtraitsAnalyzed*(1+Nperm);
		else
			Nruns = Nset*NtraitsAnalyzed;
	}
	if(AnalysisType==-1) /* leave one out CV */
	{
		Fold = N;
		A = (Astruct*) calloc (Fold, sizeof(Astruct));
		for (fold=0; fold<Fold; fold++){ A[fold].ntest = 1; A[fold].obstest = (int*) calloc (1, sizeof(int)); A[fold].obstest[0] = fold;}
		if (Nset==1) Nruns = Fold; else Nruns = (Nset * CVFoldTuning + 1) * Fold;
		Nruns *= NtraitsAnalyzed;
		TotalNtest=N;
	}
	if(AnalysisType==1) /* CV for hyperparameter tuning is performed. Parameter values estimated with the best set is returned */
	{
		Fold = 1;
		A = (Astruct*) calloc (Fold, sizeof(Astruct));
		A[0].ntest = 0; A[0].obstest = (int*) calloc (1, sizeof(int)); A[0].obstest[0] = 0;
		Nruns = Nset * CVFoldTuning + 1; Nruns *= NtraitsAnalyzed;
		TotalNtest=N;
	}
	if(AnalysisType>1) /* n-fold CV */
	{
		Fold=AnalysisType;
		if (Nset==1) Nruns = Fold; else Nruns = (Nset * CVFoldTuning + 1) * Fold;
		Nruns *= NtraitsAnalyzed;
		TotalNtest=N;
	}
	if(AnalysisType==-9) /* The assignments of individuals are specified by the Partitionfile */ 
	{
		jj=CountRowsColumns(&Maxntest, &Fold, Partitionfile);if(jj)return(0);

		fp = fopen( Partitionfile, "r" );
		A = (Astruct*) calloc (Fold, sizeof(Astruct));
		for(fold=0; fold<Fold; fold++)
		{
			A[fold].obstest = (int*) calloc (Maxntest, sizeof(int));
			A[fold].ntest = Maxntest;
		}

		for(record=0; record<Maxntest; record++)
			for (fold=0; fold<Fold; fold++)
			{
				if(fscanf(fp, "%d", &(A[fold].obstest[record]))<1) {fprintf(stderr,"Number of elements in the partition file\n"); return(0);}
				if(A[fold].obstest[record] != -9) { A[fold].obstest[record]--;}
			}
		if(fscanf(fp, "%d", &ii) >0) {fprintf(stderr,"Number of elements in the partition file\n"); return(0);}

		for (fold=0; fold<Fold; fold++)
			for(record=0; record<Maxntest; record++)
				if(A[fold].obstest[record] == -9) {A[fold].ntest = record; break;}

		for (fold=0, TotalNtest=0; fold<Fold; fold++)
			TotalNtest += A[fold].ntest;

		fclose(fp);

		if(Printinfo) { 
			printf("# CV folds     : %d\n", N);
			printf("# tested inds per fold: %d\n", Maxntest);
		}

		if (Nset==1) Nruns = Fold; else Nruns = (Nset * CVFoldTuning + 1) * Fold;
		Nruns *= NtraitsAnalyzed;
	}


	for(ntrait=0, Currentrun = 1; ntrait<NtraitsAnalyzed; ntrait++)
	{
		/* Initialize the seed */
		usedseed = (int) time ( NULL );
		init_genrand(usedseed);

		/* Target trait */
		if(K==0) {trait = ntrait;} else {trait = K-1;}

		/* Count missing records */
		for(record=0, Nmissing=0; record<N; record++)
			if(Y[trait].observations[record]==Missingvalue) Nmissing++;
		Nobserved = N - Nmissing;

		if(Printinfo) printf ("# non-missing records : %d (# missing: %d)\n\n", Nobserved, Nmissing);

		if (AnalysisType==0)
		{
			/* Perform model fitting with each hyperparameter set */
			for (set=0; set<Nset; set++)
			{
				/* Model fitting */
				FittingTime = 0.0;
				LBmonitor = (double*) calloc (Maxiteration, sizeof(double));
				Rmonitor =	(double*) calloc (Maxiteration, sizeof(double));
				if(Printinfo) printf ("====Progress:%d/%d trait:%s set:%d ====\n", Currentrun, Nruns, Y[trait].name, set+1);
				GenomeWideRegression (Algorithm, Priortype, Methodcode, CondResidual, P, F, N, Thresholdvalue, Maxiteration, Y+trait, X, Q, H+set, &expDelta2, Tau0, LBmonitor, Rmonitor, &FittingTime, Printinfo, RandomIni);
				Currentrun ++;

				/* Fitted values */
				PredictedValues = (double*) calloc(N, sizeof(double));
				BreedingValues = (double*) calloc(N, sizeof(double));
				for(record=0; record<N; record++)
				{
					PredictedValues [record] = PredictYhat (Q, F, X, P, record, Methodcode);
					BreedingValues  [record] = PredictBV (X, P, record, Methodcode);
				}

				/* Permutation test */
				if(Nperm>0)
					PermutationTestMax(Nperm, PermutationThreshold, Algorithm, Priortype, Methodcode, CondResidual, P, F, N, Thresholdvalue, Maxiteration, Yperm+trait, Xperm, Qperm, H+set, &Currentrun, Nruns, set, Printinfo, RandomIni);

				/* Output */
				ValueList = (double*) calloc (1, sizeof(double));
				WriteResults (Algorithm, Methodcode, Phenotypefile, Genotypefile, Covariatefile, Partitionfile, AnalysisType, Ntraits, P, F, N, Threshold, Maxiteration, usedseed,
					Y, H, X, Q, Rmonitor, LBmonitor, expDelta2, trait, set, CVFoldTuning, Nset, ValueList, PredictedValues, BreedingValues, FittingTime, MPP, NDCGrank, Nperm, PermutationThreshold, RandomIni);

				/* free */
				free(LBmonitor); free(Rmonitor); free(ValueList); free(PredictedValues); free(BreedingValues);

			}	/* set */
		}
		else
		{
			/* When AnalysisType is 1, CV is performed with each hyperparameter set. Then, model fitting is performed with the set that scores the leaset MSE */

			/* When AnalysisType is -9, -1, or n>1, all hyperparameter sets are tested in each Cross-validation procedure. 
			Then the set that scores the least MSE is used for prediction.*/

			/* Randomly partition samples for n-fold cross validation */
			if(AnalysisType>1)
			{
				A = (Astruct*) calloc (Fold, sizeof(Astruct));
				Order = (int*) calloc (Nobserved, sizeof(int));
				for(record=0, ii=0; record<N; record++)
					if(Y[trait].observations[record]!=Missingvalue) {Order[ii] = record; ii++;}

				/* Randomize the order */
				Shuffle (Order, Nobserved, Nobserved);

				if(Nobserved%Fold==0) { Maxntest = Nobserved/Fold; } else { Maxntest = Nobserved/Fold+1;}

				for (fold=0; fold<Fold; fold++)
				{
					if (fold==(Fold-1)) { A[fold].ntest = Nobserved-Maxntest*fold; } else {A[fold].ntest = Maxntest;}
					A[fold].obstest = (int*) calloc (A[fold].ntest, sizeof(int));
					for (record=0; record<A[fold].ntest; record++)
						A[fold].obstest[record] = Order[record + Maxntest*fold];
				}

				/* Output */
			//	sprintf ( Output, "%s_%s_%s.partition", Method, AlgorithmLabel, Y[trait].name );
				sprintf ( Output, "%s_%s.partition", Method, Y[trait].name );

				fp = fopen ( Output, "w" );
				for (record=0; record<Maxntest; record++)
				{
					for (fold=0; fold<Fold; fold++)
						if (record<A[fold].ntest){fprintf (fp, "%d ", A[fold].obstest[record]+1);} else {fprintf (fp, "-9 ");}
					fprintf (fp, "\n");
				}

				fclose (fp);
				free (Order);
			}

			PredictedValues = (double*) calloc (TotalNtest, sizeof(double));
			BreedingValues = (double*) calloc (TotalNtest, sizeof(double));
			CumNtest = 0;
			BestValue = (double*) calloc (Fold, sizeof(double));
			BestSet = (int*) calloc (Fold, sizeof(int));
			if (Fold==1) {ValueList = (double*) calloc (Nset, sizeof(double));}
		
			for (fold=0, CVstart=clock(); fold<Fold; fold++)
			{
				Ycv = (Ystruct*) calloc (1, sizeof(Ystruct));
				Ycv[0].observations = (double*) calloc (N, sizeof(double));
				memcpy (Ycv[0].observations, Y[trait].observations, sizeof(double)*N);
				for (record=0; record<A[fold].ntest; record++)
					Ycv[0].observations[A[fold].obstest[record]] = Missingvalue;

				if(Nset>1) {

					/* Count all missing values */
					for (record=0, Nmissing=0; record<N; record++)
						if (Ycv[0].observations[record]==Missingvalue) Nmissing++;
					Nobserved = N - Nmissing;

					/* Partition of records */
					Order = (int*) calloc (Nobserved, sizeof(int));
					for(record=0, ii=0; record<N; record++) 
						if(Ycv[0].observations[record]!=Missingvalue) { Order[ii] = record; ii++;}

					/* Randomize the order */
					Shuffle (Order, Nobserved, Nobserved);

					if(Nobserved%CVFoldTuning==0) { Maxntest = Nobserved/CVFoldTuning; } else { Maxntest = Nobserved/CVFoldTuning+1;}
					A2 = (Astruct*) calloc (CVFoldTuning, sizeof(Astruct));
					for (fold2=0; fold2<CVFoldTuning; fold2++)
					{
						if (fold2==(CVFoldTuning-1)) { A2[fold2].ntest = Nobserved-Maxntest*fold2; } else {A2[fold2].ntest = Maxntest;}
						A2[fold2].obstest = (int*) calloc (A2[fold2].ntest, sizeof(int));
						for (record=0; record<A2[fold2].ntest; record++)
							A2[fold2].obstest[record] = Order[record + Maxntest*fold2];
					}
					free(Order);

					switch(MPP)
					{
						case 1: BestValue[fold]= 1.0e+100;	break;
						case 2: BestValue[fold]= 0.0;		break;
						case 3: BestValue[fold]=-1.0e+100;	break;
						case 4: BestValue[fold]=-1.0e+100;	break;
						case 5: BestValue[fold]=-1.0e+100;	break;
						case 6: BestValue[fold]= 0.0;		break;
					}

					for (set=0; set<Nset; set++)
					{
						PredictedValues2 = (double*) calloc (N, sizeof(double));

						ValueAtFold=0.0;
						for (fold2=0; fold2<CVFoldTuning; fold2++)
						{
							/* Create response variables */
							Ycv2 = (Ystruct*) calloc (1, sizeof(Ystruct));

							Ycv2[0].observations = (double*) calloc (N, sizeof(double));
							memcpy (Ycv2[0].observations, Ycv[0].observations, sizeof(double)*N);
							for (record=0; record<A2[fold2].ntest; record++)
								Ycv2[0].observations[A2[fold2].obstest[record]] = Missingvalue;

							/* Model fitting */
							FittingTime = 0.0;
							LBmonitor = (double*) calloc (Maxiteration, sizeof(double));
							Rmonitor =  (double*) calloc (Maxiteration, sizeof(double));
							if(Printinfo) printf ("----Progress:%d/%d trait:%s fold:%d set:%d fold2:%d ----\n", Currentrun, Nruns, Y[trait].name, fold+1, set+1, fold2+1);
							GenomeWideRegression (Algorithm, Priortype, Methodcode, CondResidual, P, F, N, Thresholdvalue, Maxiteration, Ycv2, X, Q, H+set, &expDelta2, Tau0, LBmonitor, Rmonitor, &FittingTime, Printinfo, RandomIni);
							Currentrun ++;

							if(MPP==4)
							{
								/*LPDS*/
								/* simulate predictive log density */
								ValueAtFold += LPDS(Ycv, Q, F, X, P, A2, fold2, Methodcode, Tau0, A2[fold2].ntest, NsampleLPDS);
							}
							else
							{
								/* MSE, Pearson, NDCG, Kendall*/
								/* Predicted values */
								for ( record=0; record<A2[fold2].ntest; record++ )
									PredictedValues2 [A2[fold2].obstest[record]] = PredictYhat (Q, F, X, P, A2[fold2].obstest[record], Methodcode);
							}
							free(LBmonitor); free(Rmonitor); free(Ycv2[0].observations); free(Ycv2);
						}

						/* Measure predictive performance */
						switch(MPP){
							case 1: /* MSE */
								ValueAtFold = MSE (PredictedValues2, Ycv[0].observations, N);
								if(ValueAtFold < BestValue [fold]) { BestValue [fold] = ValueAtFold; BestSet [fold] = set;}
								if (Fold==1) ValueList[set] = ValueAtFold;
								break;
							case 2: /* Pearson */
								ValueAtFold = Pearson (PredictedValues2, Ycv[0].observations, N);
								if(ValueAtFold > BestValue [fold]) { BestValue [fold] = ValueAtFold; BestSet [fold] = set;}
								if (Fold==1) ValueList[set] = ValueAtFold;
								break;
							case 3: /* NDCG@NDCGrank */
								ValueAtFold = NDCG (NDCGrank, Ycv[0].observations, PredictedValues2, N);
								if(ValueAtFold > BestValue [fold]) { BestValue [fold] = ValueAtFold; BestSet [fold] = set;}
								if (Fold==1) ValueList[set] = ValueAtFold;
								break;
							case 4: /* LPDS */
								ValueAtFold /= (double)CVFoldTuning;
								if(ValueAtFold > BestValue [fold]) { BestValue [fold] = ValueAtFold; BestSet [fold] = set;}
								if (Fold==1) ValueList[set] = ValueAtFold;
								break;
							case 5: /* MeanNDCG@NDCGrank */
								ValueAtFold = 0.0;
								for(ii=1; ii<=NDCGrank; ii++){ValueAtFold += NDCG (ii,Ycv[0].observations, PredictedValues2, N);}
								ValueAtFold /= (double)NDCGrank;
								if(ValueAtFold > BestValue [fold]) { BestValue [fold] = ValueAtFold; BestSet [fold] = set;}
								if (Fold==1) ValueList[set] = ValueAtFold;
								break;
							case 6: /* Kendall */
								ValueAtFold = Kendall (PredictedValues2, Ycv[0].observations, N);
								if(ValueAtFold > BestValue [fold]) { BestValue [fold] = ValueAtFold; BestSet [fold] = set;}
								if (Fold==1) ValueList[set] = ValueAtFold;
								break;
						}
						free(PredictedValues2);

					}	/* set */

					for (fold2=0; fold2<CVFoldTuning; fold2++) free(A2[fold2].obstest);free (A2);
				}
				else
				{	/* if Nset==1 */
					BestSet[fold] = 0;
				}

				/* Prediction with the BestSet */
				FittingTime = 0.0;
				LBmonitor = (double*) calloc (Maxiteration, sizeof(double));
				Rmonitor =  (double*) calloc (Maxiteration, sizeof(double));
				if(Printinfo) printf ("====Progress:%d/%d trait:%s fold:%d set:%d ====\n", Currentrun, Nruns, Y[trait].name, fold+1, BestSet[fold]+1);
				GenomeWideRegression (Algorithm, Priortype, Methodcode, CondResidual, P, F, N, Thresholdvalue, Maxiteration, Ycv, X, Q, H+BestSet[fold], &expDelta2, Tau0, LBmonitor, Rmonitor, &FittingTime, Printinfo, RandomIni);
				Currentrun ++;

				/* When AnalysisType == 1, fitted values with the best set are returned. Otherwise, predictions are returned */
				if(AnalysisType==1)
				{
					for ( record=0; record<N; record++ )
					{
						PredictedValues[record] = PredictYhat (Q, F, X, P, record, Methodcode);
						BreedingValues [record] = PredictBV (X, P, record, Methodcode);
					}
				}
				else
				{	/* Predicted values are pooled across folds */
					for ( record=0; record<A[fold].ntest; record++ )
					{
						PredictedValues[record+CumNtest] = PredictYhat (Q, F, X, P, A[fold].obstest[record], Methodcode);
						BreedingValues [record+CumNtest] = PredictBV (X, P, A[fold].obstest[record], Methodcode);
					}
					CumNtest += A[fold].ntest;
				}
		
				if (Fold>1) { free(LBmonitor); free(Rmonitor);}	
				free(Ycv[0].observations); free(Ycv);

			}	/* fold */
			CVend = clock();

			/* Output */
			if(Printinfo) printf ("Total calculation time (sec.): %f\n", (double)(CVend-CVstart)/CLOCKS_PER_SEC);

			if(Fold>1)	/* When AnalysisType == -9, -1, or >1. Only predicted values are returned */
			{
				WriteResultsPrediction (Algorithm, Methodcode, Phenotypefile, Genotypefile, Covariatefile, Partitionfile,
					AnalysisType, Ntraits, P, F, TotalNtest, Threshold, Maxiteration, usedseed, 
					Y, H, X, Q, Rmonitor, LBmonitor, expDelta2, PredictedValues, BreedingValues, trait, Fold, BestSet, BestValue, Nset, (double)(CVend-CVstart)/CLOCKS_PER_SEC, CVFoldTuning, A, MPP, NDCGrank, Nperm, RandomIni);
			}
			else
			{		/* When AnalysisType == 1. Parameter values estimated with the best set are returned */

				WriteResults (Algorithm, Methodcode, Phenotypefile, Genotypefile, Covariatefile, Partitionfile, AnalysisType, Ntraits, P, F, N, Threshold, Maxiteration, usedseed,
					Y, H, X, Q, Rmonitor, LBmonitor, expDelta2, trait, BestSet[0], CVFoldTuning, Nset, ValueList, PredictedValues, BreedingValues, FittingTime, MPP, NDCGrank, Nperm, PermutationThreshold, RandomIni);

				free(LBmonitor); free(Rmonitor); free(ValueList);

			}	/* if(Fold>1) */

			/* free */
			free(PredictedValues); free(BreedingValues); free(BestValue); free(BestSet);
			if(AnalysisType>1) { for (fold=0; fold<Fold; fold++) free(A[fold].obstest); free (A);}

		}	/* AnalysisType if-else */

	}		/* ntrait */

	if(AnalysisType==-1||AnalysisType==1||AnalysisType==-9) { for (fold=0; fold<Fold; fold++) free(A[fold].obstest); free (A);}
	for (locus=0; locus<P; locus++) free(X[locus].covariates);
	free (X);
	for (trait=0; trait<Ntraits; trait++) free(Y[trait].observations);
	free (Y);
	free (H);
	for (locus=0; locus<F; locus++) free(Q[locus].covariates);
	free (Q);
	if(Nperm>0)
	{
		for (locus=0; locus<P; locus++) free(Xperm[locus].covariates);
		free (Xperm);
		for (trait=0; trait<Ntraits; trait++) free(Yperm[trait].observations);
		free (Yperm);
		for (locus=0; locus<F; locus++) free(Qperm[locus].covariates);
		free (Qperm);
	}

	return(0);
}
