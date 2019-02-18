#include <R.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double sum(double *x, int dim) {
	double s;
	int i;
	s = 0.0;
	for (i=0; i<dim; ++i) s+=x[i];
	return(s);
}

/*****************************************************************************
 * Function returns the years lived by those who died in the age interval (ax)
 * for ages 0 and 1-4 (abridged life table)
 * based on mx(0,1) and sex
 * Formulas re-estimates from the Coale/Demeny separation factors
 * on mx(0,1) insted of qx(0,1)
 * Source from Preston et al. 2001, p.48
 *****************************************************************************/
double * get_a05(double mx0, int sex) {
	static double ax[2];
	if(sex > 1) {/* female*/
		if (mx0 < 0.107) {
			ax[0] = 0.053 + 2.8 * mx0;      /*1a0*/
			ax[1] = 1.522 - 1.518 * mx0;    /*4a1*/
		} else {
			ax[0] = 0.35;
			ax[1] = 1.361;
		}
	} else { /* male */
		if (mx0 < 0.107) {
			ax[0] = 0.045 + 2.684 * mx0;
			ax[1] = 1.651 - 2.816 * mx0;
		} else {
			ax[0] = 0.33;
			ax[1] = 1.352;
		}
	}
	return(ax);
}
/*****************************************************************************
 * Function calculates an abridged life table from age-specific mortality 
 * rates
 * Input: 
 * * mx   age specific mortality rates 
 * * sex  sex
 * * nage number of age groups
 * Output: 
 * qx Probabilities of dying
 * lx Survivors to exact age
 * Lx Person years lived
 * ax proportion of years lived by those who died
 *****************************************************************************/
void doLifeTable(int sex, int nage, double *mx, 
				double *Lx, double *lx, double *qx, double *ax) {
	
	int i;
	double k;     /* correcting factor in Greville approximation */
	double *tmpa; /* pointer to estimated ax[0] and ax[1] values */
	int nage1;
	nage1 = nage -1;

	tmpa = get_a05(mx[0], sex);
	ax[0] = tmpa[0];
	ax[1] = tmpa[1];
	qx[0] = mx[0] / (1 + (1 - ax[0]) * mx[0]);                    /* 1q0 */	
	qx[1] = 4 * mx[1] / (1 + (4 - ax[1]) * mx[1]);				  /* 4q1 */	
	lx[0] = 1;                                                    /* l0 */
	lx[1] = lx[0] * (1 - qx[0]);                                  /* l1 */
	lx[2] = lx[1] * (1 - qx[1]);                                  /* l5 = l1 * (1-4q1) */
	Lx[0] =  lx[1] + ax[0] * (lx[0] - lx[1]);                     /* 1L0 */
	Lx[1] =  4 * lx[2] + ax[1] * (lx[1] - lx[2]);                 /* 4L1 */
		
    /*Rprintf("\nnage=%i, L0=%f, ax0-1=%f %f, l1-2=%f %f, mx0-1=%f %f", nage, Lx[0], ax[0], ax[1], lx[1], lx[2], mx[0], mx[1]);*/
    /* Age 5-9, .... 125-129 
	 Greville formula used in Mortpak and UN MLT (1982)*/
	/* TB: corrected for age group 3 (5-9), tentativley set to 2.5 or n/2 */
	/*  ax[2] = 2.5; */

	k     = 0.1 * log(fmax(mx[4] / fmax(mx[2], DBL_MIN), DBL_MIN));
	ax[2] = 2.5 - (25 / 12.0) * (mx[2] - k);

	for(i = 3; i < nage1; ++i) {
		k     = 0.1 * log(fmax(mx[i+1] / fmax(mx[i-1], DBL_MIN), DBL_MIN));
		ax[i] = 2.5 - (25 / 12.0) * (mx[i] - k);
	}
	/* penultimate ax calculated with k from previous age group */
	ax[nage1] = 2.5 - (25 / 12.0) * (mx[i] - k);

	/* correcting out-of (reasonable) bounds ax for older ages             */ 
	/* 0.97=1-5*exp(-5)/(1-exp(-5)), for constant mu=1, Kannisto assumption*/
	for(i = 10; i < nage; i++) {
		if(ax[i] < 0.97) {
			ax[i] = 0.97;
		}
	}
	
	/* caculate life table variables from mx and ax */
	for(i = 2; i < nage; ++i) {		
		/*Rprintf("ax%i=%f, mx%i=%f", i, ax[i], i-1, mx[i-1]);*/
		qx[i] = 5 * mx[i] / (1 + (5 - ax[i]) * mx[i]);
		lx[i+1] = fmax(lx[i] * (1-qx[i]), DBL_MIN);
		Lx[i] = 5 * lx[i+1] + ax[i] * (lx[i] - lx[i+1]);
	}
	
	/* Open ended age interval */
	Lx[nage] = lx[nage] / fmax(mx[nage], DBL_MIN); /* Assuming Mx levels off at age 130 */
	qx[nage] = 1.0;
	
	ax[nage] = Lx[nage];
}



/* Function returns collapsed Lx and lx columns of life table */
/* function calls doLifeTable first, then collapsesLx and lx  */
void LifeTableC(int sex, int nage, double *mx,
                double *Lxx, double *lxx) {
  /* life table variables returned from doLifeTable */
  /* need to be declared with nage+1 elements       */
  /* Lx[nage+1]
   * lx[nage+1]
   * qx[nage+1]
   * ax[nage+1];
   * input to doLifeTableC mx is declared outside
   * output from LifeTableC Lxx, lxx is declared 
   * with nage-1 elements outside 
   */
  
  double Lx[nage+1], lx[nage+1], qx[nage+1], ax[nage+1];
  int i;
  
  /* do life table called with nage as last index */
  doLifeTable(sex, nage, mx, Lx, lx, qx, ax);
  /* collapse 1L0 and 4L1 into 5L0 */
  Lxx[0] = Lx[0] + Lx[1];
  lxx[0] = lx[0];
  for(i = 1; i < nage; ++i) {
    Lxx[i] = Lx[i+1];
    lxx[i] = lx[i+1];
  }
}

double get_constrained_mortality(double a, double b, double k, double constraint) {
	double mx;
	
	mx = exp(a + b*k);
	if(constraint > 0 && mx < constraint) mx = constraint;
	return mx;
}

void LCEoKtC(int sex, int nage, double *ax, double *bx, 
			 double eop, double kl, double ku, double *constraints, 
			 double *LLm, double *lm, double *Mx) {
	double LTl[27], LTu[nage-1], mxm[nage], LTeo, k2;
	int i, dim;
	dim = nage-1;

	/* check if the eop lies outside of the bounds */
	for (i=0; i < nage; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], kl, constraints[i]);
	}
	LifeTableC(sex, dim, mxm, LTl, lm);
	
	if(eop < sum(LTl, dim)) {
		for (i=0; i < dim; ++i) LLm[i]=LTl[i];
		for (i=0; i < nage; ++i) Mx[i] = mxm[i];
		return;
	}
	for (i=0; i < nage; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], ku, constraints[i]);
	}
	LifeTableC(sex, dim, mxm, LTu, lm);

	if(eop > sum(LTu, dim)) {
		for (i=0; i < dim; ++i) LLm[i]=LTu[i]; 
		for (i=0; i < nage; ++i) Mx[i] = mxm[i];
		return;
	}
	/* Bi-section method */
	k2 = 0.5 * (kl + ku);
	for (i=0; i < nage; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], k2, constraints[i]);
	}
	LifeTableC(sex, dim, mxm, LLm, lm);
	LTeo = sum(LLm, dim);
	while(fabs(LTeo - eop) > 0.01) {
		if(LTeo < eop) kl = k2;
		else ku = k2;
		k2 = 0.5 * (kl + ku);
		for (i=0; i < nage; ++i) {
			mxm[i] = get_constrained_mortality(ax[i], bx[i], k2, constraints[i]);
		}
		LifeTableC(sex, dim, mxm, LLm, lm);
		LTeo = sum(LLm, dim);
	}
	for (i=0; i < nage; ++i) Mx[i] = mxm[i];
}

void get_sx(double *LLm, double *sx, int n, int Ldim) {
    /* compute survival ratios from Lx where the first age group is 0-5 (also for Lx)*/
    int i, oei;
    double sumLL;
    oei=n-1;
    /* Survival Ratios, radix of life table assumed to be 1.0  */
    sx[0] = LLm[0] / 5.0;
    for(i=1; i < oei; ++i) {
        if(LLm[i-1] == 0) sx[i] = exp(-5);
        else sx[i] = LLm[i]/LLm[i-1];
    }
    /* Last age group */
    sumLL = 0;
    for(i=oei; i < Ldim; ++i) {
        sumLL += LLm[i];
    }
    if((sumLL + LLm[oei-1]) == 0 ||  sumLL == 0) sx[oei] = exp(-5);
    else sx[oei] = sumLL/(sumLL+LLm[oei-1]);
    if(sx[oei] > sx[oei-1]) sx[oei] = sx[oei-1];
}


/*****************************************************************************
 * Lee Carter model
 * Produces a projection of age-specific mortality rates
 * 
 *****************************************************************************/
void LC(int *Npred, int *Sex, int *Nage, double *ax, double *bx, 
		double *Eop, double *Kl, double *Ku, int *constrain, double *FMx, double *FEop, 
		double *LLm, double *Sr, double *lx, double *Mx) {
	double eop, sx[*Nage-1], Lm[*Nage-1], mxm[*Nage], fmx[*Nage], lm[*Nage], locbx[*Nage], locax[*Nage];
	int i, sex, npred, pred, nage, nagem1, cage;
	
	npred = *Npred;
	sex=*Sex;
	nage=*Nage;
	nagem1 = nage-1;
	cage = -1;
	
	if(*constrain == 1) cage = 22; /* constrain old ages only */ 
	if(*constrain == 2) cage = 0;  /* constrain all ages */
	
	for (i=0; i < nage; ++i) fmx[i] = -1;
	for (pred=0; pred < npred; ++pred) {
		eop = Eop[pred];
		if(*constrain > 0 && nage > cage) {		
			if(FEop[pred] > eop) {
				for (i=cage; i < nage; ++i) {fmx[i] = FMx[i + pred*nage];}
			} else {
				for (i=cage; i < nage; ++i) {fmx[i] = -1;}
			}
		}
		for (i=0; i < nage; ++i) {
			locbx[i] = bx[i + pred*nage];
			locax[i] = ax[i + pred*nage];
		}
		/*Rprintf("\n%i: eop=%lf", pred, eop);*/
		LCEoKtC(sex, nage, locax, locbx, eop, Kl[pred], Ku[pred], fmx, Lm, lm, mxm);		
		get_sx(Lm, sx, nagem1, nagem1);
		
		for (i=0; i < nagem1; ++i) {
			Sr[i + pred*(nagem1)] = sx[i];
			/*Rprintf("\nLLm=%lf, Sr=%lf", LLm[i], Sr[i + pred*27]);*/
			Mx[i + pred*nage] = mxm[i];
			lx[i + pred*nage] = lm[i];
		}
		Mx[nagem1 + pred*nage] = mxm[nagem1];
		lx[nagem1 + pred*nage] = lm[nagem1];
		for (i=0; i < nage-2; ++i) {
			LLm[i + pred*(nagem1)] = Lm[i];
		}
	}
}

/*****************************************************************************
 * PMD model
 * Produces a projection of age-specific mortality rates
 * 
 *****************************************************************************/
void PMD(int *Npred, int *Sex, int *Nage, double *mx0, double *rho, 
        double *Eop, double *Kl, double *Ku, double *Constr, int *Nconstr,
        double *LLm, double *Sr, double *lx, double *Mx) {
    double eop, sx[*Nage-1], Lm[*Nage-1], mxm[*Nage], lm[*Nage], locrho[*Nage], locmx[*Nage], constr[*Nage];
    int i, sex, npred, pred, nage, nagem1, nconstr;
    
    npred = *Npred;
    sex=*Sex;
    nage=*Nage;
    nagem1 = nage-1;
    nconstr = *Nconstr;
    
    for (i=0; i < nage; ++i) {
        locmx[i] = log(mx0[i]);
    }
    for (pred=0; pred < npred; ++pred) {
        eop = Eop[pred];
        for (i=0; i < nage; ++i) {
            locrho[i] = -rho[i + pred*nage];
            constr[i] = -1;
        }
        if(nconstr > 0) {
            for(i=0; i < nconstr; ++i) {
                constr[i] = Constr[i + pred*nconstr];
            }
        }

        /*Rprintf("\nconstr: sex %i period %i: ", sex, pred);
        for (i=0; i < nage; ++i) 
            Rprintf("%lf, ", i+1, constr[i]);*/
        
        /*Rprintf("\n%i: eop=%lf", pred, eop);*/
        LCEoKtC(sex, nage, locmx, locrho, eop, Kl[0], Ku[0], constr, Lm, lm, mxm);		
        get_sx(Lm, sx, nagem1, nagem1);
        
        for (i=0; i < nagem1; ++i) {
            Sr[i + pred*(nagem1)] = sx[i];
            /*Rprintf("\nLLm=%lf, Sr=%lf", LLm[i], Sr[i + pred*27]);*/
            Mx[i + pred*nage] = mxm[i];
            lx[i + pred*nage] = lm[i];
            locmx[i] = log(mxm[i]);
        }
        Mx[nagem1 + pred*nage] = mxm[nagem1];
        lx[nagem1 + pred*nage] = lm[nagem1];
        locmx[nagem1] = log(mxm[nagem1]);
        for (i=0; i < nage-2; ++i) {
            LLm[i + pred*(nagem1)] = Lm[i];
        }
    }
}

