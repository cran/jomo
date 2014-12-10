#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "pdflib.h"
#include "rnglib.h"
#include "wishart.h"
#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>

SEXP jomo1con(SEXP Y, SEXP Yimp, SEXP Yimp2, SEXP X, SEXP beta, SEXP betapost, SEXP omega, SEXP omegapost, SEXP nstep, SEXP Sp, SEXP flagrng){
int i,j,k, IY,JY, IX, JX, Io, Jo, Ib, Jb, ns, nmiss=0,t, countm=0, counto=0, countmm=0, countmo=0, countoo=0, jj, tt,fl;
SEXP RdimY, RdimX, Rdimo, Rdimb;
double *betaX, *Yobs, *Ymiss, *mumiss, *omegadrawmiss, *betamiss, *betaobs, *omegaoo, *omegaom, *omegamo, *omegamm, *invomega, *invomega2, *help, *help2, *help3, *imp, *zi;
double *sumzy, *incrzz, *incrzy, *mu, *mu2, *newbeta, *newomega, *ei, *sumzi, *yi, *invomega3, *help4, *help5, *help6,*help7, *missing;

/* Protecting R objects from garbage collection and saving matrices dimensions*/ 

RdimY=PROTECT(getAttrib(Y,R_DimSymbol));
IY=INTEGER(RdimY)[0];
JY=INTEGER(RdimY)[1];
RdimX=PROTECT(getAttrib(X,R_DimSymbol));
IX=INTEGER(RdimX)[0];
JX=INTEGER(RdimX)[1];if(IX!=IY) error("Covariates and Responses matrices have different length");
Rdimb=PROTECT(getAttrib(beta,R_DimSymbol));
Ib=INTEGER(Rdimb)[0];
Jb=INTEGER(Rdimb)[1];
Rdimo=PROTECT(getAttrib(omega,R_DimSymbol));
Io=INTEGER(Rdimo)[0];
Jo=INTEGER(Rdimo)[1];
Y=PROTECT(coerceVector(Y,REALSXP));
Yimp=PROTECT(coerceVector(Yimp,REALSXP));
Yimp2=PROTECT(coerceVector(Yimp2,REALSXP));
X=PROTECT(coerceVector(X,REALSXP));
beta=PROTECT(coerceVector(beta,REALSXP));
betapost=PROTECT(coerceVector(betapost,REALSXP));
omega=PROTECT(coerceVector(omega,REALSXP));
omegapost=PROTECT(coerceVector(omegapost,REALSXP));
Sp=PROTECT(coerceVector(Sp,REALSXP));
nstep=PROTECT(coerceVector(nstep,INTSXP));
ns=INTEGER(nstep)[0];
flagrng=PROTECT(coerceVector(flagrng,INTSXP));
fl=INTEGER(flagrng)[0];

/*Allocating memory for C objects in R*/

help = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
invomega= (double * ) R_alloc ( JY * JY , sizeof ( double ) );
sumzi = ( double * ) R_alloc ( JY * JX * JY * JX , sizeof ( double ) );
sumzy = ( double * ) R_alloc ( JY * JX , sizeof ( double ) );
zi = ( double * ) R_alloc ( JY * JX * JY , sizeof ( double ) );
yi = ( double * ) R_alloc ( JY , sizeof ( double ) );
help2 = ( double * ) R_alloc ( JY *JX * JY , sizeof ( double ) );
incrzz = ( double * ) R_alloc ( JY *JX * JY *JX , sizeof ( double ) );
incrzy = ( double * ) R_alloc ( JY *JX , sizeof ( double ) );
help3 = ( double * ) R_alloc ( JY * JX * JY * JX ,sizeof ( double ) );
invomega2= (double * ) R_alloc ( JY * JX * JY * JX , sizeof ( double ) );
mu = ( double * ) R_alloc ( JY * JX, sizeof ( double ) );
newbeta = ( double * ) R_alloc ( JY * JX ,sizeof ( double ) );
mu2 = ( double * ) R_alloc ( JY * JY,sizeof ( double ) );
ei = ( double * ) R_alloc ( JY ,sizeof ( double ) );
newomega = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
betaX=( double * ) R_alloc ( JY , sizeof ( double ) );
imp=( double * ) R_alloc ( IY * JY,sizeof ( double ) );
Yobs=( double * ) R_alloc ( JY , sizeof ( double ) );
Ymiss=( double * ) R_alloc ( JY , sizeof ( double ) );
mumiss = ( double * ) R_alloc ( JY , sizeof ( double ) );
omegadrawmiss = ( double * ) R_alloc ( JY * JY ,sizeof ( double ) );
betamiss = ( double * ) R_alloc ( JY ,sizeof ( double ) );
betaobs = ( double * ) R_alloc ( JY, sizeof ( double ) );
omegaoo= ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
omegaom= ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
omegamo= ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
omegamm= ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
invomega3= ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
help4 = ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
help5 = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
help6 = ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
help7 = ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
missing = ( double * ) R_alloc ( IY , sizeof ( double ) );

/* Some initializations */

GetRNGstate();
for (j=0; j<IY; j++) {
	missing[j]=0;
	for (k=0;k<JY;k++) {
		if (ISNAN(REAL(Y)[j+k*IY])) {
			missing[j]++;
		}
	}
}
r8mat_copy_new(IY, JY, REAL(Yimp), imp);
for (i=0;i<Ib*Jb;i++) REAL(betapost)[i]=0;

/* Running ns iterations of Gibbs sampler*/

for (i=0;i<ns;i++) {

	/*Updating beta*/

	r8mat_pofac(JY,REAL(omega), help,1);
	r8mat_poinv(JY,help, invomega);
	for (j=1;j<Io;j++) for (t=0;t<j;t++) invomega[j+Io*t]=invomega[t+Io*j];

	for (j=0;j<JY*JY*JX*JX;j++) sumzi[j]=0;
	for (j=0;j<JY*JX;j++) sumzy[j]=0;
	for (j=0;j<IY;j++) {
		for (t=0;t<JY*JX*JY;t++) zi[t]=0; 
		for (t=0;t<JY;t++) {
			yi[t]=imp[j+t*IY];
			for (k=0;k<JX;k++) {
				zi[t+(k+t*JX)*JY]=REAL(X)[j+IY*k];
			}
		}
		r8mat_mtm_new(JY*JX,JY,JY,zi,invomega, help2);
		r8mat_mm_new(JY*JX,JY,JX*JY,help2,zi, incrzz);
		r8mat_mmt_new(JY*JX,JY,1,help2,yi, incrzy);
		r8mat_add(JY*JX,JY*JX,incrzz,sumzi);
		r8mat_add(JY*JX,1,incrzy,sumzy);
	}
	r8mat_pofac(JY * JX,sumzi, help3,2);
	r8mat_poinv(JY * JX, help3, invomega2);
	for (jj=1;jj<JX*JY;jj++) for (tt=0;tt<jj;tt++) invomega2[jj+JX*JY*tt]=invomega2[tt+JX*JY*jj];
	r8mat_mm_new(JY*JX,JY*JX,1,invomega2,sumzy,mu);
	r8mat_pofac(JY * JX,invomega2,help3,3);
	r8vec_multinormal_sample(JY*JX, mu,help3, REAL(beta),newbeta,fl);
	r8mat_add(Ib,Jb,REAL(beta),REAL(betapost));

	/*Updating omega*/
	
	for (j=0;j<JY*JY;j++) mu2[j]=0;
	for (j=0;j<IY;j++) {
		for (k=0;k<JY;k++) {
			ei[k]=imp[j+IY*k];
			for (t=0;t<JX;t++) ei[k]=ei[k]-REAL(beta)[t + Ib*k]*REAL(X)[j+t*IX];
		}
		r8mat_mmt_new(JY,1,JY,ei,ei,help);
		r8mat_add(JY,JY,help,mu2);
	}
	r8mat_add(JY,JY,REAL(Sp),mu2);
	r8mat_pofac(JY,mu2,help,4);
	r8mat_poinv(JY, help,invomega);
	for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
	wishart_sample(JY,IY-1,invomega,newomega,help,omegaoo,omegamm,omegaom,fl);
	r8mat_pofac(JY,newomega,help,5);
	r8mat_poinv(JY, help,invomega);
	for(k=0;k<Io;k++)  {
		for(j=0;j<k+1;j++)  {
			REAL(omega)[k+Io*j]=invomega[j+Io*k];
			REAL(omega)[j+Io*k]=REAL(omega)[k+Io*j];
		}
	}
	r8mat_add(Io,Jo,REAL(omega),REAL(omegapost));

	/*imputing missing values*/

	for (j=0; j<IY; j++) {
		for (k=0;k<JY;k++) betaX[k]=0;
		for (t=0;t<JX;t++) {
			for (k=0;k<JY;k++) {
				betaX[k]=betaX[k]+REAL(beta)[t+k*Ib]*REAL(X)[j+t*IX];
			}
		}
		nmiss=missing[j];
		if (nmiss>0) {
			for (k=0;k<JY;k++) {
				if (ISNAN(REAL(Y)[j+k*IY])) {
					betamiss[countm]=betaX[k];
					countm++;
				}
				else {
					Yobs[counto]=REAL(Y)[j+k*IY];
					betaobs[counto]=betaX[k];
					counto++;
				}
				for (t=0;t<JY;t++) {
					if (ISNAN(REAL(Y)[j+k*IY])&ISNAN(REAL(Y)[j+t*IY])) {
						omegamm[countmm]=REAL(omega)[k+t*Io];
						countmm++;
					}
					else if (ISNAN(REAL(Y)[j+t*IY])) {
						omegamo[countmo]=REAL(omega)[k+t*Io];
						countmo++;	
					}
					else if (!ISNAN(REAL(Y)[j+k*IY])&!ISNAN(REAL(Y)[j+t*IY])){
						omegaoo[countoo]=REAL(omega)[k+t*Io];
						countoo++;	
					}
				}
			}
			r8mat_pofac((JY-nmiss),omegaoo,help4,11);
			r8mat_poinv((JY-nmiss),help4,invomega3);
			for (jj=1;jj<JY-nmiss;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(JY-nmiss)*tt]=invomega3[tt+(JY-nmiss)*jj];
			r8mat_mmt_new((JY-nmiss),(JY-nmiss),nmiss,invomega3,omegamo,help5);
			r8mat_divide((JY-nmiss),1,-1,betaobs);
			r8mat_add((JY-nmiss),1,betaobs,Yobs);
			r8mat_mtm_new(1,(JY-nmiss),nmiss,Yobs,help5,mumiss);
			r8mat_add(1,nmiss,betamiss,mumiss);
			r8mat_mm_new(nmiss,(JY-nmiss),nmiss,omegamo,help5,omegadrawmiss);
			r8mat_divide(nmiss,nmiss,-1,omegadrawmiss);
			r8mat_add(nmiss,nmiss,omegamm,omegadrawmiss);
			r8mat_pofac(nmiss,omegadrawmiss,help7,12);
			r8vec_multinormal_sample(nmiss,mumiss,help7,Ymiss,help6,fl);
			countm=0;
			for (k=0;k<JY;k++) {
				if (ISNAN(REAL(Y)[j+k*IY])) {
					imp[j+k*IY]=Ymiss[countm];
					countm++;	
				}	
			}
			nmiss=0;
			countm=0;
			counto=0;
			countmm=0;
			countmo=0;
			countoo=0;
			
		}
	}
if ((i+1)%10==0) Rprintf("Iteration %d completed\n",i+1);
}
for(i=0;i<IY;i++)  {
	for(j=0;j<JY;j++)  {
		REAL(Yimp2)[i+IY*j]=imp[i+IY*j];
	}
}
r8mat_divide(Ib,Jb,ns,REAL(betapost));
r8mat_divide(Io,Jo,ns,REAL(omegapost));
PutRNGstate();
UNPROTECT(15);
return R_NilValue;
}
