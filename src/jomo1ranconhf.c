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

SEXP jomo1ranconhf(SEXP Y, SEXP Yimp, SEXP Yimp2, SEXP X, SEXP Z, SEXP clus, SEXP beta, SEXP u, SEXP betapost, SEXP upost, SEXP omega, SEXP omegapost, SEXP covu, SEXP covupost, SEXP nstep, SEXP Sp, SEXP Sup, SEXP flagrng){
int i,j,k, IY,JY, IX, JX, Io, Jo, Ib, Jb, IZ, JZ, Iu, Ju, ns,nj, nmiss=0,t, countm=0, counto=0, countmm=0, countmo=0, countoo=0,countt=0, jj, tt, c,flag=0, fl;
SEXP RdimY, RdimX, Rdimo, Rdimb, RdimZ, Rdimu;
double *betaX, *Yobs, *Ymiss, *mumiss, *omegadrawmiss, *betamiss, *betaobs, *omegaoo, *omegaom, *omegamo, *omegamm, *invomega, *help, *imp, *zi, *sumxi, *yi, *xi, *help5;
double *sumzy, *incrxx, *incrzy, *incrzz, *mu, *newbeta, *newomega, *ziu, *uj, *sumzi, *newu, *help2, *help3, *help4,*invomega2, *sumxy, *incrxy;
double *invomega3, *mu2, *mu3, *mu4, *newomega2, *help6, *help7, *help8, *help9, *invomega4, *missing, *clusnum, *cumclus, *allinvomega;

/* Protecting R objects from garbage collection and saving matrices dimensions*/ 

RdimY=PROTECT(getAttrib(Y,R_DimSymbol));
IY=INTEGER(RdimY)[0];
JY=INTEGER(RdimY)[1];
RdimX=PROTECT(getAttrib(X,R_DimSymbol));
IX=INTEGER(RdimX)[0];
JX=INTEGER(RdimX)[1];
RdimZ=PROTECT(getAttrib(Z,R_DimSymbol));
IZ=INTEGER(RdimZ)[0];
JZ=INTEGER(RdimZ)[1];if(IX!=IY) error("Covariates and Responses matrices have different length");
Rdimb=PROTECT(getAttrib(beta,R_DimSymbol));
Ib=INTEGER(Rdimb)[0];
Jb=INTEGER(Rdimb)[1];
Rdimu=PROTECT(getAttrib(u,R_DimSymbol));
Iu=INTEGER(Rdimu)[0];
Ju=INTEGER(Rdimu)[1];
Rdimo=PROTECT(getAttrib(omega,R_DimSymbol));
Io=INTEGER(Rdimo)[0];
Jo=INTEGER(Rdimo)[1];
Y=PROTECT(coerceVector(Y,REALSXP));
Yimp=PROTECT(coerceVector(Yimp,REALSXP));
Yimp2=PROTECT(coerceVector(Yimp2,REALSXP));
X=PROTECT(coerceVector(X,REALSXP));
Z=PROTECT(coerceVector(Z,REALSXP));
clus=PROTECT(coerceVector(clus,INTSXP));
beta=PROTECT(coerceVector(beta,REALSXP));
u=PROTECT(coerceVector(u,REALSXP));
upost=PROTECT(coerceVector(upost,REALSXP));
betapost=PROTECT(coerceVector(betapost,REALSXP));
omega=PROTECT(coerceVector(omega,REALSXP));
covu=PROTECT(coerceVector(covu,REALSXP));
omegapost=PROTECT(coerceVector(omegapost,REALSXP));
covupost=PROTECT(coerceVector(covupost,REALSXP));
Sp=PROTECT(coerceVector(Sp,REALSXP));
Sup=PROTECT(coerceVector(Sup,REALSXP));
nstep=PROTECT(coerceVector(nstep,INTSXP));
ns=INTEGER(nstep)[0];
flagrng=PROTECT(coerceVector(flagrng,INTSXP));
fl=INTEGER(flagrng)[0];
nj=Iu;

/*Allocating memory for C objects in R*/

imp=( double * )  R_alloc ( IY * JY , sizeof(double) );
help = ( double * )  R_alloc ( JY * JY , sizeof(double));
invomega= (double * )  R_alloc ( JY * JY , sizeof(double) );
sumxi = ( double * )  R_alloc ( JY * JX * JY * JX,  sizeof(double) );
sumxy = ( double * )  R_alloc ( JY * JX, sizeof(double));
zi = ( double * )  R_alloc ( JY * JZ * JY , sizeof(double));
xi = ( double * )  R_alloc ( JY * JX * JY , sizeof(double) );
yi = ( double * )  R_alloc ( JY ,  sizeof(double) );
uj = ( double * )  R_alloc ( JY * JZ , sizeof(double)  );
help2 = ( double * )  R_alloc ( JY *JX * JY ,  sizeof(double));
incrxx = ( double * )  R_alloc ( JY *JX * JY *JX ,sizeof(double)  );
incrxy = ( double * )  R_alloc ( JY *JX , sizeof(double) );
ziu = ( double * )  R_alloc ( JY ,sizeof(double) );
help3 = ( double * )  R_alloc ( JY * JX * JY * JX , sizeof(double) );
invomega2= (double * )  R_alloc ( JY * JX * JY * JX , sizeof(double));
mu = ( double * )  R_alloc ( JY * JX , sizeof(double)  );
newbeta = ( double * )  R_alloc ( JY * JX , sizeof(double) );
sumzi = ( double * )  R_alloc ( JY * JZ * JY * JZ ,sizeof(double)  );
sumzy = ( double * )  R_alloc ( JY * JZ,  sizeof(double) );
help4 = ( double * )  R_alloc ( JY *JZ * JY ,  sizeof(double) );
incrzz = ( double * )  R_alloc ( JY *JZ * JY *JZ , sizeof(double) );
incrzy = ( double * )  R_alloc ( JY *JZ , sizeof(double) );
help5 = ( double * )  R_alloc ( JY * JZ * JY * JZ,sizeof(double) );
invomega3= (double * )  R_alloc ( JY * JZ * JY * JZ ,sizeof(double) );
mu2 = ( double * )  R_alloc ( JY * JZ , sizeof(double) );
newu = ( double * )  R_alloc ( JY * JZ , sizeof(double) );
mu3 = ( double * )  R_alloc ( JY * JZ * JY * JZ , sizeof(double));
newomega = ( double * )  R_alloc ( JY * JZ * JY * JZ ,  sizeof(double));
mu4 = ( double * )  R_alloc ( JY * JY , sizeof(double));
newomega2 = ( double * )  R_alloc ( JY * JY , sizeof(double));
betaX=( double * )  R_alloc ( JY , sizeof(double));
help6 = ( double * )  R_alloc ( JY ,  sizeof(double) );
Yobs=( double * )  R_alloc ( JY  , sizeof(double));
Ymiss=( double * )  R_alloc ( JY ,  sizeof(double));
mumiss = ( double * )  R_alloc ( JY , sizeof(double) );
omegadrawmiss = ( double * )  R_alloc ( JY * JY , sizeof(double) );
betamiss = ( double * )  R_alloc ( JY , sizeof(double));
betaobs = ( double * )  R_alloc ( JY , sizeof(double));
omegaoo= ( double * )  R_alloc ( JY *JY,sizeof(double));
omegaom= ( double * )  R_alloc (JY *JY , sizeof(double));
omegamo= ( double * )  R_alloc ( JY *JY ,  sizeof(double));
omegamm= ( double * )  R_alloc ( JY *JY , sizeof(double));
invomega4= ( double * )  R_alloc ( JY *JY ,sizeof(double) );
help7 = ( double * )  R_alloc (JY *JY, sizeof(double) );
help8 = ( double * )  R_alloc ( JY *JY, sizeof(double));
help9 = ( double * )  R_alloc ( JY * JY , sizeof(double));
missing = ( double * )  R_alloc ( IY , sizeof(double));
clusnum = ( double * )  R_alloc ( nj , sizeof(double));
cumclus = ( double * )  R_alloc ( nj+1 , sizeof(double));
allinvomega = ( double * )  R_alloc ( nj* JY * JY , sizeof(double));

/* Some initializations */

for (j=0; j<IY; j++) {
	missing[j]=0;
	for (k=0;k<JY;k++) {
		if (ISNAN(REAL(Y)[j+k*IY])) {
			missing[j]++;
		}
	}
}
cumclus[0]=0;
for (i=0;i<nj;i++) {
	clusnum[i]=0;
	for (j=0; j<IY; j++) {		
		if (INTEGER(clus)[j]==i) {
			clusnum[i]++;
			
		}
	}
	cumclus[i+1]=cumclus[i]+clusnum[i];
}
cumclus[1]=clusnum[1];
r8mat_copy_new(IY, JY, REAL(Yimp), imp);
for (i=0;i<Ib*Jb;i++) REAL(betapost)[i]=0;
GetRNGstate();

/* Running ns iterations of Gibbs sampler*/

for (i=0;i<ns;i++) {
	for (c=0;c<nj;c++) {
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) newomega[j+t*JY]=REAL(omega)[(c*JY+j)+t*(JY*nj)];
			}
		r8mat_pofac(JY,newomega, help,1);
		r8mat_poinv(JY,help, invomega);
		for (jj=1;jj<JY;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY*tt]=invomega[tt+JY*jj];
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) allinvomega[(c*JY+j)+t*(JY*nj)]=invomega[j+t*JY];
			}
		}
	for (j=0;j<JY*JY*JX*JX;j++) sumxi[j]=0;	
	for (j=0;j<JY*JX;j++) sumxy[j]=0;
	
	for (c=0;c<nj;c++) {
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) invomega[j+t*JY]=allinvomega[(c*JY+j)+t*(JY*nj)];
			}
		for (j=cumclus[c];j<cumclus[c+1];j++) {
			for (t=0;t<JY*JZ*JY;t++) zi[t]=0;
			for (t=0;t<JY*JX*JY;t++) xi[t]=0; 
			for (t=0;t<JY;t++) {
				yi[t]=imp[j+t*IY];
				for (k=0;k<JX;k++) {
					xi[t+(k+t*JX)*JY]=REAL(X)[j+IY*k];
				}
				for (k=0;k<JZ;k++) {
					zi[t+(k+t*JZ)*JY]=REAL(Z)[j+IY*k];
					uj[k+t*JZ]=REAL(u)[(INTEGER(clus)[j])+nj*(k+t*JZ)];
				}
			}
			r8mat_mtm_new(JY*JX,JY,JY,xi,invomega,help2);
			r8mat_mm_new(JY*JX,JY,JX*JY,help2,xi,incrxx);
			r8mat_mm_new(JY,JY*JZ,1,zi,uj,ziu);
			r8mat_divide(JY,1,-1,ziu);
			r8mat_add(JY,1,ziu,yi);
			r8mat_mm_new(JY*JX,JY,1,help2,yi,incrxy);
			r8mat_add(JY*JX,JY*JX,incrxx,sumxi);
			r8mat_add(JY*JX,1,incrxy,sumxy);
		}	
	}

	r8mat_pofac(JY * JX,sumxi,help3,2);
	r8mat_poinv(JY * JX, help3,invomega2);
	for (jj=1;jj<JX*JY;jj++) for (tt=0;tt<jj;tt++) invomega2[jj+JX*JY*tt]=invomega2[tt+JX*JY*jj];
	r8mat_mm_new(JY*JX,JY*JX,1,invomega2,sumxy,mu);
	r8mat_pofac(JY * JX,invomega2,help3,3);
	r8vec_multinormal_sample(JY*JX, mu,help3, REAL(beta),newbeta, fl);
	r8mat_add(Ib,Jb,REAL(beta),REAL(betapost));
	for (c=0;c<nj;c++) {
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) invomega[j+t*JY]=allinvomega[(c*JY+j)+t*(JY*nj)];
			}
		for (j=0;j<JY*JY*JZ*JZ;j++) sumzi[j]=0;
		for (j=0;j<JY*JZ;j++) sumzy[j]=0;
		for (j=0;j<IY;j++) {
			if (INTEGER(clus)[j]==c) {
				for (t=0;t<JY*JZ*JY;t++) zi[t]=0;
				for (t=0;t<JY*JX*JY;t++) xi[t]=0; 
				for (t=0;t<JY;t++) {
					yi[t]=imp[j+t*IY];
					for (k=0;k<JX;k++) {
						xi[t+(k+t*JX)*JY]=REAL(X)[j+IX*k];
					}
					for (k=0;k<JZ;k++) {
						zi[t+(k+t*JZ)*JY]=REAL(Z)[j+IZ*k];
					}
				}
				r8mat_mtm_new(JY*JZ,JY,JY,zi,invomega,help4);
				r8mat_mm_new(JY*JZ,JY,JZ*JY,help4,zi,incrzz);
				r8mat_mm_new(JY,JY*JX,1,xi,REAL(beta),ziu);
				r8mat_divide(JY,1,-1,ziu);
				r8mat_add(JY,1,ziu,yi);
				r8mat_mm_new(JY*JZ,JY,1,help4,yi,incrzy);
				r8mat_add(JY*JZ,JY*JZ,incrzz,sumzi);
				r8mat_add(JY*JZ,1,incrzy,sumzy);
				
			}
		}
		
		r8mat_pofac(JY * JZ,REAL(covu),help5,4);
		r8mat_poinv(JY * JZ, help5,invomega3);
		for (jj=1;jj<JZ*JY;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+JZ*JY*tt]=invomega3[tt+JZ*JY*jj];
		r8mat_add(JY*JZ,JY*JZ,invomega3,sumzi);
		
		r8mat_pofac(JY * JZ,sumzi,help5,5);
		r8mat_poinv(JY * JZ, help5,invomega3);
		for (jj=1;jj<JZ*JY;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+JZ*JY*tt]=invomega3[tt+JZ*JY*jj];
		r8mat_mm_new(JY*JZ,JY*JZ,1,invomega3,sumzy,mu2);

		r8mat_pofac(JY * JZ,invomega3,help5,6); 
		r8vec_multinormal_sample(JY*JZ, mu2,help5,newu, incrzy, fl);
		for (t=0;t<JY;t++) for (k=0;k<JZ;k++) REAL(u)[c+nj*(k+t*JZ)] = newu[k+t*JZ];
		
	}
	r8mat_add(Iu,Ju,REAL(u),REAL(upost));


	for (j=0;j<JY*JY*JZ*JZ;j++) mu3[j]=0;
	for (j=0;j<nj; j++) {
		for (t=0;t<JY;t++) for (k=0;k<JZ;k++) uj[k+t*JZ]=REAL(u)[j+nj*(k+t*JZ)];
		r8mat_mmt_new(JY*JZ,1,JY*JZ,uj,uj,help5);
		r8mat_add(JY*JZ,JY*JZ,help5,mu3);
		
	}
	r8mat_add(JY*JZ,JY*JZ,REAL(Sup),mu3);
	r8mat_pofac(JY*JZ,mu3,help5,7);
	r8mat_poinv(JY*JZ, help5, invomega3);
	for (jj=1;jj<(JY*JZ);jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(JZ*JY)*tt]=invomega3[tt+(JZ*JY)*jj];
	wishart_sample(JY*JZ,nj-1,invomega3,newomega, help5,sumzi,incrzz,mu3, fl);
	
	r8mat_pofac(JY * JZ,newomega, help5,8);
	r8mat_poinv(JY * JZ, help5,invomega3);
	for (jj=1;jj<(JY*JZ);jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(JZ*JY)*tt]=invomega3[tt+(JZ*JY)*jj];
	for(k=0;k<(JY*JZ);k++)  for(j=0;j<(JY*JZ);j++)  REAL(covu)[j+(JZ*JY)*k]=invomega3[j+(JZ*JY)*k];
	r8mat_add(JY*JZ,JY*JZ,REAL(covu),REAL(covupost));
	
	for (c=0;c<nj;c++) {	
		for (j=0;j<JY*JY;j++) mu4[j]=0;
		for (j=cumclus[c];j<cumclus[c+1];j++) {
			for (t=0;t<JY*JZ*JY;t++) zi[t]=0;
			for (t=0;t<JY*JX*JY;t++) xi[t]=0; 
			for (t=0;t<JY;t++) {
				yi[t]=imp[j+t*IY];
				for (k=0;k<JX;k++) {
					xi[t+(k+t*JX)*JY]=REAL(X)[j+IY*k];
				}
				for (k=0;k<JZ;k++) {
					zi[t+(k+t*JZ)*JY]=REAL(Z)[j+IY*k];
					uj[k+t*JZ]=REAL(u)[c+nj*(k+t*JZ)];
				}
			}
			r8mat_mm_new(JY,JY*JZ,1,zi,uj,ziu);
			r8mat_divide(JY,1,-1,ziu);
			r8mat_add(JY,1,ziu,yi);
			r8mat_mm_new(JY,JY*JX,1,xi,REAL(beta),ziu);
			r8mat_divide(JY,1,-1,ziu);
			r8mat_add(JY,1,ziu,yi);
			r8mat_mmt_new(JY,1,JY,yi,yi,help);
			r8mat_add(JY,JY,help,mu4);		
		}
		r8mat_add(JY,JY,REAL(Sp),mu4);
	
		r8mat_pofac(JY,mu4,help,9);
		r8mat_poinv(JY, help,invomega);
		for (jj=1;jj<JY;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY*tt]=invomega[tt+JY*jj];
		wishart_sample(JY,clusnum[c],invomega,newomega2,help, omegaoo,omegaom,omegamm, fl);	
		r8mat_pofac(JY,newomega2,help,10);
		r8mat_poinv(JY, help,invomega);
		for (jj=1;jj<JY;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY*tt]=invomega[tt+JY*jj];
		for(k=0;k<JY;k++)  for(j=0;j<JY;j++)  REAL(omega)[(c*JY+j)+k*(JY*nj)]=invomega[j+JY*k];
	}
	r8mat_add(Io,Jo,REAL(omega),REAL(omegapost));

	for (j=0; j<IY; j++) {
		for (k=0;k<JY;k++) betaX[k]=0;		
		for (t=0;t<JY*JZ*JY;t++) zi[t]=0;
		for (t=0;t<JY*JX*JY;t++) xi[t]=0; 
		for (t=0;t<JY;t++) {
			yi[t]=imp[j+t*IY];
			for (k=0;k<JX;k++) {
				xi[t+(k+t*JX)*JY]=REAL(X)[j+IY*k];
			}
			for (k=0;k<JZ;k++) {
				zi[t+(k+t*JZ)*JY]=REAL(Z)[j+IY*k];
				uj[k+t*JZ]=REAL(u)[(INTEGER(clus)[j])+nj*(k+t*JZ)];
			}
		}
		r8mat_mm_new(JY,JY*JZ,1,zi,uj,ziu);
		r8mat_mm_new(JY,JY*JX,1,xi,REAL(beta),help6);
		r8mat_add(JY,1,ziu,betaX);
		r8mat_add(JY,1,help6,betaX);
		
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
						omegamm[countmm]=REAL(omega)[(INTEGER(clus)[j]*JY+t)+k*(JY*nj)];
						countmm++;
					}
					else if (ISNAN(REAL(Y)[j+t*IY])) {
						omegamo[countmo]=REAL(omega)[(INTEGER(clus)[j]*JY+k)+t*(JY*nj)];	
						countmo++;
						flag=1;
					}
					else if (!ISNAN(REAL(Y)[j+k*IY])&!ISNAN(REAL(Y)[j+t*IY])) {
						omegaoo[countoo]=REAL(omega)[(INTEGER(clus)[j]*JY+t)+k*(JY*nj)];
						countoo++;	
					}
				}
				if (flag==1)countt++;
				flag=0;
			}
			r8mat_pofac((JY-nmiss),omegaoo,help7,11);
			r8mat_poinv((JY-nmiss),help7,invomega4);
			for (jj=1;jj<(JY-nmiss);jj++) for (tt=0;tt<jj;tt++) invomega4[jj+(JY-nmiss)*tt]=invomega4[tt+(JY-nmiss)*jj];
			r8mat_mmt_new((JY-nmiss),(JY-nmiss),nmiss,invomega4,omegamo,help8);
			r8mat_divide((JY-nmiss),1,-1,betaobs);
			r8mat_add((JY-nmiss),1,betaobs,Yobs);
			r8mat_mtm_new(1,(JY-nmiss),nmiss,Yobs,help8,mumiss);
			r8mat_add(1,nmiss,betamiss,mumiss);
			r8mat_mm_new(nmiss,(JY-nmiss),nmiss,omegamo,help8,omegadrawmiss);
			r8mat_divide(nmiss,nmiss,-1,omegadrawmiss);
			r8mat_add(nmiss,nmiss,omegamm,omegadrawmiss);
			r8mat_pofac(nmiss,omegadrawmiss,help9,12);
			r8vec_multinormal_sample(nmiss,mumiss,help9,Ymiss,help6, fl);
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
r8mat_divide(Iu,Ju,ns,REAL(upost));
r8mat_divide(Io,Jo,ns,REAL(omegapost));
r8mat_divide(JY*JZ,JY*JZ,ns,REAL(covupost));
PutRNGstate();
UNPROTECT(24);
return R_NilValue;
}
