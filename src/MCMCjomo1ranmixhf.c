#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "pdflib.h"
#include "wishart.h"
#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>

SEXP MCMCjomo1ranmixhf(SEXP Y, SEXP Yimp, SEXP Yimp2, SEXP Yimpcat, SEXP X, SEXP Z, SEXP clus, SEXP beta, SEXP u, SEXP betapost, SEXP upost, SEXP omega,SEXP omegapost, SEXP covu, SEXP covupost, SEXP nstep, SEXP Sp, SEXP Sup, SEXP Y_numcat, SEXP num_con, SEXP a_start, SEXP flagrng, SEXP collectimp){
int indic=0,i,j,k, IY,JY, IX, JX, Io, Jo, Ib, Jb, ns, nmiss=0,t, countm=0, counto=0, countmm=0, countmo=0,countoo=0, jj, tt, kk, ncon,ncat, pos,flag=0,nmaxx,h;
int Iu, Ju, IZ, JZ, nj,c, fl,currncat;
SEXP RdimY, RdimX, Rdimo, Rdimb, RdimZ, Rdimu;
double *betaX, *Yobs, *Ymiss, *mumiss, *omegadrawmiss, *betamiss, *betaobs, *omegaoo, *omegamo, *omegamm, *invomega, *invomega2, *help, *help2, *help3, *imp, *zi;
double *sumzy, *incrzz, *incrzy, *mu, *mu2, *newbeta, *newomega, *sumzi, *yi, *invomega3, *help4, *help5, *help6, *missing, *fixomega,*resid, sdom, meanom, detom,logLH, newlogLH;
double maxx,maxim,maxim2, *sumxy, *sumxi, *uj, *xi, *ziu, *incrxx, *incrxy, *newu, *mu3,*mu4, *help7, *help8, *help9, *invomega4, *newomega2,a, *clusnum;
double *cumclus, *allinvomega;

/* Protecting R objects from garbage collection and saving matrices dimensions*/ 

RdimY=PROTECT(getAttrib(Yimp,R_DimSymbol));
IY=INTEGER(RdimY)[0];
JY=INTEGER(RdimY)[1];
RdimX=PROTECT(getAttrib(X,R_DimSymbol));
IX=INTEGER(RdimX)[0];
JX=INTEGER(RdimX)[1];
RdimZ=PROTECT(getAttrib(Z,R_DimSymbol));
IZ=INTEGER(RdimZ)[0];
JZ=INTEGER(RdimZ)[1];
Rdimb=PROTECT(getAttrib(beta,R_DimSymbol));
Ib=INTEGER(Rdimb)[0];
Jb=INTEGER(Rdimb)[1];
Rdimo=PROTECT(getAttrib(omega,R_DimSymbol));
Io=INTEGER(Rdimo)[0];
Jo=INTEGER(Rdimo)[1];
Rdimu=PROTECT(getAttrib(u,R_DimSymbol));
Iu=INTEGER(Rdimu)[0];
Ju=INTEGER(Rdimu)[1];
Y=PROTECT(coerceVector(Y,REALSXP));
Yimpcat=PROTECT(coerceVector(Yimpcat,REALSXP));
Y_numcat=PROTECT(coerceVector(Y_numcat,INTSXP));
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
num_con=PROTECT(coerceVector(num_con,INTSXP));
ncon=INTEGER(num_con)[0];
flagrng=PROTECT(coerceVector(flagrng,INTSXP));
fl=INTEGER(flagrng)[0];
ncat=length(Y_numcat);
a_start=PROTECT(coerceVector(a_start,REALSXP));
a=REAL(a_start)[0];
nj=Iu;
collectimp=PROTECT(coerceVector(collectimp,REALSXP));

/*Allocating memory for C objects in R*/

help = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
invomega= (double * ) R_alloc ( JY * JY , sizeof ( double ) );
fixomega = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
sumzi = ( double * ) R_alloc ( JY * JZ * JY * JZ , sizeof ( double ) );
sumzy = ( double * ) R_alloc ( JY * JZ , sizeof ( double ) );
sumxi = ( double * ) R_alloc ( JY * JX * JY * JX , sizeof ( double ) );
sumxy = ( double * ) R_alloc ( JY * JX , sizeof ( double ) );
zi = ( double * ) R_alloc ( JY * JZ * JY , sizeof ( double ) );
xi = ( double * ) R_alloc ( JY * JX * JY , sizeof ( double ) );
yi = ( double * ) R_alloc ( JY , sizeof ( double ) );
uj = ( double * ) R_alloc ( JY * JZ , sizeof ( double ) );
newu = ( double * ) R_alloc ( JY * JZ , sizeof ( double ) );
ziu = ( double * ) R_alloc ( JY , sizeof ( double ) );
help2 = ( double * ) R_alloc ( JY *JX * JY , sizeof ( double ) );
incrzz = ( double * ) R_alloc ( JY *JZ * JY *JZ , sizeof ( double ) );
incrzy = ( double * ) R_alloc ( JY *JZ , sizeof ( double ) );
incrxx = ( double * ) R_alloc ( JY *JX * JY *JX , sizeof ( double ) );
incrxy = ( double * ) R_alloc ( JY *JX , sizeof ( double ) );
help3 = ( double * ) R_alloc ( JY * JX * JY * JX ,sizeof ( double ) );
invomega2= (double * ) R_alloc ( JY * JX * JY * JX , sizeof ( double ) );
mu = ( double * ) R_alloc ( JY * JX, sizeof ( double ) );
newbeta = ( double * ) R_alloc ( JY * JX ,sizeof ( double ) );
mu2 = ( double * ) R_alloc ( JY * JY * JZ,sizeof ( double ) );
mu3 = ( double * ) R_alloc ( JY * JY * JZ * JZ ,sizeof ( double ) );
mu4 = ( double * ) R_alloc ( JY * JY,sizeof ( double ) );
newomega = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
newomega2 = ( double * ) R_alloc ( JY * JY *JZ * JZ , sizeof ( double ) );
betaX=( double * ) R_alloc ( JY , sizeof ( double ) );
imp=( double * ) R_alloc ( IY * JY,sizeof ( double ) );
resid=( double * ) R_alloc ( IY * JY,sizeof ( double ) );
Yobs=( double * ) R_alloc ( JY , sizeof ( double ) );
Ymiss=( double * ) R_alloc ( JY , sizeof ( double ) );
mumiss = ( double * ) R_alloc ( JY , sizeof ( double ) );
omegadrawmiss = ( double * ) R_alloc ( JY * JY ,sizeof ( double ) );
betamiss = ( double * ) R_alloc ( JY ,sizeof ( double ) );
betaobs = ( double * ) R_alloc ( JY, sizeof ( double ) );
omegaoo= ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
omegamo= ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
omegamm= ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
invomega3= ( double * ) R_alloc ( JY*JY * JZ * JZ , sizeof ( double ) );
invomega4 = ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
help4 = ( double * ) R_alloc ( JY *JY*JZ , sizeof ( double ) );
help5 = ( double * ) R_alloc ( JY * JY*JZ*JZ , sizeof ( double ) );
help6 = ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
help7 = ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
help8 = ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
help9 = ( double * ) R_alloc ( JY *JY , sizeof ( double ) );
missing = ( double * ) R_alloc ( IY , sizeof ( double ) );
clusnum = ( double * )  R_alloc ( nj , sizeof(double));
cumclus = ( double * )  R_alloc ( nj+1 , sizeof(double));
allinvomega = ( double * )  R_alloc ( nj* JY * JY , sizeof(double));




/* Some initializations */

GetRNGstate();

for (j=0; j<IY; j++) {
	missing[j]=0;
	for (k=0;k<JY;k++) {
		if (ISNAN(REAL(Yimp)[j+k*IY])) {
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

r8mat_copy_new(IY, JY, REAL(Yimp2), imp);
for (i=0;i<JY*JY;i++) fixomega[i]=0;
pos=ncon;
for (i=0;i<ncat;i++) {
	for (j=0;j<(INTEGER(Y_numcat)[i]-1);j++) {
		for (k=0;k<(INTEGER(Y_numcat)[i]-1);k++) {
			if (j==k) for (t=0;t<nj;t++) REAL(omega)[(pos+j+t*JY)+JY*nj*(pos+k)]=1;
			else for (t=0;t<nj;t++) REAL(omega)[(pos+j+t*JY)+JY*nj*(pos+k)]=0.5;
			fixomega[(pos+j)+Jo*(pos+k)]=1;
		}
	}
	pos=pos+INTEGER(Y_numcat)[i]-1;
}

/* Running ns iterations of Gibbs sampler*/

for (i=0;i<ns;i++) {
	for (c=0;c<nj;c++) {
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) newomega[j+t*JY]=REAL(omega)[(c*JY+j)+t*(JY*nj)];
			}
		r8mat_pofac(JY,newomega, help,2);
		r8mat_poinv(JY,help, invomega);
		for (jj=1;jj<JY;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY*tt]=invomega[tt+JY*jj];
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) allinvomega[(c*JY+j)+t*(JY*nj)]=invomega[j+t*JY];
			}
		}
	for (j=0;j<JY*JY*JX*JX;j++) sumxi[j]=0;	
	for (j=0;j<JY*JX;j++) sumxy[j]=0;
	
	// Rejection Sampling

	pos=ncon;
	for (j=0;j<ncat;j++) {
		currncat=INTEGER(Y_numcat)[j]-1;
		for (t=0;t<IY;t++) {
			if (!ISNAN(REAL(Y)[t+(ncon+j)*IY])) {
				for (k=0;k<currncat;k++) betaX[k]=0;
				for (tt=0;tt<JX;tt++) {
					for (k=0;k<currncat;k++) {
						betaX[k]=betaX[k]+REAL(beta)[tt+(k+pos)*Ib]*REAL(X)[t+tt*IX];
					}
				}
				for (tt=0;tt<JZ;tt++) {
					for (k=0;k<currncat;k++) {
						betaX[k]=betaX[k]+REAL(u)[(INTEGER(clus)[t])+nj*(tt+(k+pos)*JZ)]*REAL(Z)[t+tt*IZ];
					}
				}

				for (k=0;k<(JY-currncat);k++) help[k]=0;
				for (tt=0;tt<JX;tt++) {
					for (k=0;k<pos;k++) {
						help[k]=help[k]+REAL(beta)[tt+k*Ib]*REAL(X)[t+tt*IX];
					}
					for (k=(pos+currncat);k<JY;k++) {
						help[k-currncat]=help[k-currncat]+REAL(beta)[tt+k*Ib]*REAL(X)[t+tt*IX];
					}
				}
				for (tt=0;tt<JZ;tt++) {
					for (k=0;k<pos;k++) {
						help[k]=help[k]+REAL(u)[(INTEGER(clus)[t])+nj*(tt+k*JZ)]*REAL(Z)[t+tt*IZ];
					}
					for (k=(pos+currncat);k<JY;k++) {
						help[k-currncat]=help[k-currncat]+REAL(u)[(INTEGER(clus)[t])+nj*(tt+k*JZ)]*REAL(Z)[t+tt*IZ];
					}
				}

				for (k=0;k<pos;k++) {
					help[k]=imp[t+k*IY]-help[k];
				}
				
				for (k=(pos+currncat);k<JY;k++) {
					help[k-currncat]=imp[t+k*IY]-help[k-currncat];
				}
				for (k=0;k<JY;k++) {
					for (kk=0;kk<JY;kk++) {
						if (((kk<pos)||(kk>(pos+currncat-1)))&&((k<pos)||(k>(pos+currncat-1)))) {
							help4[countm]=REAL(omega)[kk+JY*INTEGER(clus)[t]+JY*nj*k];
							countm++;
						}
						else if (((kk<pos)||(kk>(pos+currncat-1)))&&((k>(pos-1))||(k<(pos+currncat)))) {
							help5[counto]=REAL(omega)[kk+JY*INTEGER(clus)[t]+JY*nj*k];
							counto++;
						}
					}
				}
		
				countm=0;
				counto=0;
				r8mat_pofac((JY-currncat),help4, help6,1);
				r8mat_poinv((JY-currncat),help6, invomega);
				for (jj=1;jj<(JY-currncat);jj++) for (tt=0;tt<jj;tt++) invomega[jj+(JY-currncat)*tt]=invomega[tt+(JY-currncat)*jj];

				r8mat_mm_new((JY-currncat),(JY-currncat),currncat,invomega,help5, help2);

				r8mat_mm_new(1,(JY-currncat),currncat,help,help2, mumiss);
				r8mat_add(currncat,1,betaX,mumiss);
				r8mat_mtm_new(currncat,(JY-currncat),currncat,help2,help5, omegadrawmiss);
				r8mat_divide(currncat,currncat,-1,omegadrawmiss);
				for (k=0;k<currncat;k++) {
					omegadrawmiss[k+k*currncat]=omegadrawmiss[k+k*currncat]+1;
					for (kk=0;kk<currncat;kk++) if (k!=kk) omegadrawmiss[k+kk*currncat]=omegadrawmiss[k+kk*currncat]+0.5;
				}
				r8mat_pofac(currncat,omegadrawmiss, omegamm,2);
				flag=0;
				kk=0;

				if (REAL(Y)[t+(ncon+j)*IY]==INTEGER(Y_numcat)[j]) {
					while ((flag==0)&(kk<10000)) {
						r8vec_multinormal_sample((INTEGER(Y_numcat)[j]-1), mumiss,omegamm, newbeta,mu4,fl);
						maxim=maxvec((INTEGER(Y_numcat)[j]-1),newbeta);
						if (maxim<0) {
						
							for (k=0;k<(INTEGER(Y_numcat)[j]-1);k++) imp[t+(k+pos)*IY]=newbeta[k];
							flag=1;
							indic++;
						}
						kk++;
					}
				}
				else {
					while ((flag==0)&(kk<10000)) {
						r8vec_multinormal_sample((INTEGER(Y_numcat)[j]-1), mumiss,omegamm, newbeta,mu4,fl);
						maxim=argmaxvec((INTEGER(Y_numcat)[j]-1),newbeta);
						maxim2=maxvec((INTEGER(Y_numcat)[j]-1),newbeta);
						if (((maxim+1)==REAL(Y)[t+(ncon+j)*IY])&(maxim2>0)) {
							for (k=0;k<(INTEGER(Y_numcat)[j]-1);k++) imp[t+(k+pos)*IY]=newbeta[k];
							flag=1;
							indic++;
						}
						kk++;
					}
				}
				flag=0;

			}
		}
		pos=pos+INTEGER(Y_numcat)[j]-1;
	}

	//Updating beta


	for (j=0;j<JY*JY*JX*JX;j++) sumxi[j]=0;	
	for (j=0;j<JY*JX;j++) sumxy[j]=0;
	for (j=0;j<IY;j++) {
		for (t=0;t<JY;t++) for (k=0;k<JY;k++) invomega[t+k*JY]=allinvomega[INTEGER(clus)[j]*JY+t+k*nj*JY];
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
	r8mat_pofac(JY * JX,sumxi,help3,5);
	r8mat_poinv(JY * JX, help3,invomega2);
	for (jj=1;jj<JX*JY;jj++) for (tt=0;tt<jj;tt++) invomega2[jj+JX*JY*tt]=invomega2[tt+JX*JY*jj];
	r8mat_mm_new(JY*JX,JY*JX,1,invomega2,sumxy,mu);
	r8mat_pofac(JY * JX,invomega2,help3,6);
	r8vec_multinormal_sample(JY*JX, mu,help3, REAL(beta),newbeta,fl);
	
	for (j=0;j<Ib;j++) {
		for (t=0;t<Jb;t++) {
			REAL(betapost)[j+Ib*t+i*Ib*Jb]=REAL(beta)[j+Ib*t];
			}
		}

	//Updating random effects

	for (c=0;c<nj;c++) {
		for (jj=0;jj<JY;jj++) for (tt=0;tt<JY;tt++) invomega[jj+JY*tt]=allinvomega[c*JY+jj+JY*nj*tt];
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
		
		r8mat_pofac(JY * JZ,REAL(covu),help5,7);
		r8mat_poinv(JY * JZ, help5,invomega3);
		for (jj=1;jj<JZ*JY;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+JZ*JY*tt]=invomega3[tt+JZ*JY*jj];
		r8mat_add(JY*JZ,JY*JZ,invomega3,sumzi);
		
		r8mat_pofac(JY * JZ,sumzi,help5,8);
		r8mat_poinv(JY * JZ, help5,invomega3);
		for (jj=1;jj<JZ*JY;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+JZ*JY*tt]=invomega3[tt+JZ*JY*jj];
		r8mat_mm_new(JY*JZ,JY*JZ,1,invomega3,sumzy,mu2);
		r8mat_pofac(JY * JZ,invomega3,help5,9); 
		r8vec_multinormal_sample(JY*JZ, mu2,help5,newu, incrzy,fl);
		for (t=0;t<JY;t++) for (k=0;k<JZ;k++) REAL(u)[c+nj*(k+t*JZ)] = newu[k+t*JZ];
		
	}
	for (j=0;j<Iu;j++) {
		for (t=0;t<Ju;t++) {
			REAL(upost)[j+Iu*t+i*Iu*Ju]=REAL(u)[j+Iu*t];
			}
		}

	//Updating level 2 covariance matrix


	for (j=0;j<JY*JY*JZ*JZ;j++) mu3[j]=0;
	for (j=0;j<nj; j++) {
		for (t=0;t<JY;t++) for (k=0;k<JZ;k++) uj[k+t*JZ]=REAL(u)[j+nj*(k+t*JZ)];
		r8mat_mmt_new(JY*JZ,1,JY*JZ,uj,uj,help5);
		r8mat_add(JY*JZ,JY*JZ,help5,mu3);
		
	}
	r8mat_add(JY*JZ,JY*JZ,REAL(Sup),mu3);
	r8mat_pofac(JY*JZ,mu3,help5,10);
	r8mat_poinv(JY*JZ, help5, invomega3);
	for (jj=1;jj<(JY*JZ);jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(JZ*JY)*tt]=invomega3[tt+(JZ*JY)*jj];
	wishart_sample(JY*JZ,(nj+JY*JZ),invomega3,newomega2, help5,sumzi,incrzz,mu3,fl);
	
	r8mat_pofac(JY * JZ,newomega2, help5,11);
	r8mat_poinv(JY * JZ, help5,invomega3);
	for (jj=1;jj<(JY*JZ);jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(JZ*JY)*tt]=invomega3[tt+(JZ*JY)*jj];
	for(k=0;k<(JY*JZ);k++)  for(j=0;j<(JY*JZ);j++)  REAL(covu)[j+(JZ*JY)*k]=invomega3[j+(JZ*JY)*k];
	
	for (j=0;j<JY*JZ;j++) {
		for (t=0;t<JY*JZ;t++) {
			REAL(covupost)[j+JY*JZ*t+i*JY*JY*JZ*JZ]=REAL(covu)[j+JY*JZ*t];
			}
		}

	//Updating residuals

	for (t=0;t<IY;t++) {
		for (j=0; j<JY; j++) {
			resid[t+j*IY]=imp[t+j*IY];
			for (k=0;k<JX;k++) resid[t+j*IY]=resid[t+j*IY]-REAL(beta)[k+j*JX]*REAL(X)[t+k*IX];
			for (k=0;k<JZ;k++) resid[t+j*IY]=resid[t+j*IY]-REAL(u)[INTEGER(clus)[t]+nj*(k+j*JZ)]*REAL(Z)[t+k*IZ];
		}
	}

	//Updating omega

	for (c=0;c<nj;c++) {	
	for (kk=0;kk<JY;kk++) {
		for (tt=0;tt<JY;tt++) {
			newomega[kk+tt*JY]=REAL(omega)[(c*JY+kk)+tt*(JY*nj)];
			}
		}
	flag=0;
	r8mat_pofac(JY,newomega,help,6);
	detom=r8mat_podet(JY, help);

	r8mat_poinv(JY, help,invomega);
	for (jj=1;jj<JY;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY*tt]=invomega[tt+JY*jj];

	logLH=0;
	for (t=cumclus[c];t<cumclus[c+1];t++) {
		for (jj=0;jj<JY;jj++) help4[jj]=resid[t+jj*IY];
		r8mat_mm_new(1,JY,JY,help4,invomega,help5);
		r8mat_mmt_new(1,JY,1,help5,help4,help6);
		logLH=logLH+help6[0];
	}
	logLH=logLH*(-0.5)-clusnum[c]*log(detom)/2;
	for (j=0;j<JY;j++) {
		for (k=0;k<JY;k++) {
			if ((fixomega[j+JY*k]==0)&(j<=k)) {
				if (j==k) {
					sdom=REAL(omega)[(c*JY+j)+k*(JY*nj)]*sqrt(11.6/clusnum[c]);
					meanom=REAL(omega)[(c*JY+j)+k*(JY*nj)];
				}
				else {
					sdom=0.1*sqrt(REAL(omega)[(c*JY+j)+j*(JY*nj)]*REAL(omega)[(c*JY+k)+k*(JY*nj)]);
					meanom=REAL(omega)[(c*JY+j)+k*(JY*nj)];
				}
				kk=0;
				while ((flag==0)&(kk<100)) {
					newomega[j+JY*k]=r8_normal_sample(meanom,sdom,fl);
					newomega[k+JY*j]=newomega[j+JY*k];
					flag=checkposdef(JY,newomega,help,help5);
					kk++;
				}
					
				r8mat_pofac(JY,newomega,help,7);
				detom=r8mat_podet(JY, help);
				r8mat_poinv(JY, help,invomega);
				for (jj=1;jj<JY;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY*tt]=invomega[tt+JY*jj];
				newlogLH=0;
				for (t=cumclus[c];t<cumclus[c+1];t++) {
					for (jj=0;jj<JY;jj++) help4[jj]=resid[t+jj*IY];
					r8mat_mm_new(1,JY,JY,help4,invomega,help5);
					r8mat_mmt_new(1,JY,1,help5,help4,help6);
					newlogLH=newlogLH+help6[0];
				}
				newlogLH=newlogLH*(-0.5)-clusnum[c]*log(detom)/2;
				if (((( double ) unif_rand ( ))<exp(newlogLH-logLH))&(flag==1)) {
					REAL(omega)[(c*JY+j)+k*(JY*nj)]=newomega[j+JY*k];
					REAL(omega)[(c*JY+k)+j*(JY*nj)]=newomega[k+JY*j];
					logLH=newlogLH;
				}
				else {
					newomega[j+JY*k]=REAL(omega)[(c*JY+j)+k*(JY*nj)];
					newomega[k+JY*j]=REAL(omega)[(c*JY+k)+j*(JY*nj)];
				}
				flag=0;
			}
		}
	}

	}

	for (j=0;j<Io;j++) {
		for (t=0;t<Jo;t++) {
			REAL(omegapost)[j+Io*t+i*Jo*Io]=REAL(omega)[j+Io*t];
			}
		}

	//imputing missing values

	flag=0;
	for (j=0; j<IY; j++) {
		for (t=0;t<JY;t++) for (k=0;k<JY;k++) invomega[t+k*JY]=REAL(omega)[INTEGER(clus)[j]*JY+t+k*nj*JY];
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
				if (ISNAN(REAL(Yimp)[j+k*IY])) {
					betamiss[countm]=betaX[k];
					countm++;
				}
				else {
					Yobs[counto]=imp[j+k*IY];
					betaobs[counto]=betaX[k];
					counto++;
				}
				for (t=0;t<JY;t++) {
					if (ISNAN(REAL(Yimp)[j+k*IY])&ISNAN(REAL(Yimp)[j+t*IY])) {
						omegamm[countmm]=REAL(omega)[(INTEGER(clus)[j]*JY+t)+k*(JY*nj)];
						countmm++;
					}
					else if (ISNAN(REAL(Yimp)[j+t*IY])) {
						omegamo[countmo]=REAL(omega)[(INTEGER(clus)[j]*JY+k)+t*(JY*nj)];	
						countmo++;	
					}
					else if (!ISNAN(REAL(Yimp)[j+k*IY])&!ISNAN(REAL(Yimp)[j+t*IY])){
						omegaoo[countoo]=REAL(omega)[(INTEGER(clus)[j]*JY+t)+k*(JY*nj)];
						countoo++;	
					}
				}
			}
			r8mat_pofac((JY-nmiss),omegaoo,help7,15);
			r8mat_poinv((JY-nmiss),help7,invomega4);
			for (jj=1;jj<JY-nmiss;jj++) for (tt=0;tt<jj;tt++) invomega4[jj+(JY-nmiss)*tt]=invomega4[tt+(JY-nmiss)*jj];
			r8mat_mmt_new((JY-nmiss),(JY-nmiss),nmiss,invomega4,omegamo,help8);
			r8mat_divide((JY-nmiss),1,-1,betaobs);
			r8mat_add((JY-nmiss),1,betaobs,Yobs);
			r8mat_mtm_new(1,(JY-nmiss),nmiss,Yobs,help8,mumiss);
			r8mat_add(1,nmiss,betamiss,mumiss);
			r8mat_mm_new(nmiss,(JY-nmiss),nmiss,omegamo,help8,omegadrawmiss);
			r8mat_divide(nmiss,nmiss,-1,omegadrawmiss);
			r8mat_add(nmiss,nmiss,omegamm,omegadrawmiss);
			r8mat_pofac(nmiss,omegadrawmiss,help9,16+1000*j);
			r8vec_multinormal_sample(nmiss,mumiss,help9,Ymiss,help6,fl);
			countm=0;
			for (k=0;k<JY;k++) {
				if (ISNAN(REAL(Yimp)[j+k*IY])) {
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
	
	for (j=0;j<IY;j++) {
			for (t=0;t<JY;t++) {
				REAL(collectimp)[j+IY*t+i*IY*JY]=imp[j+IY*t];
			}
		}
	if ((i+1)%10==0) Rprintf("Iteration %d completed\n",i+1);
}
for(i=0;i<IY;i++)  {
	for(j=0;j<JY;j++)  {
		REAL(Yimp2)[i+IY*j]=imp[i+IY*j];
	}
	h=0;
	for (j=0;j<ncat;j++) {
		maxx=imp[i+(ncon+h)*IY];
		nmaxx=0;
		for (k=1;k<(INTEGER(Y_numcat)[j]-1);k++) {
			if (imp[i+(ncon+h+k)*IY]>maxx) {
				maxx=imp[i+(ncon+h+k)*IY];
				nmaxx=k;
			}
		}
		if (maxx>0) REAL(Yimpcat)[i+IY*j]=nmaxx+1;
		else REAL(Yimpcat)[i+IY*j]=INTEGER(Y_numcat)[j];
		h=h+INTEGER(Y_numcat)[j]-1;
	}

}
REAL(a_start)[0]=a;
PutRNGstate();
UNPROTECT(29);
return R_NilValue;
}
