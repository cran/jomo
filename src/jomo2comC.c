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

SEXP jomo2comC(SEXP Y, SEXP Yimp, SEXP Yimp2, SEXP Yimpcat, SEXP Y2, SEXP Y2imp, SEXP Y2imp2, SEXP Y2impcat, SEXP X, SEXP X2, SEXP Z, SEXP clus, SEXP beta, SEXP beta2, SEXP u, SEXP betapost, SEXP beta2post, SEXP upost, SEXP omega, SEXP omegapost, SEXP covu, SEXP covupost, SEXP nstep, SEXP Sp, SEXP Sup, SEXP Y_numcat, SEXP Y2_numcat, SEXP num_con, SEXP num_con2, SEXP flagrng){
int indic=0,i,j,k, IY,JY, IX, JX, Io, Jo, Ib, Jb, ns, nmiss=0,t, countm=0, counto=0, countmm=0, countmo=0,countoo=0, jj, tt, kk, ncon,ncat, pos,flag=0,nmaxx,h;
int Iu, Ju, IZ, JZ, nj,c,fl, currncat, IY2, JY2, IX2,JX2,Ib2, Jb2, ncon2, ncat2, JYm, JXm;
SEXP RdimY, RdimX, Rdimo, Rdimb, RdimZ, Rdimu, RdimY2, RdimX2, Rdimb2;
double *betaX, *Yobs, *Ymiss, *mumiss, *omegadrawmiss, *betamiss, *betaobs, *omegaoo, *omegamo, *omegamm, *invomega, *invomega2, *help, *help2, *help3, *imp,*imp2, *zi;
double *sumzy, *incrzz, *incrzy, *mu, *mu2, *newbeta, *newomega, *sumzi, *yi, *invomega3, *help4, *help5, *help6, *missing, *fixomega,meanom,sdom, *resid, logLH, newlogLH,detom;
double maxx,maxim,maxim2, *sumxy, *sumxi, *uj, *xi, *ziu, *incrxx, *incrxy, *newu, *mu3,  *mu4, *help7, *help8, *help9, *invomega4, *newomega2,*missing2, *Y2red, *X2red;
double *covu1, *covu2, *covu12, *fixomega2, *Y2impred;

/* Protecting R objects from garbage collection and saving matrices dimensions*/ 

RdimY=PROTECT(getAttrib(Yimp,R_DimSymbol));
IY=INTEGER(RdimY)[0];
JY=INTEGER(RdimY)[1];
RdimY2=PROTECT(getAttrib(Y2imp,R_DimSymbol));
IY2=INTEGER(RdimY2)[0];
JY2=INTEGER(RdimY2)[1];
RdimX=PROTECT(getAttrib(X,R_DimSymbol));
IX=INTEGER(RdimX)[0];
JX=INTEGER(RdimX)[1];
RdimX2=PROTECT(getAttrib(X2,R_DimSymbol));
IX2=INTEGER(RdimX2)[0];
JX2=INTEGER(RdimX2)[1];
RdimZ=PROTECT(getAttrib(Z,R_DimSymbol));
IZ=INTEGER(RdimZ)[0];
JZ=INTEGER(RdimZ)[1];
Rdimb=PROTECT(getAttrib(beta,R_DimSymbol));
Ib=INTEGER(Rdimb)[0];
Jb=INTEGER(Rdimb)[1];
Rdimb2=PROTECT(getAttrib(beta2,R_DimSymbol));
Ib2=INTEGER(Rdimb2)[0];
Jb2=INTEGER(Rdimb2)[1];
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
Y2=PROTECT(coerceVector(Y2,REALSXP));
Y2impcat=PROTECT(coerceVector(Y2impcat,REALSXP));
Y2_numcat=PROTECT(coerceVector(Y2_numcat,INTSXP));
Y2imp=PROTECT(coerceVector(Y2imp,REALSXP));
Y2imp2=PROTECT(coerceVector(Y2imp2,REALSXP));
X=PROTECT(coerceVector(X,REALSXP));
X2=PROTECT(coerceVector(X2,REALSXP));
Z=PROTECT(coerceVector(Z,REALSXP));
clus=PROTECT(coerceVector(clus,INTSXP));
beta=PROTECT(coerceVector(beta,REALSXP));
beta2=PROTECT(coerceVector(beta2,REALSXP));
u=PROTECT(coerceVector(u,REALSXP));
upost=PROTECT(coerceVector(upost,REALSXP));
betapost=PROTECT(coerceVector(betapost,REALSXP));
beta2post=PROTECT(coerceVector(beta2post,REALSXP));
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
num_con2=PROTECT(coerceVector(num_con2,INTSXP));
ncon2=INTEGER(num_con2)[0];
flagrng=PROTECT(coerceVector(flagrng,INTSXP));
fl=INTEGER(flagrng)[0];
if (REAL(Yimpcat)[0]==(-999)) ncat=0;
else ncat=length(Y_numcat);
if (REAL(Y2impcat)[0]==(-999)) ncat2=0;
else ncat2=length(Y2_numcat);
nj=Iu;
JXm=JX;
if (JX2>JX) JXm=JX2;
JYm=JY;
if (JY2>JY) JYm=JY2;

/*Allocating memory for C objects in R*/

help = ( double * ) R_alloc ( (Ju*Ju) , sizeof ( double ) );
invomega= (double * ) R_alloc ( (JYm*JYm) , sizeof ( double ) );
fixomega = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
fixomega2 = ( double * ) R_alloc ( Ju * Ju , sizeof ( double ) );
sumzi = ( double * ) R_alloc ( Ju * Ju , sizeof ( double ) );
sumzy = ( double * ) R_alloc ( JY * JZ , sizeof ( double ) );
sumxi = ( double * ) R_alloc ( JYm * JXm * JYm * JXm , sizeof ( double ) );
sumxy = ( double * ) R_alloc ( JYm * JXm , sizeof ( double ) );
zi = ( double * ) R_alloc ( JY * JZ * JY , sizeof ( double ) );
xi = ( double * ) R_alloc ( JYm * JXm * JYm , sizeof ( double ) );
yi = ( double * ) R_alloc ( JYm , sizeof ( double ) );
uj = ( double * ) R_alloc ( Ju * Ju , sizeof ( double ) );
newu = ( double * ) R_alloc ( JY * JZ , sizeof ( double ) );
ziu = ( double * ) R_alloc ( JY , sizeof ( double ) );
help2 = ( double * ) R_alloc ( JYm *JXm * Ju , sizeof ( double ) );
incrzz = ( double * ) R_alloc ( Ju * Ju, sizeof ( double ) );
incrzy = ( double * ) R_alloc ( JY *JZ , sizeof ( double ) );
incrxx = ( double * ) R_alloc ( JYm *JXm * JYm *JXm , sizeof ( double ) );
incrxy = ( double * ) R_alloc ( JYm *JXm , sizeof ( double ) );
help3 = ( double * ) R_alloc ( JYm * JXm * JYm * JXm ,sizeof ( double ) );
invomega2= (double * ) R_alloc ( Ju * JXm * Ju * JXm, sizeof ( double ) );
mu = ( double * ) R_alloc ( Ju * JXm, sizeof ( double ) );
newbeta = ( double * ) R_alloc ( JYm*JXm  ,sizeof ( double ) );
mu2 = ( double * ) R_alloc ( JY * JZ, sizeof ( double ) );
mu3 = ( double * ) R_alloc ( Ju * Ju ,sizeof ( double ) );
mu4 = ( double * ) R_alloc ( JYm * JYm,sizeof ( double ) );
newomega = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
newomega2 = ( double * ) R_alloc ( Ju*Ju , sizeof ( double ) );
betaX=( double * ) R_alloc ( JYm, sizeof ( double ) );
imp=( double * ) R_alloc ( IY * JY,sizeof ( double ) );
imp2=( double * ) R_alloc ( Iu * JY2,sizeof ( double ) );
Y2impred=( double * ) R_alloc ( Iu * JY2,sizeof ( double ) );
Y2red=( double * ) R_alloc ( Iu * JY2,sizeof ( double ) );
X2red=( double * ) R_alloc ( Iu * JX2,sizeof ( double ) );
resid=( double * ) R_alloc ( IY * JY,sizeof ( double ) );
Yobs=( double * ) R_alloc ( Ju, sizeof ( double ) );
Ymiss=( double * ) R_alloc ( JYm, sizeof ( double ) );
mumiss = ( double * ) R_alloc ( JYm, sizeof ( double ) );
omegadrawmiss = ( double * ) R_alloc ( JYm*JYm ,sizeof ( double ) );
betamiss = ( double * ) R_alloc ( JYm ,sizeof ( double ) );
betaobs = ( double * ) R_alloc ( Ju, sizeof ( double ) );
omegaoo= ( double * ) R_alloc ( Ju*Ju , sizeof ( double ) );
omegamo= ( double * ) R_alloc (  Ju*JYm  , sizeof ( double ) );
omegamm= ( double * ) R_alloc (  JYm*JYm  , sizeof ( double ) );
invomega3= ( double * ) R_alloc ( Ju * Ju , sizeof ( double ) );
invomega4 = ( double * ) R_alloc ( Ju *Ju , sizeof ( double ) );
help4 = ( double * ) R_alloc ( Ju*Ju , sizeof ( double ) );
help5 = ( double * ) R_alloc ( Ju*Ju , sizeof ( double ) );
help6 = ( double * ) R_alloc ( Ju*Ju , sizeof ( double ) );
help7 = ( double * ) R_alloc ( Ju*Ju , sizeof ( double ) );
help8 = ( double * ) R_alloc ( JYm *Ju , sizeof ( double ) );
help9 = ( double * ) R_alloc ( JYm *JYm , sizeof ( double ) );
missing = ( double * ) R_alloc ( IY , sizeof ( double ) );
missing2 = ( double * ) R_alloc ( Iu , sizeof ( double ) );
covu1= (double * ) R_alloc ( JY * JZ * JY * JZ , sizeof ( double ) );
covu2= (double * ) R_alloc ( JY2 * JY2, sizeof ( double ) );
covu12= (double * ) R_alloc ( JY * JZ * JY2 , sizeof ( double ) );

/* Some initializations */

for (j=0; j<IY; j++) {
	missing[j]=0;
	for (k=0;k<JY;k++) {
		if (ISNAN(REAL(Yimp)[j+k*IY])) {
			missing[j]++;
		}
	}
}
r8mat_copy_new(IY, JY, REAL(Yimp2), imp);
for (i=0;i<Iu;i++) {
	indic=0;
	j=0;
	while (indic==0) {
		if (INTEGER(clus)[j]==i) {
				for (k=0;k<JY2;k++) {
					imp2[i+Iu*k]=REAL(Y2imp2)[j+k*IY2];
					Y2impred[i+Iu*k]=REAL(Y2imp)[j+k*IY2];
				}			
				for (k=0;k<(ncon2+ncat2);k++) Y2red[i+Iu*k]=REAL(Y2)[j+k*IY2];
				for (k=0;k<JX2;k++) X2red[i+Iu*k]=REAL(X2)[j+k*IX2];
				indic++;
		}
		else j++;
	}
}
for (j=0; j<Iu; j++) {
	missing2[j]=0;
	for (k=0;k<JY2;k++) {
		if (ISNAN(Y2impred[j+k*Iu])) {
			missing2[j]++;
		}
	}
}
for (i=0;i<Ib*Jb;i++) REAL(betapost)[i]=0;
for (i=0;i<Ib2*Jb2;i++) REAL(beta2post)[i]=0;
for (i=0;i<JY*JY;i++) fixomega[i]=0;
for (i=0;i<Ju*Ju;i++) fixomega2[i]=0;
pos=ncon;
if (ncat>0) {
	for (i=0;i<ncat;i++) {
		for (j=0;j<(INTEGER(Y_numcat)[i]-1);j++) {
			for (k=0;k<(INTEGER(Y_numcat)[i]-1);k++) {
				if (j==k) REAL(omega)[(pos+j)+Jo*(pos+k)]=1;
				else REAL(omega)[(pos+j)+Jo*(pos+k)]=0.5;
				fixomega[(pos+j)+Jo*(pos+k)]=1;
			}
		}
		pos=pos+INTEGER(Y_numcat)[i]-1;
	}
}
pos=JY*JZ+ncon2;
if (ncat2>0) {
	for (i=0;i<ncat2;i++) {
		for (j=0;j<(INTEGER(Y2_numcat)[i]-1);j++) {
			for (k=0;k<(INTEGER(Y2_numcat)[i]-1);k++) {
				if (j==k) REAL(covu)[(pos+j)+Ju*(pos+k)]=1;
				else REAL(covu)[(pos+j)+Ju*(pos+k)]=0.5;
				fixomega2[(pos+j)+Ju*(pos+k)]=1;
			}
		}
		pos=pos+INTEGER(Y2_numcat)[i]-1;
	}
}
for (i=0;i<JY*JY;i++) newomega[i]=REAL(omega)[i];
for (i=0;i<Ju*Ju;i++) newomega2[i]=REAL(covu)[i];

GetRNGstate();

/* Running ns iterations of Gibbs sampler*/

for (i=0;i<ns;i++) {

	// Rejection Sampling
	if (ncat>0) {
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
									help4[countm]=REAL(omega)[kk+JY*k];
									countm++;
								}
								else if (((kk<pos)||(kk>(pos+currncat-1)))&&((k>(pos-1))||(k<(pos+currncat)))) {
									help5[counto]=REAL(omega)[kk+JY*k];
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
								r8vec_multinormal_sample((INTEGER(Y_numcat)[j]-1), mumiss,omegamm, newbeta,mu4,0);
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
								r8vec_multinormal_sample((INTEGER(Y_numcat)[j]-1), mumiss,omegamm, newbeta,mu4,0);
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

		}
		// Rejection sampling for level 2 variables
		
			if (ncat2>0) {
			pos=ncon2;
			for (j=0;j<ncat2;j++) {
				currncat=INTEGER(Y2_numcat)[j]-1;
				for (t=0;t<Iu;t++) {
					if (!ISNAN(Y2red[t+(ncon2+j)*Iu])) {
						for (k=0;k<currncat;k++) betaX[k]=0;
						for (tt=0;tt<JX2;tt++) {
							for (k=0;k<currncat;k++) {
								betaX[k]=betaX[k]+REAL(beta2)[tt+(k+pos)*Ib2]*X2red[t+tt*Iu];
							}
						}
						for (k=0;k<(pos+JY*JZ);k++) {
							help[k]=REAL(u)[t+k*Iu];
						}
						for (k=(pos+JY*JZ+currncat);k<Ju;k++) {
							help[k-currncat]=REAL(u)[t+k*Iu];
						}
						for (k=0;k<Ju;k++) {
							for (kk=0;kk<Ju;kk++) {
								if (((kk<(pos+JY*JZ))||(kk>(pos+JY*JZ+currncat-1)))&&((k<(pos+JY*JZ))||(k>(pos+JY*JZ+currncat-1)))) {
									help4[countm]=REAL(covu)[kk+k*Ju];
									countm++;
								}
								else if (((kk<(pos+JY*JZ))||(kk>(pos+JY*JZ+currncat-1)))&&((k>(pos+JY*JZ-1))||(k<(pos+JY*JZ+currncat)))) {
									help5[counto]=REAL(covu)[kk+k*Ju];
									counto++;
								}
							}
						}
						countm=0;
						counto=0;
						r8mat_pofac((Ju-currncat),help4, help6,1);
						r8mat_poinv((Ju-currncat),help6, invomega4);			
						for (jj=1;jj<(Ju-currncat);jj++) for (tt=0;tt<jj;tt++) invomega4[jj+(Ju-currncat)*tt]=invomega4[tt+(Ju-currncat)*jj];
						r8mat_mtm_new(currncat,(Ju-currncat),(Ju-currncat),help5,invomega4, help2);
						r8mat_mm_new(currncat,(Ju-currncat),1,help2,help, mumiss);
						r8mat_add(currncat,1,betaX,mumiss);
						r8mat_mm_new(currncat,(Ju-currncat),currncat,help2,help5, omegadrawmiss);
						r8mat_divide(currncat,currncat,-1,omegadrawmiss);
						for (k=0;k<currncat;k++) {
							omegadrawmiss[k+k*currncat]=omegadrawmiss[k+k*currncat]+1;
							for (kk=0;kk<currncat;kk++) if (k!=kk) omegadrawmiss[k+kk*currncat]=omegadrawmiss[k+kk*currncat]+0.5;
						}
						r8mat_pofac(currncat,omegadrawmiss, omegamm,2);
						flag=0;
						kk=0;
						if (Y2red[t+(ncon2+j)*Iu]==INTEGER(Y2_numcat)[j]) {
							while ((flag==0)&(kk<10000)) {
								r8vec_multinormal_sample((INTEGER(Y2_numcat)[j]-1), mumiss,omegamm, newbeta,mu4,0);
								maxim=maxvec((INTEGER(Y2_numcat)[j]-1),newbeta);
								if (maxim<0) {
						
									for (k=0;k<currncat;k++) {
										imp2[t+(k+pos)*Iu]=newbeta[k];
										REAL(u)[t+(JY*JZ+k+pos)*Iu]=newbeta[k]-betaX[k];
									}
									flag=1;
									indic++;
								}
							kk++;
							}
						}
						else {
							while ((flag==0)&(kk<10000)) {
								r8vec_multinormal_sample((INTEGER(Y2_numcat)[j]-1), mumiss,omegamm, newbeta,mu4,0);
								maxim=argmaxvec((INTEGER(Y2_numcat)[j]-1),newbeta);
								maxim2=maxvec((INTEGER(Y2_numcat)[j]-1),newbeta);
								if (((maxim+1)==Y2red[t+(ncon2+j)*Iu])&(maxim2>0)) {
									for (k=0;k<currncat;k++) {
										imp2[t+(k+pos)*Iu]=newbeta[k];
										REAL(u)[t+(JY*JZ+k+pos)*Iu]=newbeta[k]-betaX[k];
									}
									flag=1;
									indic++;
								}
								kk++;
							}
						}
						flag=0;

					}
				}
				pos=pos+INTEGER(Y2_numcat)[j]-1;
			}

		}

	//Updating beta

	r8mat_pofac(JY,REAL(omega), help,3);
	r8mat_poinv(JY,help, invomega);

	for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];

	for (j=0;j<JY*JY*JX*JX;j++) sumxi[j]=0;	
	for (j=0;j<JY*JX;j++) sumxy[j]=0;
	for (j=0;j<IY;j++) {
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
	r8mat_pofac(JY * JX,sumxi,help3,4);
	r8mat_poinv(JY * JX, help3,invomega2);
	for (jj=1;jj<JX*JY;jj++) for (tt=0;tt<jj;tt++) invomega2[jj+JX*JY*tt]=invomega2[tt+JX*JY*jj];
	r8mat_mm_new(JY*JX,JY*JX,1,invomega2,sumxy,mu);
	r8mat_pofac(JY * JX,invomega2,help3,5);
	r8vec_multinormal_sample(JY*JX, mu,help3, REAL(beta),newbeta,0);
	r8mat_add(Ib,Jb,REAL(beta),REAL(betapost));

	// Partitioning covu
	
	for (j=0;j<JY*JZ;j++) {
		for (k=0;k<JY*JZ;k++) covu1[j+JY*JZ*k]=REAL(covu)[j+Ju*k];
		for (k=0;k<JY2;k++) covu12[j+JY*JZ*k]=REAL(covu)[j+Ju*(k+JY*JZ)];
	}
	for (j=0;j<JY2;j++) for (k=0;k<JY2;k++) covu2[j+JY2*k]=REAL(covu)[(j+JY*JZ)+Ju*(k+JY*JZ)];

	// Updating level 2 beta
	
	
	r8mat_pofac(JY2,covu2, help,3);
	r8mat_poinv(JY2,help, invomega);
	for (jj=1;jj<JY2;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY2*tt]=invomega[tt+JY2*jj];

	for (j=0;j<JY2*JY2*JX2*JX2;j++) sumxi[j]=0;	
	for (j=0;j<JY2*JX2;j++) sumxy[j]=0;
	for (j=0;j<Iu;j++) {
		for (t=0;t<JY2*JX2*JY2;t++) xi[t]=0; 
		for (t=0;t<JY2;t++) {
			yi[t]=imp2[j+t*Iu];
			for (k=0;k<JX2;k++) {
				xi[t+(k+t*JX2)*JY2]=X2red[j+Iu*k];
			}
		}
		r8mat_mtm_new(JY2*JX2,JY2,JY2,xi,invomega,help2);
		r8mat_mm_new(JY2*JX2,JY2,JX2*JY2,help2,xi,incrxx);
		r8mat_mm_new(JY2*JX2,JY2,1,help2,yi,incrxy);
		r8mat_add(JY2*JX2,JY2*JX2,incrxx,sumxi);
		r8mat_add(JY2*JX2,1,incrxy,sumxy);	
	}
	r8mat_pofac(JY2 * JX2,sumxi,help3,4);
	r8mat_poinv(JY2 * JX2, help3,invomega2);
	for (jj=1;jj<JX2*JY2;jj++) for (tt=0;tt<jj;tt++) invomega2[jj+JX2*JY2*tt]=invomega2[tt+JX2*JY2*jj];
	r8mat_mm_new(JY2*JX2,JY2*JX2,1,invomega2,sumxy,mu);
	r8mat_pofac(JY2 * JX2,invomega2,help3,5);
	r8vec_multinormal_sample(JY2*JX2, mu,help3, REAL(beta2),newbeta,0);
	r8mat_add(Ib2,Jb2,REAL(beta2),REAL(beta2post));
	
	// Calculating level 2 residuals 
	
	for (j=0;j<Iu;j++) {
		for (t=0;t<JY2*JX2*JY2;t++) xi[t]=0; 
		for (t=0;t<JY2;t++) {
			yi[t]=imp2[j+t*Iu];
			for (k=0;k<JX2;k++) {
				xi[t+(k+t*JX2)*JY2]=X2red[j+Iu*k];
			}
		}
		r8mat_mm_new(JY2,JY2*JX2,1,xi,REAL(beta2),help2);
		r8mat_divide(JY2,1,-1,help2);
		r8mat_add(JY2,1,help2,yi);	
		for (t=0;t<JY2;t++) REAL(u)[j+Iu*(JY*JZ+t)]=yi[t];
	}
	
	//Updating random effects
	
	r8mat_pofac(JY2,covu2, help,1);
	r8mat_poinv(JY2,help, invomega);
	for (jj=1;jj<JY2;jj++) for (tt=0;tt<jj;tt++) invomega[jj+JY2*tt]=invomega[tt+JY2*jj];
	r8mat_mm_new(JY*JZ,JY2,JY2,covu12,invomega,help6);
	r8mat_mmt_new(JY*JZ,JY2,JY*JZ,help6,covu12,help);
	r8mat_divide(JY*JZ,JY*JZ,-1,help);
	r8mat_add(JY*JZ,JY*JZ,covu1,help);	
	r8mat_pofac(JY*JZ,help, help5,1);
	r8mat_poinv(JY*JZ,help5, invomega2);
	for (jj=1;jj<JY*JZ;jj++) for (tt=0;tt<jj;tt++) invomega2[jj+JY*JZ*tt]=invomega2[tt+JY*JZ*jj];
	
	for (c=0;c<nj;c++) {
		for (j=0;j<JY2;j++) yi[j]=REAL(u)[c+Iu*(JY*JZ+j)];
		r8mat_mm_new(JY*JZ,JY2,1,help6,yi,mu);
		r8mat_pofac(JY,REAL(omega),help,6);
		r8mat_poinv(JY,help,invomega);
		for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
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
		
		r8mat_add(JY*JZ,JY*JZ,invomega2,sumzi);
		
		r8mat_pofac(JY * JZ,sumzi,help5,8);
		r8mat_poinv(JY * JZ, help5,invomega3);
		for (jj=1;jj<JZ*JY;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+JZ*JY*tt]=invomega3[tt+JZ*JY*jj];
		r8mat_mm_new(JY*JZ,JY*JZ,1,invomega2,mu,help5);
		r8mat_add(JY*JZ,1,help5,sumzy);
		r8mat_mm_new(JY*JZ,JY*JZ,1,invomega3,sumzy,mu2);
		r8mat_pofac(JY * JZ,invomega3,help5,9); 
		r8vec_multinormal_sample(JY*JZ, mu2,help5,newu, incrzy,0);
		for (t=0;t<JY;t++) for (k=0;k<JZ;k++) REAL(u)[c+nj*(k+t*JZ)] = newu[k+t*JZ];
		
	}
	r8mat_add(Iu,Ju,REAL(u),REAL(upost));
	
	//Updating level 2 covariance matrix

	if (ncat2==0) {
		for (j=0;j<Ju*Ju;j++) mu3[j]=0;
		for (j=0;j<nj; j++) {
			for (t=0;t<Ju;t++) uj[t]=REAL(u)[j+nj*t];
			r8mat_mmt_new(Ju,1,Ju,uj,uj,help4);
			r8mat_add(Ju,Ju,help4,mu3);	
		}
		r8mat_add(Ju,Ju,REAL(Sup),mu3);
		r8mat_pofac(Ju,mu3,help4,10);
		r8mat_poinv(Ju, help4, invomega3);
		for (jj=1;jj<(Ju);jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(Ju)*tt]=invomega3[tt+(Ju)*jj];

		wishart_sample(Ju,(nj+Ju),invomega3,newomega2, help4,sumzi,incrzz,mu3,0);
	
		r8mat_pofac(Ju,newomega2, help4,11);
		r8mat_poinv(Ju, help4,invomega3);
		for (jj=1;jj<(Ju);jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(Ju)*tt]=invomega3[tt+(Ju)*jj];
		for(k=0;k<(Ju);k++)  for(j=0;j<(Ju);j++)  REAL(covu)[j+(Ju)*k]=invomega3[j+(Ju)*k];
	
		r8mat_add(Ju,Ju,REAL(covu),REAL(covupost));
	}
	else {
			flag=0;
			r8mat_pofac(Ju,REAL(covu),help,12);
			detom=r8mat_podet(Ju, help);
			r8mat_poinv(Ju, help,invomega3);
			for (jj=1;jj<Ju;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+Ju*tt]=invomega3[tt+Ju*jj];
			logLH=0;
			for (t=0;t<Iu;t++) {
				for (jj=0;jj<Ju;jj++) help4[jj]=REAL(u)[t+jj*Iu];
				r8mat_mm_new(1,Ju,Ju,help4,invomega3,help5);
				r8mat_mmt_new(1,Ju,1,help5,help4,help6);
				logLH=logLH+help6[0];
			}
			logLH=logLH*(-0.5)-Iu*log(detom)/2;
			for (j=0;j<Ju;j++) {
				for (k=0;k<Ju;k++) {
					if ((fixomega2[j+Ju*k]==0)&(j<=k)) {
						if (j==k) {
						sdom=REAL(covu)[j+Ju*k]*sqrt(11.6/Iu);
						meanom=REAL(covu)[j+Ju*k];
						}
						else {
							sdom=0.1*sqrt(REAL(covu)[j+Ju*j]*REAL(covu)[k+Ju*k]);
							meanom=REAL(covu)[j+Ju*k];
						}
						kk=0;
						while ((flag==0)&(kk<100)) {
							newomega2[j+Ju*k]=r8_normal_sample(meanom,sdom,0);
							newomega2[k+Ju*j]=newomega2[j+Ju*k];
							flag=checkposdef(Ju,newomega2,help,help5);
							kk++;
						}
						r8mat_pofac(Ju,newomega2,help,13);
						detom=r8mat_podet(Ju, help);
						r8mat_poinv(Ju, help,invomega3);
						for (jj=1;jj<Ju;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+Ju*tt]=invomega3[tt+Ju*jj];
						newlogLH=0;
						for (t=0;t<Iu;t++) {
							for (jj=0;jj<Ju;jj++) help4[jj]=REAL(u)[t+jj*Iu];
							r8mat_mm_new(1,Ju,Ju,help4,invomega3,help5);
							r8mat_mmt_new(1,Ju,1,help5,help4,help6);
							newlogLH=newlogLH+help6[0];
						}
						newlogLH=newlogLH*(-0.5)-Iu*log(detom)/2;

						if (((( double ) unif_rand ( ) )<exp(newlogLH-logLH))&(flag==1)) {
							REAL(covu)[j+Ju*k]=newomega2[j+Ju*k];
							REAL(covu)[k+Ju*j]=newomega2[k+Ju*j];
							logLH=newlogLH;
						}
						else {
							newomega2[j+Ju*k]=REAL(covu)[j+Ju*k];
							newomega2[k+Ju*j]=REAL(covu)[k+Ju*j];
						}
						flag=0;
					}
				}
			}
						
		//	r8mat_print(Ju,Ju,REAL(covu),"covu");
			r8mat_add(Ju,Ju,REAL(covu),REAL(covupost));
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
	
	if (ncat==0) {
		for (j=0;j<JY*JY;j++) mu4[j]=0;
		for (j=0;j<IY;j++) {
			for (t=0;t<JY;t++) {
				yi[t]=resid[j+t*IY];
			}
			r8mat_mmt_new(JY,1,JY,yi,yi,help);
			r8mat_add(JY,JY,help,mu4);		
		}
		r8mat_add(JY,JY,REAL(Sp),mu4);

		r8mat_pofac(JY,mu4,help,9);
		r8mat_poinv(JY, help,invomega);
		for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
		wishart_sample(JY,IY-1,invomega,newomega,help, omegaoo,omegamo,omegamm,0);	
		r8mat_pofac(JY,newomega,help,10);
		r8mat_poinv(JY, help,invomega);
		for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
		for(k=0;k<Io;k++)  for(j=0;j<Io;j++)  REAL(omega)[k+Io*j]=invomega[k+Io*j];
		r8mat_add(JY,JY,REAL(omega),REAL(omegapost));
	}
	else {
		flag=0;
		r8mat_pofac(JY,newomega,help,12);
		detom=r8mat_podet(JY, help);
		r8mat_poinv(JY, help,invomega);
		for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
		logLH=0;
		for (t=0;t<IY;t++) {
			for (jj=0;jj<JY;jj++) help4[jj]=resid[t+jj*IY];
			r8mat_mm_new(1,JY,JY,help4,invomega,help5);
			r8mat_mmt_new(1,JY,1,help5,help4,help6);
			logLH=logLH+help6[0];
		}
		logLH=logLH*(-0.5)-IY*log(detom)/2;
		for (j=0;j<JY;j++) {
			for (k=0;k<JY;k++) {
				if ((fixomega[j+JY*k]==0)&(j<=k)) {
					if (j==k) {
						sdom=REAL(omega)[j+JY*k]*sqrt(11.6/IY);
						meanom=REAL(omega)[j+JY*k];
					}
					else {
						sdom=0.1*sqrt(REAL(omega)[j+JY*j]*REAL(omega)[k+JY*k]);
						meanom=REAL(omega)[j+JY*k];
					}
					kk=0;
					while ((flag==0)&(kk<100)) {
						newomega[j+JY*k]=r8_normal_sample(meanom,sdom,0);
						newomega[k+JY*j]=newomega[j+JY*k];
						flag=checkposdef(JY,newomega,help,help5);
						kk++;
					}
					r8mat_pofac(JY,newomega,help,13);
					detom=r8mat_podet(JY, help);
					r8mat_poinv(JY, help,invomega);
					for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
					newlogLH=0;
					for (t=0;t<IY;t++) {
						for (jj=0;jj<JY;jj++) help4[jj]=resid[t+jj*IY];
						r8mat_mm_new(1,JY,JY,help4,invomega,help5);
						r8mat_mmt_new(1,JY,1,help5,help4,help6);
						newlogLH=newlogLH+help6[0];
					}
					newlogLH=newlogLH*(-0.5)-IY*log(detom)/2;
					if (((( double ) unif_rand ( ) )<exp(newlogLH-logLH))&(flag==1)) {
						REAL(omega)[j+JY*k]=newomega[j+JY*k];
						REAL(omega)[k+JY*j]=newomega[k+JY*j];
						logLH=newlogLH;
					}
					else {
						newomega[j+JY*k]=REAL(omega)[j+JY*k];
						newomega[k+JY*j]=REAL(omega)[k+JY*j];
					}
					flag=0;
				}
			}
		}
		r8mat_add(JY,JY,REAL(omega),REAL(omegapost));

	}
		

	//imputing missing values

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
						omegamm[countmm]=REAL(omega)[k+t*Io];
						countmm++;
					}
					else if (ISNAN(REAL(Yimp)[j+t*IY])) {
						omegamo[countmo]=REAL(omega)[k+t*Io];
						countmo++;	
					}
					else if (!ISNAN(REAL(Yimp)[j+k*IY])&!ISNAN(REAL(Yimp)[j+t*IY])){
						omegaoo[countoo]=REAL(omega)[k+t*Io];
						countoo++;	
					}
				}
			}
			r8mat_pofac((JY-nmiss),omegaoo,help7,14);
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
			r8mat_pofac(nmiss,omegadrawmiss,help9,15);
			r8vec_multinormal_sample(nmiss,mumiss,help9,Ymiss,help6,0);
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
	
	//Imputing level 2 variables
	
for (j=0; j<Iu; j++) {
		for (k=0;k<JY2;k++) betaX[k]=0;		
		for (t=0;t<JY2*JX2*JY2;t++) xi[t]=0; 
		for (t=0;t<JY2;t++) {
			yi[t]=imp2[j+t*Iu];
			for (k=0;k<JX2;k++) {
				xi[t+(k+t*JX2)*JY2]=X2red[j+Iu*k];
			}
		}
		r8mat_mm_new(JY2,JY2*JX2,1,xi,REAL(beta2),help6);
		r8mat_add(JY2,1,help6,betaX);
		nmiss=missing2[j];
		if (nmiss>0) {
			for (k=0;k<Ju;k++) {
				if (k<JY*JZ) Yobs[k]=REAL(u)[j+k*Iu];
				else if ((ISNAN(Y2impred[j+(k-JY*JZ)*Iu]))) {
					betamiss[countm]=betaX[(k-JY*JZ)];
					countm++;
				}
				else {
					Yobs[JY*JZ+counto]=REAL(u)[j+k*Iu];
					betaobs[counto]=betaX[k];
					counto++;
				}
				for (t=0;t<Ju;t++) {
					if (t<(JY*JZ)) {
						if ((k<JY*JZ)||((k>(JY*JZ-1))&&(!ISNAN(Y2impred[j+(k-JY*JZ)*Iu])))) {
							omegaoo[countoo]=REAL(covu)[k+t*Ju];
							countoo++;
						}
					} 
					else {
						if (((k<JY*JZ)&&(!ISNAN(Y2impred[j+(t-JY*JZ)*Iu])))||((k>(JY*JZ-1))&&(!ISNAN(Y2impred[j+(k-JY*JZ)*Iu]))&&(!ISNAN(Y2impred[j+(t-JY*JZ)*Iu])))) {
							omegaoo[countoo]=REAL(covu)[k+t*Ju];
							countoo++;
						}
						if (((k<JY*JZ)&&(ISNAN(Y2impred[j+(t-JY*JZ)*Iu])))||((k>(JY*JZ-1))&&(!ISNAN(Y2impred[j+(k-JY*JZ)*Iu]))&&(ISNAN(Y2impred[j+(t-JY*JZ)*Iu])))) {
							omegamo[countmo]=REAL(covu)[k+t*Ju];
							countmo++;
						}
						if ((k>(JY*JZ-1))&&(ISNAN(Y2impred[j+(k-JY*JZ)*Iu]))&&(ISNAN(Y2impred[j+(t-JY*JZ)*Iu]))) {
							omegamm[countmm]=REAL(covu)[k+t*Ju];
							countmm++;
						}
					}
					
				}
			}
			r8mat_pofac((Ju-nmiss),omegaoo,help7,14);
			r8mat_poinv((Ju-nmiss),help7,invomega4);
			for (jj=1;jj<Ju-nmiss;jj++) for (tt=0;tt<jj;tt++) invomega4[jj+(Ju-nmiss)*tt]=invomega4[tt+(Ju-nmiss)*jj];
			r8mat_mm_new(nmiss,(Ju-nmiss),(Ju-nmiss),omegamo,invomega4,help8);
			r8mat_mm_new(nmiss,(Ju-nmiss),1,help8,Yobs,mumiss);
			r8mat_add(1,nmiss,betamiss,mumiss);
			r8mat_mmt_new(nmiss,(Ju-nmiss),nmiss,help8,omegamo,omegadrawmiss);
			r8mat_divide(nmiss,nmiss,-1,omegadrawmiss);
			r8mat_add(nmiss,nmiss,omegamm,omegadrawmiss);
			r8mat_pofac(nmiss,omegadrawmiss,help9,15);
			r8vec_multinormal_sample(nmiss,mumiss,help9,Ymiss,help6,0);
			countm=0;
			for (k=0;k<JY2;k++) {
				if (ISNAN(Y2impred[j+k*Iu])) {
					imp2[j+k*Iu]=Ymiss[countm];
					REAL(u)[j+(JY*JZ+k)*Iu]=Ymiss[countm]-betamiss[countm];
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
	
	if ((i+1)%fl==0) Rprintf("Iteration %d completed\n",i+1);
}


for(i=0;i<IY;i++)  {
	for(j=0;j<JY;j++)  {
		REAL(Yimp2)[i+IY*j]=imp[i+IY*j];
	}
	for(j=0;j<JY2;j++)  {
		REAL(Y2imp2)[i+IY2*j]=imp2[INTEGER(clus)[i]+Iu*j];
	}
	if (ncat>0) {
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
}
if (ncat2>0) {
	for (i=0;i<Iu;i++) {
	h=0;
	for (j=0;j<ncat2;j++) {
		maxx=imp2[i+(ncon2+h)*Iu];
		nmaxx=0;
		for (k=1;k<(INTEGER(Y2_numcat)[j]-1);k++) {
			if (imp2[i+(ncon2+h+k)*Iu]>maxx) {
				maxx=imp2[i+(ncon2+h+k)*Iu];
				nmaxx=k;
			}
		}
		if (maxx>0) {
			for (jj=0;jj<IY;jj++) if (INTEGER(clus)[jj]==i) REAL(Y2impcat)[jj+IY*j]=nmaxx+1;
		}
		else {
			for (jj=0;jj<IY;jj++) if (INTEGER(clus)[jj]==i) REAL(Y2impcat)[jj+IY*j]=INTEGER(Y2_numcat)[j];
		}
		h=h+INTEGER(Y2_numcat)[j]-1;
		}
	}

}
	
r8mat_divide(Ib,Jb,ns,REAL(betapost));
r8mat_divide(Ib2,Jb2,ns,REAL(beta2post));
r8mat_divide(Iu,Ju,ns,REAL(upost));
r8mat_divide(JY,JY,ns,REAL(omegapost));
r8mat_divide(Ju,Ju,ns,REAL(covupost));
PutRNGstate();
UNPROTECT(39);
return R_NilValue;
}
