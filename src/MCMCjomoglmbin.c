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

SEXP MCMCjomoglmbin(SEXP Ysub, SEXP Ysubimp, SEXP Ysubcat, SEXP submod, SEXP ordersub, SEXP Y, SEXP Yimp, SEXP Yimp2, SEXP Yimpcat, SEXP X,SEXP betaY, SEXP betaYpost, SEXP beta, SEXP betapost, SEXP varY, SEXP varYpost, SEXP omega, SEXP omegapost, SEXP nstep, SEXP varYprior, SEXP Sp, SEXP Y_numcat, SEXP num_con, SEXP flagrng){
int indic=0,i,j,k, IY,JY, IX, JX, Io, Jo, Ib, Jb, ns, nmiss=0,t, countm=0, counto=0,countoo=0, jj, tt, kk, ncon,ncat, pos,flag=0,nmaxx,h=0;
int fl, currncat, Is,  Il=0, JXm, accratio=0, totprop=0, nconnoaux, nconcat, ncatnoaux;
SEXP RdimY, RdimX, Rdimo, Rdimb, Rdims;
double *betaX, *Yobs, *Ymiss, *mumiss, *omegadrawmiss, *betamiss, *betaobs, *omegaoo, *omegamo, *omegamm, *invomega, *invomega2, *help, *help2, *help3, *imp, *yicategorized, *impsub;
double *mu, *mu2, *newbeta, *newomega, *yi, *invomega3, *help4, *help5, *help6, *missing, *fixomega,meanom,sdom, *resid, logLH, newlogLH,detom, *residsub;
double maxx,maxim,maxim2, *sumxy, *sumxi, *xi,  *incrxx, *incrxy, *Xsub, *Xsubprop;

/* Protecting R objects from garbage collection and saving matrices dimensions*/ 

RdimY=PROTECT(getAttrib(Yimp,R_DimSymbol));
IY=INTEGER(RdimY)[0];
JY=INTEGER(RdimY)[1];
Rdims=PROTECT(getAttrib(submod,R_DimSymbol));
Is=INTEGER(Rdims)[0];
RdimX=PROTECT(getAttrib(X,R_DimSymbol));
IX=INTEGER(RdimX)[0];
JX=INTEGER(RdimX)[1];
Rdimb=PROTECT(getAttrib(beta,R_DimSymbol));
Ib=INTEGER(Rdimb)[0];
Jb=INTEGER(Rdimb)[1];
Rdimo=PROTECT(getAttrib(omega,R_DimSymbol));
Io=INTEGER(Rdimo)[0];
Jo=INTEGER(Rdimo)[1];
Ysub=PROTECT(coerceVector(Ysub,REALSXP));
submod=PROTECT(coerceVector(submod,INTSXP));
ordersub=PROTECT(coerceVector(ordersub,INTSXP));
Ysubimp=PROTECT(coerceVector(Ysubimp,REALSXP));
Ysubcat=PROTECT(coerceVector(Ysubcat,REALSXP));
Y=PROTECT(coerceVector(Y,REALSXP));
Yimpcat=PROTECT(coerceVector(Yimpcat,REALSXP));
Y_numcat=PROTECT(coerceVector(Y_numcat,INTSXP));
Yimp=PROTECT(coerceVector(Yimp,REALSXP));
Yimp2=PROTECT(coerceVector(Yimp2,REALSXP));
X=PROTECT(coerceVector(X,REALSXP));
beta=PROTECT(coerceVector(beta,REALSXP));
betaY=PROTECT(coerceVector(betaY,REALSXP));
betapost=PROTECT(coerceVector(betapost,REALSXP));
betaYpost=PROTECT(coerceVector(betaYpost,REALSXP));
varY=PROTECT(coerceVector(varY,REALSXP));
varYpost=PROTECT(coerceVector(varYpost,REALSXP));
omega=PROTECT(coerceVector(omega,REALSXP));
omegapost=PROTECT(coerceVector(omegapost,REALSXP));
varYprior=PROTECT(coerceVector(varYprior,REALSXP));
Sp=PROTECT(coerceVector(Sp,REALSXP));
nstep=PROTECT(coerceVector(nstep,INTSXP));
ns=INTEGER(nstep)[0];
num_con=PROTECT(coerceVector(num_con,INTSXP));
ncon=INTEGER(num_con)[0];
nconnoaux=INTEGER(num_con)[1];
nconcat=INTEGER(num_con)[2];
ncatnoaux=INTEGER(num_con)[3];
flagrng=PROTECT(coerceVector(flagrng,INTSXP));
fl=INTEGER(flagrng)[0];
if (REAL(Yimpcat)[0]==(-999)) ncat=0;
else ncat=length(Y_numcat);




for (i=0;i<XLENGTH(ordersub);i++) {
	for (j=0;j<INTEGER(ordersub)[i];j++) {
		maxx=maxx*(INTEGER(submod)[3+4*h]);
		h++;
	}
	Il=Il+maxx;
	maxx=1;
}
Il=Il+1;
JXm=JY;
if (Il>JY) JXm=Il;
if (JX>JXm) JXm=JX;

/*Allocating memory for C objects in R*/

help = ( double * ) R_alloc ( (JY*JY) , sizeof ( double ) );
invomega= (double * ) R_alloc ( (JY*JY) , sizeof ( double ) );
fixomega = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
sumxi = ( double * ) R_alloc ( JY * JXm * JY * JXm , sizeof ( double ) );
sumxy = ( double * ) R_alloc ( JY * JXm , sizeof ( double ) );
xi = ( double * ) R_alloc ( JY * JXm * JY , sizeof ( double ) );
yi = ( double * ) R_alloc ( JY , sizeof ( double ) );
help2 = ( double * ) R_alloc ( JY *JX * JY , sizeof ( double ) );
incrxx = ( double * ) R_alloc ( JY *JXm * JY *JXm , sizeof ( double ) );
incrxy = ( double * ) R_alloc ( JY *JXm , sizeof ( double ) );
help3 = ( double * ) R_alloc ( JY * JXm * JY * JXm ,sizeof ( double ) );
invomega2= (double * ) R_alloc ( JY * JXm * JY * JXm, sizeof ( double ) );
mu = ( double * ) R_alloc ( JY * JXm, sizeof ( double ) );
newbeta = ( double * ) R_alloc ( JY*JXm  ,sizeof ( double ) );
mu2 = ( double * ) R_alloc ( JY * JY, sizeof ( double ) );
newomega = ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
betaX=( double * ) R_alloc ( JY, sizeof ( double ) );
imp=( double * ) R_alloc ( IY * JY,sizeof ( double ) );
impsub=( double * ) R_alloc ( IY ,sizeof ( double ) );
resid=( double * ) R_alloc ( IY * JY,sizeof ( double ) );
residsub=( double * ) R_alloc ( IY,sizeof ( double ) );
Yobs=( double * ) R_alloc ( JY, sizeof ( double ) );
Ymiss=( double * ) R_alloc ( JY, sizeof ( double ) );
mumiss = ( double * ) R_alloc ( JY, sizeof ( double ) );
omegadrawmiss = ( double * ) R_alloc ( JY*JY ,sizeof ( double ) );
betamiss = ( double * ) R_alloc ( JY ,sizeof ( double ) );
betaobs = ( double * ) R_alloc ( JY, sizeof ( double ) );
omegaoo= ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
omegamo= ( double * ) R_alloc (  JY*JY  , sizeof ( double ) );
omegamm= ( double * ) R_alloc (  JY*JY  , sizeof ( double ) );
invomega3= ( double * ) R_alloc ( JY * JY , sizeof ( double ) );
help4 = ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
help5 = ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
help6 = ( double * ) R_alloc ( JY*JY , sizeof ( double ) );
missing = ( double * ) R_alloc ( IY , sizeof ( double ) );
Xsub = ( double * ) R_alloc ( IY* Il , sizeof ( double ) );
Xsubprop = ( double * ) R_alloc ( Il , sizeof ( double ) );
yicategorized=( double * ) R_alloc (  JY,sizeof ( double ) );

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
r8mat_copy_new(IY, 1, REAL(Ysubimp), impsub);

for (i=0;i<JY*JY;i++) fixomega[i]=0;
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

for (i=0;i<JY*JY;i++) newomega[i]=REAL(omega)[i];

for (i=0;i<(IY*Il);i++) Xsub[i]=1;
h=0;
for (i=0;i<IY;i++) {
	for (j=0;j<XLENGTH(ordersub);j++) {
		pos=1;
		currncat=1;
		for (k=0;k<INTEGER(ordersub)[j];k++) {
			pos=pos*(INTEGER(submod)[3+(h+k)*4]);
		}
		for (k=0;k<INTEGER(ordersub)[j];k++) {
			if (INTEGER(submod)[1+h*4]==1) {
				for (jj=0;jj<pos;jj++) {
					Xsub[i+IY*(1+jj+indic)]=Xsub[i+IY*(1+jj+indic)]*pow(imp[i+IY*(INTEGER(submod)[h*Is]-1)],INTEGER(submod)[2+h*Is]);
				}	
			}
			else {
				for (jj=0;jj<pos;jj++) {
					kk=(jj*currncat)%INTEGER(submod)[3+h*4]+2;
					Xsub[i+IY*(1+jj+indic)]=Xsub[i+IY*(1+jj+indic)]*(REAL(Yimpcat)[i+IY*(INTEGER(submod)[h*Is]-1)]==kk);
				}
			}
			currncat=currncat*INTEGER(submod)[3+h*4];
			h=h+1;
		}
		currncat=1;
		indic=indic+pos;
	}
	h=0;
	indic=0;
}

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
						
						for (k=0;k<(JY-currncat);k++) help[k]=0;
						for (tt=0;tt<JX;tt++) {
							for (k=0;k<pos;k++) {
								help[k]=help[k]+REAL(beta)[tt+k*Ib]*REAL(X)[t+tt*IX];
							}
							for (k=(pos+currncat);k<JY;k++) {
								help[k-currncat]=help[k-currncat]+REAL(beta)[tt+k*Ib]*REAL(X)[t+tt*IX];
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
								r8vec_multinormal_sample((INTEGER(Y_numcat)[j]-1), mumiss,omegamm, newbeta,mu2,0);
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
								r8vec_multinormal_sample((INTEGER(Y_numcat)[j]-1), mumiss,omegamm, newbeta,mu2,0);
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
		

	//Updating beta

	r8mat_pofac(JY,REAL(omega), help,3);
	r8mat_poinv(JY,help, invomega);

	for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];

	for (j=0;j<JY*JY*JX*JX;j++) sumxi[j]=0;	
	for (j=0;j<JY*JX;j++) sumxy[j]=0;
	for (j=0;j<IY;j++) {
		for (t=0;t<JY*JX*JY;t++) xi[t]=0; 
		for (t=0;t<JY;t++) {
			yi[t]=imp[j+t*IY];
			for (k=0;k<JX;k++) {
				xi[t+(k+t*JX)*JY]=REAL(X)[j+IY*k];
			}
		}
		r8mat_mtm_new(JY*JX,JY,JY,xi,invomega,help2);
		r8mat_mm_new(JY*JX,JY,JX*JY,help2,xi,incrxx);
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
	for (j=0;j<Ib;j++) {
		for (t=0;t<Jb;t++) {
			REAL(betapost)[j+Ib*t+i*Ib*Jb]=REAL(beta)[j+Ib*t];
			}
		}	
	
	//Updating residuals

	for (t=0;t<IY;t++) {
		for (j=0; j<JY; j++) {
			resid[t+j*IY]=imp[t+j*IY];
			for (k=0;k<JX;k++) resid[t+j*IY]=resid[t+j*IY]-REAL(beta)[k+j*JX]*REAL(X)[t+k*IX];
		}
	}

	//Updating omega
	
	if (ncat==0) {
		for (j=0;j<JY*JY;j++) mu2[j]=0;
		for (j=0;j<IY;j++) {
			for (t=0;t<JY;t++) {
				yi[t]=resid[j+t*IY];
			}
			r8mat_mmt_new(JY,1,JY,yi,yi,help);
			r8mat_add(JY,JY,help,mu2);		
		}
		r8mat_add(JY,JY,REAL(Sp),mu2);

		r8mat_pofac(JY,mu2,help,9);
		r8mat_poinv(JY, help,invomega);
		for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
		wishart_sample(JY,IY-1,invomega,newomega,help, omegaoo,omegamo,omegamm,0);	
		r8mat_pofac(JY,newomega,help,10);
		r8mat_poinv(JY, help,invomega);
		for (jj=1;jj<Io;jj++) for (tt=0;tt<jj;tt++) invomega[jj+Io*tt]=invomega[tt+Io*jj];
		for(k=0;k<Io;k++)  for(j=0;j<Io;j++)  REAL(omega)[k+Io*j]=invomega[k+Io*j];
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) {
				REAL(omegapost)[j+JY*t+i*JY*JY]=REAL(omega)[j+JY*t];
				}
			}	}
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
		for (j=0;j<JY;j++) {
			for (t=0;t<JY;t++) {
				REAL(omegapost)[j+JY*t+i*JY*JY]=REAL(omega)[j+JY*t];
				}
			}
	}
		

	//imputing missing values compatibly with substantive model
	
	r8mat_pofac(JY,REAL(omega),help2,14);
	r8mat_poinv(JY,help2,invomega);
	for (jj=1;jj<JY;jj++) for (tt=0;tt<jj;tt++) invomega[jj+(JY)*tt]=invomega[tt+(JY)*jj];
	for (j=0; j<IY; j++) {
		nmiss=missing[j];
		k=0;
		if (nmiss>0) {
			for (t=0;t<JY;t++) betaX[t]=0;		
			for (t=0;t<JY*JX*JY;t++) xi[t]=0; 
			for (t=0;t<Il;t++) Xsubprop[t]=1;
			for (t=0;t<JY;t++) {
				yi[t]=imp[j+t*IY];
				for (tt=0;tt<JX;tt++) {
					xi[t+(tt+t*JX)*JY]=REAL(X)[j+IY*tt];
				}
			}

			r8mat_mm_new(JY,JY*JX,1,xi,REAL(beta),help6);
			r8mat_add(JY,1,help6,betaX);
			
			for (t=0;t<JY;t++) help4[t]=yi[t]-betaX[t];
			r8mat_mm_new(1,JY,JY,help4,invomega,help5);
			r8mat_mmt_new(1,JY,1,help5,help4,help6);
			mu2[0]=0;
			for (t=0;t<Il;t++) mu2[0]=mu2[0]+REAL(betaY)[t]*Xsub[j+IY*t];
			if (REAL(Ysubcat)[j]==1) mu2[0]=-mu2[0];
			logLH=log(normal_cdf(mu2[0]))-help6[0]/2;

		}
		while (nmiss>0) {
			if (ISNAN(REAL(Yimp)[j+k*IY])) {
				betamiss[0]=betaX[k];
				omegamm[0]=REAL(omega)[k+k*Io];
				for (t=0;t<JY;t++) {
					if (t!=k) {
						betaobs[counto]=betaX[t];
						omegamo[counto]=REAL(omega)[k+t*Io];
						Yobs[counto]=imp[j+t*IY];
						for (tt=0;tt<JY;tt++) {
							if (tt!=k) {
								omegaoo[countoo]=REAL(omega)[tt+t*Io];
								countoo++;
							}
						}
						counto++;
					}
				}

				r8mat_pofac((JY-1),omegaoo,help2,14);
				r8mat_poinv((JY-1),help2,invomega3);
				for (jj=1;jj<JY-1;jj++) for (tt=0;tt<jj;tt++) invomega3[jj+(JY-1)*tt]=invomega3[tt+(JY-1)*jj];
				r8mat_mmt_new((JY-1),(JY-1),1,invomega3,omegamo,help3);
				r8mat_divide((JY-1),1,-1,betaobs);
				r8mat_add((JY-1),1,betaobs,Yobs);
				r8mat_mtm_new(1,(JY-1),1,Yobs,help3,mumiss);

				mumiss[0]=mumiss[0]+betamiss[0];
				r8mat_mm_new(1,(JY-1),1,omegamo,help3,omegadrawmiss);
				omegadrawmiss[0]=omegamm[0]-omegadrawmiss[0];
				
				if ((k<nconnoaux)||((k>=ncon)&(k<nconcat))) {
					Ymiss[0]=r8_normal_sample(yi[k],sqrt(omegamm[0]/10),0);
					mu2[0]=0;
					for (t=0;t<Il;t++) mu2[0]=mu2[0]+REAL(betaY)[t]*Xsub[j+IY*t];
					if (REAL(Ysubcat)[j]==1) mu2[0]=-mu2[0];
					logLH=log(normal_cdf(mu2[0]))-help6[0]/2;

					// Controllare se accettabile
					help4[k]=help4[k]+Ymiss[0]-yi[k];
					r8mat_mm_new(1,JY,JY,help4,invomega,help5);
					r8mat_mmt_new(1,JY,1,help5,help4,help);
					yi[k]=Ymiss[0];								
					if ((ncatnoaux>0)) {
						h=0;
						for (jj=0;jj<(ncatnoaux);jj++) {
							maxx=yi[(ncon+h)];
							nmaxx=0;
							if (INTEGER(Y_numcat)[jj]>2) {
								for (kk=1;kk<(INTEGER(Y_numcat)[jj]-1);kk++) {
									if (yi[(ncon+h+kk)]>maxx) {
										maxx=yi[(ncon+h+kk)];
										nmaxx=kk;
									}
								}
							}			
							if (maxx>0) yicategorized[jj]=nmaxx;
							else yicategorized[jj]=INTEGER(Y_numcat)[jj]-1;
							h=h+INTEGER(Y_numcat)[jj]-1;
						}
					}
					
						//Update Xsubprop
					h=0;
					indic=0;
					for (t=0;t<Il;t++) Xsubprop[t]=1;

					for (jj=0;jj<XLENGTH(ordersub);jj++) {
						pos=1;
						currncat=1;
						for (t=0;t<INTEGER(ordersub)[jj];t++) {
							pos=pos*(INTEGER(submod)[3+(h+t)*4]);
						}
						for (t=0;t<INTEGER(ordersub)[jj];t++) {
							if (INTEGER(submod)[1+h*4]==1) {
								for (tt=0;tt<pos;tt++) {
									Xsubprop[(1+tt+indic)]=Xsubprop[(1+tt+indic)]*pow(yi[(INTEGER(submod)[h*Is]-1)],INTEGER(submod)[2+h*Is]);

								}	
							} else {
								for (tt=0;tt<pos;tt++) {
									kk=(tt*currncat)%INTEGER(submod)[3+h*4]+1;
									Xsubprop[(1+tt+indic)]=Xsubprop[(1+tt+indic)]*(yicategorized[(INTEGER(submod)[h*Is]-1)]==kk);
								}
							}	
							currncat=currncat*INTEGER(submod)[3+h*4];
							h=h+1;
						}
						currncat=1;
						indic=indic+pos;
					}		
 			
  					mu2[0]=0;
					for (t=0;t<Il;t++) mu2[0]=mu2[0]+REAL(betaY)[t]*Xsubprop[t];
					if (REAL(Ysubcat)[j]==1) mu2[0]=-mu2[0];	
					newlogLH=log(normal_cdf(mu2[0]))-help[0]/2;
					
					if ((( double ) unif_rand ( ) )<exp(newlogLH-logLH)) {	

						imp[j+k*IY]=Ymiss[0];
						help6[0]=help[0];
						for (t=0;t<Il;t++) Xsub[j+t*IY]=Xsubprop[t];
						logLH=newlogLH;
						accratio=accratio+1;
					} else {
						help4[k]=help4[k]-Ymiss[0]+imp[j+k*IY];
						yi[k]=imp[j+k*IY];
					}
					totprop=totprop+1;				
				} else {
					Ymiss[0]=r8_normal_sample(mumiss[0],sqrt(omegadrawmiss[0]),0);
					imp[j+k*IY]=Ymiss[0];
					help4[k]=help4[k]+Ymiss[0]-yi[k];
					r8mat_mm_new(1,JY,JY,help4,invomega,help5);
					r8mat_mmt_new(1,JY,1,help5,help4,help6);
				}
				nmiss--;
				counto=0;
				countoo=0;
			}				
			k++;			
		}			
		counto=0;
		countoo=0;			
	}

		// Rejection sampling for latent normal outcome
	
		for (t=0;t<IY;t++) {
			if (!ISNAN(REAL(Ysub)[t])) {
				kk=0;
				flag=0;
				mu2[0]=0;
				for (k=0;k<Il;k++) mu2[0]=mu2[0]+REAL(betaY)[k]*Xsub[t+k*IX];
				while (flag==0&&kk<10000) {
					yi[0]=r8_normal_sample(mu2[0],sqrt(REAL(varY)[0]),0);
					if ((REAL(Ysub)[t]==2&&yi[0]>0)||(REAL(Ysub)[t]==1&&yi[0]<0)) {
						impsub[t]=yi[0];
						flag=1;
					} else {
						kk++;
					}
				}
			}
		}
		
	// Update beta of substantive model
	
	for (j=0;j<Il*Il;j++) sumxi[j]=0;	
	for (j=0;j<Il;j++) sumxy[j]=0;
	for (j=0;j<IY;j++) {
		yi[0]=impsub[j];
		for (k=0;k<Il;k++) {
			xi[k]=Xsub[j+IY*k];
		}
		r8mat_mtm_new(Il,1,Il,xi,xi,incrxx);
		r8mat_divide(Il,Il,REAL(varY)[0],incrxx);
		r8mat_mtm_new(Il,1,1,xi,yi,incrxy);
		r8mat_divide(Il,1,REAL(varY)[0],incrxy);
		r8mat_add(Il,Il,incrxx,sumxi);
		r8mat_add(Il,1,incrxy,sumxy);	
	}
	r8mat_pofac(Il,sumxi,help3,4);
	r8mat_poinv(Il, help3,invomega2);
	for (jj=1;jj<Il;jj++) for (tt=0;tt<jj;tt++) invomega2[jj+Il*tt]=invomega2[tt+Il*jj];
	r8mat_mm_new(Il,Il,1,invomega2,sumxy,mu);
	r8mat_pofac(Il,invomega2,help3,5);
	r8vec_multinormal_sample(Il, mu,help3, REAL(betaY),newbeta,0);	
	for (j=0;j<Il;j++) REAL(betaYpost)[j+i*Il]=REAL(betaY)[j];

			//Updating residuals

	for (t=0;t<IY;t++) {
		residsub[t]=impsub[t];
		for (k=0;k<Il;k++) {
			residsub[t]=residsub[t]-REAL(betaY)[k]*Xsub[t+k*IX];
		}
	}
	
//Not updating omega (fixed to 1)

	REAL(varYpost)[i]=REAL(varY)[0];

		// Imputing missing outcomes
	for (t=0;t<IY;t++) {
		if (ISNAN(REAL(Ysub)[t])) {
			mu2[0]=0;
			for (k=0;k<Il;k++) mu2[0]=mu2[0]+REAL(betaY)[k]*Xsub[t+k*IX];
			impsub[t]=r8_normal_sample(mu2[0],sqrt(REAL(varY)[0]),0);
			REAL(Ysubcat)[t]=(impsub[t]>0)+1;
		}
	}

	if ((i+1)%fl==0) Rprintf(".");
}
	if (fl==1) Rprintf("\n");
if (((double)accratio/((double)totprop))<0.15) Rprintf("Warning: acceptance ratio = %f. This might be a sign that the chain did not mix well. \n" , ((double)accratio/((double)totprop)));
for(i=0;i<IY;i++)  {
	for(j=0;j<JY;j++)  {
		REAL(Yimp2)[i+IY*j]=imp[i+IY*j];
	}
	
	REAL(Ysubimp)[i]=impsub[i];

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

PutRNGstate();
UNPROTECT(29);
return R_NilValue;
}
