# include <stdlib.h>
# include <R.h>
# include <Rinternals.h>


# include <math.h>

# include <time.h>

# include <string.h>


# include "pdflib.h"
# include "rnglib.h"




int checkposdef(int dim, double matr[], double matrh[],double matrh2[]) 

{
int hg,hj,hi,flag=1;
if (matr[0]<=0) flag=0;
for (hg=2;hg<(dim+1);hg++) {
	for (hj=0;hj<hg;hj++) {
		for (hi=0;hi<hg;hi++) {
			matrh[hj+hg*hi]=matr[hj+hi*dim];
		}
	}
	r8mat_pofac(hg,matrh,matrh2,17);
	if(ISNAN(r8mat_podet(hg,matrh2))) flag=0;	
}

return flag;

}

double argmaxvec (int card, double vec[])

{
int gf, argmaxx;
double maxv;

maxv=vec[0];
argmaxx=0;

for (gf=1; gf<card; gf++) {
	if (vec[gf]>maxv) {
		maxv=vec[gf];
		argmaxx=gf;
	}
}
return argmaxx;
}


double maxvec (int card, double vec[])

{
int gf;
double maxv;

maxv=vec[0];
for (gf=1; gf<card; gf++) {
	if (vec[gf]>maxv) maxv=vec[gf];
}
return maxv;
}




/******************************************************************************/



double r8_chi_sample ( double df, int fl)



/******************************************************************************/


/*
  Purpose:

    R8_CHI_SAMPLE generates a Chi-Square random deviate.

  Discussion:

    This procedure generates a random deviate from the chi square distribution
    with DF degrees of freedom random variable.

    The algorithm exploits the relation between chisquare and gamma.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2013

  Author:

    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.

  Parameters:

    Input, double DF, the degrees of freedom.
    0.0 < DF.

    Output, double R8_CHI_SAMPLE, a random deviate from the distribution.
*/


{
  

double arg1;
  
double arg2;
  
double value;

  
if ( df <= 0.0 )
  {
    
	Rprintf ( "\n" );
    
	Rprintf ( "R8_CHI_SAMPLE - Fatal error!\n" );
    
	Rprintf ( "  DF <= 0.\n" );
    
	Rprintf ( "  Value of DF: %g\n", df );
    
 }

  	
arg1 = 1.0;
  
arg2 = df / 2.0;

  
value = 2.0 * r8_gamma_sample ( arg1, arg2, fl );

  
return value;

}


/******************************************************************************/




double r8_epsilon ( void )



/******************************************************************************/


/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/


{
  
const double value = 2.220446049250313E-016;

  
return value;

}


/******************************************************************************/




double r8_exponential_01_sample ( int fl )



/******************************************************************************/


/*
  Purpose:

    R8_EXPONENTIAL_01_SAMPLE samples the standard exponential PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EXPONENTIAL_01_SAMPLE, a sample of the PDF.
*/


{
  
double r;
  
double value;

  
r = r8_uniform_01_sample ( fl );

  
value = - log ( r );

  
return value;

}


/******************************************************************************/

double r8_gamma_sample ( double a, double r, int fl )



/******************************************************************************/


/*
  Purpose:

    R8_GAMMA_SAMPLE generates a Gamma random deviate.

  Discussion:

    This procedure generates random deviates from the gamma distribution whose
    density is (A^R)/Gamma(R) * X^(R-1) * Exp(-A*X)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2013

  Author:

    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.

  Reference:

    Joachim Ahrens, Ulrich Dieter,
    Generating Gamma Variates by a Modified Rejection Technique,
    Communications of the ACM,
    Volume 25, Number 1, January 1982, pages 47-54.

    Joachim Ahrens, Ulrich Dieter,
    Computer Methods for Sampling from Gamma, Beta, Poisson and
    Binomial Distributions,
    Computing,
    Volume 12, Number 3, September 1974, pages 223-246.

  Parameters:

    Input, double A, the location parameter.

    Input, double R, the shape parameter.

    Output, double R8_GAMMA_SAMPLE, a random deviate from the distribution.
*/


{
  
double value;

  
value = r8_gamma_01_sample ( r, fl ) / a;

  
return value;

}


/******************************************************************************/





double r8_gamma_01_sample ( double a, int fl  )



/******************************************************************************/


/*
  Purpose:

    R8_GAMMA_01_SAMPLE samples the standard Gamma distribution.

  Discussion:

    This procedure corresponds to algorithm GD in the reference.

    pdf ( a; x ) = 1/gamma(a) * x^(a-1) * exp ( - x )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2013

  Author:

    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.

  Reference:

    Joachim Ahrens, Ulrich Dieter,
    Generating Gamma Variates by a Modified Rejection Technique,
    Communications of the ACM,
    Volume 25, Number 1, January 1982, pages 47-54.

  Parameters:

    Input, double A, the shape parameter of the standard gamma
    distribution.  0.0 < A.

    Output, double R8_GAMMA_01_SAMPLE, a random deviate from the distribution.
*/


{
  
double a1 =  0.3333333;
  
double a2 = -0.2500030;
  
double a3 =  0.2000062;
  
double a4 = -0.1662921;
  
double a5 =  0.1423657;
  
double a6 = -0.1367177;
  
double a7 =  0.1233795;
  
double b;
  
double c;
  
double d;
  
double e;
  
double e1 = 1.0;
  
double e2 = 0.4999897;
  
double e3 = 0.1668290;
  
double e4 = 0.0407753;
  
double e5 = 0.0102930;
  
double p;
  
double q;
  
double q0;
  
double q1 =  0.04166669;
  
double q2 =  0.02083148;
  
double q3 =  0.00801191;
  
double q4 =  0.00144121;
  
double q5 = -0.00007388;
  
double q6 =  0.00024511;
  
double q7 =  0.00024240;
  
double r;
  
double s;
  
double s2;
  
double si;
  
double sqrt32 = 5.6568542494923801952;
  
double t;
  
double u;
  
double v;
  
double value=0;
  
double w;
  
double x;

  
if ( 1.0 <= a )
  {
    
	s2 = a - 0.5;
    
	s = sqrt ( s2 );
    
	d = sqrt32 - 12.0 * s;
	
/*
  Immediate acceptance.
*/
    
	t = r8_normal_01_sample ( fl );
    
	x = s + 0.5 * t;
    
	value = x * x;

    
	if ( 0.0 <= t )
    {
      
		return value;
    
		}

	/*
  Squeeze acceptance.
*/
    
	u = r8_uniform_01_sample ( fl );
    
	if ( d * u <= t * t * t )
    {
      
		return value;
    
		}

    
	r = 1.0 / a;
    
	q0 = (((((( 
q7 * r + q6 ) * r 
+ q5 ) * r 
+ q4 ) * r
 + q3 ) * r
 + q2 ) * r
 + q1 ) * r;


	/*
  Approximation depending on size of parameter A.
*/
    
	if ( 13.022 < a )
    {
      
		b = 1.77;
      
		si = 0.75;
      
		c = 0.1515 / s;
    
		}
    
	else if ( 3.686 < a )
    {
      
		b = 1.654 + 0.0076 * s2;
      
		si = 1.68 / s + 0.275;
      
		c = 0.062 / s + 0.024;
    
		}
    
	else
    {
      
		b = 0.463 + s + 0.178 * s2;
      
		si = 1.235;
      
		c = 0.195 / s - 0.079 + 0.16 * s;
    
		}


	/*
  Quotient test.
*/
    
	if ( 0.0 < x )
    {
      
		v = 0.5 * t / s;

      
		if ( 0.25 < fabs ( v ) )
      {
        
			q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
      
			}
      
		else
      {
        
			q = q0 + 0.5 * t * t * (((((( 
a7 * v 
+ a6 ) * v 
+ a5 ) * v 
+ a4 ) * v 
+ a3 ) * v 
+ a2 ) * v 
+ a1 ) * v;
      
			}

      
		if ( log ( 1.0 - u ) <= q )
      {
        
			return value;
      
			}
    
		}

    
	for ( ; ; )
    {
      
		e = r8_exponential_01_sample ( fl );
      
		u = 2.0 * r8_uniform_01_sample ( fl ) - 1.0;
 
      
		if ( 0.0 <= u )
      {
        
			t = b + fabs ( si * e );
      
			}
      
		else
      {
        
			t = b - fabs ( si * e );
      
			}


	/*
  Possible rejection.
*/
      

		if ( t < -0.7187449 )
      {
        
			continue;
      
			}


	/*
  Calculate V and quotient Q.
*/
      
		v = 0.5 * t / s;

      
		if ( 0.25 < fabs ( v ) )
      {
        
			q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
      
			}
      
		else
      {
        
			q = q0 + 0.5 * t * t * (((((( 
a7 * v 
+ a6 ) * v 
+ a5 ) * v 
+ a4 ) * v 
+ a3 ) * v 
+ a2 ) * v 
+ a1 ) * v;
      
			}


	/*
  Hat acceptance.
*/
      

		if ( q <= 0.0 )
      {
        
			continue;
      
			}

      
		if ( 0.5 < q )
      {
        
			w = exp ( q ) - 1.0;
      
			}	
      
		else
      {
        
			w = (((( 
e5 * q 
+ e4 ) * q 
+ e3 ) * q 
+ e2 ) * q 
+ e1 ) * q;
      
			}


	/*
  May have to sample again.
*/
      

		if ( c * fabs ( u ) <= w * exp ( e - 0.5 * t * t ) )
      
			{
        
			break;
      
			}
    
		}

    

	x = s + 0.5 * t;
    
	value = x * x;
  
	}


/*
  Method for A < 1.
*/
  

else if ( a < 1.0 )
  {
    
	b = 1.0 + 0.3678794 * a;

    	
	for ( ; ; )
    {
      
		p = b * r8_uniform_01_sample ( fl );

      
		if ( p < 1.0 )
      {
        
			value = exp ( log ( p ) / a );
        
			if ( value <= r8_exponential_01_sample ( fl ) )
        {
          
				break;
        
				}
      
			}
      
		else
      {
        
			value = - log ( ( b - p ) / a );
        
			if ( ( 1.0 - a ) * log ( value ) <= r8_exponential_01_sample ( fl ) )
        {
          
				break;
        
				}
      
			}
    
		}	
  
	}
  
return value;

}



/******************************************************************************/




double r8_max ( double x, double y )



/******************************************************************************/


/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/


{
  
double value;

  
if ( y < x )
  {
    
	value = x;
  
	}
  
else
  {
    
	value = y;
  
	}
  
return value;

}


/******************************************************************************/



double r8_min ( double x, double y )



/******************************************************************************/


/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/


{
  

double value;

  
if ( y < x )
  {
    
	value = y;
  
	}
  
else
  {
    
	value = x;
  
	}
  
return value;


}




/******************************************************************************/



double r8_normal_sample ( double av, double sd, int fl )



/******************************************************************************/


/*
  Purpose:

    R8_NORMAL_SAMPLE generates a normal random deviate.

  Discussion:

    This procedure generates a single random deviate from a normal distribution
    with mean AV, and standard deviation SD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2013

  Author:

    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.

  Reference:

    Joachim Ahrens, Ulrich Dieter,
    Extensions of Forsythe's Method for Random
    Sampling from the Normal Distribution,
    Mathematics of Computation,
    Volume 27, Number 124, October 1973, page 927-937.

  Parameters:

    Input, double AV, the mean.

    Input, double SD, the standard deviation.

    Output, double R8_NORMAL_SAMPLE, a random deviate from the distribution.
*/


{
  
double value;

  
value = sd * r8_normal_01_sample ( fl ) + av;

  
return value;

}


/******************************************************************************/





double r8_normal_01_sample ( int fl )



/******************************************************************************/


/*
  Purpose:

    R8_NORMAL_01_SAMPLE samples the standard normal probability distribution.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    The Box-Muller method is used, which is efficient, but
    generates two values at a time.

    Typically, we would use one value and save the other for the next call.
    However, the fact that this function has saved memory makes it difficult
    to correctly handle cases where we want to re-initialize the code,
    or to run in parallel.  Therefore, we will instead use the first value
    and DISCARD the second.

    EFFICIENCY must defer to SIMPLICITY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R8_NORMAL_01_SAMPLE, a normally distributed random value.
*/


{
  
const double pi = 3.14159265358979323;
  
double r1;
  
double r2;
  
double x;

  
r1 = r8_uniform_01_sample ( fl );
  
r2 = r8_uniform_01_sample ( fl );

  
x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );

  
return x;

}



/******************************************************************************/



double r8_uniform_sample ( double low, double high, int fl )



/******************************************************************************/


/*
  Purpose:

    R8_UNIFORM_SAMPLE generates a uniform random deviate.

  Discussion:

    This procedure generates a real deviate uniformly distributed between
    LOW and HIGH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2013

  Author:

    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.

  Parameters:

    Input, double LOW, HIGH, the lower and upper bounds.

    Output, double R8_UNIFORM_SAMPLE, a random deviate from the distribution.
*/


{
  
double value;

  
value = low + ( high - low ) * r8_uniform_01_sample ( fl );

  
return value;

}




/******************************************************************************/



double r8_uniform_01_sample ( int fl )



/******************************************************************************/


/*
  Purpose:

    R8_UNIFORM_01_SAMPLE generates a uniform random deviate from [0,1].

  Discussion:

    This function should be the only way that the package accesses random
    numbers.

    Setting OPTION to 0 accesses the R8_UNI_01() function in RNGLIB,
    for which there are versions in various languages, which should result
    in the same values being returned.  This should be the only place in
    this library that accesses a function in RNGLIB.

    Setting OPTION to 1 in the C version calls the system random number
    generator function "rand()".

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 August 2013

  Author:

    John Burkardt.

  Parameters:

    Output, double R8_UNIFORM_01_SAMPLE, a random deviate.
*/


{
  
  
  
double value;

  
if ( fl == 0 )
  {
    
	value = r8_uni_01 ( );
  
	}
  
else
  {
    
	value = ( double ) unif_rand ( ) / ( double ) RAND_MAX;
  
	}

  
return value;


}



/******************************************************************************/



double r8mat_podet ( int n, double r[] )



/******************************************************************************/


/*
  Purpose:

    R8MAT_PODET computes the determinant of a factored positive definite matrix.

  Discussion:

    This routine expects to receive R, the upper triangular factor of A,
    computed by R8MAT_POFAC, with the property that A = R' * R.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2013

  Author:

    C version by John Burkardt.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double R[N*N], the Cholesky factor of A.

    Output, double R8MAT_PODET, the determinant of A.
*/


{
  
double det;
  
int i;

  
det = 1.0;
  
for ( i = 0; i < n; i++ )
  {
    
	det = det * r[i+i*n] * r[i+i*n];
  
}
  
return det;

}


/******************************************************************************/


void r8mat_pofac ( int n, double a[], double r[] , int indica)



/******************************************************************************/


/*
  Purpose:

    R8MAT_POFAC factors a real symmetric positive definite matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 April 2013

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
    Pete Stewart,
    C version by John Burkardt.

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
    LINPACK User's Guide,
    SIAM, (Society for Industrial and Applied Mathematics),
    3600 University City Science Center,
    Philadelphia, PA, 19104-2688.
    ISBN 0-89871-172-X

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the matrix to be  factored.

    Output, double R8MAT_POFAC[N*N], an upper triangular matrix such that
    A = R'*R.
*/


{
  

double dot;
  
int i;
  
int j;
  
int k;
  
double s;
  
double t;

  
for ( j = 0; j < n; j++ )
  {
    
	for ( i = 0; i <= j; i++ )
    {
      
		r[i+j*n] = a[i+j*n];
    
		}
    
	for ( i = j + 1; i < n; i++ )
    {
      
		r[i+j*n] = 0.0;
    
		}
  
	}

  
for ( j = 0; j < n; j++ )
  {
    
	s = 0.0;

    
	for ( k = 0; k < j; k++ )
    {
      
		dot = 0.0;
      
		for ( i = 0; i < k; i++ )
      {
        
			dot = dot + r[i+k*n] * r[i+j*n];
      
			}
      t = r[k+j*n] - dot;
      
		t = t / r[k+k*n];
      
		r[k+j*n] = t;
      
		s = s + t * t;
    
		}

    
	s = r[j+j*n] - s;

   
	r[j+j*n] = sqrt ( s );
  
	}

  
return;


}


/******************************************************************************/



void r8mat_poinv ( int n, double r[], double b[] )



/******************************************************************************/


/*
  Purpose:

    R8MAT_POINV computes the inverse of a factored positive definite matrix.

  Discussion:

    This routine expects to receive R, the upper triangular factor of A,
    computed by R8MAT_POFAC, with the property that A = R' * R.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2013

  Author:

    C version by John Burkardt.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix A.

    Input, double R[N*N], the Cholesky factor of A.

    Input, double R8MAT_POINV[N*N], the inverse of A.
*/


{
  
int i;
  
int j;
  
int k;
  
double t;

  
for ( j = 0; j < n; j++ )
  {
    
	for ( i = 0; i < n; i++ )
    {
      
		b[i+j*n] = r[i+j*n];
    
		}
  
	}

  
for ( k = 0; k < n; k++ )
  {
    
	b[k+k*n] = 1.0 / b[k+k*n];
    
	for ( i = 0; i < k; i++ )
    {
      
		b[i+k*n] = - b[i+k*n] * b[k+k*n];
    
		}
    
	for ( j = k + 1; j < n; j++ )
    {
      
		t = b[k+j*n];
      
		b[k+j*n] = 0.0;
      
		for ( i = 0; i <= k; i++ )
      {
        
			b[i+j*n] = b[i+j*n] + t * b[i+k*n];
      
			}
    
		}
  
	}


/*
  Form inverse(R) * (inverse(R))'.
*/
  

for ( j = 0; j < n; j++ )
  {
    
	for ( k = 0; k < j; k++ )
    {
      
		t = b[k+j*n];
      
		for ( i = 0; i <= k; i++ )
      {
        
			b[i+k*n] = b[i+k*n] + t * b[i+j*n];

		        }
    
		}
    
	t = b[j+j*n];
    
	for ( i = 0; i <= j; i++ )
    {
      
		b[i+j*n] = b[i+j*n] * t;
    
		}
  
	}
  
return;


}


/******************************************************************************/



void r8vec_multinormal_sample ( int n, double mu[], double r[], double x[], double z[], int fl )



/******************************************************************************/


/*
  Purpose:

    R8VEC_MULTINORMAL_SAMPLE samples a multivariate normal PDF.

  Discussion:

    PDF ( MU(1:N), C(1:N,1:N); X(1:N) ) = 
      1 / ( 2 * pi ) ^ ( N / 2 ) * 1 / det ( C )
      * exp ( - ( X - MU )' * inverse ( C ) * ( X - MU ) / 2 )

    Here,

      X is the argument vector of length N,
      MU is the mean vector of length N,
      C is an N by N positive definite symmetric covariance matrix.

    The properties of C guarantee that it has an upper triangular
    matrix R, the Cholesky factor, such that C = R' * R.  It is the
    matrix R that is required by this routine.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension.

    Input, double MU[N], the mean vector.

    Input, double R[N*N], the upper triangular Cholesky
    factor of the covariance matrix C.

    Output, double R8VEC_MULTINORMAL_SAMPLE[N], a sample of the distribution.
*/


{
  

int i;
  
int j;


/*
  Compute X = MU + R' * Z
  where Z is a vector of standard normal variates.
*/
 

for ( j = 0; j < n; j++ )
  {
    
	z[j] = r8_normal_01_sample ( fl );
  
	}

  
for ( i = 0; i < n; i++ )
  {
    
	x[i] = mu[i];
    
	for ( j = 0; j <= i; j++ )
    {
      
		x[i] = x[i] + r[j+i*n] * z[j];
    
		}
  
	}

   
return;

}

/******************************************************************************/

double r8_gamma_log ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA_LOG evaluates the logarithm of the gamma function.

  Discussion:

    This routine calculates the LOG(GAMMA) function for a positive real
    argument X.  Computation is based on an algorithm outlined in
    references 1 and 2.  The program uses rational functions that
    theoretically approximate LOG(GAMMA) to at least 18 significant
    decimal digits.  The approximation for X > 12 is from reference
    3, while approximations for X < 12.0 are similar to those in
    reference 1, but are unpublished.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2013

  Author:

    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.

  Reference:

    William Cody, Kenneth Hillstrom,
    Chebyshev Approximations for the Natural Logarithm of the
    Gamma Function,
    Mathematics of Computation,
    Volume 21, Number 98, April 1967, pages 198-203.

    Kenneth Hillstrom,
    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
    May 1969.

    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_GAMMA_LOG, the value of the function.
*/
{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04,
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double corr;
  const double d1 = -5.772156649015328605195174E-01;
  const double d2 = 4.227843350984671393993777E-01;
  const double d4 = 1.791759469228055000094023;
  const double frtbig = 2.25E+76;
  int i5;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = { 
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+10, 
    1.702665737765398868392998E+11, 
    4.926125793377430887588120E+11, 
    5.606251856223951465078242E+11 };
  double q1[8] = { 
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = { 
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = { 
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+10, 
    1.016803586272438228077304E+11, 
    3.417476345507377132798597E+11, 
    4.463158187419713286462081E+11 };
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  const double xbig = 2.55E+305;
  double xden;
  const double xinf = 1.79E+308;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double y;
  double ysq;

  y = x;

  if ( 0.0 < y && y <= xbig )
  {
    if ( y <= r8_epsilon ( ) )
    {
      res = - log ( y );
    }
/*
  EPS < X <= 1.5.
*/
    else if ( y <= 1.5 )
    {
      if ( y < 0.6796875 )
      {
        corr = -log ( y );
        xm1 = y;
      }
      else
      {
        corr = 0.0;
        xm1 = ( y - 0.5 ) - 0.5;
      }

      if ( y <= 0.5 || 0.6796875 <= y )
      {
        xden = 1.0;
        xnum = 0.0;
        for ( i5 = 0; i5 < 8; i5++ )
        {
          xnum = xnum * xm1 + p1[i5];
          xden = xden * xm1 + q1[i5];
        }
        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
      }
      else
      {
        xm2 = ( y - 0.5 ) - 0.5;
        xden = 1.0;
        xnum = 0.0;
        for ( i5 = 0; i5 < 8; i5++ )
        {
          xnum = xnum * xm2 + p2[i5];
          xden = xden * xm2 + q2[i5];
        }
        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );
      }
    }
/*
  1.5 < X <= 4.0.
*/
    else if ( y <= 4.0 )
    {
      xm2 = y - 2.0;
      xden = 1.0;
      xnum = 0.0;
      for ( i5 = 0; i5 < 8; i5++ )
      {
        xnum = xnum * xm2 + p2[i5];
        xden = xden * xm2 + q2[i5];
      }
      res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
    }
/*
  4.0 < X <= 12.0.
*/
    else if ( y <= 12.0 )
    {
      xm4 = y - 4.0;
      xden = -1.0;
      xnum = 0.0;
      for ( i5 = 0; i5 < 8; i5++ )
      {
        xnum = xnum * xm4 + p4[i5];
        xden = xden * xm4 + q4[i5];
      }
      res = d4 + xm4 * ( xnum / xden );
    }
/*
  Evaluate for 12 <= argument.
*/
    else
    {
      res = 0.0;

      if ( y <= frtbig )
      {
        res = c[6];
        ysq = y * y;
        for ( i5 = 0; i5 < 6; i5++ )
        {
          res = res / ysq + c[i5];
        }
      }
      res = res / y;
      corr = log ( y );
      res = res + sqrtpi - 0.5 * corr;
      res = res + y * ( corr - 1.0 );
    }
  }
/*
  Return for bad arguments.
*/
  else
  {
    res = xinf;
  }
/*
  Final adjustments and return.
*/
  return res;
}


/******************************************************************************/

double log_mul_gamma ( int p, double a )



/******************************************************************************/


{
  

int j2;
  
double g=0;




for ( j2 = 0; j2 < p; j2++ )
  {
    
	g = g+(lgamma(a + (1-(double)j2)/2));
	}

  
return g;

}


/******************************************************************************/

double r8_chi_pdf ( double df, double rval)

/******************************************************************************/

/*
  Purpose:

    R8_CHI_PDF evaluates the PDF of a chi-squared distribution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2013

  Author:

    Original FORTRAN90 version by Guannan Zhang.
    C  by John Burkardt.

  Parameters:

    Input, double DF, the degrees of freedom.
    0.0 < DF.

    Input, double RVAL, the point where the PDF is evaluated.

    Output, double R8_CHI_PDF, the value of the PDF at RVAL.
*/
{
  double temp1;
  double temp2;
  double value;

  if ( df <= 0.0 )
  {
    Rprintf ( "\n" );
    Rprintf ( "R8_CHI_PDF - Fatal error!\n" );
    Rprintf ( "  Degrees of freedom must be positive.\n" );
  }
      
  if ( rval <= 0.0 )

  {
    value = 0.0;
  }
  else
  {
    temp2 = df * 0.5;

    temp1 = ( temp2 - 1.0 ) * log ( rval ) - 0.5 * rval 
      - temp2 * log ( 2.0 ) - r8_gamma_log ( temp2 );
    value = exp ( temp1 );
  }
  return value;
}

/******************************************************************************/


double wishart_dens ( double df, int dim, double X[], double invA[], double help[], double help2[])



/******************************************************************************/


{
  

double d1,d2,res;



r8mat_pofac(dim,X,help,18);
d1=
r8mat_podet ( dim, help )

;
r8mat_pofac(dim,invA,help,19);
d2=
r8mat_podet ( dim, help )

;
res=(-df*dim/2)*log(2)-(df/2)*log(d2)-log_mul_gamma(dim,df/2)+((df-dim-1)/2)*log(d1);
return res;

}




/******************************************************************************/


double log_f_u ( double eta, double a, int dim, int nclus, double allinvomega[], double omega[], double invA[], double help[], double help2[], double gamma, double Gammastar[])



/******************************************************************************/


{
  
double res;
int i7;
  
int j4;


int t1;

res=log(r8_chi_pdf(eta, a));
for (i7=0;i7<nclus;i7++) {
	for (j4=0;j4<dim;j4++) {
		for (t1=0;t1<dim;t1++) omega[j4+t1*dim]=allinvomega[(i7*dim+j4)+t1*(dim*nclus)];
	}
	res=res+wishart_dens(a,dim,omega,invA,help,help2);
}
return res;


}


/******************************************************************************/


double derive_log_f_u ( double dx, double eta, double u, int dim, int nclus, double allomega[], double omega[], double invA[], double help[], double help2[], double gamma, double Gammastar[])



/******************************************************************************/


{
  

double res;
res=(log_f_u(eta, (u+dx), dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar)-log_f_u(eta, (u-dx), dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar))/(2*dx);
return res;

}


/******************************************************************************/


double derive2_log_f_u ( double dx, double eta, double u, int dim, int nclus, double allomega[], double omega[], double invA[], double help[], double help2[], double gamma, double Gammastar[])



/******************************************************************************/


{
  

double res;
res=(log_f_u(eta, (u+dx), dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar)-2*log_f_u(eta, u, dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar)+log_f_u(eta, (u-dx), dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar))/(dx*dx);
return res;

}


/******************************************************************************/


double derive2_f_u ( double dx, double eta, double u, int dim, int nclus, double allomega[], double omega[], double invA[], double help[], double help2[], double gamma, double Gammastar[], double K)



/******************************************************************************/


{
  

double res;
res=(exp(K+log_f_u(eta, (u+dx), dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar))-2*exp(K+log_f_u(eta, u, dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar))+exp(K+log_f_u(eta, (u-dx), dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar)))/(dx*dx);
return res;

}


/******************************************************************************/

double newton_raphson ( double x, double precision, double dx, double eta, int dim, int nclus, double allomega[], double omega[], double invA[], double help[], double help2[], double gamma, double Gammastar[])



/******************************************************************************/


{
  

double res=0,d1,d2;
int NMAX=50;
int i8;
int flag=0;

for (i8=0;i8<NMAX;i8++) {
	if (flag==0) {
		d1=derive_log_f_u(dx, eta, x, dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar);
		d2=derive2_log_f_u(dx,eta, x, dim, nclus, allomega, omega,  invA,  help,  help2, gamma, Gammastar);
		res=x-(d1/d2);
		if (fabs((res-x)/res)<precision) flag=1;
		x=res;
	}
}
if (flag==0) {
	res=-9999;
}
return res;

}


/******************************************************************************/


double h_u ( double u, double u_m, double lambda )



/******************************************************************************/


{
  

double res;
res=pow((1+((u-u_m)*(u-u_m))/(4*lambda*lambda)),(-5/2));
return res;

}


/******************************************************************************/



double t_sample ( double df , int fl)



/******************************************************************************/


{
  

double arg1;
  
double arg2;
  
double value;

  
if ( df <= 0.0 )
  {
    
	Rprintf ( "\n" );
    
	Rprintf ( "R8_CHI_SAMPLE - Fatal error!\n" );
    
	Rprintf ( "  DF <= 0.\n" );
    
	Rprintf ( "  Value of DF: %g\n", df );
    
 }

  	
arg1=r8_normal_01_sample( fl );
arg2=r8_chi_sample(df, fl)/df;
value=arg1/sqrt(arg2);
return value;

}
