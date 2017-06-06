# include <stdlib.h>

# include <Rinternals.h>
# include <R.h>

# include <math.h>

# include <time.h>

# include <string.h>


# include "wishart.h"

# include "pdflib.h"





void r8mat_add ( int m, int n, double a[], double b[] )



/******************************************************************************/

/*
  Purpose:

 
 
  R8MAT_ADD adds one R8MAT to another.



    Discussion:

  

  An R8MAT is a doubly dimensioned array of R8 values, stored as a vector

    in column-major order.

  

    Licensing:

   

  This code is distributed under the GNU LGPL license.

  

    Modified:

 
   
  31 July 2013

 
  
    Author:

    
  
  John Burkardt

 
 
    Parameters:

    Input, int M, N, the number of rows and columns.

    
Input, double A[M*N], the matrix to add.

    Input/output, double B[M*N], the matrix to be incremented.
*/


{
  
	int i;
  
	int j;

  
	for ( j = 0; j < n; j++ )
  {
    
		for ( i = 0; i < m; i++ )
    {
      
			b[i+j*m] = b[i+j*m] + a[i+j*m];
    
			}
  
		}
  
	return;

}



/******************************************************************************/


void r8mat_cholesky_factor_upper ( int n, double a[], double c[], int *flag )



/******************************************************************************/


/*
  Purpose:

    R8MAT_CHOLESKY_FACTOR_UPPER: upper Cholesky factor of a symmetric R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix must be symmetric and positive semidefinite.

    For a positive semidefinite symmetric matrix A, the Cholesky factorization
    is an upper triangular matrix R such that:

      A = R' * R

    Note that the usual Cholesky factor is a LOWER triangular matrix L
    such that

      A = L * L'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix A.

    Input, double A[N*N], the N by N matrix.

    Output, int *FLAG, an error flag.
    0, no error was detected.
    1, the matrix was not positive definite.  A NULL factor was returned.

    Output, double R8MAT_CHOLESKY_FACTOR_UPPER[N*N], the N by N upper triangular
    "Choresky" factor.
*/


{
  
 
int i;
  
int j;
  
int k;
  
double sum2;

  
*flag = 0;

  
r8mat_copy_new ( n, n, a, c );

  
for ( j = 0; j < n; j++ )
  {
    
	for ( i = 0; i < j; i++ )
    {
      
		c[j+i*n] = 0.0;
    
		}
    
	for ( i = j; i < n; i++ )
    {
      
		sum2 = c[i+j*n];
      
		for ( k = 0; k < j; k++ )
      {
        
			sum2 = sum2 - c[k+j*n] * c[k+i*n];
      
			}
      
		if ( i == j )
      {
        
			if ( sum2 <= 0.0 )
        {
          
				*flag = 1;
          
				return;
        
				}
        
			c[j+i*n] = sqrt ( sum2 );
      
			}
      
		else
      {
        
			if ( c[j+j*n] != 0.0 )
        
				{
          
				c[j+i*n] = sum2 / c[j+j*n];
        
				}
        
			else
        
				{
          
				c[j+i*n] = 0.0;
        
				}
      
			}
    
		}
  
	}

  
return;

}


/******************************************************************************/


void r8mat_copy_new ( int m, int n, double a1[], double a2[] )



/******************************************************************************/


/*
  
Purpose:

    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
*/


{
  

int i;
  
int j;

  
for ( j = 0; j < n; j++ )
  {
    
	for ( i = 0; i < m; i++ )
    {
     
		a2[i+j*m] = a1[i+j*m];
    
		}
  
	}

  
return;

}



/******************************************************************************/



void r8mat_divide ( int m, int n, double s, double a[] )



/******************************************************************************/


/*
  Purpose:

    R8MAT_DIVIDE divides an R8MAT by a scalar.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double S, the divisor

    Input/output, double A[M*N], the matrix to be scaled.
*/


{
  

int i;
  
int j;

  
for ( j = 0; j < n; j++ )
  {
    
	for ( i = 0; i < m; i++ )
    {
      
		a[i+j*m] = a[i+j*m] / s;
    
		}
  
	}
  
return;

}


/******************************************************************************/

void r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[], double c[] )



/******************************************************************************/


/*
  Purpose:

    R8MAT_MM_NEW multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
*/


{
  

int i;
  
int j;
  
int k;

  
for ( i = 0; i < n1; i++ )
  {
    
	for ( j = 0; j < n3; j++ )
    {
      
		c[i+j*n1] = 0.0;
      
		for ( k = 0; k < n2; k++ )
      {
        
			c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      
			}
    
		}
  
	}

  
return;

}


/******************************************************************************/



void r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[], double c[] )



/******************************************************************************/


/*
  Purpose:

    R8MAT_MMT_NEW computes C = A * B'.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N3*N2], the matrices to multiply.

    Output, double R8MAT_MMT[N1*N3], the product matrix C = A * B'.
*/


{
   

int i;
  
int j;
  
int k;

  
for ( i = 0; i < n1; i++ )
  {
    
	for ( j = 0; j < n3; j++ )
    {
      
		c[i+j*n1] = 0.0;
      
		for ( k = 0; k < n2; k++ )
      {
        
			c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[j+k*n3];
      
			}
    
		}
  
	}

  
return;

}


/******************************************************************************/



void r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[], double c[] )



/******************************************************************************/


/*
  Purpose:

    R8MAT_MTM_NEW computes C = A' * B.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N2*N1], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MTM_NEW[N1*N3], the product matrix C = A' * B.
*/


{
  

int i;
  
int j;
  
int k;

  
for ( i = 0; i < n1; i++ )
  {
    
	for ( j = 0; j < n3; j++ )
    {
      
		c[i+j*n1] = 0.0;
      
		for ( k = 0; k < n2; k++ )
      {
        
			c[i+j*n1] = c[i+j*n1] + a[k+i*n2] * b[k+j*n2];
      
			}
    
		}
  
	}

  
return;

}


/******************************************************************************/



void r8mat_print ( int m, int n, double a[], char *title )



/******************************************************************************/


/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/


{
  
r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  
return;

}


/******************************************************************************/



void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )



/******************************************************************************/


/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/


{

# define INCX 5

  
int i;
  
int i2hi;
  
int i2lo;
  
int j;
  
int j2hi;
  
int j2lo;

  
Rprintf ( "\n" );
  
Rprintf ("%s\n", title );

  
if ( m <= 0 || n <= 0 )
  {
    
	Rprintf (  "\n" );
    
	Rprintf (  "  (None)\n" );
    
	return;
  
	}

/*
  Print the columns of the matrix, in strips of 5.
*/
  

for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    
	j2hi = j2lo + INCX - 1;
    
	if ( n < j2hi )
    {
      
		j2hi = n;
    
		}
    
	if ( jhi < j2hi )
    {
      
		j2hi = jhi;
    
		}

    
	Rprintf ( "\n" );


/*
  For each column J in the current range...

  Write the header.
*/


	Rprintf ( "  Col:  ");
    
	for ( j = j2lo; j <= j2hi; j++ )
    {
      
		Rprintf (  "  %7d     ", j - 1 );
    
		}
    
	Rprintf ( "\n" );
    
	Rprintf (  "  Row\n" );
    
	Rprintf (  "\n" );


/*
  Determine the range of the rows in this strip.
*/
    
	if ( 1 < ilo )
    {
      
		i2lo = ilo;
    
		}
    
	else
    {
      
		i2lo = 1;
    
		}
    
	if ( m < ihi )
    {
      
		i2hi = m;
    
		}
    
	else
    {
      
		i2hi = ihi;
    
		}

    
	for ( i = i2lo; i <= i2hi; i++ )
    {

/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
  
		Rprintf ( "%5d:", i - 1 );
      
		for ( j = j2lo; j <= j2hi; j++ )
      {
        
			Rprintf (  "  %14f", a[i-1+(j-1)*m] );
      
			}
      
		Rprintf (  "\n" );
    
		}
  
	}

  
return;

# undef INCX

}



/******************************************************************************/



void wishart_sample ( int m, int df, double sigma[], double a[], double au[], double aur[], double r[], double h[], int fl )



/******************************************************************************/


/*
  Purpose:

    WISHART_SAMPLE samples the Wishart distribution.

  Discussion:

    This function requires functions from the PDFLIB and RNGLIB libraries.

    The "initialize()" function from RNGLIB must be called before using
    this function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt

  Reference:

    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.

  Parameters:

    Input, int M, the order of the matrix.

    Input, int DF, the number of degrees of freedom.
    M <= DF.

    Input, double SIGMA[M*M], the covariance matrix, which should be 
    a symmetric positive definite matrix.

    Output, double WISHART_SAMPLE[M*M], the sample matrix from 
    the Wishart distribution.
*/


{
  

 
//int flag;
  
if ( df < m )
  {
    
	Rprintf ( "\n" );
    
	Rprintf ( "WISHART_SAMPLE - Error!\n" );
    
	Rprintf ( "  DF = %d < M = %d.\n Setting df=m instead.\n", df, m );
	
	df=m;
 
	}


/*
  Get R, the upper triangular Cholesky factor of SIGMA.
*/
  

r8mat_pofac ( m, sigma, r ,22);

  
/*if ( flag != 0 )
  {
    
	Rprintf ( "\n" );
    
	Rprintf ( "WISHART_SAMPLE - Fatal error!\n" );
    
	Rprintf ( 
      "  Unexpected error return from R8MAT_CHOLESKY_FACTOR_UPPER.\n" );
    
	Rprintf ( "  FLAG = %d\n", flag );
  
  	
	}


*/
/*
  Get AU, a sample from the unit Wishart distribution.
*/
  

wishart_unit_sample ( m, df, au, h, fl);


/*
  Construct the matrix A = R' * AU * R.
*/
  

r8mat_mm_new ( m, m, m, au, r,aur );
  
r8mat_mtm_new ( m, m, m, r, aur,a );


return;


}



/******************************************************************************/


void wishart_unit_sample ( int m, int df, double a[], double c[] , int fl)



/******************************************************************************/


/*
  Purpose:

    WISHART_UNIT_SAMPLE samples the unit Wishart distribution.

  Discussion:

    This function requires functions from the PDFLIB and RNGLIB libraries.

    The "initialize()" function from RNGLIB must be called before using
    this function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 October 2013

  Author:

    John Burkardt

  Reference:

    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.

  Parameters:

    Input, int M, the order of the matrix.

    Input, int DF, the number of degrees of freedom.
    M <= DF.

    Output, double WISHART_UNIT_SAMPLE[M*M], the sample matrix from the 
    unit Wishart distribution.
*/


{
  

double df_chi;
  
int i;
  
int j;

  
if ( df < m )
  {
    
	Rprintf ( "\n" );
    
	Rprintf ( "  DF = %d < M = %d.\n Setting df=m instead.", df, m );
	
	df=m;
 
  
	}

  
for ( i = 0; i < m; i++ )
  {
    
	for ( j = 0; j < i; j++ )
    {
      
		c[i+j*m] = 0.0;
    
		}
    
	df_chi = ( double ) ( df - i );
    
	c[i+i*m] = sqrt ( r8_chi_sample ( df_chi, fl ) );
    
	for ( j = i + 1; j < m; j++ )
    {
      
		c[i+j*m] = r8_normal_01_sample ( fl );
    
		}
  
	}

  
r8mat_mtm_new ( m, m, m, c, c,a );

   

return;


}
