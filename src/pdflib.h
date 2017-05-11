double normal_cdf(double value);
int checkposdef(int dim, double matr[], double matrh[],double matrh2[]); 
double argmaxvec(int card, double vec[]);
double maxvec(int card, double vec[]); 
double r8_chi_sample ( double df , int fl);

double r8_epsilon ( void );


double r8_exponential_01_sample ( int fl );

double r8_gamma_sample ( double a, double r, int fl );


double r8_gamma_01_sample ( double a, int fl );

double r8_max ( double x, double y );
double r8_min ( double x, double y );

double r8_normal_sample ( double av, double sd, int fl );
double r8_normal_01_sample ( int fl );



double r8_uniform_sample ( double low, double high, int fl );

double r8_uniform_01_sample ( int fl );


double r8mat_podet ( int n, double r[] );

void r8mat_pofac ( int n, double a[], double r[], int indica );
void r8mat_poinv ( int n, double r[], double b[] );

void r8vec_multinormal_sample ( int n, double mu[], double r[],double x[], double z[] , int fl);
double r8_gamma_log ( double x ); 
double log_mul_gamma ( int p, double a ); 
double r8_chi_pdf ( double df, double rval ); 
double wishart_dens ( double df, int dim, double X[], double invA[], double help[], double help2[]); 
double log_f_u ( double eta, double u, int dim, int nclus, double allinvomega[], double omega[], double invgamma[], double help[], double help2[]);
double derive_log_f_u ( double dx, double eta, double u, int dim, int nclus, double allinvomega[], double omega[], double invgamma[], double help[], double help2[]);
double derive2_log_f_u ( double dx, double eta, double u, int dim, int nclus, double allinvomega[], double omega[], double invgamma[], double help[], double help2[]);
double derive2_f_u ( double dx, double eta, double u, int dim, int nclus, double allinvomega[], double omega[], double invgamma[], double help[], double help2[], double K);
double newton_raphson ( double x, double precision, double dx, double eta, int dim, int nclus, double allinvomega[], double omega[], double invgamma[], double help[], double help2[]); 
double h_u ( double u, double u_m, double lambda ); 
double t_sample ( double df , int fl);
