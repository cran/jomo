
void r8mat_add ( int m, int n, double a[], double b[] );

void r8mat_cholesky_factor_upper ( int n, double a[], double c[], int *flag );
void r8mat_copy_new ( int m, int n, double a1[], double a2[] );

void r8mat_divide ( int m, int n, double s, double a[] );
void r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[], double c[]);
void r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[], double c[] );
 
void r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[], double c[] );

void r8mat_print ( int m, int n, double a[], char *title );

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void wishart_sample ( int m, int df, double sigma[], double a[], double au[], double aur[], double r[], double h[], int fl );
void wishart_unit_sample ( int m, int df, double a[], double c[] , int fl);
