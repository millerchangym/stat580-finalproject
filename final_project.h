#ifndef final_project_h
#define final_project_h
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef long double (*function)(long double, double,
			unsigned long **, int,
				unsigned long **, int);
int stringNumericCheck(char *str);
int max(int a, int b);
int min(int a, int b);
unsigned long binom(int n, int k);
long double generate_poly(long double theta, double alpha,
			unsigned long **numer, int numer_rows,
				unsigned long **denom, int denom_rows);
long double generate_poly_deriv(long double theta, double alpha,
				unsigned long **numer, int numer_rows,
					unsigned long **denom, int denom_rows);
long double newton_raphson(function f, function f_prime,
			double theta_0, double alpha, 
				unsigned long **numer, int numer_rows,
					unsigned long **denom, int denom_rows);

#endif
