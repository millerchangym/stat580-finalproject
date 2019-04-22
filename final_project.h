#ifndef final_project_h
#define final_project_h
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef double (*function)(double, double,
			unsigned long **, int,
				unsigned long **, int);
int stringNumericCheck(char *str);
int max(int a, int b);
int min(int a, int b);
unsigned long binom(int n, int k);
double generate_poly(double theta, double alpha,
			unsigned long **numer, int numer_rows,
				unsigned long **denom, int denom_rows);
double generate_poly_deriv(double theta, double alpha,
				unsigned long **numer, int numer_rows,
					unsigned long **denom, int denom_rows);

#endif
