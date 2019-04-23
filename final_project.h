#ifndef final_project_h
#define final_project_h
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef double (*function)(double, double,
			double **, int,
				double **, int);
int stringNumericCheck(char *str);
int maxInt(int a, int b);
int minInt(int a, int b);
double max(double a, double b);
double min(double a, double b);
double binom(int n, int k);
double generate_poly(double theta, double alpha,
			double **numer, int numer_rows,
				double **denom, int denom_rows);
double generate_poly_deriv(double theta, double alpha,
				double **numer, int numer_rows,
					double **denom, int denom_rows);
double newton_raphson(function f, function f_prime,
			double theta_0, double alpha,
				double **numer, int numer_rows,
					double **denom, int denom_rows,
						int verbose);

#endif
