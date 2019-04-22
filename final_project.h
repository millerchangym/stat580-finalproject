#ifndef final_project_h
#define final_project_h
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int stringNumericCheck(char *str);
int max(int a, int b);
int min(int a, int b);
unsigned long binom(int n, int k);
double generate_poly(double theta,
			unsigned long **numer, int numer_rows,
				unsigned long **denom, int denom_rows);

#endif
