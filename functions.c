#include "final_project.h"

// check if string has only numeric values
int stringNumericCheck(char *str) {
        int i = 0;
        int notNumericCount = 0;

        for (i = 0; i < strlen(str); i++) {
                notNumericCount += !(isdigit(str[i]));
        }

        return !(notNumericCount > 0);
}

// https://en.wikipedia.org/wiki/Binomial_coefficient#Binomial_coefficient_in_programming_languages 
unsigned long binom(int n, int k) {
	unsigned long c;
	int i = 0;

	if (k < 0 || k > n)
		return 0;
	if (k == 0 || k == n)
		return 1;

	// k = min(k, n-k)
	k = min(k, (n - k));

	c = 1;

	for (i = 0; i < k; i++)
		c = c * (n - i) / (i + 1);

	return c;

}

int max(int a, int b) {
	return a * (a >= b) + b * (a < b);
}

int min(int a, int b) {
	return a * (a < b) + b * (a >= b);
}

// generate the polynomial of interest by taking 
// the numerator coefficients and powers and
// the denominator coefficients and powers
double generate_poly(double theta, double alpha,
			unsigned long **numer, int numer_rows,
				unsigned long **denom, int denom_rows) {
	int i;
	double num_val, denom_val, out;

	num_val = 0;
	denom_val = 0;

	// numerator coefficients and powers
	for (i = 0; i < numer_rows; i++)
		num_val = num_val + numer[i][0] * pow(theta, numer[i][1]);

	// denominator coefficients and powers
	for (i = 0; i < denom_rows; i++)
		denom_val = denom_val + denom[i][0] * pow(theta, denom[i][1]);

	out = num_val + (alpha/2) * denom_val;

	return out;
}

