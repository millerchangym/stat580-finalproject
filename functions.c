#include "final_project.h"

// check if string has only numeric values or decimal places
int stringNumericCheck(char *str) {
        int i = 0;
        int notNumericCount = 0;
	int decimalCount = 0;
	int minusCount = 0;

        for (i = 0; i < strlen(str); i++) {
		// count of non-numeric characters, minus periods
		// and minuses
                notNumericCount += !(isdigit(str[i])) - (str[i] == '.') -
					(str[i] == '-');
		decimalCount += (str[i] == '.');
		minusCount += (str[i] == '-');
        }

        return !(decimalCount >= 2 || notNumericCount > 0 || minusCount >= 2);
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

	out = num_val - (alpha/2) * denom_val;

	return out;
}

// derivative of the generated polynomial, for Newton-Raphson
double generate_poly_deriv(double theta, double alpha,
				unsigned long **numer, int numer_rows,
					unsigned long **denom, int denom_rows) {

	int i;
	double num_val, denom_val, out;

	num_val = 0;
	denom_val = 0;

	// numerator derivative coefficients and powers
	for (i = 0; i < numer_rows; i++)
		num_val = num_val + numer[i][0] * numer[i][1] *
			pow(theta, numer[i][1] - 1);

	// denominator derivative coefficients and powers
	for (i = 0; i < denom_rows; i++)
		denom_val = denom_val + denom[i][0] * denom[i][1] *
			pow(theta, denom[i][1] - 1);

	out = num_val - (alpha/2) * denom_val;

	return out;
}

// implement Newton-Raphson
double newton_raphson(function f, function f_prime, 
			double theta_0,
				double alpha,
					unsigned long **numer, int numer_rows,
						unsigned long **denom, int denom_rows) {

	double threshold, theta_next, theta_prev;
	double func_value, deriv_value;
	int n_iter, n_iter_MAX;
	int stop;

	threshold = 0.001;
	n_iter = 0;
	n_iter_MAX = 5000;
	theta_prev = theta_0;

	do {
		n_iter++;
		printf("Iteration %d ", n_iter);
		func_value = f(theta_prev, alpha, numer, numer_rows,
				denom, denom_rows);
		deriv_value = f_prime(theta_prev, alpha, numer, numer_rows,
				denom, denom_rows);
		printf("function value: %f ", func_value);
		printf("derivative value: %f ", deriv_value);
		if (deriv_value == 0) {
			printf("ZERO DERIVATIVE ERROR! STOPPING...\n");
			exit(1);
		}

		stop = (
			fabs(theta_next - theta_prev) <= threshold ||
				n_iter >= n_iter_MAX ||
					deriv_value == 0
			);
						
		theta_next = theta_prev - (func_value/deriv_value);
		printf("New value: %f\n", theta_next);
		// fabs is the floating-point absolute value
		stop = (fabs(theta_next - theta_prev) <= threshold ||
				n_iter >= n_iter_MAX ||
				isnan(theta_next));
		if (!stop)
			theta_prev = theta_next;
		else
			printf("STOPPING...\n");
	} while (!stop);

	return theta_next;
}
