#include "final_project.h"

// documentation function
// simply prints documentation
void documentation() {
	printf("REQUIRED:\n");
	printf(" a b c d  -------");
	printf(" four numeric inputs in row-major order\n");
	printf("                  ");
	printf("These must precede all other inputs.\n");
	printf("-theta e  -------");
	printf(" e is an initial guess for the odds ratio\n\n");
	printf("OPTIONAL:\n");
	printf("-alpha f  -------");
	printf(" f is a number to the nearest thousandth\n");
	printf("                  ");
	printf("with 0 < f < 1, with 1 - alpha\n");
	printf("                  ");
	printf("being the level of the confidence interval.\n");
	printf("                  ");
	printf("f is set to 0.05 by default.\n");
	printf("-v g      -------");
	printf(" g is either 0 or 1. Setting g to 1 shows\n");
	printf("                  ");
	printf("all Newton-Raphson output, including function\n");
	printf("                  ");
	printf("outputs, derivative outputs, and updated theta\n");
	printf("                  ");
	printf("values.\n");
	printf("                  ");
	printf("g is set to 0 by default.\n");
}

// check if string has only numeric values,
// a decimal point, or a minus sign
int stringNumericCheck(char *str)
{
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

// Binomial coefficient
// https://en.wikipedia.org/wiki/Binomial_coefficient#Binomial_coefficient_in_programming_languages
double binom(int n, int k)
{
	double c;
	int i = 0;

	if (k < 0 || k > n)
		return 0;
	if (k == 0 || k == n)
		return 1;

	// k = min(k, n-k)
	k = minInt(k, (n - k));

	c = 1;

	for (i = 0; i < k; i++)
		c = c * (n - i) / (i + 1);

	return c;

}

// max for integers
int maxInt(int a, int b)
{
	return a * (a >= b) + b * (a < b);
}

// min for integers
int minInt(int a, int b)
{
	return a * (a < b) + b * (a >= b);
}

// max for doubles
double max(double a, double b)
{
	return a * (a >= b) + b * (a < b);
}

// min for doubles
double min(double a, double b)
{
	return a * (a < b) + b * (a >= b);
}


// generate the polynomial of interest by taking
// the numerator coefficients and powers and
// the denominator coefficients and powers
double generate_poly(double theta, double alpha,
			double **numer, int numer_rows,
				double **denom, int denom_rows)
{
	int i;
	double num_val, denom_val, out;

	num_val = 0;
	denom_val = 0;

	// numerator coefficients and powers
	for (i = 0; i < numer_rows; i++)
		num_val += numer[i][0] * pow(theta, numer[i][1]);

	// denominator coefficients and powers
	for (i = 0; i < denom_rows; i++)
		denom_val += denom[i][0] * pow(theta, denom[i][1]);

	out = num_val - (alpha/2) * denom_val;

	return out;
}

// derivative of the generated polynomial, for Newton-Raphson
double generate_poly_deriv(double theta, double alpha,
				double **numer, int numer_rows,
					double **denom, int denom_rows)
{

	int i;
	double num_val, denom_val, out;

	num_val = 0;
	denom_val = 0;

	// numerator derivative coefficients and powers
	for (i = 0; i < numer_rows; i++) {
		if (numer[i][1] == 0)
			num_val = num_val;
		else
			num_val += numer[i][0] * numer[i][1] *
					pow(theta, numer[i][1] - 1);
	}

	// denominator derivative coefficients and powers
	for (i = 0; i < denom_rows; i++) {
		if (denom[i][1] == 0)
			denom_val = denom_val;
		else
			denom_val += denom[i][0] * denom[i][1] *
				pow(theta, denom[i][1] - 1);
	}

	out = num_val - (alpha/2) * denom_val;

	return out;
}

// implement Newton-Raphson
double newton_raphson(function f, function f_prime,
			double theta_0,
				double alpha,
					double **numer, int numer_rows,
						double **denom, int denom_rows,
							int verbose)
{

	double threshold;
	double theta_next, theta_prev;
	double func_value, deriv_value;
	int n_iter, n_iter_MAX;
	int stop;

	// threshold for stopping (i.e., differences
	// between two iterations of the parameter
	// theta)
	threshold = 0.001;

	// iteration counter
	n_iter = 0;

	// max number of iterations allowed
	n_iter_MAX = 5000;
	theta_prev = theta_0;

	// verbose prints additional information about
	// the Newton-Raphson procedure when it is 1
	if (verbose)
		printf("Initial theta: %f\n", theta_0);

	do {
		n_iter++;
		if (verbose)
			printf("Iteration %d ", n_iter);

		// compute the function and its derivative evaluated
		// at the given value of theta (i.e., theta_prev)
		func_value = f(theta_prev, alpha, numer, numer_rows,
				denom, denom_rows);
		deriv_value = f_prime(theta_prev, alpha, numer, numer_rows,
				denom, denom_rows);
		if (verbose) {
			printf("function value: %f ", func_value);
			printf("derivative value: %f ", deriv_value);
		}

		// avoids division by zero
		if (deriv_value == 0) {
			printf("ZERO DERIVATIVE ERROR! STOPPING...\n");
			exit(1);
		}

		// Newton update equation
		theta_next = theta_prev - (func_value/deriv_value);

		if (verbose)
			printf("New theta: %f\n", theta_next);

		// fabs is the floating-point absolute value
		// stop is a boolean used to tell the program to stop
		// executing Newton-Raphson when the threshold is met,
		// or you've exceeded the number of iterations allowed,
		// or you've hit a zero-derivative point
		stop = (fabs(theta_next - theta_prev) <= threshold ||
				n_iter >= n_iter_MAX ||
				isnan(theta_next) ||
				deriv_value == 0);

		if (!stop) // keep updating if stop is not 1
			theta_prev = theta_next;
		else if (verbose) // otherwise, print if verbose
			printf("STOPPING...\n");

	} while (!stop);

	// print error for negative odds ratio
	if (theta_next < 0) {
		printf("ERROR: Newton-Raphson generated a negative ");
		printf("estimate for theta. Use a different value ");
		printf("for the theta switch.\n");
		exit(1);
	}

	// print error indicating the threshold
	// was not met, or the number of iterations
	// exceeded the number allowable
	if (fabs(theta_next - theta_prev) > threshold ||
		n_iter >= n_iter_MAX) {
		printf("ERROR: Newton-Raphson did not converge. ");
		printf("Use a different value for the theta switch.\n");
		exit(1);
	}

	// print error indicating a zero derivative
	if (deriv_value == 0) {
		printf("ERROR: Newton-Raphson did not converge ");
		printf("due to a zero-derivative point.");
		printf("Use a different value for the theta switch.\n");
		exit(1);
	}

	return theta_next;
}
