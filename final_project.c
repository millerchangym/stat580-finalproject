#include "final_project.h"

// compile as
// gcc -o final_project final_project.c functions.c -pedantic -Wall -lm
// run as ./final_project a b c d -alpha e, where a, b, c, and d are
// integers, and e is a float

int main(int argc, char *argv[])
{
	/////////////////// ERROR CHECKING ///////////////////

	// if the number of parameters is less than 5 (i.e., less
	// than just the file name of the .exe and 4 numbers),
	// return error
	if (argc <= 6 || argc == 8 || argc >= 10) {
		printf("Please input a 2-by-2 matrix of counts in the form ");
		printf("a b c d in row major order,\n");
		printf("with required switch -theta e ");
		printf("(e an initial guess for the odds ratios)\n");
		printf("and optional switch -alpha f (f between 0 and 1 ");
		printf("inclusive.)\n");
		printf("Note that f is set to 0.05 by default.\n");
		exit(1);
	}

	// check for numeric matrix counts
	if (!stringNumericCheck(argv[1]) ||
		!stringNumericCheck(argv[2]) ||
			!stringNumericCheck(argv[3]) ||
				!stringNumericCheck(argv[4])) {
		printf("Please use numeric inputs for ");
		printf("the four matrix components.\n");
		exit(1);

	};

	if (argc >= 6) {
		// check for use of the correct switch after the numeric counts
		if (strcmp(argv[5], "-alpha") != 0 && strcmp(argv[5], "-theta") != 0) {
			printf("Your first switch after the matrix entries ");
			printf("is not a valid switch. ");
			printf("Use either -alpha or -theta as switches.\n");
			exit(1);
		}
	}

	// check for a numeric 7th component, after the first switch
	if (argc >= 7) {
		if (!stringNumericCheck(argv[6])) {
			printf("Please use a numeric value for the first switch.\n");
			exit(1);
		}
	}

	// check for use of the correct second switch
	if (argc >= 8) {
		if (strcmp(argv[7], "-alpha") != 0 && strcmp(argv[7], "-theta") != 0) {
			printf("Your second switch after the matrix entries ");
			printf("is not a valid switch. ");
			printf("Use either -alpha or -theta as switches.\n");
			exit(1);
		}

		if ((strcmp(argv[5], "-alpha") == 0 &&
			strcmp(argv[7], "-alpha") == 0) ||
				(strcmp(argv[5], "-theta") == 0 &&
					strcmp(argv[7], "-theta") == 0)) {
			printf("You cannot have two of the same switch.\n");
			exit(1);
		}
	}

	// check for numeric value after second switch
	if (argc >= 9) {
		if (!stringNumericCheck(argv[8])) {
			printf("Please use a numeric value for the second switch.\n");
			exit(1);
		}
	}


	/////////////////// VARIABLE DECLARATION ///////////////////

	// alpha and level for confidence interval
	double alpha, level, theta_0;
	// input: four components of a 2 x 2 contingency
	// table, in row-major order
	int a, b, c, d;
	int n_11;
	int n_1plus, n_2plus, n_plus1, n_plus2;
	int n;
	int m_plus, m_minus;
	unsigned long **numer_coeff_upper, **numer_coeff_lower;
	unsigned long **denom_coeff;
	int n_numer_coeff_upper, n_numer_coeff_lower;
	int n_denom_coeff;
	int i;
	int t, u;
	function poly, poly_deriv;
	long double estimate;

	// Use alpha = 5% by default
	alpha = 0.05;
	if (argc >= 7 && strcmp(argv[5], "-alpha") == 0)
		alpha = atof(argv[6]);
	if (argc >= 9 && strcmp(argv[7], "-alpha") == 0)
		alpha = atof(argv[8]);

	// use nonsensical value for theta_0 by default
	// for error checking
	theta_0 = -9999;
	if (argc >= 7 && strcmp(argv[5], "-theta") == 0)
		theta_0 = atof(argv[6]);
	if (argc >= 9 && strcmp(argv[7], "-theta") == 0)
		theta_0 = atof(argv[8]);
	if (theta_0 == -9999) {
		printf("No value for theta has been provided. ");
		printf("Provide a value for theta.\n");
		exit(1);
	}

	if (theta_0 < 0) {
		printf("Provide a non-negative value for theta.\n");
		exit(1);
	}

	// check for alpha in (0, 1)
	if (!(alpha > 0 && alpha < 1)) {
		printf("Please input a value of alpha between ");
		printf("0 and 1 exclusive.\n");
		exit(1);
	}

	level = 1 - alpha;

	a = atoi(argv[1]);
	b = atoi(argv[2]);
	c = atoi(argv[3]);
	d = atoi(argv[4]);

	if (a < 0 || b < 0 || c < 0 || d < 0) {
		printf("The integer inputs must all be positive.\n");
		exit(1);
	}

	// necessary variables
	n_11 = a;

	// Row and column totals
	n_1plus = a + b;
	n_2plus = c + d;
	n_plus1 = a + c;
	n_plus2 = b + d;

	printf("Assuming fixed row and column totals...\n");

	n = a + b + c + d;

	// variables for PMF
	m_minus = max(0, n_1plus + n_plus1 - n);
	m_plus = min(n_1plus, n_plus1);

	/////// POLYNOMIAL DETERMINATION (NUMERATOR) ///////

	// Note that the denominator is not dependent on
	// t, which is the variable to be summed for
	// computation of the 100(1-alpha)% confidence
	// interval

	// note m_minus <= t <= m_plus

	// t >= n_11 for upper
	n_numer_coeff_upper = m_plus - n_11 + 1;
	// t <= n_11 for lower
	n_numer_coeff_lower = n_11 - m_minus + 1;

	// memory allocation of numerator coefficients
	numer_coeff_lower = (unsigned long **)
				calloc(n_numer_coeff_lower,
					sizeof(unsigned long *));
	if (!numer_coeff_lower) {
		printf("Calloc for numer_coeff_lower failed !\n");
		exit(1);
	}

	numer_coeff_upper = (unsigned long **)
				calloc(n_numer_coeff_upper,
					sizeof(unsigned long *));

	if (!numer_coeff_upper) {
		printf("Calloc for numer_coeff_upper failed !\n");
		exit(1);
	}

	// determine lower-bound confidence interval coefficients
	for (i = 0; i < n_numer_coeff_lower; i++) {
		numer_coeff_lower[i] = (unsigned long *)
			calloc(2, sizeof(unsigned long));
		if (!numer_coeff_lower[i]) {
			printf("Calloc for numer_coeff_lower[i] failed !\n");
			exit(1);
		}
	}

	for (i = 0; i < n_numer_coeff_lower; i++) {
		t = m_minus + i;
		numer_coeff_lower[i][0] = binom(n_1plus, t) *
						binom(n - n_1plus, n_plus1 - t);
		numer_coeff_lower[i][1] = t;
	}

	// determine upper-bound confidence interval coefficients
	for (i = 0; i < n_numer_coeff_upper; i++) {
		numer_coeff_upper[i] = (unsigned long *)
			calloc(2, sizeof(unsigned long));
		if (!numer_coeff_upper[i]) {
			printf("Calloc for numer_coeff_upper[i] failed !\n");
			exit(1);
		}
	}

	for (i = 0; i < n_numer_coeff_upper; i++) {
		t = n_11 + i;
		numer_coeff_upper[i][0] = binom(n_1plus, t) *
						binom(n - n_1plus, n_plus1 - t);
		numer_coeff_upper[i][1] = t;
	}

	/////// POLYNOMIAL DETERMINATION (DENOMINATOR) ///////

	n_denom_coeff = m_plus - m_minus + 1;
	denom_coeff = (unsigned long **)
			calloc(n_denom_coeff, sizeof(unsigned long *));
	if (!denom_coeff) {
		printf("Calloc for denom_coeff failed !\n");
		exit(1);
	}

	for (i = 0; i < n_denom_coeff; i++) {
		denom_coeff[i] = (unsigned long *)
			calloc(2, sizeof(unsigned long));
		if (!denom_coeff[i]) {
			printf("Calloc for denom_coeff[i] failed !\n");
			exit(1);
		}
	}

	for (i = 0; i < n_denom_coeff; i++) {
		u = m_minus + i;
		denom_coeff[i][0] = binom(n_1plus, u) *
					binom(n - n_1plus, n_plus1 - u);
		denom_coeff[i][1] = u;
	}

	////////// GENERATE POLYNOMIALS ///////////

	poly = generate_poly;
	poly_deriv = generate_poly_deriv;

	//// IMPLEMENT NEWTON-RAPHSON FOR LOWER BOUND ////
	printf("****** BEGIN NEWTON-RAPHSON ******\n");
	estimate = newton_raphson(poly, poly_deriv, theta_0, alpha,
			numer_coeff_upper, n_numer_coeff_upper,
			denom_coeff, n_denom_coeff);

	////////// FREE MEMORY SPACE //////////

	for (i = 0; i < n_numer_coeff_lower; i++) {
		free(numer_coeff_lower[i]);
		numer_coeff_lower[i] = NULL;
	}
	free(numer_coeff_lower);
	numer_coeff_lower = NULL;

	for (i = 0; i < n_numer_coeff_upper; i++) {
		free(numer_coeff_upper[i]);
		numer_coeff_upper[i] = NULL;
	}
	free(numer_coeff_upper);
	numer_coeff_upper = NULL;

	for (i = 0; i < n_denom_coeff; i++) {
		free(denom_coeff[i]);
		denom_coeff[i] = NULL;
	}
	free(denom_coeff);
	denom_coeff = NULL;

	return 0;
}

