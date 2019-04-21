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
	if (argc < 5 || argc > 7) {
		printf("Please input a 2-by-2 matrix of counts in the form ");
		printf("a b c d in row major order,\n");
		printf("with optional switch -alpha e,");
		printf(" e a number between 0 and 1 exclusive.\n");
		printf("Note that e is set to 0.05 by default.");
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

	if (argc == 7) {
		// check for use of the correct switch after the numeric counts
		if (strcmp(argv[5], "-alpha") != 0) {
			printf("The switch you used is incorrect.\n");
			printf("Use the correct switch -alpha e, with ");
			printf("e a number between 0 and 1 exclusive.\n");
			exit(1);
		}

		// check for a numeric alpha value
		if (!stringNumericCheck(argv[6])) {
			printf("Please use a numeric value for the alpha.\n");
			exit(1);
		}
	}

	/////////////////// VARIABLE DECLARATION ///////////////////

	// alpha and level for confidence interval
	float alpha, level;
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

	// Use alpha = 5% by default
	if (argc == 7)
		alpha = atof(argv[6]);
	else
		alpha = 0.05;

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


	// FREE MEMORY SPACE //
	FREE_2ARRAY(denom_coeff);
	FREE_2ARRAY(numer_coeff_lower);
	FREE_2ARRAY(numer_coeff_upper);
	return 0;
}

