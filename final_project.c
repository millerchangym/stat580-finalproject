#include "final_project.h"

//compile as
//gcc -o final_project final_project.c functions.c -pedantic -Wall -lm
//run as ./final_project a b c d -alpha e, where a, b, c, and d are integers, and e is a float

int main(int argc, char *argv[])
{
	// if the number of parameters is less than 5 (i.e., less than just the file name of the .exe and 4 numbers),
	// return error
	if (argc < 5 || argc > 7) {
		printf("Please input a 2-by-2 matrix of counts in the form a b c d in row-major order, \n");
		printf("with optional switch -alpha e, e a number between 0 and 1 exclusive.\n");
		printf("Note that e is set to 0.05 by default.");
		exit(1);
	}

	// check for numeric matrix counts
        if (!stringNumericCheck(argv[1]) ||
                !stringNumericCheck(argv[2]) ||
                        !stringNumericCheck(argv[3]) ||
                                !stringNumericCheck(argv[4])) {
                printf("Please use numeric inputs for the four matrix components.\n");
                exit(1);

        };

	if (argc == 7) {
		// check for use of the correct switch after the numeric counts
		if (strcmp(argv[5], "-alpha") != 0) {
			printf("The switch you used is incorrect.\n");
			printf("Use the correct switch -alpha e, with e a number between 0 and 1 exclusive.\n");
			exit(1);
		}

		// check for a numeric alpha value
		if(!stringNumericCheck(argv[6])) {
			printf("Please use a numeric value for the alpha.\n");
			exit(1);
		}
	}

	// alpha and level for confidence interval
	float alpha, level;
	// input: four components of a 2 x 2 contingency table, in row-major order
	int a, b, c, d;
	int n_1plus, n_2plus, n_plus1, n_plus2;
	int n;
	int m_plus, m_minus;
	unsigned long **coeff;
	int n_coeff;
	int i;
	int u;

	// Use alpha = 5% by default
	if (argc == 7) {
		alpha = atof(argv[6]);
	} else {
		alpha = 0.05;
	}

	// check for alpha in (0, 1)
	if (!(alpha > 0 && alpha < 1)) {
		printf("Please input a value of alpha between 0 and 1 exclusive.\n");
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

	// Row and column totals
	n_1plus = a + b;
	n_2plus = c + d;
	n_plus1 = a + c;
	n_plus2 = b + d;

	printf("Assuming fixed totals...\n");
	printf("Row totals: %d, %d\n", n_1plus, n_2plus);
	printf("Column totals: %d, %d\n", n_plus1, n_plus2);

	n = a + b + c + d;

	// variables for PMF
	m_minus = max(0, n_1plus + n_plus1 - n);
	m_plus = min(n_1plus, n_plus1);

	// determine coefficients
	n_coeff = m_plus - m_minus + 1;
	printf("\nNumber of coefficients: %d\n", n_coeff);
	coeff = (unsigned long **) calloc(n_coeff, sizeof(unsigned long *));
	if (!coeff) {
		printf("Calloc for coeff failed!\n");
		exit(1);
	}

	for (i = 0; i < n_coeff; i++) {
		coeff[i] = (unsigned long *) calloc(2, sizeof(unsigned long));
		if (!coeff[i]) {
			printf("Calloc for coeff[%d] failed!\n", i);
			exit(1);
		}
	}

	for (i = 0; i < n_coeff; i++) {
		u = m_minus + i;
		coeff[i][1] = binom(n_1plus, u) * binom(n - n_1plus, n_plus1 - u);
		coeff[i][2] = u;
	}

	FREE_2ARRAY(coeff);

	return 0;
}

