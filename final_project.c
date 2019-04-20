#include "final_project.h"

//compile as
//gcc -o final_project final_project.c -pedantic -Wall -lm
//run as ./final_project a b c d -alpha e, where a, b, c, and d are integers, and e is a float

int main(int argc, char *argv[])
{
	// if the number of parameters is less than 5 (i.e., less than just the file name of the .exe and 4 numbers),
	// return error
	if (argc < 5 || argc == 6 || argc > 7) {
		printf("Please input a 2-by-2 matrix of counts in the form a b c d in row-major order, \n");
		printf("with optional switch -alpha e, e a number between 0 and 1 exclusive.\n");
		exit(1);
	}

	if (strcmp(argv[5], "-alpha") != 0) {
		printf("The switch you used is incorrect.\n");
		printf("Use the correct switch -alpha e, with e a number between 0 and 1 exclusive.\n");
		exit(1);
	}

	// alpha and level for confidence interval
	float alpha, level;
	// input: four components of a 2 x 2 contingency table, in row-major order
	int a, b, c, d;

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

	if (!stringNumericCheck(argv[1]) ||
		!stringNumericCheck(argv[2]) ||
			!stringNumericCheck(argv[3]) || 
				!stringNumericCheck(argv[4])) {
		printf("Please use numeric inputs for the four matrix components.\n");
		exit(1);

	};

	a = atoi(argv[1]);
	b = atoi(argv[2]);
	c = atoi(argv[3]);
	d = atoi(argv[4]);

	if (a < 0 || b < 0 || c < 0 || d < 0) {
		printf("The integer inputs must all be positive.\n");
		exit(1);
	}

	return 0;
}

