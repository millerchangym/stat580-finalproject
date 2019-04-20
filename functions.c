#include "final_project.h"

int stringNumericCheck(char *str) {
        int i = 0;
        int notNumericCount = 0;

        for (i = 0; i < strlen(str); i++) {
                notNumericCount += !(isdigit(str[i]));
        }

        return !(notNumericCount > 0);
}

// Recursive formula: https://en.wikipedia.org/wiki/Binomial_coefficient#Recursive_formula
int binom(int n, int k) {
	if (n == k || k == 0) {
		return 1;
	}
	return binom(n-1, k-1) + binom(n-1, k);
}
