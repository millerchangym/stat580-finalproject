#include "final_project.h"

int stringNumericCheck(char *str) {
        int i = 0;
        int notNumericCount = 0;

        for (i = 0; i < strlen(str); i++) {
                notNumericCount += !(isdigit(str[i]));
        }

        return !(notNumericCount > 0);
}

// https://en.wikipedia.org/wiki/Binomial_coefficient#Binomial_coefficient_in_programming_languages 
unsigned long binom(unsigned long n, unsigned long k) {
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
		c *= (n - i) / (i + 1);

	return c;

}

int max(int a, int b) {
	return a * (a >= b) + b * (a < b);
}

int min(int a, int b) {
	return a * (a < b) + b * (a >= b);
}
