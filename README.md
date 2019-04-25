# Final Project for STAT 580 - Yeng Miller-Chang

Please see the file final_project_writeup.pdf for more information.

The project should be compiled as

```
gcc -o final_project final_project.c functions.c -pedantic -Wall -lm
```
and should be run as

```
./final_project a b c d -theta e
```
where a, b, c, and d are integers, and e is a double. There are also
optional switches -v (1 for verbose, 0 otherwise; 0 default), and
-alpha (for which 1 - alpha is the confidence level) which must
be to the nearest thousandth.

## Examples

```
./final_project 8 38 48 89 -theta 0.05

Input matrix:
8 38
48 89

Assuming fixed row and column totals...
Row totals: 46, 137
Column totals: 56, 127

ERROR: Newton-Raphson generated a negative estimate for theta. Use a different value for the theta switch.
```
```
./final_project 8 38 48 89 -theta 2

Input matrix:
8 38
48 89

Assuming fixed row and column totals...
Row totals: 46, 137
Column totals: 56, 127

The 95.0% confidence interval for the odds ratio is [0.146, 0.943].
```
```
./final_project 8 38 48 89 -theta 2 -alpha 0.02

Input matrix:
8 38
48 89

Assuming fixed row and column totals...
Row totals: 46, 137
Column totals: 56, 127

The 98.0% confidence interval for the odds ratio is [0.121, 1.083].
```
```
./final_project 2 2 3 4 -theta 50 -alpha 0.02 -v 1

Input matrix:
2 2
3 4

Assuming fixed row and column totals...
Row totals: 4, 7
Column totals: 5, 6

****** BEGIN NEWTON-RAPHSON ******
Initial theta: 50.000000
Iteration 1 function value: -15799.210000 derivative value: -20371.400000 New theta: 49.224442
Iteration 2 function value: -575.405841 derivative value: -18896.411616 New theta: 49.193991
Iteration 3 function value: -0.865462 derivative value: -18839.581357 New theta: 49.193945
STOPPING...
Initial theta: 50.000000
Iteration 1 function value: 54227179.790000 derivative value: 4109488.600000 New theta: 36.804397
Iteration 2 function value: 17142871.895470 derivative value: 1735188.464870 New theta: 26.924852
Iteration 3 function value: 5415943.363325 derivative value: 733123.609588 New theta: 19.537362
Iteration 4 function value: 1709212.528601 derivative value: 310074.937754 New theta: 14.025106
Iteration 5 function value: 538433.401650 derivative value: 131377.481132 New theta: 9.926737
Iteration 6 function value: 169109.335502 derivative value: 55825.121585 New theta: 6.897468
Iteration 7 function value: 52854.995553 derivative value: 23831.839604 New theta: 4.679636
Iteration 8 function value: 16391.615741 derivative value: 10248.482299 New theta: 3.080218
Iteration 9 function value: 5022.091638 derivative value: 4456.453803 New theta: 1.953292
Iteration 10 function value: 1510.895250 derivative value: 1969.215468 New theta: 1.186035
Iteration 11 function value: 443.032197 derivative value: 888.938453 New theta: 0.687651
Iteration 12 function value: 125.726074 derivative value: 411.509114 New theta: 0.382127
Iteration 13 function value: 34.400731 derivative value: 195.464351 New theta: 0.206132
Iteration 14 function value: 9.076039 derivative value: 95.152923 New theta: 0.110748
Iteration 15 function value: 2.298882 derivative value: 47.746659 New theta: 0.062601
Iteration 16 function value: 0.537595 derivative value: 25.613843 New theta: 0.041612
Iteration 17 function value: 0.097751 derivative value: 16.336365 New theta: 0.035629
Iteration 18 function value: 0.007800 derivative value: 13.732315 New theta: 0.035061
STOPPING...
The 98.0% confidence interval for the odds ratio is [0.035, 49.194].
```
