# Final Project for STAT 580 - Yeng Miller-Chang

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
-alpha (for which 1 - alpha is the confidence level).

Please see the file final_project_writeup.pdf for more information.
