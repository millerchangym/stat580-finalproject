#ifndef final_project_h
#define final_project_h
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAKE_2ARRAY(a, m, n) do {                                            \
	size_t ARRAY_H2RESERVED;                                             \
	MAKE_1ARRAY(a, (m)+1);                                               \
	if (a == NULL)                                                       \
		break;                                                       \
	(a)[m] = NULL;                                                       \
	for (ARRAY_H2RESERVED = 0; ARRAY_H2RESERVED < (size_t)(m);           \
	ARRAY_H2RESERVED++) {                                                \
	MAKE_1ARRAY((a)[ARRAY_H2RESERVED], (n));                             \
	if ((a)[ARRAY_H2RESERVED] == NULL) {                                 \
		FREE_2ARRAY(a);                                              \
		break;                                                       \
	}                                                                    \
}                                                                            \
} while (0)

#define FREE_2ARRAY(a) do {                                                           \
	size_t ARRAY_H2RESERVED;                                                      \
	if (a == NULL)                                                                \
		break;                                                                \
	for (ARRAY_H2RESERVED = 0; (a)[ARRAY_H2RESERVED] != NULL; ARRAY_H2RESERVED++) \
	FREE_1ARRAY((a)[ARRAY_H2RESERVED]);                                           \
	FREE_1ARRAY(a);                                                               \
} while (0)


#endif 
