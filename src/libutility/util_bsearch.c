#include "util_bsearch.h"
#include <stdlib.h>
#include <stdio.h>


extern bool
util_bsearch(const void *key, const void *base,
             size_t nmemb, size_t size,
             int (*compar)(const void *, const void *),
             void **result)
{
	size_t pos;
	size_t step;
	bool done = false;
	int comparResult;

	/* Setting first try position */
	pos = nmemb / 2;
	step = nmemb / 4;
	if (step == 0)
		step = 1;

	/* See if there is anything to look at */
	if (nmemb == (size_t)0) {
		if (result != NULL)
			*result = (void *)(((char *)base)+0);
		return false;
	}

	/* Do a binary search to find the matching object or the object
	 * closest to it position
	 */
	while (!done) {
		comparResult = compar(key, (void *)((char *)base+(pos*size)));
		switch (comparResult) {
			case 0:
				*result = (void *)((char *)base+(pos*size));
				return true;
				break;
			case -1:
				if (pos == 0) {
					/* Closet position found */
					done = true;
				} else {
					pos -= step;
					if (step > 1) {
						step /= 2;
					} else {
						comparResult =
						   compar(key, (void *)((char *)base+(pos*size)));
						if (comparResult == 1) {
							pos++;
							done = true;
						}
					}
				}
				break;
			case 1:
				if ((pos == nmemb-1) || (step == 0)) {
					/* Closet position found */
					done = true;
				} else {
					pos += step;
					if (step > 1) {
						step /=2;
					} else {
						comparResult =
						   compar(key, (void *)((char *)base+(pos*size)));
						if (comparResult == -1) {
							done = true;
						}
					}
				}
				break;
			default:
				fprintf(stderr,
				        "Wrong compare function in %s!\n",
				        __func__);
				exit(EXIT_FAILURE);
		}
	} /* End of while loop */

	/* The key was not found, but this is the closet we got */
	*result = (void *)((char *)base+(pos*size));

	/* The object was not found */
	return false;
}
