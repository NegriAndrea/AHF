#ifndef UTIL_BSEARCH_H
#define UTIL_BSEARCH_H

#include <stdlib.h>
#include <stdbool.h>


extern bool
util_bsearch(const void *key, const void *base,
             size_t nmemb, size_t size,
             int (*compar)(const void *, const void *),
             void **result);


#endif /* UTIL_BSEARCH_H */
