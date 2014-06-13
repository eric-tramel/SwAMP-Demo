#include "../swamp.h"

/* Random order */
void sort_rand( int n, int *seq ) {
    int i, key, exchg;

    for (key = 0; key < n - 1; key++) {
        exchg = key + rand_int(n - key);
        i = seq[key]; seq[key] = seq[exchg]; seq[exchg] = i;
    }
}
