#ifndef __TTQ__HH__
#define __TTQ__HH__

#include <cstdio>

#ifdef DEBUG
#define ttqDebug(fmt, ...) printf("DEBUG %s:%d -> "fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#else
#define ttqDebug(_msg) (_msg)
#endif

#endif
