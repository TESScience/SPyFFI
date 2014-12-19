// From twister.c

//
// uint32 must be an unsigned integer type capable of holding at least 32
// bits; exactly 32 should be fastest, but 64 is better on an Alpha with
// GCC at -O3 optimization so try your options and see what's best for you
//

typedef unsigned long uint32;
extern void seedMT(uint32 seed);
extern uint32 reloadMT(void);
extern uint32 randomMT(void);
extern double ranMT(void);
