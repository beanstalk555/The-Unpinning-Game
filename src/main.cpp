#include "app.h"
#include "permutation.h"

int main()
{
#ifdef RUN_PERMUTATION_TESTS
    RunPermutationTests();
#endif
    RunApp();
}
