#include "INCLUDE/LK.h"

void DisExclude(Node * ta, Node * tb)
{
    if (ta == tb->OldPred)
        tb->OldPredExcluded = 0;
    else if (ta == tb->OldSuc)
        tb->OldSucExcluded = 0;
    if (tb == ta->OldPred)
        ta->OldPredExcluded = 0;
    else if (tb == ta->OldSuc)
        ta->OldSucExcluded = 0;
}