#include "INCLUDE/Segment.h"
#include "INCLUDE/LK.h"


SwapMove *SwNeighborhoodBreadth(SwapMove *currentMove) {
    // Variables
    Node *t1, *t2;
    Node *t3, *t4, *T1, *T2, *T3;
    Candidate *Nt2;
    SwapMove *newSwapMove;
    long G1, G2, moveGain = LONG_MIN;
    int randomIndex, reelRandomIndex, lambdaCandidates;

    // Variables de la file
    int rear = 0, front = 0;

    t1 = currentMove->t1;
    t2 = currentMove->t2;

    Queue[rear++] = currentMove;

    while (rear != front) {
        currentMove = Queue[front++];

        if (SUC(currentMove->t1) != currentMove->t2)
            Reversed ^= 1;

        // if not root
        if (currentMove->parent != NULL) {
            RestoreTour();
            ReInitSwaps(currentMove);
            t1 = currentMove->t1;
            t2 = currentMove->t4;
            if (SUC(t1) != t2)
                Reversed ^= 1;
        }

        lambdaCandidates = Lambda;
        int Count = 0;
        for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++)
            Count++;

        // if lambda is bigger than the number of reel candidates
        if (Count < Lambda)
            lambdaCandidates = Count;

        int PossibleIndexes[Count];
        for (int i = 0; i < Count; i++)
            PossibleIndexes[i] = i;

        int l = 0;
        while (l < lambdaCandidates) {
            // Sampling bias : Giving higher sampling probability to candidates order
            if (!SamplingBiasUsed)
                reelRandomIndex = l;
            else {
                if (SamplingBiasUsed == 2)
                    randomIndex = (int) trunc(-log((double) GenerateRandomNonZero()) / log(2));
                else if (SamplingBiasUsed == 1)
                    randomIndex = (int) Random() % Count;

                if (randomIndex > Count - l - 1)
                    randomIndex = Count - l - 1;

                reelRandomIndex = PossibleIndexes[randomIndex];
            }

            /* Choose (t2,t3) as a candidate edge emanating from t2 */
            Nt2 = &(t2->CandidateSet[reelRandomIndex]);
            t3 = Nt2->To;

            // if the Candidat is not feasible
            if (t3 == t2->Pred || t3 == t2->Suc ||
                ((G1 = currentMove->G2Gain - Nt2->Cost) <= 0 && GainCriterionUsed && ProblemType != HCP &&
                 ProblemType != HPP)) {
                if (SamplingBiasUsed != 0){
                    ShiftToEnd(PossibleIndexes, Count, randomIndex);
                }
                l++;
                continue;
            } else {
                // Choose t4 (only one choice gives a closed tour)
                t4 = PRED(t3);
                // Feasible move !
                G2 = G1 + C(t3, t4);
                moveGain = G2 - C(t4, t1);

                // Depth reached
                if (currentMove->firstChild == NULL) {
                    RestoreTour();
                    return BestSwapTree;
                }

                if (currentMove->firstChild->IsActive == 0) {
                    newSwapMove = currentMove->firstChild;
                    Queue[rear++] = currentMove->firstChild;

                    SetSwapMove(&(currentMove->firstChild), t1, t2, t3, t4, moveGain, G2);
                    coNodes++;
                } else {
                    SwapMove *lastChild = currentMove->firstChild;
                    while (lastChild->nextSibling->IsActive == 1) {
                        lastChild = lastChild->nextSibling;
                    }
                    newSwapMove = lastChild->nextSibling;
                    Queue[rear++] = lastChild->nextSibling;

                    SetSwapMove(&(lastChild->nextSibling), t1, t2, t3, t4, moveGain, G2);
                    coNodes++;
                }
                currentMove->NumChildren++;

                // Record the best tree node or return it if the gain is positive
                if (moveGain > 0) {
                    BestSwapTree = newSwapMove;
                    return BestSwapTree;
                } else {
                    if (moveGain > BestSwapTree->Gain)
                        BestSwapTree = newSwapMove;
                }
            }
            if (SamplingBiasUsed != 0){
                ShiftToEnd(PossibleIndexes, Count, randomIndex);
            }
            l++;
        }
    }

    RestoreTour();
    return BestSwapTree;
}


void ReInitSwaps(SwapMove *swapMove) {
    int pathLength;
    SwapMove **pathToRoot = FindPathToRoot(swapMove, &pathLength);
    ApplyMoves(pathToRoot, pathLength - 1);
    free(pathToRoot);
}


