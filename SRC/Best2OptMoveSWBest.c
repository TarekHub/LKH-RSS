#include "INCLUDE/Segment.h"
#include "INCLUDE/LK.h"

/*
   The Best2OptMove function makes sequential edge exchanges. If possible, it makes a
   2-opt move that improves the tour. Otherwise, it makes the most promising
   2-opt move that fulfils the positive gain criterion. To prevent an infinity chain
   of moves the last edge in a 2-opt move must not previously have been included in
   the chain.

   The edge (t1,t2) is the first edge to be exchanged. G0 is a pointer to the
   accumulated gain.

   In case a 2-opt move is found that improves the tour, the improvement of the cost
   is made available to the caller through the parameter Gain. If *Gain > 0, an
   improvement of the current tour has been found.

   Otherwise, the best 2-opt move is made, and a pointer to the node that was connected
   to t1 (in order to close the tour) is returned. The new accumulated gain is made
   available to the caller through the parameter G0.

   If no move can be made, the function returns 0.

   The function is called from the LinKernighan function.
*/
void SaveImprovingMove(Node *t1, Node *t2, Node *t3, Node *t4, long PositiveGain);
Node* ApplyBestSwap();

Node *Best2OptMoveSWBest(Node *t1, Node *t2, long *G0, long *Gain) {
    //printf("\nSwaps number : %ld", Swaps);

    Node *t3, *t4, *T3, *T4 = 0;
    Candidate *Nt2, *NNa;
    long G1, G2, BestG2 = LONG_MIN;
    BestImprovingMoves->Gain = 0;
    int randomIndex, reelRandomIndex, lambdaCandidates;

    if (SUC(t1) != t2)
        Reversed ^= 1;

    lambdaCandidates = Lambda;
    int Count = 0;
    int PossibleIndexes[MaxCandidates];

    for (NNa = t2->CandidateSet; NNa->To; NNa++) {
        PossibleIndexes[Count] = Count;
        Count++;
    }
    // test
    // if lambda is bigger than the number of reel candidates
    if (Count < Lambda)
        lambdaCandidates = Count;

    int l = 0;
    while (l < lambdaCandidates) {
        // Sampling bias : Giving higher sampling probability to candidates order
        if (SamplingBiasUsed == 2)
            randomIndex = (int) trunc(-log((double) GenerateRandomNonZero()) / log(2));
        else if(SamplingBiasUsed == 1)
            randomIndex = (int) Random() % Count;

        if (randomIndex > Count - l - 1)
            randomIndex = Count - l - 1;

        reelRandomIndex = PossibleIndexes[randomIndex];

        /* Choose (t2,t3) as a candidate edge emanating from t2 */
        Nt2 = &(t2->CandidateSet[reelRandomIndex]);
        t3 = Nt2->To;

        // if the Candidat is not feasible
        if (t3 == t2->Pred || t3 == t2->Suc ||
            ((G1 = *G0 - Nt2->Cost) <= 0 && GainCriterionUsed && ProblemType != HCP && ProblemType != HPP)) {
            ShiftToEnd(PossibleIndexes, Count, randomIndex);
            l++;
            continue;
        } else {
            t4 = PRED(t3);
            G2 = G1 + C(t3, t4);

            ShiftToEnd(PossibleIndexes, Count, randomIndex);
            l++;

            if (!Forbidden(t4, t1) &&
                (!c || G2 - c(t4, t1) > 0) && (*Gain = G2 - C(t4, t1)) > 0) {
                SaveImprovingMove(t1, t2, t3, t4, *Gain);
                continue;
            }
            if (GainCriterionUsed && G2 - Precision < t4->Cost) {
                continue;
            } else if (G2 > BestG2 &&
                       Swaps < MaxSwaps &&
                       Excludable(t3, t4) && !InOptimumTour(t3, t4)) {
                T3 = t3;
                T4 = t4;
                BestG2 = G2;
            }
        }
    }

    if(BestImprovingMoves->Gain != 0){
        *Gain = BestImprovingMoves->Gain;
        return ApplyBestSwap();
    }

    *Gain = 0;
    if (T4) {
        /* Make the best 2-opt move */
        Swap1(t1, t2, T3);
        Exclude(t1, t2);
        Exclude(T3, T4);
        *G0 = BestG2;
    }
    return T4;
}

void SaveImprovingMove(Node *t1, Node *t2, Node *t3, Node *t4, long PositiveGain){
    if(PositiveGain > BestImprovingMoves->Gain){
        BestImprovingMoves->t1 = t1;
        BestImprovingMoves->t2 = t2;
        BestImprovingMoves->t3 = t3;
        BestImprovingMoves->t4 = t4;
        BestImprovingMoves->Gain = PositiveGain;
    }
}

Node* ApplyBestSwap(){
    Swap1(BestImprovingMoves->t1, BestImprovingMoves->t2, BestImprovingMoves->t3);
    return BestImprovingMoves->t4;
}