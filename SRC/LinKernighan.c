#include <time.h>
#include "INCLUDE/Segment.h"
#include "INCLUDE/LK.h"

/*
   The LinKenighan function seeks to improve a tour by sequential and nonsequential
   edge exchanges.

   The function returns the cost of the resulting tour. 
*/
int ApplyAnyFeasibleMove(Node *t1, Node *t2, long G0, long *Gain);

int ApplyBestFeasibleMove(Node *t1, Node *t2, long G0, long *Gain);

Node *ChooseRandomNode(Node *firstNode);

double LinKernighan() {
    Node *t1, *t2, *SUCt1;
    long Gain, G0, i;
    double Cost, minimumCost;
    Candidate *Nt1;
    Segment *S;
    int X2, it = 0;
    Reversed = 0;
    S = FirstSegment;
    i = 0;
    do {
        S->Size = 0;
        S->Rank = ++i;
        S->Reversed = 0;
        S->First = S->Last = 0;
    } while ((S = S->Suc) != FirstSegment);
    i = 0;
    Hash = 0;
    Swaps = 0;
    FirstActive = LastActive = 0;

    /* Compute the cost of the initial tour, Cost.
       Compute the corresponding hash value, Hash.
       Initialize the segment list.
       Make all nodes "active" (so that they can be used as t1). */
    Cost = 0;
    t1 = FirstNode;
    do {
        t2 = t1->OldSuc = t1->Next = t1->Suc;
        t1->OldPred = t1->Pred;
        t1->Rank = ++i;
        Cost += C(t1, t2) - t1->Pi - t2->Pi;
        Hash ^= Rand[t1->Id] * Rand[t2->Id];
        t1->Cost = LONG_MAX;
        for (Nt1 = t1->CandidateSet; t2 = Nt1->To; Nt1++)
            if (t2 != t1->Pred && t2 != t1->Suc && Nt1->Cost < t1->Cost)
                t1->Cost = Nt1->Cost;
        t1->Parent = S;
        S->Size++;
        if (S->Size == 1)
            S->First = t1;
        S->Last = t1;
        if (S->Size == GroupSize)
            S = S->Suc;
        t1->OldPredExcluded = t1->OldSucExcluded = 0;
        t1->Next = 0;
        Activate(t1);
    } while ((t1 = t1->Suc) != FirstNode);

    minimumCost = Cost;
    printf("Initial Cost = %0.0f", Cost / Precision);
    printf(" (Gap = %0.2f%%)\n",
           100.0 * (minimumCost / Precision - Optimum) / Optimum);

    if (HashSearch(HTable, Hash, Cost))
        return Cost / Precision;

    double LastTime = GetTime();
    double LastImporvementTime = GetTime();
    /* Loop as long as improvements are found */
    /* Choose t1 as the first "active" node */

    while ((t1 = RemoveFirstActive()) && (GetTime() - LastTime < TrialTimeBudget)) {
        SUCt1 = SUC(t1);

        /* Choose t2 as one of t1's two neighbor nodes on the tour */
        for (X2 = 1; X2 <= 2; X2++) {
            t2 = X2 == 1 ? PRED(t1) : SUCt1;
            if ((RestrictedSearch && Near(t1, t2)) || Fixed(t1, t2))
                continue;
            G0 = C(t1, t2);
            /* Make sequential moves */
            while ((t2 = BacktrackMove ?
                        BacktrackMove(t1, t2, &G0, &Gain) :
                        BestMove(t1, t2, &G0, &Gain)) && (GetTime() - LastTime < TrialTimeBudget)) {
                if (Gain > 0) {
                    /* An improvement has been found */
                    Cost -= Gain;
                    if (minimumCost > Cost){
                      minimumCost = Cost;
                      if (TraceLevel >= 3 || (TraceLevel == 2 && Cost / Precision < BetterCost)) {
                        printf("Cost = %0.0f, Time = %0.3f sec.\n",
                               Cost / Precision, GetTime() - LastTime);
                        fflush(stdout);
                      }
                    }
                    StoreTour(1);
                    /* Make t1 "active" again */
                    Activate(t1);
                    goto Next_t1;
                } else {
                    Activate(t1);
                }
            }
            RestoreTour();
            Activate(t1);
        }
        Next_t1:;
    }

    printf("Cost = %0.0f", minimumCost / Precision);
    if (Optimum != -DBL_MAX && Optimum != 0)
        printf(", Gap = %0.2f%%",
               100.0 * (minimumCost / Precision - Optimum) / Optimum);
    printf(", Time = %0.0f sec.\n", fabs(GetTime() - LastTime));

    End_LinKernighan:
    NormalizeNodeList();
    return minimumCost / Precision;
}

int ApplyBestFeasibleMove(Node *t1, Node *t2, long G0, long *Gain) {
    *Gain = 0;

    Candidate *Nt2, *NNa;
    Node *t3, *t4 = 0;
    long G1, G2 = LONG_MIN;

    for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
        if (t3 == t2->Pred || t3 == t2->Suc)
            continue;
        else {
            G1 = G0 - Nt2->Cost;
            t4 = PRED(t3);
            G2 = G1 + C(t3, t4);
            *Gain = G2 - C(t4, t1);

            Swap1(t1, t2, t3);
            return 1;
        }
    }
    // no feasible candidat was found
    return 0;
}

int ApplyAnyFeasibleMove(Node *t1, Node *t2, long G0, long *Gain) {
    *Gain = 0;

    Candidate *Nt2, *NNa;
    Node *t3, *t4 = 0;
    int randomIndex, reelRandomIndex;
    long G1, G2 = LONG_MIN;
    int Count = 0;
    int PossibleIndexes[50];

    for (NNa = t2->CandidateSet; NNa->To; NNa++) {
        PossibleIndexes[Count] = Count;
        Count++;
    }

    for (int i = 0; i < Count; i++) {
        randomIndex = (int) Random() % Count;

        if (randomIndex > Count - i - 1)
            randomIndex = Count - i - 1;

        reelRandomIndex = PossibleIndexes[randomIndex];

        Nt2 = &(t2->CandidateSet[reelRandomIndex]);
        t3 = Nt2->To;

        // if the Candidat is not feasible
        if (t3 == t2->Pred || t3 == t2->Suc) {
            ShiftToEnd(PossibleIndexes, Count, randomIndex);
            i++;
            continue;
        }

        G1 = G0 - Nt2->Cost;
        t4 = PRED(t3);
        G2 = G1 + C(t3, t4);
        *Gain = G2 - C(t4, t1);

        Swap1(t1, t2, t3);
        return 1;
    }
    // no feasible candidat was found
    return 0;
}

Node *ChooseRandomNode(Node *firstNode) {
    Node *current = firstNode;
    Node *selected = NULL;
    int count = 0;

    do {
        count++;
        if (Random() % count == 0)
            selected = current;

        current = SUC(current->Suc);

    } while (current != firstNode);

    return selected;
}