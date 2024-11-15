#include "INCLUDE/Segment.h"
#include "INCLUDE/LK.h"

void printExcludedEdges(ExcludedEdges *head);

void SwNeighborhood2Opt(Node *t1, Node *t2, long *G0, SwapMove *parent, int depth, ExcludedEdges *edges) {

    if (depth == 0) {
        printf("Limit depth ... \n");
        return;
    }

    if (SUC(t1) != t2)
        Reversed ^= 1;

    if (GetTime() - SwNeighborhoodStartTime > SwNeighborhoodLimitTimeExceed) {
        ExitSwNeighborhood = 1;
        return;
    }

    Node *t3, *t4, *T1, *T2, *T3;
    SwapMove *CurrentMove;
    Candidate *Nt2;
    long G1, G2, MoveGain = LONG_MIN;
    //SRandom(++Seed);

    int lambdaCandidates = Lambda;
    int Count = 0;
    for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
        Count++;
    }

    // Lambda is bigger than the number of reel candidates
    if (Count < Lambda) {
        lambdaCandidates = Count;
    }

    int PossibleIndexes[Count];
    for (int i = 0; i < Count; i++) {
        PossibleIndexes[i] = i;
    }

    for (int i = 0; i < lambdaCandidates; ++i) {
        int randomIndex;
        if (SamplingBiasUsed) {
            // Sampling bias : Giving higher sampling probability to candidates order
            double r = GenerateRandomNonZero();
            randomIndex = (int) trunc(-log((double) r) / log(2));
        } else {
            int r = (int) Random();
            randomIndex = (int) r % Count;
        }
        if (randomIndex > Count - i - 1) {
            randomIndex = Count - i - 1;
        }
        int reelRandomIndex = PossibleIndexes[randomIndex];
        int currentIndex = 0;

        for (Nt2 = t2->CandidateSet; t3 = Nt2->To; Nt2++) {
            if (t3 == t2->Pred || t3 == t2->Suc ||
                ((G1 = *G0 - Nt2->Cost) <= 0 && GainCriterionUsed && ProblemType != HCP && ProblemType != HPP)) {
                if (currentIndex >= reelRandomIndex)
                    break;
                else {
                    currentIndex++;
                    continue;
                }
            }
            if (currentIndex != reelRandomIndex)
                currentIndex++;
            else {
                /* Choose t4 (only one choice gives a closed tour) */
                t4 = PRED(t3);
                if (Fixed(t3, t4) || t4 == t1)
                    break;
                else if ((NotExcludable(edges, t1->Id, t2->Id)) ||
                           (NotExcludable(edges, t3->Id, t4->Id)))
                    break;

                // save the tree node
                G2 = G1 + C(t3, t4);
                MoveGain = G2 - C(t4, t1);

                if (parent->firstChild->IsActive == 0) {
                    CurrentMove = parent->firstChild;
                    SetSwapMove(&(parent->firstChild), t1 , t2, t3, t4, MoveGain, G2);
                    coNodes++;
                } else {
                    SwapMove *lastChild = parent->firstChild;
                    while (lastChild->nextSibling->IsActive == 1) {
                        lastChild = lastChild->nextSibling;
                    }
                    CurrentMove = lastChild->nextSibling;
                    SetSwapMove(&(lastChild->nextSibling), t1 , t2, t3, t4, MoveGain, G2);
                    coNodes++;
                }
                parent->NumChildren++;

                // Record the best tree node or return it if it's gain is positive
                if (MoveGain > 0) {
                    BestSwapTree = CurrentMove;
                    return;
                } else {
                    if (BestSwapTree->Gain < MoveGain) {
                        BestSwapTree = CurrentMove;
                    }
                }

                /*if (t4->Cost > G2 - Precision){
                    break;
                }*/

                // exclude
                ExcludeEdge(&edges, t1->Id, t2->Id);
                ExcludeEdge(&edges, t3->Id, t4->Id);
                // apply the movement
                Swap1(t1, t2, t3);

                // recursive call
                SwNeighborhood2Opt(t1, t4, &G2, CurrentMove, depth - 1, edges);

                if(ExitSwNeighborhood || BestSwapTree->Gain > 0){
                    RestoreMovement(&T1, &T2, &T3, t1, t2);
                    return;
                }

                DisExcludLastEdges(&edges);
                RestoreMovement(&T1, &T2, &T3, t1, t2);
                break;
            }
        }
        ShiftToEnd(PossibleIndexes, Count, randomIndex);
    }
}

void RestoreMovement(Node **T1, Node **T2, Node **T3, Node *t1, Node *t2) {
    Swaps--;
    // get modified nodes of the last swap
    *T1 = SwapStack[Swaps].t1;
    *T2 = SwapStack[Swaps].t2;
    *T3 = SwapStack[Swaps].t3;
    Swap1(*T3, *T2, *T1);
    Swaps--;
    if (t2 != SUC(t1))
        Reversed ^= 1;
}

void ShiftToEnd(int *tab, int n, int index) {
    if (index < 0 || index >= n) {
        printf("Index invalide\n");
        return;
    }
    int selected = tab[index];

    for (int i = index; i < n - 1; i++) {
        tab[i] = tab[i + 1];
    }

    tab[n - 1] = selected;
}

void ExcludeEdge(ExcludedEdges **head, int first, int second) {
    ExcludedEdges *newEdg = (ExcludedEdges *) malloc(sizeof(ExcludedEdges));
    newEdg->firstNode = first;
    newEdg->secondNode = second;
    newEdg->next = *head;
    *head = newEdg;
}

int NotExcludable(ExcludedEdges *head, long first, long second) {
    ExcludedEdges *current = head;
    while (current != NULL) {
        if (current->firstNode == first && current->secondNode == second) {
            return 1;
        }
        current = current->next;
    }
    return 0;
}

void DisExcludLastEdges(ExcludedEdges **head) {
    if (*head == NULL || (*head)->next == NULL) {
        return;
    }

    ExcludedEdges *temp = *head;
    *head = (*head)->next->next;

    free(temp->next);
    free(temp);
}

void printExcludedEdges(ExcludedEdges *head) {
    ExcludedEdges *current = head;
    while (current != NULL) {
        printf("(%d, %d) -> ", current->firstNode, current->secondNode);
        current = current->next;
    }
    printf("NULL\n");
}

double GenerateRandomNonZero() {
    double r;
    do {
        r = (double) Random() / RAND_MAX; // Generate a random number between 0 and 1
    } while (r == 0.0);
    return r;
}

void SetSwapMove(SwapMove **move, Node *t1, Node *t2, Node *t3, Node *t4, long gain, long g2Gain){
    (*move)->t1 = t1;
    (*move)->t2 = t2;
    (*move)->t3 = t3;
    (*move)->t4 = t4;
    (*move)->Gain = gain;
    (*move)->G2Gain = g2Gain;
    (*move)->IsActive = 1;
}
