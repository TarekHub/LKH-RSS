#include <time.h>
#include "INCLUDE/Segment.h"
#include "INCLUDE/LK.h"

void FreeSwapMove(SwapMove **root);

Neighborhood *CreateTree(SwapMove *tree);

void ApplyMovement(SwapMove *swapMove);

void RecordBestNegative(SwapMove **bestSwapMove);

void AddTreeToNeighborhood(Neighborhood **head, SwapMove *tree);

void FreeInverseSwapMove(SwapMove **move);

SwapMove *GetRandomSwapMove(Neighborhood *listHead);

void FreeNeighborhood(Neighborhood **hood);

int CountSwapMoveTreeNodes(SwapMove *root);

void ResetIsActive(SwapMove **root);

SwapMove **StoreNegativeMove(SwapMove *move, int *size);

void FreeNegativeMove(SwapMove **pathArray, int size);

double LinKernighanSW() {
    Node *t1, *t2, *SUCt1;
    long Gain, G0, i, BestNegativeGain;
    double Cost, minimumCost;
    SwapMove *pathArray[MaxDepth];
    int pathNegativeMoveSize;
    Candidate *Nt1;
    Segment *S;
    int X2, it = 0;
    double LastTime = GetTime();
    double startTime = GetTime();
    int positiveGainAppliyed;

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
    if (HashSearch(HTable, Hash, Cost))
        return Cost / Precision;

    SwapMove *initialBestSwapTree = BestSwapTree;
    //printf("\n## Start PNLS ##\n");
    int nbrNegativeBestMove;
    // Début de la recherche locale
    while (GetTime() - LastTime < TrialTimeBudget) {
        //idMove = -1;
        positiveGainAppliyed = 0;
        nbrNegativeBestMove = 0;
        BestNegativeGain = LONG_MIN;

        /*NegativeSwapList->id = -1;
        NegativeSwapList->Gain = LONG_MIN;
        NegativeSwapList->parent = NULL;
*/
        // Pour tous les nœuds du tour
        while ((t1 = RemoveFirstActive())) {

            if (TraceLevel >= 3 && fabs(GetTime() - startTime) > 1) {
                printf("Cost = %0.0f", minimumCost / Precision);
                if (Optimum != -DBL_MAX && Optimum != 0)
                    printf(", Gap = %0.2f%%",
                           100.0 * (minimumCost / Precision - Optimum) / Optimum);
                printf(", Time = %0.2f sec.\n",
                       it, fabs(GetTime() - LastTime));

                startTime = GetTime();
            }

            SUCt1 = SUC(t1);
            for (X2 = 1; X2 <= 2; X2++) {
                t2 = X2 == 1 ? PRED(t1) : SUCt1;
                if ((RestrictedSearch && Near(t1, t2)) || Fixed(t1, t2))
                    continue;

                ExcludedEdges *edges = NULL;
                G0 = C(t1, t2);

                /// Recherche échantillonnée [SW sur les candidats]
                ExitSwNeighborhood = 0;
                SwNeighborhoodStartTime = GetTime();
                SwapTree->t1 = t1;
                SwapTree->t2 = t2;
                SwapTree->IsActive = 1;
                SwapTree->G2Gain = G0;
                coNodes = 0;
                SwapMove *result = SwNeighborhoodBreadth(SwapTree);
                //SwNeighborhood2Opt(t1, t2, &G0, SwapTree, MaxDepth, edges);
                if (coNodes > maxCoNodes) {
                    //printf("Nodes Generated = %d\n", coNodes);
                    maxCoNodes = coNodes;
                }

                // À ce point, nous devrions obtenir le meilleur gain, mais avant cela :
                // 1- Si aucun mouvement ne peut être effectué : nous libérons l'arbre
                if (SwapTree->NumChildren == 0) {
                    ResetIsActive(&SwapTree);
                    BestSwapTree = initialBestSwapTree;
                    continue;
                }

                // 2- Si le meilleur gain est positif : on applique le mouvement
                if (BestSwapTree->Gain > 0) {
                    //ApplyMovement(BestSwapTree);
                    Swap1(BestSwapTree->t1, BestSwapTree->t2, BestSwapTree->t3);
                    positiveGainAppliyed = 1;
                    Cost -= BestSwapTree->Gain;
                    if (minimumCost > Cost)
                        minimumCost = Cost;
                    StoreTour(1);
                }

                    // 3- Si le meilleur gain est négatif : on l'enregistre
                /*else if (BestSwapTree->Gain < 0) {

                    // si on emploie ID_Best
                    if (IntensificationDiversificationBest) {
                        if (BestSwapTree->Gain > BestNegativeGain) {
                            if (BestNegativeGain != LONG_MIN)
                                FreeNegativeMove(NegativeSwapList, pathNegativeMoveSize);

                            pathNegativeMoveSize = 0;
                            NegativeSwapList = StoreNegativeMove(BestSwapTree, &pathNegativeMoveSize);
                            BestNegativeGain = (*NegativeSwapList)->Gain;
                        }
                    }

                    // si on emploie ID_Any
                    else if (IntensificationDiversificationAny) {
                        nbrNegativeBestMove++;
                        if ((double) Random() / RAND_MAX <= 1.0 / nbrNegativeBestMove) {
                            if (nbrNegativeBestMove > 1)
                                free(NegativeSwapList);

                            pathNegativeMoveSize = 0;
                            NegativeSwapList = StoreNegativeMove(BestSwapTree, &pathNegativeMoveSize);
                            if((*NegativeSwapList)->id == 167)
                                printf("HELLO!");
                        }
                    }
                }*/

                // Réinitialiser l'arbre et 'BestSwapTree'
                ResetIsActive(&SwapTree);
                BestSwapTree = initialBestSwapTree;
                if (HashSearch(HTable, Hash, Cost))
                    goto End_LinKernighan;
            }
        }

        // Si aucun gain positif n'a été appliqué pour l'ensemble des nœuds du tour
        if (positiveGainAppliyed == 0) {
            // Si le mouvement non séquentiel apporte de gain positif, on l'applique
            if (Gain23Used && (Gain = Gain23()) > 0) {
                Cost -= (double) Gain;
                if (minimumCost > Cost)
                    minimumCost = Cost;

                // Sinon, appliquer un mouvement aléatoire détériorant enregistré précédemment (Id Any or Best)
                //else {
                //ApplyMoves(NegativeSwapList, pathNegativeMoveSize - 1);
                //Cost -= (double) (*NegativeSwapList)->Gain;
                //}

                StoreTour(1);
                if (HashSearch(HTable, Hash, Cost))
                    goto End_LinKernighan;
            }else { goto ForceExit; }
        }

        //free(NegativeSwapList);

        // Réactiver tous les nœuds

        /*t1 = FirstNode;
        do {
            Activate(t1);
        } while ((t1 = t1->Suc) != FirstNode);*/
    }

    ForceExit:
    //printf("## End PNLS ## \n");
    if (!HashSearch(HTable, Hash, Cost))
        HashInsert(HTable, Hash, Cost);

    End_LinKernighan:
    NormalizeNodeList();
    return minimumCost / Precision;
}

SwapMove **StoreNegativeMove(SwapMove *move, int *size) {
    SwapMove **pathArray = NULL;
    int pathSize = 0;
    SwapMove *current = move;

    while (current != NULL) {
        pathSize++;
        current = current->parent;
    }

    pathArray = (SwapMove **) malloc(pathSize * sizeof(SwapMove *));

    current = move;
    for (int i = 0; i < pathSize; i++) {
        pathArray[i] = (SwapMove *) malloc(sizeof(SwapMove));
        memcpy(pathArray[i], current, sizeof(SwapMove));
        current = current->parent;
    }

    *size = pathSize;
    return pathArray;
}

void FreeNegativeMove(SwapMove **pathArray, int size) {
    free(pathArray);
}


void ApplyMovement(SwapMove *swapMove) {
    int pathLength;
    SwapMove **pathToRoot = FindPathToRoot(swapMove, &pathLength);
    ApplyMoves(pathToRoot, pathLength - 1);
}

SwapMove *CreateSwapMoveAlloc(int id, Node *t1, Node *t2, Node *t3, Node *t4, long gain, SwapMove *parent) {
    SwapMove *swapMove = (SwapMove *) malloc(sizeof(SwapMove));
    swapMove->id = id;
    swapMove->t1 = t1;
    swapMove->t2 = t2;
    swapMove->t3 = t3;
    swapMove->t4 = t4;
    swapMove->parent = parent;
    swapMove->Gain = gain;
    swapMove->NumChildren = 0;
    swapMove->firstChild = NULL;
    swapMove->nextSibling = NULL;
    swapMove->IsBest = 0;

    return swapMove;
}

SwapMove **FindPathToRoot(SwapMove *move, int *pathLength) {
    if (move == NULL) {
        return NULL;
    }

    SwapMove **path = NULL;
    *pathLength = 0;

    SwapMove *current = move;
    while (current != NULL) {
        (*pathLength)++;
        current = current->parent;
    }

    path = (SwapMove **) malloc(*pathLength * sizeof(SwapMove *));

    current = move;
    for (int i = 0; i < *pathLength; i++) {
        path[i] = current;
        current = current->parent;
    }

    return path;
}

void ApplyMoves(SwapMove **pathToRoot, int pathLength) {
    for (int i = pathLength - 1; i >= 0; i--) {
        SwapMove *currentMove = pathToRoot[i];
        // apply swap
        Swap1(currentMove->t1, currentMove->t2, currentMove->t3);
    }
}

void FreeSwapMove(SwapMove **root) {
    if ((root == NULL) || (*root == NULL)) return;

    FreeSwapMove(&(*root)->firstChild);
    FreeSwapMove(&(*root)->nextSibling);

    // Free the root itself
    if ((*root)->IsBest == 0) {
        free(*root);
        *root = NULL;
    }
}

void FreeInverseSwapMove(SwapMove **move) {

    while (*move != NULL) {
        SwapMove *parent = (*move)->parent;
        free(*move);
        *move = NULL;
        *move = parent;
    }
}

void RecordBestNegative(SwapMove **bestSwapMove) {
    (*bestSwapMove)->IsBest = 1;
    SwapMove *parent = (*bestSwapMove)->parent;
    while (parent != NULL) {
        parent->IsBest = 1;
        parent = parent->parent;
    }
}

void AddTreeToNeighborhood(Neighborhood **head, SwapMove *tree) {
    Neighborhood *newNode = CreateTree(tree);
    if (*head == NULL) {
        *head = newNode;
    } else {
        Neighborhood *current = *head;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = newNode;
    }
}

Neighborhood *CreateTree(SwapMove *tree) {
    Neighborhood *newNode = (Neighborhood *) malloc(sizeof(Neighborhood));
    newNode->tree = tree;
    newNode->next = NULL;
    return newNode;
}

SwapMove *GetRandomSwapMove(Neighborhood *listHead) {
    // Step 1: Count the number of nodes in the list
    int listLength = 0;
    Neighborhood *current = listHead;
    while (current != NULL) {
        listLength++;
        current = current->next;
    }

    // Step 2: Generate a random index
    //srand(time(NULL)); // Seed the random number generator with the current time
    int randd = (int) Random();
    int randomIndex = randd % listLength;

    // Step 3: Traverse the list again and select the node at the random index
    current = listHead;
    for (int i = 0; i < randomIndex; i++) {
        current = current->next;
    }

    // Return the SwapMove object at the selected node
    return current->tree;
}

void FreeNeighborhood(Neighborhood **hood) {
    if (hood == NULL || *hood == NULL) return;

    // Free the tree
    FreeInverseSwapMove(&(*hood)->tree);

    // Free the next neighborhood recursively
    FreeNeighborhood(&(*hood)->next);

    // Free the neighborhood itself
    free(*hood);
    *hood = NULL;
}

int CountSwapMoveTreeNodes(SwapMove *root) {
    if (root == NULL) {
        return 0;
    }

    int count = 1; // Compte le nœud courant
    SwapMove *currentChild = root->firstChild;

    while (currentChild != NULL) {
        count += CountSwapMoveTreeNodes(currentChild);
        currentChild = currentChild->nextSibling;
    }

    return count;
}

void ResetIsActive(SwapMove **root) {
    if (root == NULL) {
        return;
    }

    (*root)->IsActive = 0; // Réinitialiser IsActive à 0
    (*root)->NumChildren = 0;
    SwapMove *currentChild = (*root)->firstChild;

    while (currentChild != NULL) {
        ResetIsActive(&currentChild);
        currentChild = currentChild->nextSibling;
    }
}


