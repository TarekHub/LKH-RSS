#include "INCLUDE/LK.h"
#include "INCLUDE/Heap.h"

/*
   The MinimumSpanningTree function determines a minimum spanning tree using 
   Prim's algorithm.

   At return the Dad field of each node contains the father of the node, and the
   Cost field contains cost of the corresponding edge. The nodes are placed in a
   topological ordered list, i.e., for any node its father precedes the node in the
   list. The fields Pred and Suc of a node are pointers to the predecessor and
   successor node in this list.

   The function can be used to determine a minimum spanning tree in a dense graph,
   or in a sparse graph (a graph determined by a candidate set).

   When the graph is sparse a priority queue, implemented as a binary heap, is used 
   to speed up the determination of which edge to include next into the tree.
   The Rank field of a node is used to contain its priority (usually equal to the 
   smallest distance (Cost) to nodes of the tree).        
*/

void MinimumSpanningTree(const int Sparse)
{
    Node *Blue,                 /* Points to the last node included in the tree */
    *NextBlue,                  /* Points to the provisional next node to be included */
    *N;
    Candidate *NBlue;
    long Min, d;

    Blue = N = FirstNode;
    Blue->Dad = 0;              /* The root of the tree has no father */
    Blue->Loc = 0;              /* A blue node is not in the heap */
    if (Sparse && Blue->CandidateSet) {
        /* The graph is sparse */
        /* Insert all nodes in the heap */
        while ((N = N->Suc) != FirstNode) {
            N->Dad = Blue;
            N->Cost = N->Rank = LONG_MAX;
            Insert(N);
        }
        /* Update all neighbors to the blue node */
        for (NBlue = Blue->CandidateSet; N = NBlue->To; NBlue++) {
            if (Fixed(N, Blue)) {
                N->Dad = Blue;
                N->Cost = NBlue->Cost + N->Pi + Blue->Pi;
                N->Rank = LONG_MIN;
                SiftUp(N);
            } else if (!Blue->FixedTo2 && !N->FixedTo2) {
                N->Dad = Blue;
                N->Cost = N->Rank = NBlue->Cost + N->Pi + Blue->Pi;
                SiftUp(N);
            }
        }
        /* Loop as long as there a more nodes to include in the tree */
        while (NextBlue = DeleteMin()) {
            Follow(NextBlue, Blue);
            Blue = NextBlue;
            /* Update all neighbors to the blue node */
            for (NBlue = Blue->CandidateSet; N = NBlue->To; NBlue++) {
                if (!N->Loc)
                    continue;
                if (Fixed(N, Blue)) {
                    N->Dad = Blue;
                    N->Cost = NBlue->Cost + N->Pi + Blue->Pi;
                    N->Rank = LONG_MIN;
                    SiftUp(N);
                } else if (!Blue->FixedTo2 && !N->FixedTo2 &&
                           (d =
                            NBlue->Cost + N->Pi + Blue->Pi) < N->Cost) {
                    N->Dad = Blue;
                    N->Cost = N->Rank = d;
                    SiftUp(N);
                }
            }
        }
    } else {
        /* The graph is dense */
        while ((N = N->Suc) != FirstNode)
            N->Cost = LONG_MAX;
        /* Loop as long as there a more nodes to include in the tree */
        while ((N = Blue->Suc) != FirstNode) {
            Min = LONG_MAX;
            /* Update all non-blue nodes (the successors of Blue in the list) */
            do {
                if (Fixed(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = D(N, Blue);
                    NextBlue = N;
                    Min = LONG_MIN;
                } else {
                    if (!Blue->FixedTo2 && !N->FixedTo2 &&
                        !Forbidden(N, Blue) &&
                        (!c || c(Blue, N) < N->Cost) &&
                        (d = D(Blue, N)) < N->Cost) {
                        N->Cost = d;
                        N->Dad = Blue;
                    }
                    if (N->Cost < Min) {
                        Min = N->Cost;
                        NextBlue = N;
                    }
                }
            }
            while ((N = N->Suc) != FirstNode);
            Follow(NextBlue, Blue);
            Blue = NextBlue;
        }
    }
}
