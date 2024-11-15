#include "INCLUDE/LK.h"

/*
   This file contains the main function of the program.
*/

/* Declarations of all global variables */

long *BestTour, Dimension, MaxCandidates, AscentCandidates, InitialPeriod,
        InitialStepSize, Precision, Runs, MaxTrials, MaxSwaps;
double BestCost, WorstCost, Excess, Optimum, SwNeighborhoodStartTime, SwNeighborhoodLimitTimeExceed, PrintGapTime;
unsigned int Seed;
int Subgradient, TraceLevel, MoveType, BacktrackMoveType, RestrictedSearch, SW, Lambda, TrialTimeBudget, MaxDepth,
        GainCriterionUsed, idMove, Gain23Used, SamplingBiasUsed, ExitSwNeighborhood, IntensificationDiversificationBest, IntensificationDiversificationAny,
        coNodes, maxCoNodes;

ImprovingMove *BestImprovingMoves;
Node *NodeSet, *FirstNode, *FirstActive, *LastActive, **Heap;
SwapRecord *SwapStack;
SwapMove *SwapTree;
SwapMove *BestSwapTree;
SwapMove **NegativeSwapList;
SwapMove **Queue;

Neighborhood *NegativeNeighbors;
long Swaps, Norm, M, GroupSize, Groups, Trial, *BetterTour, *CacheVal,
        *CacheSig, *CostMatrix;
double BetterCost, LowerBound;
unsigned long Hash;
int *Rand, Reversed;
Segment *FirstSegment;
HashTable *HTable;

FILE *ParameterFile, *ProblemFile, *PiFile, *TourFile, *OutputFile,
        *InputTourFile, *CandidateFile, *InitialTourFile, *MergeTourFile[2];
char *ParameterFileName, *ProblemFileName, *PiFileName, *TourFileName, *OutputFileName,
        *InputTourFileName, *CandidateFileName, *InitialTourFileName,
        *MergeTourFileName[2];
char *Name, *Type, *EdgeWeightType, *EdgeWeightFormat, *EdgeDataFormat,
        *NodeCoordType, *DisplayDataType;
int ProblemType = -1, WeightType = -1, WeightFormat = -1, CoordType =
        NO_COORDS, CandidateSetSymmetric = 0;

long (*Distance)(Node *Na, Node *Nb);

long (*C)(Node *Na, Node *Nb);

long (*D)(Node *Na, Node *Nb);

long (*c)(Node *Na, Node *Nb);

Node *(*BestMove)(Node *t1, Node *t2, long *G0, long *Gain);

Node *(*BacktrackMove)(Node *t1, Node *t2, long *G0, long *Gain);

/* 
   The main function: 
*/

int main(int argc, char *argv[]) {
    long TrialSum, MinTrial, Successes, Run;
    double Cost, CostSum, Time, TimeSum, MinTime;
    double LastTime = GetTime();

    TrialSum = Successes = 0;
    CostSum = TimeSum = 0.0;
    MinTrial = LONG_MAX;
    MinTime = DBL_MAX;

    /* Read the specification of the problem */
    if (argc >= 2)
        ParameterFileName = argv[1];
    ReadParameters();
    ReadProblem();
    CreateCandidateSet();
    printf("Preprocessing time = %0.0f sec.\n\n", GetTime() - LastTime);
    fflush(stdout);
    if (Norm != 0) {
        BestCost = DBL_MAX;
        WorstCost = -DBL_MAX;
        Successes = 0;
    } else {
        /* The ascent has solved the problem! */
        Successes = 1;
        Runs = 0;
        RecordBetterTour();
        RecordBestTour();
        BestCost = WorstCost = Cost = CostSum = LowerBound;
        PrintBestTour();
    }
    TimeSum = 0;
    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        if (Run != 1) {
            if (InitialTourFileName){
                // Trouver la position de la partie dynamique de la chaîne de caractères
                char *dynamic_part = strstr(InitialTourFileName, "_input_");
                int length_to_end = strlen(dynamic_part) - strlen("_input_") + 1;
                sprintf(dynamic_part, "_input_%d.inputtourfile", Run);
                ReadTour(InitialTourFileName, &InitialTourFile);
            }
        }
        Cost = FindTour();      /* using the Lin-Kerninghan heuristics */
        if (Cost < BestCost) {
            RecordBestTour();
            BestCost = Cost;
            PrintBestTour();
        }
        /* Update statistics */
        if (Cost > WorstCost)
            WorstCost = Cost;
        if (Cost <= Optimum)
            Successes++;
        Time = GetTime() - LastTime;
        if (TraceLevel >= 1) {
            printf("#RUN %d : Cost = %0.0f, Seed = %u, Time = %0.3f sec.\n\n",
                   Run, Cost, Seed, Time);
            fflush(stdout);
        }
        CostSum += Cost;
        TrialSum += Trial;
        if (Trial < MinTrial)
            MinTrial = Trial;
        TimeSum += Time;
        if (Time < MinTime)
            MinTime = Time;
        srand(++Seed);
        if (Cost < Optimum || (Cost == Optimum && Successes == 1)) {
            if (Cost < Optimum) {
                Node *N;
                N = FirstNode;
                while ((N = N->OptimumSuc = N->Suc) != FirstNode);
                printf("New optimum = %f, Old optimum = %f\n", Cost,
                       Optimum);
                fflush(stdout);
                Optimum = Cost;
            }
            PrintBestTour();
        }
    }
    /* Report the resuls */
    printf("\nLAMBDA = %d MAX_CANDIDATES = %ld", Lambda, MaxCandidates);
    //printf("\nMax Nodes Generated = %d", maxCoNodes);
    printf("\nSuccesses/Runs = %ld/%ld \n", Successes, Runs);
    if (Runs == 0) {
        Runs = 1;
        MinTrial = 0;
        MinTime = 0;
    }
    printf("Cost.min = %0.0f, Cost.avg = %0.1f, Cost.max = %0.0f\n",
           BestCost, CostSum / Runs, WorstCost);
    if (Optimum == -DBL_MAX)
        Optimum = BestCost;
    printf("Gap.min = %0.3f%%, Gap.avg = %0.3f%%, Gap.max = %0.3f%%\n",
           (BestCost - Optimum) / Optimum * 100.0,
           (CostSum / Runs - Optimum) / Optimum * 100.0,
           (WorstCost - Optimum) / Optimum * 100.0);
    printf("MinTrials = %ld, Trials.avg. = %0.1f\n", MinTrial,
           1.0 * TrialSum / Runs);
    printf("Time.min = %0.3f sec., Time.avg. = %0.3f sec.\n\n", MinTime,
           TimeSum / Runs);
    fflush(stdout);
    return 0;
}
