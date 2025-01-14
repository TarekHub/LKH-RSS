#include "INCLUDE/LK.h"
#include "INCLUDE/Heap.h"

/*      
   The ReadProblem function reads the problem data in TSPLIB format from the file 
   specified in the parameter file (PROBLEM_FILE).

   The following description of the file format is extracted from the TSPLIB 
   documentation.  

   The file consists of a specification part and a data part. The specification part
   contains information on the file format and on its contents. The data part contains
   explicit data.

   (1) The Specification part

   All entries in this section are of the form <keyword> : <value>, where <keyword> 
   denotes an alphanumerical keyword and <value> denotes alphanumerical or numerical
   data. The terms <string>, <integer> and <real> denote character string, integer 
   or real data, respectively. The order of specification of the keywords in the data
   file is arbitrary (in principle), but must be consistent, i.e., whenever a keyword
   is specified, all necessary information for the correct interpretation of the
   keyword has to be known. Below is given a list of all available keywords.

   NAME : <string>e
   Identifies the data file.

   TYPE : <string>
   Specifies the type of data. Possible types are
   TSP          Data for a symmetric traveling salesman problem
   ATSP         Data for an asymmetric traveling salesman problem
   HCP          Hamiltonian cycle problem data.
   HPP          Hamiltonian path problem data (not available in TSPLIB)

   COMMENT : <string>
   Additional comments (usually the name of the contributor or the creator of the
   problem instance is given here).

   DIMENSION : < integer>
   The number of nodes.

   EDGE_WEIGHT_TYPE : <string>
   Specifies how the edge weights (or distances) are given. The values are:
   EXPLICIT     Weights are listed explicitly in the corresponding section
   EUC_2D       Weights are Euclidean distances in 2-D
   EUC_3D       Weights are Euclidean distances in 3-D
   MAX_2D       Weights are maximum distances in 2-D
   MAX_3D       Weights are maximum distances in 3-D
   MAN_2D       Weights are Manhattan distances in 2-D
   MAN_3D       Weights are Manhattan distances in 3-D
   CEIL_2D      Weights are Euclidean distances in 2-D rounded up
   CEIL_3D      Weights are Euclidean distances in 3-D rounded up
   GEO          Weights are geographical distances (TSPLIB)
   GEOM         Weights are geographical distances (used for the world TSP) 
   ATT          Special distance function for problem att48 and att532

   EDGE-WEIGHT_FORMAT : <string>
   Describes the format of the edge weights if they ar given explicitely. The values
   are:
   FUNCTION             Weights are given by a function (see above)
   FULL_MATRIX          Weights are given by a full matrix
   UPPER_ROW            Upper triangular matrix (row-wise without diagonal entries)
   LOWER_ROW            Lower triangular matrix (row-wise without diagonal entries)     
   UPPER_DIAG_ROW       Upper triangular matrix (row-wise including diagonal entries)
   LOWER_DIAG_ROW       Lower triangular matrix (row-wise including diagonal entries)
   UPPER_COL            Upper triangular matrix (column-wise without diagonal entries)
   LOWER_COL            Lower triangular matrix (column-wise without diagonal entries)  
   UPPER_DIAG_COL       Upper triangular matrix (column-wise including diagonal entries)
   LOWER_DIAG_COL       Lower triangular matrix (colum-wise including diagonal entries)

   EDGE_DATA_FORMAT : <string>
   Describes the format in which the edges of a graph are given, if the graph is
   not complete. The values are
   EDGE_LIST    The graph is given by an edge list
   ADJ_LIST     The graph is given by an adjacency list

   NODE_COORD_TYPE : <string>
   Specifies whether the coordinates are associated with each node (which, for
   example may be used for either graphical display or distance computations.
   The values are
   TWOD_COORDS          Nodes are specified by coordinates in 2-D
   THREE_COORDS         Nodes are specified by coordinates in 3-D
   NO_COORDS            The nodes do not have associated coordinates
   The default value is NO_COORDS. In the current implementation, however, the value 
   has no significance.

   DISPLAY_DATA_TYPE : <string>
   Specifies how a graphical display of the nodes can be obtained. The values are
   COORD_DISPLAY        Display is generated from the node coordinates
   TWOD_DISPLAY         Explicit coordinates in 2-D are given
   BO_DISPLAY           No graphical display is possible

   The default value is COORD_DISPLAY if node coordinates are specifies and 
   NO_DISPLAY otherwise. In the current implementation, however, the value has no 
   significance.

   EOF
   Terminates input data. The entry is optional.

   (2) The data part

   Depending on the choice of specifications some additional data may be required. 
   These data are given corresponding data sections following the specification part.
   Each data section begins with the corresponding keyword. The length of the section
   is either explicitly known form the format specification, or the section is
   terminated by an appropriate end-of-section identifier.

   NODE_COORD_SECTION :
   Node coordinates are given in this section. Each line is of the form
   <integer> <real> <real>
   if NODE_COORD_TYPE is TWOD_COORDS, or
   <integer> <real> <real> <real>
   if NODE_COORD_TYPE is THREED_COORDS. The integers give the number of the 
   respective nodes. The real numbers are the associated coordinates.

   EDGE_DATA_SECTION :
   Edges of the graph are specified in either of the two formats allowed in the
   EDGE_DATA_FORAT entry. If a type is EDGE_LIST, then the edges are given as a 
   sequence of lines of the form
   <integer> <integer>
   each entry giving the terminal nodes of some edge. The list is terminated by
   a -1. If the type is ADJ_LIST, the section consists of adjacency list for nodes.
   The adjacency list of a node x is specified as
   <integer> <integer> ... <integer> -1
   where the first integer gives the number of node x and the following integers
   (terminated by -1) the numbers of the nodes adjacent to x. The list of adjacency
   lists are terminated by an additional -1.

   FIXED_EDGES_SECTION :
   In this section, edges are listed that are required to appear in each solution
   to the problem. The edges to be fixed are given in the form (per line)
   <integer> <integer>
   meaning that the edge (arc) from the first node to the second node has to be
   contained in a solution. This section is terminated by a -1.

   DISPLAY_DATA_SECTION :
   If DISPLAY_DATA_TYPE is TWOD_DISPLAY, the 2-dimensional coordinates from which
   a display can be generated are given in the form (per line)
   <integer> <real> <real>
   The integers specify the respective nodes and the real numbers give the 
   associated coordinates. The contents of this section, however, has no 
   significance in the current implementation.

   TOUR_SECTION :
   A tour is specified in this section. The tour is given by a list of integers
   giving the sequence in which the nodes are visited in the tour. The tour is
   terminated by a -1. Note: In contrast to the TSPLIB format, only one tour can 
   be given in this section. The tour is used to limit the search (the last edge 
   to be removed in a non-gainful move must not belong to the tour). In addition, 
   the Alpha field of its edges is set to zero.

   EDGE_WEIGHT_SECTION :
   The edge weights are given in the format specifies by the EDGE_WEIGHT_FORMAT 
   entry. At present, all explicit data are integral and is given in one of the
   (self-explanatory) matrix formats, with explicitly known lengths.
*/

static const int MaxMatrixDimension = 2000;
static const char Delimiters[] = " :=\n\t\r\f\v";

static void CheckSpecificationPart();

static char *Copy(char *S);

static void CreateNodes();

static void Read_DIMENSION();

static void Read_DISPLAY_DATA_SECTION();

static void Read_DISPLAY_DATA_TYPE();

static void Read_EDGE_DATA_FORMAT();

static void Read_EDGE_DATA_SECTION();

static void Read_EDGE_WEIGHT_FORMAT();

static void Read_EDGE_WEIGHT_SECTION();

static void Read_EDGE_WEIGHT_TYPE();

static void Read_FIXED_EDGES_SECTION();

static void Read_NAME();

static void Read_NODE_COORD_SECTION();

static void Read_NODE_COORD_TYPE();

static void Read_TOUR_SECTION(FILE **File);

static void Read_TYPE();

void ReadProblem() {
    long i;
    Segment *S, *SPrev;
    char *Line, *Keyword;

    if (!(ProblemFile = fopen(ProblemFileName, "r")))
        eprintf("Cannot open %s", ProblemFileName);
    FreeStructures();
    WeightType = WeightFormat = -1;
    CoordType = NO_COORDS;
    Name = Type = EdgeWeightType = EdgeWeightFormat = 0;
    EdgeDataFormat = NodeCoordType = DisplayDataType = 0;
    Distance = 0;
    C = 0;
    c = 0;

    while (Line = ReadLine(ProblemFile)) {
        if (!(Keyword = strtok(Line, Delimiters)))
            continue;
        for (i = 0; i < strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "COMMENT"));
        else if (!strcmp(Keyword, "DEMAND_SECTION"))
            eprintf("Not implemented: %s", Keyword);
        else if (!strcmp(Keyword, "DEPOT_SECTION"))
            eprintf("Not implemented: %s", Keyword);
        else if (!strcmp(Keyword, "DIMENSION"))
            Read_DIMENSION();
        else if (!strcmp(Keyword, "DISPLAY_DATA_SECTION"))
            Read_DISPLAY_DATA_SECTION();
        else if (!strcmp(Keyword, "DISPLAY_DATA_TYPE"))
            Read_DISPLAY_DATA_TYPE();
        else if (!strcmp(Keyword, "EDGE_DATA_FORMAT"))
            Read_EDGE_DATA_FORMAT();
        else if (!strcmp(Keyword, "EDGE_DATA_SECTION"))
            Read_EDGE_DATA_SECTION();
        else if (!strcmp(Keyword, "EDGE_WEIGHT_FORMAT"))
            Read_EDGE_WEIGHT_FORMAT();
        else if (!strcmp(Keyword, "EDGE_WEIGHT_SECTION"))
            Read_EDGE_WEIGHT_SECTION();
        else if (!strcmp(Keyword, "EDGE_WEIGHT_TYPE"))
            Read_EDGE_WEIGHT_TYPE();
        else if (!strcmp(Keyword, "EOF"))
            break;
        else if (!strcmp(Keyword, "FIXED_EDGES_SECTION"))
            Read_FIXED_EDGES_SECTION();
        else if (!strcmp(Keyword, "NAME"))
            Read_NAME();
        else if (!strcmp(Keyword, "NODE_COORD_SECTION"))
            Read_NODE_COORD_SECTION();
        else if (!strcmp(Keyword, "NODE_COORD_TYPE"))
            Read_NODE_COORD_TYPE();
        else if (!strcmp(Keyword, "TOUR_SECTION"))
            Read_TOUR_SECTION(&ProblemFile);
        else if (!strcmp(Keyword, "TYPE"))
            Read_TYPE();
        else
            eprintf("Unknown Keyword: %s", Keyword);
    }
    assert(BestTour = (long *) calloc((Dimension + 1), sizeof(long)));
    assert(BetterTour = (long *) calloc((Dimension + 1), sizeof(long)));
    assert(SwapStack =
           (SwapRecord *) malloc((Dimension + 10) * sizeof(SwapRecord)));
    assert(HTable = (HashTable *) malloc(sizeof(HashTable)));
    assert(Rand = (int *) malloc((Dimension + 1) * sizeof(int)));

    assert(BestImprovingMoves  = (ImprovingMove *) malloc(sizeof(ImprovingMove)));

    if (SW) {

        assert(BestSwapTree = (SwapMove *) malloc(sizeof(SwapMove)));
        BestSwapTree->id = -1;
        BestSwapTree->Gain = LONG_MIN;
        BestSwapTree->parent = NULL;

        idMove = 0;
        printf("Creating Tree & Queue for SW moves ... \n");
        SwapTree = CreateSwapMoveTree(MaxDepth + 1, Lambda, 0);
        assert(Queue = (SwapMove **) malloc(sizeof(SwapMove *) * (pow(Lambda, MaxDepth))));
        printf("End\n");

        assert(NegativeSwapList = (SwapMove *) malloc(sizeof(SwapMove)));
    }


    if (Seed == 0)
        Seed = 1;
    srand(Seed);
    for (i = 1; i <= Dimension; i++)
        Rand[i] = rand();
    HashInitialize(HTable);
    srand(Seed);
    Swaps = 0;
    MakeHeap(Dimension);
    if (MaxCandidates < 0)
        MaxCandidates = 5;
    else if (MaxCandidates > Dimension)
        MaxCandidates = Dimension;
    if (AscentCandidates > Dimension)
        AscentCandidates = Dimension;
    if (InitialPeriod == 0) {
        InitialPeriod = Dimension / 2;
        if (InitialPeriod < 100)
            InitialPeriod = 100;
    }
    if (Precision == 0)
        Precision = 100;
    if (InitialStepSize == 0)
        InitialStepSize = 1;
    if (Excess == 0.0)
        Excess = 1.0 / Dimension;
    if (MaxTrials == 0)
        MaxTrials = Dimension;
    if (MaxSwaps < 0)
        MaxSwaps = Dimension;
    if (CostMatrix == 0 && Dimension <= MaxMatrixDimension && Distance != 0
        && Distance != Distance_1 && Distance != Distance_ATSP) {
        Node *Ni, *Nj;
        assert(CostMatrix =
                       (long *) calloc(Dimension * (Dimension - 1) / 2,
                                       sizeof(long)));
        Ni = FirstNode->Suc;
        do {
            Ni->C = &CostMatrix[(Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
            if (ProblemType != HPP || Ni->Id < Dimension)
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = Fixed(Ni, Nj) ? 0 : Distance(Ni, Nj);
            else
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = 0;
        } while ((Ni = Ni->Suc) != FirstNode);
        WeightType = EXPLICIT;
        c = 0;
    }
    C = WeightType == EXPLICIT ? C_EXPLICIT : C_FUNCTION;
    D = WeightType == EXPLICIT ? D_EXPLICIT : D_FUNCTION;
    if (WeightType != EXPLICIT && !FirstNode->CandidateSet) {
        for (i = 1; i <= Dimension; i *= 2);
        assert(CacheSig = (long *) calloc(i, sizeof(long)));
        assert(CacheVal = (long *) calloc(i, sizeof(long)));
    }
    if (MoveType == 0)
        MoveType = 5;
    if (InputTourFileName)
        ReadTour(InputTourFileName, &InputTourFile);
    if (InitialTourFileName)
        ReadTour(InitialTourFileName, &InitialTourFile);
    for (i = 0; i <= 1; i++)
        if (MergeTourFileName[i])
            ReadTour(MergeTourFileName[i], &MergeTourFile[i]);
    switch (MoveType) {
        case 2:
            BestMove = Best2OptMove;
            break;
        case 3:
            BestMove = Best3OptMove;
            break;
        case 4:
            BestMove = Best4OptMove;
            break;
        case 5:
            BestMove = Best5OptMove;
            break;
        case 6:
            BestMove = Best2OptMoveSW;
            break;
        case 7:
            BestMove = Best2OptMoveSWBest;
    }
    switch (BacktrackMoveType) {
        case 2:
            BacktrackMove = Backtrack2OptMove;
            break;
        case 3:
            BacktrackMove = Backtrack3OptMove;
            break;
        case 4:
            BacktrackMove = Backtrack4OptMove;
            break;
        case 5:
            BacktrackMove = Backtrack5OptMove;
            break;
    }
    GroupSize = sqrt(1.0 * Dimension);
    Groups = 0;
    for (i = Dimension, SPrev = 0; i > 0; i -= GroupSize, SPrev = S) {
        S = (Segment *) malloc(sizeof(Segment));
        S->Rank = ++Groups;
        if (!SPrev)
            FirstSegment = S;
        else
            Link(SPrev, S);
    }
    Link(S, FirstSegment);
    if (TraceLevel >= 1)
        PrintParameters();
    else
        printf("PROBLEM_FILE = %s\n",
               ProblemFileName ? ProblemFileName : "");
    fclose(ProblemFile);
}

void CheckSpecificationPart() {
    if (ProblemType == -1)
        eprintf("TYPE is missing");
    if (Dimension <= 0)
        eprintf("DIMENSION is not positive (or not specified)");
    if (WeightType == -1 && ProblemType != HCP)
        eprintf("EDGE_WEIGHT_TYPE is missing");
    if (WeightType == EXPLICIT && WeightFormat == -1 && !EdgeWeightFormat)
        eprintf("EDGE_WEIGHT_FORMAT is missing");
    if (WeightType == EXPLICIT && WeightFormat == FUNCTION)
        eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (WeightType != EXPLICIT
        && (WeightType != SPECIAL || CoordType != NO_COORDS)
        && WeightType != -1 && WeightFormat != -1
        && WeightFormat != FUNCTION)
        eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (ProblemType == ATSP && WeightType != EXPLICIT)
        eprintf("Conflicting TYPE and EDGE_WEIGHT_TYPE");
    if (ProblemType == ATSP && WeightFormat != FULL_MATRIX)
        eprintf("Conflicting TYPE and EDGE_WEIGHT_FORMAT");
}

static char *Copy(char *S) {
    char *Buffer;

    if (strlen(S) == 0)
        return 0;
    assert(Buffer = (char *) malloc(strlen(S) + 1));
    strcpy(Buffer, S);
    return Buffer;
}

static void CreateNodes() {
    Node *Prev, *N;
    long i;

    if (Dimension <= 0)
        eprintf("DIMENSION is not positive (or not specified)");
    if (ProblemType == ATSP)
        Dimension *= 2;
    else if (ProblemType == HPP) {
        Dimension++;
        if (Dimension > MaxMatrixDimension)
            eprintf("Dimension too large in HPP problem");
    }
    assert(NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node)));
    for (i = 1; i <= Dimension; i++, Prev = N) {
        N = &NodeSet[i];
        if (i == 1)
            FirstNode = N;
        else
            Link(Prev, N);
        N->Id = i;
    }
    Link(N, FirstNode);
}

static void Read_NAME() {
    if (!(Name = Copy(strtok(0, Delimiters))))
        eprintf("(NAME): string expected");
}

static void Read_DIMENSION() {
    if (!sscanf(strtok(0, Delimiters), "%ld", &Dimension))
        eprintf("(DIMENSION): integer expected");
}

static void Read_DISPLAY_DATA_SECTION() {
    Node *N;
    long Id, i;

    CheckSpecificationPart();
    if (ProblemType == HPP)
        Dimension--;
    if (strcmp(DisplayDataType, "TWOD_DISPLAY"))
        eprintf
                ("DISPLAY_DATA_SECTION conflicts with DISPLAY_DATA_TYPE: %s",
                 DisplayDataType);
    if (!FirstNode)
        CreateNodes();
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);
    N = FirstNode;
    for (i = 1; i <= Dimension; i++) {
        if (!fscanf(ProblemFile, "%ld", &Id))
            eprintf("Missing nodes in DIPLAY_DATA_SECTION");
        if (Id <= 0 || Id > Dimension)
            eprintf("(DIPLAY_DATA_SECTION) Node number out of range: %ld",
                    Id);
        N = &NodeSet[Id];
        if (N->V == 1)
            eprintf("(DIPLAY_DATA_SECTION) Node number occours twice: %ld",
                    N->Id);
        N->V = 1;
        if (!fscanf(ProblemFile, "%lf", &N->X))
            eprintf("Missing X-coordinate in DIPLAY_DATA_SECTION");
        if (!fscanf(ProblemFile, "%lf", &N->Y))
            eprintf("Missing Y-coordinate in DIPLAY_DATA_SECTION");
    }
    N = FirstNode;
    do
        if (!N->V)
            break;
    while ((N = N->Suc) != FirstNode);
    if (!N->V)
        eprintf("(DIPLAY_DATA_SECTION) No coordinates given for node %ld",
                N->Id);
    if (ProblemType == HPP)
        Dimension++;
}

static void Read_DISPLAY_DATA_TYPE() {
    long i;

    if (!(DisplayDataType = Copy(strtok(0, Delimiters))))
        eprintf("(DISPLAY_DATA_TYPE): string expected");
    for (i = 0; i < strlen(DisplayDataType); i++)
        DisplayDataType[i] = (char) toupper(DisplayDataType[i]);
    if (strcmp(DisplayDataType, "COORD_DISPLAY") != 0 &&
        strcmp(DisplayDataType, "TWOD_DISPLAY") != 0 &&
        strcmp(DisplayDataType, "NO_DISPLAY") != 0)
        eprintf("Unknown DISPLAY_DATA_TYPE: %s", DisplayDataType);
}

static void Read_EDGE_DATA_FORMAT() {
    long i;

    if (!(EdgeDataFormat = Copy(strtok(0, Delimiters))))
        eprintf("(EDGE_DATA_FORMAT): string expected");
    for (i = 0; i < strlen(EdgeDataFormat); i++)
        EdgeDataFormat[i] = (char) toupper(EdgeDataFormat[i]);
    if (strcmp(EdgeDataFormat, "EDGE_LIST") != 0 &&
        strcmp(EdgeDataFormat, "ADJ_LIST") != 0)
        eprintf("Unknown EDGE_DATA_FORMAT: %s", EdgeDataFormat);
}

static void Read_EDGE_DATA_SECTION() {
    Node *Ni, *Nj;
    long i, j;

    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes();
    if (ProblemType == HPP)
        Dimension--;
    if (!strcmp(EdgeDataFormat, "EDGE_LIST")) {
        if (!fscanf(ProblemFile, "%ld\n", &i))
            i = -1;
        while (i != -1) {
            if (i <= 0 || i > Dimension)
                eprintf
                        ("(EDGE_DATA_SECTION) Node number out of range: %ld",
                         i);
            fscanf(ProblemFile, "%ld\n", &j);
            if (j <= 0 || j > Dimension)
                eprintf
                        ("(EDGE_DATA_SECTION) Node number out of range: %ld",
                         j);
            if (i == j)
                eprintf("(EDGE_DATA_SECTION) Illgal edge: %ld to %ld", i,
                        j);
            Ni = &NodeSet[i];
            Nj = &NodeSet[j];
            if (!Ni->CandidateSet) {
                Ni->V = 1;
                assert(Ni->CandidateSet =
                               (Candidate *) calloc(2, sizeof(Candidate)));
                Ni->CandidateSet[0].To = Nj;
                Ni->CandidateSet[0].Cost = 0;
            } else {
                Ni->CandidateSet[Ni->V].To = Nj;
                Ni->CandidateSet[Ni->V].Cost = 0;
                assert(Ni->CandidateSet =
                               (Candidate *) realloc(Ni->CandidateSet,
                                                     (++Ni->V +
                                                      1) * sizeof(Candidate)));
                Ni->CandidateSet[Ni->V].To = 0;
            }
            if (!Nj->CandidateSet) {
                Nj->V = 1;
                assert(Nj->CandidateSet =
                               (Candidate *) calloc(2, sizeof(Candidate)));
                Nj->CandidateSet[0].To = Ni;
                Nj->CandidateSet[0].Cost = 0;
            } else {
                Nj->CandidateSet[Nj->V].To = Ni;
                Nj->CandidateSet[Nj->V].Cost = 0;
                assert(Nj->CandidateSet =
                               (Candidate *) realloc(Nj->CandidateSet,
                                                     (++Nj->V +
                                                      1) * sizeof(Candidate)));
                Nj->CandidateSet[Nj->V].To = 0;
            }
            fscanf(ProblemFile, "%ld\n", &i);
        }
    } else if (!strcmp(EdgeDataFormat, "ADJ_LIST")) {
        Ni = FirstNode;
        do
            Ni->V = 0;
        while ((Ni = Ni->Suc) != FirstNode);
        if (!fscanf(ProblemFile, "%ld\n", &i))
            i = -1;
        while (i != -1) {
            if (i <= 0 || i > Dimension)
                eprintf
                        ("(EDGE_DATA_SECTION) Node number out of range: %ld",
                         i);
            Ni = &NodeSet[i];
            fscanf(ProblemFile, "%ld\n", &j);
            while (j != -1) {
                if (j <= 0 || j > Dimension)
                    eprintf
                            ("(EDGE_DATA_SECTION) Node number out of range: %ld",
                             j);
                if (i == j)
                    eprintf("(EDGE_DATA_SECTION) Illgal edge: %ld to %ld",
                            i, j);
                Nj = &NodeSet[j];
                if (!Ni->CandidateSet) {
                    Ni->V = 1;
                    assert(Ni->CandidateSet =
                                   (Candidate *) calloc(2, sizeof(Candidate)));
                    Ni->CandidateSet[0].To = Nj;
                    Ni->CandidateSet[0].Cost = 0;
                } else {
                    Ni->CandidateSet[Ni->V].To = Nj;
                    Ni->CandidateSet[Ni->V].Cost = 0;
                    assert(Ni->CandidateSet =
                                   (Candidate *) realloc(Ni->CandidateSet,
                                                         (++Ni->V +
                                                          1) * sizeof(Candidate)));
                    Ni->CandidateSet[Ni->V].To = 0;
                }
                if (!Nj->CandidateSet) {
                    Nj->V = 1;
                    assert(Nj->CandidateSet =
                                   (Candidate *) calloc(2, sizeof(Candidate)));
                    Nj->CandidateSet[0].To = Ni;
                    Nj->CandidateSet[0].Cost = 0;
                } else {
                    Nj->CandidateSet[Nj->V].To = Ni;
                    Nj->CandidateSet[Nj->V].Cost = 0;
                    assert(Nj->CandidateSet =
                                   (Candidate *) realloc(Nj->CandidateSet,
                                                         (++Nj->V +
                                                          1) * sizeof(Candidate)));
                    Nj->CandidateSet[Nj->V].To = 0;
                }
                fscanf(ProblemFile, "%ld\n", &j);
            }
            fscanf(ProblemFile, "%ld\n", &i);
        }
    } else
        eprintf("(EDGE_DATA_SECTION) No EDGE_DATA_FORMAT specified");
    if (ProblemType == HPP)
        Dimension++;
    Distance = Distance_1;
}

static void Read_EDGE_WEIGHT_FORMAT() {
    long i;

    if (!(EdgeWeightFormat = Copy(strtok(0, Delimiters))))
        eprintf("(EDGE_WEIGHT_FORMAT): string expected");
    for (i = 0; i < strlen(EdgeWeightFormat); i++)
        EdgeWeightFormat[i] = (char) toupper(EdgeWeightFormat[i]);
    if (!strcmp(EdgeWeightFormat, "FUNCTION"))
        WeightFormat = FUNCTION;
    else if (!strcmp(EdgeWeightFormat, "FULL_MATRIX"))
        WeightFormat = FULL_MATRIX;
    else if (!strcmp(EdgeWeightFormat, "UPPER_ROW"))
        WeightFormat = UPPER_ROW;
    else if (!strcmp(EdgeWeightFormat, "LOWER_ROW"))
        WeightFormat = LOWER_ROW;
    else if (!strcmp(EdgeWeightFormat, "UPPER_DIAG_ROW"))
        WeightFormat = UPPER_DIAG_ROW;
    else if (!strcmp(EdgeWeightFormat, "LOWER_DIAG_ROW"))
        WeightFormat = LOWER_DIAG_ROW;
    else if (!strcmp(EdgeWeightFormat, "UPPER_COL"))
        WeightFormat = UPPER_COL;
    else if (!strcmp(EdgeWeightFormat, "LOWER_COL"))
        WeightFormat = LOWER_COL;
    else if (!strcmp(EdgeWeightFormat, "UPPER_DIAG_COL"))
        WeightFormat = UPPER_DIAG_COL;
    else if (!strcmp(EdgeWeightFormat, "LOWER_DIAG_COL"))
        WeightFormat = LOWER_DIAG_COL;
    else
        eprintf("Unknown EDGE_WEIGHT_FORMAT: %s", EdgeWeightFormat);
}

static void Read_EDGE_WEIGHT_SECTION() {
    Node *Ni, *Nj;
    long i, j, n, W;

    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes();
    if (ProblemType != ATSP) {
        assert(CostMatrix =
                       (long *) calloc(Dimension * (Dimension - 1) / 2,
                                       sizeof(long)));
        Ni = FirstNode->Suc;
        do {
            Ni->C = &CostMatrix[(Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
        } while ((Ni = Ni->Suc) != FirstNode);
    } else {
        n = Dimension / 2;
        assert(CostMatrix = (long *) calloc(n * n, sizeof(long)));
        for (Ni = FirstNode; Ni->Id <= n; Ni = Ni->Suc)
            Ni->C = &CostMatrix[(Ni->Id - 1) * n] - 1;
    }
    if (ProblemType == HPP)
        Dimension--;
    switch (WeightFormat) {
        case FULL_MATRIX:
            if (ProblemType == ATSP) {
                long n = Dimension / 2;
                for (i = 1; i <= n; i++) {
                    Ni = &NodeSet[i];
                    for (j = 1; j <= n; j++) {
                        if (!fscanf(ProblemFile, "%ld", &W))
                            eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                        Ni->C[j] = W;
                        if (i != j && W > M)
                            M = W;
                    }
                    Nj = &NodeSet[i + n];
                    if (!Ni->FixedTo1)
                        Ni->FixedTo1 = Nj;
                    else if (!Ni->FixedTo2)
                        Ni->FixedTo2 = Nj;
                    if (!Nj->FixedTo1)
                        Nj->FixedTo1 = Ni;
                    else if (!Nj->FixedTo2)
                        Nj->FixedTo2 = Ni;
                }
                Distance = Distance_ATSP;
                WeightType = -1;
            } else
                for (i = 1, Ni = FirstNode; i <= Dimension; i++, Ni = Ni->Suc) {
                    for (j = 1; j <= Dimension; j++) {
                        if (!fscanf(ProblemFile, "%ld", &W))
                            eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                        if (j < i)
                            Ni->C[j] = W;
                    }
                }
            break;
        case UPPER_ROW:
            for (i = 1, Ni = FirstNode; i < Dimension; i++, Ni = Ni->Suc) {
                for (j = i + 1, Nj = Ni->Suc; j <= Dimension;
                     j++, Nj = Nj->Suc) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    Nj->C[i] = W;
                }
            }
            break;
        case LOWER_ROW:
            for (i = 2, Ni = FirstNode->Suc; i <= Dimension; i++, Ni = Ni->Suc) {
                for (j = 1; j < i; j++) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    Ni->C[j] = W;
                }
            }
            break;
        case UPPER_DIAG_ROW:
            for (i = 1, Ni = FirstNode; i <= Dimension; i++, Ni = Ni->Suc) {
                for (j = i, Nj = Ni; j <= Dimension; j++, Nj = Nj->Suc) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    if (i != j)
                        Nj->C[i] = W;
                }
            }
            break;
        case LOWER_DIAG_ROW:
            for (i = 1, Ni = FirstNode; i <= Dimension; i++, Ni = Ni->Suc) {
                for (j = 1; j <= i; j++) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    if (j != i)
                        Ni->C[j] = W;
                }
            }
            break;
        case UPPER_COL:
            for (j = 2, Nj = FirstNode->Suc; j <= Dimension; j++, Nj = Nj->Suc) {
                for (i = 1; i < j; i++) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    Nj->C[i] = W;
                }
            }
            break;
        case LOWER_COL:
            for (j = 1, Nj = FirstNode; j < Dimension; j++, Nj = Nj->Suc) {
                for (i = j + 1, Ni = Nj->Suc; i <= Dimension;
                     i++, Ni = Ni->Suc) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    Ni->C[j] = W;
                }
            }
            break;
        case UPPER_DIAG_COL:
            for (j = 1, Nj = FirstNode; j <= Dimension; j++, Nj = Nj->Suc) {
                for (i = 1; i <= j; i++) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    if (i != j)
                        Nj->C[i] = W;
                }
            }
            break;
        case LOWER_DIAG_COL:
            for (j = 1, Nj = FirstNode; j <= Dimension; j++, Ni = Ni->Suc) {
                for (i = j, Ni = Nj; i <= Dimension; i++, Ni = Ni->Suc) {
                    if (!fscanf(ProblemFile, "%ld", &W))
                        eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    if (i != j)
                        Ni->C[j] = W;
                }
            }
            break;
    }
    if (ProblemType == HPP)
        Dimension++;
}

static void Read_EDGE_WEIGHT_TYPE() {
    long i;

    if (!(EdgeWeightType = Copy(strtok(0, Delimiters))))
        eprintf("(EDGE_WEIGHT_TYPE): string expected");
    for (i = 0; i < strlen(EdgeWeightType); i++)
        EdgeWeightType[i] = (char) toupper(EdgeWeightType[i]);
    if (!strcmp(EdgeWeightType, "EXPLICIT")) {
        WeightType = EXPLICIT;
        Distance = Distance_EXPLICIT;
    } else if (!strcmp(EdgeWeightType, "EUC_2D")) {
        WeightType = EUC_2D;
        Distance = Distance_EUC_2D;
        c = c_EUC_2D;
    } else if (!strcmp(EdgeWeightType, "EUC_3D")) {
        WeightType = EUC_3D;
        Distance = Distance_EUC_3D;
        c = c_EUC_3D;
    } else if (!strcmp(EdgeWeightType, "MAX_2D")) {
        WeightType = MAX_2D;
        Distance = Distance_MAX_2D;
    } else if (!strcmp(EdgeWeightType, "MAX_3D")) {
        WeightType = MAX_3D;
        Distance = Distance_MAX_3D;
    } else if (!strcmp(EdgeWeightType, "MAN_2D")) {
        WeightType = MAN_2D;
        Distance = Distance_MAN_2D;
    } else if (!strcmp(EdgeWeightType, "MAN_3D")) {
        WeightType = MAN_3D;
        Distance = Distance_MAN_3D;
    } else if (!strcmp(EdgeWeightType, "CEIL_2D")) {
        WeightType = CEIL_2D;
        Distance = Distance_CEIL_2D;
        c = c_CEIL_2D;
    } else if (!strcmp(EdgeWeightType, "CEIL_3D")) {
        WeightType = CEIL_3D;
        Distance = Distance_CEIL_3D;
        c = c_CEIL_3D;
    } else if (!strcmp(EdgeWeightType, "GEO")) {
        WeightType = GEO;
        Distance = Distance_GEO;
        c = c_GEO;
    } else if (!strcmp(EdgeWeightType, "GEOM")) {
        WeightType = GEOM;
        Distance = Distance_GEOM;
        c = c_GEOM;
    } else if (!strcmp(EdgeWeightType, "ATT")) {
        WeightType = ATT;
        Distance = Distance_ATT;
    } else if (!strcmp(EdgeWeightType, "XRAY1") ||
               !strcmp(EdgeWeightType, "XRAY2") ||
               !strcmp(EdgeWeightType, "SPECIAL"))
        eprintf("EDGE_WEIGHT_TYPE not implemented: %s", EdgeWeightType);
    else
        eprintf("Unknown EDGE_WEIGHT_TYPE: %s", EdgeWeightType);
}

static void Read_FIXED_EDGES_SECTION() {
    Node *Ni, *Nj;
    long i, j;
    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes();
    if (ProblemType == HPP)
        Dimension--;
    if (!fscanf(ProblemFile, "%ld\n", &i))
        i = -1;
    while (i != -1) {
        if (i <= 0 || i > Dimension)
            eprintf("(FIXED_EDGES_SECTION) Node number out of range: %ld",
                    i);
        fscanf(ProblemFile, "%ld\n", &j);
        if (j <= 0 || j > Dimension)
            eprintf("(FIXED_EDGES_SECTION) Node number out of range: %ld",
                    j);
        if (i == j)
            eprintf("(FIXED_EDGES_SECTION) Illgal edge: %ld to %ld", i, j);
        if (ProblemType == ATSP)
            i += Dimension / 2;
        Ni = &NodeSet[i];
        Nj = &NodeSet[j];
        if (!Ni->FixedTo1)
            Ni->FixedTo1 = Nj;
        else if (!Ni->FixedTo2)
            Ni->FixedTo2 = Nj;
        else
            eprintf("(FIXED_EDGES_SECTION) Illegal fix: %ld to %ld", i, j);
        if (!Nj->FixedTo1)
            Nj->FixedTo1 = Ni;
        else if (!Nj->FixedTo2)
            Nj->FixedTo2 = Ni;
        else
            eprintf("(FIXED_EDGES_SECTION) Illegal fix: %ld to %ld", i, j);
        fscanf(ProblemFile, "%ld\n", &i);
    }
    if (ProblemType == HPP)
        Dimension++;
}

static void Read_NODE_COORD_SECTION() {
    Node *N;
    long Id, i;

    CheckSpecificationPart();
    if (CoordType == TWOD_COORDS || CoordType == THREED_COORDS)
        eprintf("NODE_COORD_SECTION conflicts with NODE_COORDS_TYPE: %s",
                NodeCoordType);
    if (!FirstNode)
        CreateNodes();
    N = FirstNode;
    if (ProblemType == HPP)
        Dimension--;
    for (i = 1; i <= Dimension; i++) {
        if (!fscanf(ProblemFile, "%ld", &Id))
            eprintf("Missing nodes in NODE_COORD_SECTION");
        if (Id <= 0 || Id > Dimension)
            eprintf("(NODE_COORD_SECTION) Node number out of range: %ld",
                    Id);
        N = &NodeSet[Id];
        if (N->V == 1)
            eprintf("(NODE_COORD_SECTION) Node number occours twice: %ld",
                    N->Id);
        N->V = 1;
        if (!fscanf(ProblemFile, "%lf", &N->X))
            eprintf("Missing X-coordinate in NODE_COORD_SECTION");
        if (!fscanf(ProblemFile, "%lf", &N->Y))
            eprintf("Missing Y-coordinate in NODE_COORD_SECTION");
        if (CoordType == THREED_COORDS
            && !fscanf(ProblemFile, "%lf", &N->Z))
            eprintf("Missing Z-coordinate in NODE_COORD_SECTION");
        if (strcmp(Name, "tsp225") && strcmp(Name, "d657"))
            continue;
        N->X = (float) N->X;
        N->Y = (float) N->Y;
        N->Z = (float) N->Z;
    }
    N = FirstNode;
    do
        if (!N->V && N->Id <= Dimension)
            break;
    while ((N = N->Suc) != FirstNode);
    if (!N->V)
        eprintf("(NODE_COORD_SECTION) No coordinates given for node %ld",
                N->Id);
    if (ProblemType == HPP)
        Dimension++;
}

static void Read_NODE_COORD_TYPE() {
    long i;

    if (!(NodeCoordType = Copy(strtok(0, Delimiters))))
        eprintf("(NODE_COORD_TYPE): string expected");
    for (i = 0; i < strlen(NodeCoordType); i++)
        NodeCoordType[i] = (char) toupper(NodeCoordType[i]);
    if (!strcmp(NodeCoordType, "TWOD_COORDS"))
        CoordType = TWOD_COORDS;
    else if (!strcmp(NodeCoordType, "THREED_COORDS"))
        CoordType = THREED_COORDS;
    else if (!strcmp(NodeCoordType, "NO_COORDS"))
        CoordType = NO_COORDS;
    else
        eprintf("Unknown NODE_COORD_TYPE: %s", NodeCoordType);
}

static void Read_TOUR_SECTION(FILE **File) {
    Node *First, *Last, *N;
    long i, k;

    if (!FirstNode)
        CreateNodes();
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);
    if (ProblemType == HPP)
        Dimension--;
    if (!fscanf(*File, "%ld\n", &i))
        i = -1;
    for (k = 0; k < Dimension && i != -1; k++) {
        if (i <= 0 || i > Dimension)
            eprintf("(TOUR_SECTION) Node number out of range: %ld", i);
        N = &NodeSet[i];
        if (N->V == 1)
            eprintf("(TOUR_SECTION) Node number occours twice: %ld",
                    N->Id);
        N->V = 1;
        if (k == 0)
            First = Last = N;
        else {
            if (File == &InputTourFile)
                Last->OptimumSuc = N;
            else if (File == &InitialTourFile)
                Last->InitialSuc = N;
            else if (File == &MergeTourFile[0])
                Last->MergeSuc[0] = N;
            else if (File == &MergeTourFile[1])
                Last->MergeSuc[1] = N;
            Last = N;
        }
        if (ProblemType == ATSP) {
            N = &NodeSet[i + Dimension / 2];
            N->V = 1;
            if (File == &InputTourFile)
                Last->OptimumSuc = N;
            else if (File == &InitialTourFile)
                Last->InitialSuc = N;
            else if (File == &MergeTourFile[0])
                Last->MergeSuc[0] = N;
            else if (File == &MergeTourFile[1])
                Last->MergeSuc[1] = N;
            Last = N;
        }
        fscanf(*File, "%ld\n", &i);
    }
    if (Last) {
        if (File == &InputTourFile)
            Last->OptimumSuc = First;
        else if (File == &InitialTourFile)
            Last->InitialSuc = First;
        else if (File == &MergeTourFile[0])
            Last->MergeSuc[0] = First;
        else if (File == &MergeTourFile[1])
            Last->MergeSuc[1] = First;
    }
    N = FirstNode;
    do
        if (!N->V)
            break;
    while ((N = N->Suc) != FirstNode);
    if (!N->V)
        eprintf("(TOUR_SECTION) Node is missing: %ld", N->Id);
    if (ProblemType == HPP)
        Dimension++;
}

static void Read_TYPE() {
    long i;

    if (!(Type = Copy(strtok(0, Delimiters))))
        eprintf("(TYPE): string expected");
    for (i = 0; i < strlen(Type); i++)
        Type[i] = (char) toupper(Type[i]);
    if (!strcmp(Type, "TSP"))
        ProblemType = TSP;
    else if (!strcmp(Type, "ATSP"))
        ProblemType = ATSP;
    else if (!strcmp(Type, "SOP")) {
        ProblemType = SOP;
        eprintf("(TYPE) Type not implemented: %s", Type);
    } else if (!strcmp(Type, "HCP"))
        ProblemType = HCP;
    else if (!strcmp(Type, "CVRP")) {
        ProblemType = CVRP;
        eprintf("(TYPE) Type not implemented: %s", Type);
    } else if (!strcmp(Type, "TOUR")) {
        ProblemType = TOUR;
        eprintf("(TYPE) Type not implemented: %s", Type);
    } else if (!strcmp(Type, "HPP"))
        ProblemType = HPP;
    else
        eprintf("Unknown TYPE: %s", Type);
}

/*
   The ReadTour function reads a tour from a file, if requested in the parameter file
   (by using the INPUT_TOUR keyword).  The tour is used to limit the search 
   (the last edge to be removed in a non-gainful move must not belong to the tour). 

   The format is as follows: 

   OPTIMUM = <real>
   Known optimal tour length. A run will be terminated as soon as a tour 
   length less than or equal to optimum is achieved.
   Default: DBL_MAX.

   TOUR_SECTION :
   A tour is specified in this section. The tour is given by a list of integers
   giving the sequence in which the nodes are visited in the tour. The tour is
   terminated by a -1. 

   EOF
   Terminates the input data. The entry is optional.

   Other keywords in TSPLIB format may be included in the file, but they are ignored.
*/

void ReadTour(char *FileName, FILE **File) {
    char *Line, *Keyword;
    int i;
    if (!(*File = fopen(FileName, "r")))
        eprintf("Cannot open %s", FileName);
    while (Line = ReadLine(*File)) {
        if (!(Keyword = strtok(Line, Delimiters)))
            continue;
        for (i = 0; i < strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
            if (!strcmp(Keyword, "COMMENT") ||
            !strcmp(Keyword, "DEMAND_SECTION") ||
            !strcmp(Keyword, "DEPOT_SECTION") ||
            !strcmp(Keyword, "DIMENSION") ||
            !strcmp(Keyword, "DISPLAY_DATA_SECTION") ||
            !strcmp(Keyword, "DISPLAY_DATA_TYPE") ||
            !strcmp(Keyword, "EDGE_DATA_FORMAT") ||
            !strcmp(Keyword, "EDGE_DATA_SECTION") ||
            !strcmp(Keyword, "EDGE_WEIGHT_FORMAT") ||
            !strcmp(Keyword, "EDGE_WEIGHT_SECTION") ||
            !strcmp(Keyword, "EDGE_WEIGHT_TYPE") ||
            !strcmp(Keyword, "FIXED_EDGES_SECTION") ||
            !strcmp(Keyword, "NAME") ||
            !strcmp(Keyword, "NODE_COORD_SECTION") ||
            !strcmp(Keyword, "NODE_COORD_TYPE") ||
            !strcmp(Keyword, "TYPE"));
        else if (strcmp(Keyword, "OPTIMUM") == 0) {
            if (sscanf(strtok(0, Delimiters), "%lf", &Optimum) == 0)
                eprintf("(OPTIMUM): real expected");
            Optimum = floor(Optimum + 0.5);
        } else if (!strcmp(Keyword, "TOUR_SECTION"))
            Read_TOUR_SECTION(File);
        else if (strcmp(Keyword, "EOF") == 0)
            break;
        else
            eprintf("Unknown Keyword: %s", Keyword);
    }
    fclose(*File);
}

SwapMove *CreateSwapMoveTree(int maxDepth, int maxChildren, long gain) {
    if (maxDepth == 0) {
        return NULL;
    }

    SwapMove *root = CreateSwapMove(idMove, gain);
    idMove++;

    for (int i = 0; i < maxChildren; i++) {
        SwapMove *child = CreateSwapMoveTree(maxDepth - 1, maxChildren, gain);
        if (child != NULL) {
            child->parent = root;
            if (root->firstChild == NULL) {
                root->firstChild = child;

            } else {
                SwapMove *lastSibling = root->firstChild;
                while (lastSibling->nextSibling != NULL) {
                    lastSibling = lastSibling->nextSibling;
                }
                lastSibling->nextSibling = child;
            }
        }
    }

    return root;
}

SwapMove *CreateSwapMove(int id, long gain) {
    SwapMove *newNode = (SwapMove *) malloc(sizeof(SwapMove));
    newNode->id = id;
    newNode->Gain = gain;
    newNode->t1 = newNode->t2 = newNode->t3 = newNode->t4 = NULL;
    newNode->parent = NULL;
    newNode->firstChild = NULL;
    newNode->nextSibling = NULL;
    newNode->NumChildren = 0;
    newNode->IsBest = 0;
    newNode->IsActive = 0;

    return newNode;
}
