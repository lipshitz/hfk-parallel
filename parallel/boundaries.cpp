#include <sys/time.h>
#include <iostream>
#include <stdlib.h>
#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <queue>
#include <list>
#include <assert.h>
#include <algorithm>
#include "matrix-z2z3.h"
#include "reduction-functions.h"

#include <unistd.h>

using std::list;
using std::vector;

int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

// Globals

//const int gridsize = 12; // arc-index
//const int gridsize = 12;

int gridsize = 12; // arc-index
int *white;
int *black;
int default_white[12] = {9,5,11,7,8,1,10,4,0,3,2,6};
int default_black[12] = {1,0,4,3,2,6,5,9,8,11,7,10};
// Trefoil
//int white[5] = {1, 2, 3, 4, 0};
//int black[5] = {4, 0, 1, 2, 3};

//int white[10] = {8,7,6,5,4,3,2,9,1,0};
//int black[10] = {1,3,9,0,7,5,8,4,6,2};

// Kinoshita-Terasaka KT_{2,1}
//int white[11]={5,10,9,4,8,0,1,6,7,2,3};
//int black[11]={0,6,1,7,10,2,5,9,3,4,8};



// Don't waste time computing factorials.  Look them up.
// Fill in the big ones later since g++ doesn't seem to like big constants
long long Factorial[16] = {
  1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,
  0,0,0};

//braid: (xy)^{-5}
//int white[10] = {7,6,5,3,4,1,2,0,9,8};
//int black[10] = {0,9,8,6,7,4,5,3,2,1};

//another braid
//int white[11] = {5,2,3,1,4,6,7,8,9,0};
//int black[11] = {1,9,0,5,2,3,4,6,7,8};

//int white[2] = {0,1};
//int black[2] = {1,0};

// Function Prototypes

void getPerm( long long n, int *P);
void NextPerm(short counter[], int h[]);
bool RectDotFree(int xll, int yll, int xur, int yur, int which); 
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
// Decides whether one of the four rectangles on the torus with corners at(xll,yll) and (xur,yur) contains no white or black dots
/*
 1 | 2 | 1
---+---+---
 3 | 0 | 3
---+---+---
 1 | 2 | 1 

*/
long long getIndex( int *P);
int WindingNumber(int x, int y); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y []);
int NumComp(); //Returns the number of components of the link.
bool ValidGrid();
int Find(vector<long long> & V, long long x); // Returns i if V[i]=x or -1 if x isn't V[i] for any i

// Main

int main(int argc, char *argv[]){

 // set up mpi
 int n_proc, rank;
 MPI_Init( &argc, &argv );
 MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
 MPI_Comm_rank( MPI_COMM_WORLD, &rank );
 
 char *saveDir = read_string( argc, argv, "-i", NULL );
 bool printMatrices = read_int( argc, argv, "-p", 0 );
 bool justTime = read_int( argc, argv, "-t", 0 );
 int aGrading = read_int( argc, argv, "-a", -100 );
 int mGrading = read_int( argc, argv, "-m", -100 );
 if( aGrading == -100 || mGrading == -100 ) {
   printf("Must specify in which grading to calculate the boundary image and kernel dimensions with -a and -g\n");
   MPI_Finalize();
   exit(0);
 }
 if( n_proc <= 1 ) {
   printf("Run with at least 2 threads\n");
   MPI_Finalize();
   exit(0);
 }
 if( !saveDir ) {
   printf("Use -i to specify where the output of gradings was written\n");
   MPI_Finalize();
   exit(0);
 }

 char knotFile[30];
 sprintf(knotFile, "%s/knot", saveDir);
 FILE *f = fopen(knotFile, "r");
 if( !f ) {
   printf("Error opening file %s\n", knotFile);
   MPI_Finalize();
   exit(-1);
 }
 fscanf(f, "%d\n", &gridsize);
 white = (int*) malloc( gridsize*sizeof(int) );
 black = (int*) malloc( gridsize*sizeof(int) );
 for( int i = 0; i < gridsize; i++ )
   fscanf(f, "%d ", white+i);
 fscanf(f, "\n");
 for( int i = 0; i < gridsize; i++ )
   fscanf(f, "%d ", black+i);
 fclose(f);

 Factorial[13] = 13*Factorial[12];
 Factorial[14] = 14*Factorial[13];
 Factorial[15] = 15*Factorial[14];

 int amin=0;
 int amax=20;

 int  numcomp = NumComp();
 if( rank == 0 && !justTime )
   printf("Number of components: %d\n", numcomp);

 if( rank == 0 )
   if(!ValidGrid()) {printf("Invalid grid!!\n"); return 0;} // Check that the grid is valid
 double starttime;
 if( rank == 0 )
   starttime = read_timer(); // Used to record how long this takes

 // Record for later use whether every possible rectangle has a black or white dot in it
 // This will speed boundary computations.
 if( rank == 0 && !justTime )
   printf("Computing which rectangles on the torus have no black or white dots inside.\n");
 bool Rectangles[gridsize][gridsize][gridsize][gridsize][4];
 for(int xll=0; xll < gridsize; xll++) {
  for(int xur=xll+1; xur < gridsize; xur++) {
   for(int yll=0; yll < gridsize; yll++) {
    for(int yur=yll+1; yur < gridsize; yur++) {
     Rectangles[xll][yll][xur][yur][0] = RectDotFree(xll,yll,xur,yur,0);
     Rectangles[xll][yll][xur][yur][1] = RectDotFree(xll,yll,xur,yur,1);
     Rectangles[xll][yll][xur][yur][2] = RectDotFree(xll,yll,xur,yur,2);
     Rectangles[xll][yll][xur][yur][3] = RectDotFree(xll,yll,xur,yur,3);     
    }
   }
  }
 }

 vector<long long> rows, cols;
 int g[gridsize];
#define NGEN 5
#define GEN 6

 // Read in the sorted list of generators for both rows and columns
 if( rank == 0 ) {
   printf("reading in rows\n");
   char rowFile[30];
   sprintf(rowFile, "%s/gen%d,%d.dat", saveDir, aGrading+30, mGrading+30);
   FILE *f;
   f = fopen(rowFile, "r");
   if( f ) {
     int num;
     fscanf(f, "%d\n", &num);
     rows.reserve(num);
     for( int i = 0; i < num; i++ ) {
       long long g;
       fscanf(f, "%lld\n", &g);
       rows.push_back(g);
     }
     fclose(f);
   } else {
     printf("waring, no such file %s\n", rowFile);
   }
   printf("reading in columns\n");
   char colFile[30];
   sprintf(colFile, "%s/gen%d,%d.dat", saveDir, aGrading+30, mGrading+31);
   f = fopen(colFile, "r");
   if( f ) {
     int num;
     fscanf(f, "%d\n", &num);
     cols.reserve(num);
     for( int i = 0; i < num; i++ ) {
       long long g;
       fscanf(f, "%lld\n", &g);
       cols.push_back(g);
     }
     fclose(f);
   } else {
     printf("waring, no such file %s\n", colFile);
   }
 }


 // Send rows and cols to everyone
 int num_generators;
 if( rank == 0 )
   num_generators = rows.size();
 MPI_Bcast( &num_generators, 1, MPI_INT, 0, MPI_COMM_WORLD );
 long long *buffer;
 if( rank == 0 )
   buffer = &(rows.front());
 else
   buffer = (long long*) malloc( num_generators * sizeof(long long) );
 MPI_Bcast( buffer, num_generators, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD );
 if( rank != 0 )
   rows.assign(buffer, buffer+num_generators);

 if( rank == 0 )
   num_generators = cols.size();
 MPI_Bcast( &num_generators, 1, MPI_INT, 0, MPI_COMM_WORLD );

 if( rank == 0 )
   buffer = &(cols.front());
 else
   buffer = (long long*) malloc( num_generators * sizeof(long long) );
 MPI_Bcast( buffer, num_generators, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD );
 if( rank != 0 )
   cols.assign(buffer, buffer+num_generators);
 
 // if j is 0, we are done, because there is no calculation with those cols
 //if (j == 0)
 //continue;

 int imageDimension, kernelDimension;

 if( cols.size() == 0 || rows.size() == 0 ) {
   imageDimension = 0;
   kernelDimension = cols.size();
 } else {
   // Decide who has ownership of what subset of the generators (without communication).
   int colsPerProcMin = num_generators/n_proc;
   int extraCols = num_generators - n_proc * colsPerProcMin;
   int firstCol, numCols;
   if( rank < extraCols ) {
     firstCol = (rank)*(colsPerProcMin+1);
     numCols = colsPerProcMin+1;
   } else {
     firstCol = extraCols + (rank)*colsPerProcMin;
     numCols = colsPerProcMin;
   }
   printf("(%d) ownership determined, my block is %lu %d starting at %d\n", rank, rows.size(), numCols, firstCol);
   //printf("Here %d %d %d\n", i, j, rank);
   // Each proc calculates the part of the matrix that it owns, storing both rows and columns as in the serial code.  It holds entire columns, but only certain rows
   vector<Generator> GraphIn( rows.size() ); // Will hold boundary data.
   vector<Generator> GraphOut( numCols ); // Will hold boundary data.
   //printf("Filling %d %d\n", I, J);
   for(int index=0; index < numCols; index++) {
     getPerm(cols[firstCol+index],g);
     bool firstrect, secondrect;
     for(int i=0; i<gridsize; i++) {
       for(int j=i+1; j<gridsize; j++) {
	 if(g[i]<g[j]) {
	   firstrect = Rectangles[i][g[i]][j][g[j]][0];
	   for(int k=i+1; k<j && firstrect; k++) {
	     if(g[i] < g[k] && g[k] < g[j]) firstrect=0;
	   }
	   secondrect = Rectangles[i][g[i]][j][g[j]][1];
	   for(int k=0; k<i && secondrect; k++) {
	     if(g[k]<g[i] || g[k] > g[j]) secondrect=0;
	   }
	   for(int k=j+1; k<gridsize && secondrect; k++) {
	     if(g[k]<g[i] || g[k] > g[j]) secondrect=0;
	   }
	 }
	 if(g[j]<g[i]) {
	   firstrect = Rectangles[i][g[j]][j][g[i]][2];
	   for(int k=i+1; k<j && firstrect; k++) {
	     if(g[k]<g[j] || g[k] > g[i]) firstrect=0;
	   }
	   secondrect = Rectangles[i][g[j]][j][g[i]][3];
	   for(int k=0; k<i && secondrect; k++) {
	     if(g[k]>g[j] && g[k]<g[i]) secondrect=0;
	   }
	   for(int k=j+1; k<gridsize && secondrect; k++) {
	     if(g[k]>g[j] && g[k]<g[i]) secondrect=0;
	   }
	 }
	 if(firstrect != secondrect) { // Exactly one rectangle is a boundary
	   int gij [gridsize];
	   for(int k=0; k<i; k++) {
	     gij[k] = g[k];
	   }
	   gij[i]=g[j];
	   for(int k=i+1; k<j; k++) {
	     gij[k]=g[k];
	   } 
	   gij[j] = g[i];
	   for(int k=j+1; k<gridsize; k++) {
	     gij[k] = g[k];
	   }
	   long long Indexgij = getIndex(gij);
	   int indexgij = Find(rows,Indexgij);
	   if(indexgij==-1) {printf("Error with Alexander grading: %lld->%lld\n", cols[index], Indexgij); return 0; }
#ifdef FIELD_Z3
	   printf("Not implemented\n");
	   exit(-1);
	   //GraphOut[index].ones.push_back( indexgij );
	   //GraphIn[indexgij].ones.push_back( index );
	   // This is, of course, wrong
#else
	   GraphOut[index].ones.push_back( indexgij );
	   GraphIn[indexgij].ones.push_back( index+firstCol );
#endif
	 }
       }
     }
   }
   printf("(%d) Matrix filled\n", rank);
   
   // this skips over the first matrix in 10_125
   //if( i == 30 && j == 26 )
   //  continue;
   
   /* Now do the reduction.  Algorithm is:
      Each proc stores 'turn' which says from whom it next expects to recieve a column to subtract from its matrix (below where it is)
      Asynchronously wait to recieve a column from processor 'turn', on recepit.  Send that column on to the next processor (unless it came from them), immediatly subtract from my current top column (if appropriate), then enqueue for future processing.  increment 'turn'.  If it is now our turn, enqueue our top column and send it to the next processor.  Take a new top column, subtract everything in queue from it (if it has a 1 where the column in queue has its first 1).
      Meanwhile, while not receiving anything, take the first column off the queue and subtract it from all of our columns that have weight at its pivot
   */
#define TAG_COLUMN_SIZE 10
#define TAG_COLUMN 11
#define TAG_FINISHED 12

   int turn = 0; // the id of the matrix whose turn it is
   int globalTurn = 0; // This just keeps increasing
   std::deque<int*> colsToReduceBy;
   std::deque<MPI_Request*> outRequests; // keep these so we can make sure we have sent a given column when we dequeue it befor we delete it
   bool waitingForSize = false;
   bool waitingForCol = false;
   int waitingSource = -1; // the processor from which we are currently waiting
   bool colOfMineOnStack = false; // we don't want to ever have multiple columns of our own on the stack (it means the one we are about to add hasn't been reduced sufficiently yet)
   int nextProc = rank + 1;
   if( nextProc == n_proc )
     nextProc = 0;
   int prevProc = rank - 1;
   if( prevProc < 0 )
     prevProc = n_proc - 1;
   MPI_Request *inRequest;
   MPI_Request *finishedInRequest;
   MPI_Request *finishedOutRequest;
   int finishedOutSignal[2];
   int finishedInSignal[2];
   finishedInRequest = new MPI_Request();
   finishedOutRequest = NULL;
   MPI_Irecv( finishedInSignal, 2, MPI_INT, prevProc, TAG_FINISHED, MPI_COMM_WORLD, finishedInRequest );
   int inSize;
   int *inCol=0;
   //bool starting = true;
   list<int>::iterator indexInCol;
   list<int>::iterator indexInColEndValue = indexInCol;
   int *currentCol;
   MPI_Request *currentColOutRequest;
   int c = 0;
   int numProcReportingFinished = 0;
   int procsReportingFinished[n_proc]; // this stores the globalTurn at which each proc has reported it is finished, or -1 if it has not reported finished
   for( int p = 0; p < n_proc; p++ )
     procsReportingFinished[p] = -1;
   bool selfFinished = false;
   int zero = 0;
   
   while( numProcReportingFinished < n_proc ) { // main loop for the reduction.  The key principle is that we only do one thing per loop iteration
     bool allDone = true;
     for( int p = 0; p < n_proc; p++ )
       if( procsReportingFinished[p] == -1 || procsReportingFinished[p] > globalTurn ) {
	 allDone = false;
	 break;
       }
     if( allDone )
       break;
     // make sure turn is someone who is actually still working
     while( procsReportingFinished[turn] != -1 && procsReportingFinished[turn] <= globalTurn ) {
       increment(turn, 0, n_proc);
       globalTurn++;
     }
     //printf("(%d) in reduction, c = %d/%d, turn = %d(%d), numfinished = %d, selfFinished = %d(%d), colsToReduceBy.size = %lu, outRequests.size() = %lu, waiting-for-size=%d waiting-for-col=%d, colOfMine = %d\n", rank, c, numCols, turn, globalTurn, numProcReportingFinished, selfFinished, procsReportingFinished[rank], colsToReduceBy.size(), outRequests.size(), waitingForSize, waitingForCol, colOfMineOnStack);
     
     // check if we have received a notification that another processor is done
     if( finishedOutRequest ) {
       int outDone = 0;
       MPI_Test( finishedOutRequest, &outDone, MPI_STATUS_IGNORE );
       if( outDone ) {
	 delete finishedOutRequest;
	 finishedOutRequest = NULL;
	 if( numProcReportingFinished < n_proc-1 ||
	     numProcReportingFinished == n_proc-1 && selfFinished ) {
	   finishedInRequest = new MPI_Request();
	   //printf("(%d) prepared to receive another done signal\n", rank);
	   MPI_Irecv( finishedInSignal, 2, MPI_INT, prevProc, TAG_FINISHED, MPI_COMM_WORLD, finishedInRequest );
	 }
	 continue;
	 }
     }
     
     if( finishedInRequest ) {
       int someoneDone = 0;
       MPI_Test( finishedInRequest, &someoneDone, MPI_STATUS_IGNORE );
       if( someoneDone ) {
	 delete finishedInRequest;
	 finishedInRequest = NULL;
	 //printf("(%d) got finished signal from %d (nextProc=%d)\n", rank, finishedInSignal[0], nextProc);
	 procsReportingFinished[finishedInSignal[0]] = finishedInSignal[1];
	 numProcReportingFinished++;
	 if( finishedInSignal[0] != nextProc ) {
	   finishedOutRequest = new MPI_Request();
	   //printf("(%d) sending on finished signal from %d\n", rank, finishedInSignal[0]);
	   MPI_Isend( finishedInSignal, 2, MPI_INT, nextProc, TAG_FINISHED, MPI_COMM_WORLD, finishedOutRequest );
	 } else {
	   // put out a new request if there are any procs left not reporting done
	   if( numProcReportingFinished < n_proc-1 ||
	       numProcReportingFinished == n_proc-1 && selfFinished ) {
	     finishedInRequest = new MPI_Request();
	     //printf("(%d) prepared to receive another done signal\n", rank);
	     MPI_Irecv( &finishedInSignal, 2, MPI_INT, prevProc, TAG_FINISHED, MPI_COMM_WORLD, finishedInRequest );
	   }	   
	 }
	 continue;
       }
     }
     // check if we ourselves are done
     if( !selfFinished && c >= numCols ) {
       selfFinished = true;
       procsReportingFinished[rank] = globalTurn;
       numProcReportingFinished++;
       // send a message along the chain
       MPI_Request req;
       //printf("(%d) sending out my done signal\n", rank);
       finishedOutSignal[0] = rank;
       finishedOutSignal[1] = globalTurn;
       MPI_Isend( finishedOutSignal, 2, MPI_INT, nextProc, TAG_FINISHED, MPI_COMM_WORLD, &req);
       continue;
     }
     
     // check if we should prepare to receive
     if( turn != rank && !waitingForSize && !waitingForCol ) {
       inRequest = new MPI_Request();
       MPI_Irecv( &inSize, 1, MPI_INT, prevProc, TAG_COLUMN_SIZE, MPI_COMM_WORLD, inRequest );
       waitingSource = turn;
       waitingForSize = true;
       continue;
     }
     
     // check if we are waiting to receive a size (only continue if we actually received something)
     if( waitingForSize ) {
       // check if we have received the size, if so toggle waitingForSize, toggle waitingFolCol, issue the issue Irecv for the column
       //printf("(%d) Here 1\n", rank);
       int ready = 0;
       MPI_Status stat;
       MPI_Test( inRequest, &ready, &stat );
       //printf("(%d) Here 2\n", rank);
       if( ready ) {
	 delete inRequest;
	 waitingForSize = false;
	 if( inSize ) {
	   inCol = (int*) malloc( (inSize+2)*sizeof(int) );
	   inRequest = new MPI_Request();
	   MPI_Irecv( inCol, inSize+2, MPI_INT, prevProc, TAG_COLUMN, MPI_COMM_WORLD, inRequest );
	   waitingForCol = true;
	 } else {
	   waitingSource = -1;
	   globalTurn += rank-turn;
	   if( rank < turn )
	     globalTurn += n_proc;
	   turn = rank;
	   // a special case: if this process is finished, it needs to pass on a zero size so that the calculation can continue
	   if( selfFinished ) {
	     MPI_Request req;
	     MPI_Isend( &zero, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, MPI_COMM_WORLD, &req);
	   }
	 }
	 continue;
       }
     }
     
     // check if we are waiting to receive a column (only continue if we actually receive one)
     if( waitingForCol ) {
       int ready = 0;
       MPI_Test( inRequest, &ready, MPI_STATUS_IGNORE );
       if( ready ) {
	 delete inRequest;
	 for( int e = 2; e < 2+inSize; e++ )
	   if( inCol[e] >= rows.size() ) {
	     printf("(%d) FATAL receive error e=%d inCol[e]=%d numRows=%lu (%d)\n", rank, e, inCol[e], rows.size(), c ); 
	     printf("(%d) FATAL ", rank);
	     for( int f = 0; f < 2+inSize; f++ )
	       printf("%d ", inCol[f]);
	     printf("\n");
	     exit(0);
	   }
	 
	 waitingForCol = false;
	 waitingSource = -1;
	 colsToReduceBy.push_back(inCol);
	 globalTurn += inCol[0]-turn;
	 if( inCol[0] < turn )
	   globalTurn += n_proc;
	 //printf("(%d) received column %d %d %d\n", rank, inCol[0], inCol[1], inCol[2]);
	 turn = inCol[0]; // this should take care of dropped blank columns
	 increment(turn, 0, n_proc);
	 globalTurn++;
	 if( inCol[1] == 0 ) {
	   printf("(%d) Warning: just received a column of zero size, size should have been %d\n", rank, inSize);
	 }
	 if( inCol[0] != nextProc ) { // if it is, this column has made it all the way around
	   MPI_Request request1, *request2;
	   request2 = new MPI_Request();
	   //printf("(%d) passing along column %d %d %d\n", rank, inCol[0], inCol[1], inCol[2]);
	   MPI_Isend( inCol+1, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, MPI_COMM_WORLD, &request1 ); 
	   MPI_Isend( colsToReduceBy.back(), inCol[1]+2, MPI_INT, nextProc, TAG_COLUMN, MPI_COMM_WORLD, request2 );
	   outRequests.push_back(request2);
	 } else {
	   // need to keep outRequests in sync with colsToReduceBy
	   MPI_Request rnull = MPI_REQUEST_NULL;
	   outRequests.push_back(&rnull);
	 }
	 if( c < numCols ) {
	   resolveCols(GraphOut, GraphIn, inCol+2, inCol[1], c, firstCol);
	 }
	 continue;
       }
     }
     
     if( turn == rank && selfFinished ) {
       assert(false);
       /*
	 MPI_Request request1;
	 MPI_Isend( &zero, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, MPI_COMM_WORLD, &request1 );
	 increment(turn, 0, n_proc);
	 globalTurn++;
	 */
     }
     
     if( turn == rank && !selfFinished && // it is our turn to send 
	 !colOfMineOnStack // but only send if we don't already have a column on the stack (perhaps redundant with the next condition
	 ){
       int numEntries, *column;
       
       // Should we just check if all the columns have been applied at least this far?, perhaps outsize this if statement
       
       numEntries = GraphOut[c].ones.size();
       
       column = (int*) malloc((numEntries+2)*sizeof(int));
       column[0] = rank;
       column[1] = numEntries;
       int e = 2;
       for( list<int>::iterator it = GraphOut[c].ones.begin(); e < numEntries+2; e++, it++ ) {
	 column[e] = *it;
	 if( column[e] >= rows.size() ) {
	   printf("(%d) FATAL error %d %d %lu\n", rank, e, column[e], rows.size() ); 
	 }
       }
       //printf("(%d) sending column %d %d %d\n", rank, column[0], column[1], column[2]);
       increment(turn, 0, n_proc);
       globalTurn++;
       // send this column to the next processor 
       // first send the size
       MPI_Request request1, *request2;
       request2 = new MPI_Request();
       MPI_Isend( column+1, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, MPI_COMM_WORLD, &request1 );
       if( numEntries ) {
	 MPI_Isend( column, numEntries+2, MPI_INT, nextProc, TAG_COLUMN, MPI_COMM_WORLD, request2 );
	 outRequests.push_back(request2); // shouldn't be a need to save request1
	 colsToReduceBy.push_back(column);
	 colOfMineOnStack = true;
       }
       
       
       do {
	 c++;
	 if( c < numCols && GraphOut[c].ones.size() ) {
	   // all the columns on the queue
	   for( std::deque<int*>::iterator it = colsToReduceBy.begin(); it != colsToReduceBy.end(); it++ ) {
	     int *columnFromQueue = *it;
	     resolveCols(GraphOut, GraphIn, columnFromQueue+2, columnFromQueue[1], c, firstCol);
	   }
	   // also the column we are currently working on
	   if( indexInCol != indexInColEndValue ) {
	     // this should just deal with it, whether or not there is anything to resolve
	     resolveCols(GraphOut, GraphIn, currentCol+2, currentCol[1], c, firstCol);
	   }
	 }
       } while( c < numCols && (!GraphOut[c].alive || GraphOut[c].ones.size()==0) );
       
       continue;
     }
     
     
     // do some work on the column we are working on
     if( indexInCol != indexInColEndValue || colsToReduceBy.size() ) {
       if( indexInCol == indexInColEndValue ) {
	 // pull a new column out of the queue
	 currentCol = colsToReduceBy.front();
	 currentColOutRequest = outRequests.front();
	 colsToReduceBy.pop_front();
	 outRequests.pop_front();
	 indexInCol = GraphIn[currentCol[2]].ones.begin();
	 indexInColEndValue = GraphIn[currentCol[2]].ones.end();
       }
       const int STEPS_PER_COM = 10; // number of steps to take between communication, unless we finish a column
       for( int count = 0; count < STEPS_PER_COM; count++, indexInCol++ ) {
	 if( indexInCol == indexInColEndValue ) {
	   break;
	   }
	 // indexInCol is what was called j in the serial code
	 if( *indexInCol - firstCol <= c ) {
	   continue; // these columns have already been dealt with
	 }
	 for( int k = 2; k < currentCol[1]+2; k++ ) {
	   list<int>::iterator search = find( GraphOut[*indexInCol - firstCol].ones.begin(), GraphOut[*indexInCol - firstCol].ones.end(), currentCol[k] );
	   if( search != GraphOut[*indexInCol - firstCol].ones.end() ) {
	     GraphOut[*indexInCol - firstCol].ones.erase(search);
	     if( k != 2 )
	       GraphIn[currentCol[k]].ones.remove(*indexInCol);
	   } else {
	     GraphOut[*indexInCol - firstCol].ones.push_back(currentCol[k]);
	     if( k != 2 )
	       GraphIn[currentCol[k]].ones.push_back(*indexInCol);	       
	   }
	 }
       }
       if( indexInCol == indexInColEndValue ) {
	 GraphIn[currentCol[2]].ones.clear();
	 // make sure we have sent on this column before deleting it
	 MPI_Wait(currentColOutRequest, MPI_STATUS_IGNORE);
	 if( currentCol[0] == rank )
	   colOfMineOnStack = false; // because there could be only one and we just finished with it
	 
	 free(currentCol);
       }
       continue;
     }
     
     if (c >= numCols && colsToReduceBy.size() ) {
       // there's nothing to do, but we should clear out colsToReduceBy
       while( colsToReduceBy.size() ) {
	 currentCol = colsToReduceBy.front();
	 currentColOutRequest = outRequests.front();
	 colsToReduceBy.pop_front();
	 outRequests.pop_front();
	 MPI_Wait(currentColOutRequest, MPI_STATUS_IGNORE);
	 delete currentColOutRequest;
	 free(currentCol);
       }
       continue;
     } 
     
     // Nothing useful to do, apparently
     
   }
   printf("(%d) Outside main reduction loop\n", rank);
   
   
   //printf("(%d) imposing barrier\n", rank);
   // Probably don't need this
   //MPI_Barrier(MPI_COMM_WORLD);


   // result of the reduction
   unsigned int localKD = 0;
   for( int col = 0; col < numCols; col++ )
     if( GraphOut[col].ones.size() == 0 || !GraphOut[col].alive )
       localKD++;
     else {
       /*
	 printf("(%d) TAG %d:", rank, col+firstCol);
	 for( list<int>::iterator it = GraphOut[col].ones.begin(); it != GraphOut[col].ones.end(); it++ )
	 printf("%d ", *it);
	 printf("\n");
       */
     }
   unsigned int localID = numCols - localKD;
   
   MPI_Reduce( &localKD, &kernelDimension, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &localID, &imageDimension, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
   if( rank == 0 )
     printf("(%d) TAG gradings %d %d\n", rank, kernelDimension, imageDimension );
   
 }
     
     
  if( rank == 0 ) {
   printf("Time to compute and reduce the matrix %f\n", read_timer()-starttime); 
  }

 // Should now output the results

 if( rank == 0 ) {
   printf("Result: %d %d\n", kernelDimension, imageDimension);
   char outFile[30];
   sprintf(outFile, "%s/bnd%d,%d.dat", saveDir, aGrading, mGrading);
   FILE *f;
   f = fopen(outFile, "w");
   fprintf(f, "%d %d\n", kernelDimension, imageDimension);
   fclose(f);
 }
 MPI_Finalize();
 return 0;
}

// Actual Functions

int NumComp(){
  int nc = 0;
  int c[gridsize];
  int d=0;
  int k;
  for (int i=0; i<gridsize; i++) c[i]=i;
  int dblack[gridsize];
  //int dwhite[gridsize];
  for (int i =0; i<gridsize; i++){
    dblack[black[i]]=i;
    //dwhite[white[i]]=i;
  }
  bool t=0;
  
  while(!t){
    d=0;
    k=0;
    while(c[k]==-1){
      d++;
      k++;
    }
    c[d]=-1;
    int l = dblack[white[d]];
    while(l!=d){
      c[l]=-1;
      l=dblack[white[l]];
    }
    nc++;
    t=1;    
    for (int j=0; j<gridsize; j++) t = t&& (c[j] ==-1);
  }
  
  return nc;
}


int Find(vector<long long> & V, long long x) {
 int above=V.size()-1;
 int below=0;
 while(above - below > 1) {
  if (x >= V[below+(above-below)/2]) below += (above-below)/2;
  else above = below+(above-below)/2;
 }
 if (V[below] == x) return below;
 if (V[above] == x) return above;
 return -1;
}

int WindingNumber(int x, int y){ // Return winding number around (x,y)
 int ret=0;
 for(int i=0; i<x; i++) {
  if ((black[i] >= y) && (white[i] < y)) ret++;
  if ((white[i] >= y) && (black[i] < y)) ret--;
 }
 return ret;
}

int MaslovGrading(int y []) {

 // Use the formula:
 // 4M(y) = 4M(white)+4P_y(R_{y, white})+4P_{white}(R_{y, white})-8W(R_{y, white})
 //       = 4-4*gridsize+4P_y(R_{y, white})+4P_{x_0}(R_{y, white})-8W(R_{y, white})

 int P=4-4*gridsize; // Four times the Maslov grading
 for(int i=0; i<gridsize; i++) {

  // Calculate incidence number R_{y x_0}.S for each of the four
  // squares S having (i,white[i]) as a corner and each of the 
  // four squares having (i,y[i]) as a corner and shift P appropriately

  for(int j=0; j<=i; j++) { // Squares whose BL corners are (i,white[i]) and (i,y[i])
   if ((white[j] > white[i]) && (y[j] <= white[i])) P-=7; // because of the -8W(R_{y, white}) contribution
   if ((y[j] > white[i]) && (white[j] <= white[i])) P+=7; 
   if ((white[j] > y[i]) && (y[j] <= y[i])) P++;
   if ((y[j] > y[i]) && (white[j] <= y[i])) P--; 
  }
  for(int j=0; j<=((i-1)% gridsize); j++) { // Squares whose BR corners are (i,white[i]) and (i,y[i])  (mod gridsize)
   if ((white[j] > white[i]) && (y[j] <= white[i])) P++;
   if ((y[j] > white[i]) && (white[j] <= white[i])) P--; 
   if ((white[j] > y[i]) && (y[j] <= y[i])) P++;
   if ((y[j] > y[i]) && (white[j] <= y[i])) P--; 
  }
  for(int j=0; j<=((i-1) % gridsize); j++) { // Squares whose TR corners are...
   if ((white[j] > ((white[i]-1) % gridsize)) && (y[j] <= ((white[i]-1) % gridsize))) P++;
   if ((y[j] > ((white[i]-1) % gridsize) ) && (white[j] <= ((white[i]-1) % gridsize))) P--; 
   if ((white[j] > ((y[i]-1) % gridsize)) && (y[j] <= ((y[i]-1) % gridsize))) P++;
   if ((y[j] > ((y[i]-1) % gridsize) ) && (white[j] <= ((y[i]-1) % gridsize))) P--; 
  }
  for(int j=0; j<=i; j++) { // Squares whose TL corners are...
   if ((white[j] > ((white[i]-1) % gridsize)) && (y[j] <= ((white[i]-1) % gridsize))) P++;
   if ((y[j] > ((white[i]-1) % gridsize) ) && (white[j] <= ((white[i]-1) % gridsize))) P--; 
   if ((white[j] > ((y[i]-1) % gridsize)) && (y[j] <= ((y[i]-1) % gridsize))) P++;
   if ((y[j] > ((y[i]-1) % gridsize) ) && (white[j] <= ((y[i]-1) % gridsize))) P--; 
  }
 }
 return (P/4);
}

bool RectDotFree(int xll, int yll, int xur, int yur, int which) {
 bool dotfree = 1;
 switch (which) {
  case 0: 
   for(int x=xll; x<xur && dotfree; x++) {
    if (white[x] >= yll && white[x] < yur) dotfree = 0;
    if (black[x] >= yll && black[x] < yur) dotfree = 0;
   }
   return dotfree;
  case 1:
   for(int x=0; x<xll && dotfree; x++) {
    if (white[x] < yll || white[x] >= yur) dotfree = 0;
    if (black[x] < yll || black[x] >= yur) dotfree = 0;
   }
   for(int x=xur; x<gridsize && dotfree; x++) {
    if (white[x] < yll || white[x] >= yur) dotfree = 0;
    if (black[x] < yll || black[x] >= yur) dotfree = 0;
   }
   return dotfree;
  case 2:
   for(int x=xll; x<xur && dotfree; x++) {
    if (white[x] < yll || white[x] >= yur) dotfree = 0;
    if (black[x] < yll || black[x] >= yur) dotfree = 0;
   }
   return dotfree;
  case 3:
   for(int x=0; x<xll && dotfree; x++) {
    if (white[x] >= yll && white[x] < yur) dotfree = 0;
    if (black[x] >= yll && black[x] < yur) dotfree = 0;
   }
   for(int x=xur; x<gridsize && dotfree; x++) {
    if (white[x] >= yll && white[x] < yur) dotfree = 0;
    if (black[x] >= yll && black[x] < yur) dotfree = 0;
   }
   return dotfree;
 }
 return 0; //Error!
}

bool ValidGrid() {
 int numwhite=0;
 int numblack=0;
 for(int i=0; i<gridsize; i++) {
  for(int j=0; j<gridsize; j++) {
   if (white[j]==i) numwhite++;
   if (black[j]==i) numblack++;
  }
  if (numwhite != 1 || numblack != 1) {
   std::cout << "\nInvalid Grid!\n";
   return 0;
  }
  numwhite=0;
  numblack=0;
 }
 return 1;
}

long long getIndex( int *P ) {
  long long index = 0;
  for( int i = gridsize-2; i >= 0; i-- ) {
    int r = P[i];
    int m = 0;
    for( int j = 0; j < i; j++ )
      if( P[j] < r )
	m++;
    index += Factorial[gridsize-1-i]*(r-m);
  }
  return index;
}

// Inverse mapping, from integers < gridsize! to permutations of size n
// Writes the permutation corresponding to n into the array P.
void getPerm( long long n, int *P ) {
  int taken[gridsize];
  int offset;
  for( int i = 0; i < gridsize; i++ )
    taken[i] = 0;
  for( int i = 0; i < gridsize; i++ ) {
    offset = n / Factorial[gridsize-1-i];
    n -= offset*Factorial[gridsize-1-i];
    for( int j = 0; j <= offset; j++ )
      if( taken[j] )
	offset++;
    P[i] = offset;
    taken[P[i]] = 1;
  }
}
