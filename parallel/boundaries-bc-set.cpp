#include <sys/time.h>
#include <iostream>
#include <stdlib.h>
#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <queue>
#include <set>
#include <assert.h>
#include <algorithm>
#include "matrix-z2z3-set.h"
#include "reduction-functions-set.h"

#include <unistd.h>

using std::set;
using std::pair;
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

 int  numcomp = NumComp();
 if( rank == 0 && !justTime )
   printf("Number of components: %d\n", numcomp);

 if( rank == 0 )
   if(!ValidGrid()) {printf("Invalid grid!!\n"); return 0;} // Check that the grid is valid

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
   if( !justTime )
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
     if( !justTime )
       printf("waring, no such file %s\n", colFile);
   }
 }

 double starttime;
 if( rank == 0 )
   starttime = read_timer(); // Used to record how long this takes


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
 

 int imageDimension, kernelDimension;

 int BLOCKSIZE = 10;

 if( cols.size() == 0 || rows.size() == 0 ) {
   imageDimension = 0;
   kernelDimension = cols.size();
 } else {
   // Decide who has ownership of what subset of the generators (without communication).
   // The division will be block cyclic, with size set by BLOCKSIZE
   int blocksPerProcMin = num_generators/BLOCKSIZE/n_proc;
   int extraBlocks = (num_generators - n_proc * blocksPerProcMin * BLOCKSIZE)/BLOCKSIZE;
   int extraCols = num_generators - n_proc * blocksPerProcMin * BLOCKSIZE - extraBlocks*BLOCKSIZE;
   int numFullBlocks, tailSize, numBlocks; // the number of complete blocks, the size of the last block (usually zero), and the total number of blocks including that one
   if( rank < extraBlocks ) {
     numFullBlocks = blocksPerProcMin+1;
     tailSize = 0;
   } else if( rank == extraBlocks ) {
     numFullBlocks = blocksPerProcMin;
     tailSize = extraCols;
   } else {
     numFullBlocks = blocksPerProcMin;
     tailSize = 0;
   }
   if( tailSize )
     numBlocks = numFullBlocks+1;
   else
     numBlocks = numFullBlocks;
   int numCols = numFullBlocks*BLOCKSIZE + tailSize;
   
   if( !justTime )
     printf("(%d) ownership determined, I have %d blocks of size %d and %d extra columns at the end for a total of %d columns\n", rank, numFullBlocks, BLOCKSIZE, tailSize, numCols);

   // Each proc calculates the part of the matrix that it owns, storing both rows and columns as in the serial code.  It holds entire columns, but only certain rows
   vector<Generator> GraphIn( rows.size() ); // Will hold boundary data.
   vector<Generator> GraphOut( numCols ); // Will hold boundary data.
   //printf("Filling %d %d\n", I, J);
   int globalIndex = BLOCKSIZE*rank-BLOCKSIZE*(n_proc-1)-1;
   for(int index=0; index < numCols; index++) {
     // update globalIndex
     if( index % BLOCKSIZE == 0 )
       globalIndex += BLOCKSIZE*(n_proc-1);
     globalIndex++;
     //printf("(%d) li=%d gi=%d\n", rank, index, globalIndex);
     getPerm(cols[globalIndex],g);
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
	   GraphOut[index].ones.insert( indexgij );
	   GraphIn[indexgij].ones.insert( GraphIn[indexgij].ones.end(), index);
	   //printf("filling entry %d %d\n", globalIndex+1, indexgij+1);
#endif
	 }
       }
     }
   }
   if( !justTime )
     printf("(%d) Matrix filled\n", rank);

   //char fileName[10];
   //sprintf(fileName, "outmat%d", rank); 
   //printMatrix(fileName, GraphOut);

   // this skips over the first matrix in 10_125
   //if( i == 30 && j == 26 )
   //  continue;
   
   /* 
      Unlike previous versions, send a block at a time.  The format of a message is: <origProc> <col1Size> <col1Entry1> ... <col1lastEntry> <col2size> ...
      Everything below this line needs to be updated
   */
#define TAG_COLUMN_SIZE 10
#define TAG_COLUMN 11
#define TAG_FINISHED 12

   int turn = 0; // the id of the matrix whose turn it is
   int globalTurn = 0; // This just keeps increasing
   std::deque<int*> colsToReduceBy;
   std::deque<MPI_Request*> outRequests; // keep these so we can make sure we have sent a given column when we dequeue it befor we delete it
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
   //bool starting = true;
   int indexInCol=-1;
   int colInCol;
   int *currentCol;
   MPI_Request *currentColOutRequest;
   int emptyCol[2];
   emptyCol[0] = rank;
   emptyCol[1] = 0;
   int block = 0;
   int currentBlockSize = BLOCKSIZE;
   if( numFullBlocks == 0 )
     currentBlockSize = tailSize;
   int numProcReportingFinished = 0;
   int procsReportingFinished[n_proc]; // this stores the globalTurn at which each proc has reported it is finished, or -1 if it has not reported finished
   for( int p = 0; p < n_proc; p++ )
     procsReportingFinished[p] = -1;
   bool selfFinished = false;
   
   while( numProcReportingFinished < n_proc ) { // main loop for the reduction.  The key principle is that we only do one thing per loop iteration

     // make sure turn is someone who is actually still working
     while( procsReportingFinished[turn] != -1 && procsReportingFinished[turn] <= globalTurn ) {
       increment(turn, 0, n_proc);
       globalTurn++;
     }
     //printf("(%d) in reduction, block = %d/%d, turn = %d(%d), numfinished = %d, selfFinished = %d(%d), colsToReduceBy.size = %lu, outRequests.size() = %lu, colOfMine = %d\n", rank, block, numBlocks, turn, globalTurn, numProcReportingFinished, selfFinished, procsReportingFinished[rank], colsToReduceBy.size(), outRequests.size(), colOfMineOnStack);

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
	 procsReportingFinished[finishedInSignal[0]] = finishedInSignal[1];
	 numProcReportingFinished++;
	 if( finishedInSignal[0] != nextProc ) {
	   finishedOutRequest = new MPI_Request();
	   MPI_Isend( finishedInSignal, 2, MPI_INT, nextProc, TAG_FINISHED, MPI_COMM_WORLD, finishedOutRequest );
	 } else {
	   // put out a new request if there are any procs left not reporting done
	   if( numProcReportingFinished < n_proc-1 ||
	       numProcReportingFinished == n_proc-1 && selfFinished ) {
	     finishedInRequest = new MPI_Request();
	     MPI_Irecv( &finishedInSignal, 2, MPI_INT, prevProc, TAG_FINISHED, MPI_COMM_WORLD, finishedInRequest );
	   }	   
	 }
	 continue;
       }
     }
     // check if we ourselves are done
     if( !selfFinished && (block >= numBlocks) ) {
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
     
     // if it is not our turn, see if there is a message to pick up
     if( turn != rank ) {
       MPI_Status stat;
       int flag;
       MPI_Iprobe( prevProc, TAG_COLUMN, MPI_COMM_WORLD, &flag, &stat);
       if( flag ) {
	 int inSize;
	 MPI_Get_count(&stat, MPI_INT, &inSize);
	 int *inCols;
	 inCols = (int*) malloc( inSize*sizeof(int) );
	 MPI_Recv( inCols, inSize, MPI_INT, prevProc, TAG_COLUMN, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	 /*
	 printf("(%d) received a block: ", rank);
	 for( int l = 0; l < inSize; l++ )
	   printf("%d ", inCols[l]);
	 printf("\n");
	 printf("(%d) received while block*BLOCKSIZE=%d\n", rank, block*BLOCKSIZE);
	 */
	 globalTurn += inCols[0] - turn;
	 if( inCols[0] < turn )
	   globalTurn += n_proc;
	 turn = inCols[0];
	 increment(turn, 0, n_proc);
	 globalTurn++;

	 if( inCols[1] ) { // the new columns are not zero, add it to the queue
	   colsToReduceBy.push_back(inCols);
	   if( inCols[0] != nextProc ) { //pass on, unless it has made it all the way around
	     MPI_Request *outReq;
	     outReq = new MPI_Request();
	     MPI_Isend( colsToReduceBy.back(), inSize, MPI_INT, nextProc, TAG_COLUMN, MPI_COMM_WORLD, outReq );
	     outRequests.push_back(outReq);
	   } else {
	     // need to keep outRequests in sync with colsToReduceBy
	     MPI_Request rnull = MPI_REQUEST_NULL;
	     outRequests.push_back(&rnull);
	   }
	   
	   if( !selfFinished ) {
	     // resolve the newly received columns with our next set
	     int colStartIndex = 2;
	     for( int col = 0; col < inCols[1]; col++, colStartIndex += inCols[colStartIndex]+1 )
	       for( int target = block*BLOCKSIZE; target < block*BLOCKSIZE+currentBlockSize; target++ )
		 resolveCols(GraphOut, GraphIn, inCols+colStartIndex+1, inCols[colStartIndex], target);
	   }
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
       // resolve a block of columns with eachother, fully (there is no need to do it fully, unless we are going to have multiple threads work on different ones at once later
       /*
       //debug code
       printf("(%d) about to fully reduce a block starting at %d\n", rank, block*BLOCKSIZE);
       printf("(%d) starting at: ", rank);
       for( int source=block*BLOCKSIZE; source < block*BLOCKSIZE+currentBlockSize; source++ ) {
	 int size = GraphOut[source].ones.size();
	 if( size ) {
	   printf("%d ", size);
	     for( list<int>::iterator it = GraphOut[source].ones.begin(); it != GraphOut[source].ones.end(); it++ )
	       printf("%d ", *it);
	   }
       }
       printf("\n");
       // end debug code
       */

       for( int source = block*BLOCKSIZE; source < block*BLOCKSIZE+currentBlockSize; source++ )
	 for( int target = block*BLOCKSIZE; target < block*BLOCKSIZE+currentBlockSize; target++ )
	   if( source != target ) // for full reduction
	   //if( source < target ) // for partial reduction
	     resolveColsInternal(GraphOut, GraphIn, source, target);
       /*
       //debug code
       printf("(%d) finished with: ", rank);
       for( int source=block*BLOCKSIZE; source < block*BLOCKSIZE+currentBlockSize; source++ ) {
	 int size = GraphOut[source].ones.size();
	 if( size ) {
	   printf("%d ", size);
	     for( list<int>::iterator it = GraphOut[source].ones.begin(); it != GraphOut[source].ones.end(); it++ )
	       printf("%d ", *it);
	   }
       }
       printf("\n");
       // end debug code
       */

       // assemble them to send       
       int numColumns=0, *columns, totalEntries=0;
       for( int source = block*BLOCKSIZE; source < block*BLOCKSIZE+currentBlockSize; source++ ) {
	 int size = GraphOut[source].ones.size();
	 if( size ) {
	   totalEntries += size;
	   numColumns++;
	 }
       }

       if( numColumns ) { // if the message isn't all zeros
	 totalEntries += 2 + numColumns; // this is the size of the message
	 columns = (int*) malloc( totalEntries*sizeof(int) );
	 columns[0] = rank;
	 columns[1] = numColumns;
	 for( int source=block*BLOCKSIZE, e=2; source < block*BLOCKSIZE+currentBlockSize; source++ ) {
	   int size = GraphOut[source].ones.size();
	   if( size ) {
	     columns[e] = size;
	     e++;
	     for( set<int>::iterator it = GraphOut[source].ones.begin(); it != GraphOut[source].ones.end(); it++, e++ )
	       columns[e] = *it;
	   }
	 }
	 
	 //printf("(%d) sending column %d %d %d\n", rank, column[0], column[1], column[2]);
	 // send this column to the next processor 
	 MPI_Request *req;
	 req = new MPI_Request();
	 MPI_Isend( columns, totalEntries, MPI_INT, nextProc, TAG_COLUMN, MPI_COMM_WORLD, req );
	 outRequests.push_back(req); 
	 colsToReduceBy.push_back(columns);
	 colOfMineOnStack = true;
       } else { // just send out the empty message, but don't put in on the queue
	 MPI_Request req;
	 MPI_Isend( emptyCol, 2, MPI_INT, nextProc, TAG_COLUMN, MPI_COMM_WORLD, &req );
       }
       increment(turn, 0, n_proc);
       globalTurn++;

       block++;
       if( block == numFullBlocks )
	 currentBlockSize = tailSize;
       
       // resolve the new block with all the columns on the queue
       if( block < numBlocks ) {
	 for( std::deque<int*>::iterator it = colsToReduceBy.begin(); it != colsToReduceBy.end(); it++ ) {
	   int *columnsFromQueue = *it;
	   int colStartIndex = 2;
	   for( int col = 0; col < columnsFromQueue[1]; col++, colStartIndex += columnsFromQueue[colStartIndex]+1 )
	     for( int target = block*BLOCKSIZE; target < block*BLOCKSIZE+currentBlockSize; target++ )
	       resolveCols(GraphOut, GraphIn, columnsFromQueue+colStartIndex+1, columnsFromQueue[colStartIndex], target);
	 }
       }
       continue;
 }
     
     
     // do some work on the columns we are working on
     if( colsToReduceBy.size() ) {
       //printf("(%d) beginning work cycle\n", rank);
       if( indexInCol == -1 ) {
	 //printf("(%d) getting a new set of columns out of the queue\n", rank);
	 // pull a new column out of the queue
	 currentCol = colsToReduceBy.front();
	 currentColOutRequest = outRequests.front();
	 colInCol = 0;
	 indexInCol = 2; // this is very different from what it used to mean
	 //printf("(%d) done getting a new set of columns out of the queue\n", rank);
       }
       //printf("(%d) currentCol[0]=%d currentCol[1]=%d\n", rank, currentCol[0], currentCol[1]);
       for( set<int>::iterator it = GraphIn[currentCol[indexInCol+1]].ones.lower_bound(block*BLOCKSIZE+currentBlockSize); it != GraphIn[currentCol[indexInCol+1]].ones.end(); it++ ) {
	 //printf("(%d) %d %d %lu\n", rank, *it, currentCol[indexInCol+1], GraphIn[currentCol[indexInCol+1]].ones.size());
	 // *it is what was called j in the serial code
	 if( *it < block*BLOCKSIZE+currentBlockSize ) {
	   assert(false); // the lower_bound above should take care of these
	 }
	 for( int k = indexInCol+1; k < indexInCol+1+currentCol[indexInCol]; k++ ) {
	   //printf("(%d) k=%d indexInCol=%d currentCol[indexInCol]=%d, currentCol[1]=%d\n", rank, k, indexInCol, currentCol[indexInCol], currentCol[1]);
	   pair<set<int>::iterator,set<int>::iterator> search = GraphOut[*it].ones.equal_range(currentCol[k]);
	   if( search.first != search.second ) {
	     GraphOut[*it].ones.erase(search.first);
	     if( k != indexInCol+1 )
	       GraphIn[currentCol[k]].ones.erase(*it);
	   } else {
	     GraphOut[*it].ones.insert(search.first, currentCol[k]);
	     if( k != indexInCol+1 )
	       GraphIn[currentCol[k]].ones.insert(*it);	       
	   }
	 }
       }
       GraphIn[currentCol[indexInCol+1]].ones.clear();
       colInCol++;
       indexInCol += currentCol[indexInCol]+1;
       if( colInCol >= currentCol[1] ) { // this is the end condition
	 // make sure we have sent on this column before deleting it
	 MPI_Wait(currentColOutRequest, MPI_STATUS_IGNORE);
	 if( currentCol[0] == rank )
	   colOfMineOnStack = false; // because there could be only one and we just finished with it
	 
	 free(currentCol);
	 indexInCol = -1;
	 colsToReduceBy.pop_front();
	 outRequests.pop_front();
       }
       //printf("(%d) work cycle finished\n", rank);
       continue;
     }

     /*
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
	 colsToReduceBy.pop_front();
	 outRequests.pop_front();
       }
       continue;
     } 
     */
     // Nothing useful to do, apparently
     
   }
   if( !justTime )
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
	 printf("(%d) TAG %d:", rank, col);
	 for( set<int>::iterator it = GraphOut[col].ones.begin(); it != GraphOut[col].ones.end(); it++ )
	 printf("%d ", *it);
	 printf("\n");
       */
     }
   unsigned int localID = numCols - localKD;
   
   MPI_Reduce( &localKD, &kernelDimension, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &localID, &imageDimension, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
   if( rank == 0 )
     printf("(%d) TAG gradings kdim=%d rank=%d\n", rank, kernelDimension, imageDimension );
   
 }
     
     
  if( rank == 0 ) {
   printf("Time to compute and reduce the matrix %f\n", read_timer()-starttime); 
  }

 // Should now output the results

 if( rank == 0 && !justTime ) {
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
