// essentially a serial version of boundaries-bc
#include <sys/time.h>
#include <iostream>
#include <stdlib.h>
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

//#define PROFILE

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

//Globals 
int gridsize = 12; // arc-index
int *white;
int *black;
int default_white[12] = {9,5,11,7,8,1,10,4,0,3,2,6};
int default_black[12] = {1,0,4,3,2,6,5,9,8,11,7,10};

long long *Factorial;

// Function Prototypes
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
void getPerm( long long n, int *P);
void NextPerm(short counter[], int h[]);
bool RectDotFree(int xll, int yll, int xur, int yur, int which); 
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
int Find(vector<long long> & V, long long x); // Returns i if V[i]=x or -1 if x isn't V[i] for any i, assuming V is sorted

// Main

int main(int argc, char *argv[]){

 // set up mpi
 int n_proc=1, rank=0;
 
 char *saveDir = read_string( argc, argv, "-i", NULL );
 bool justTime = read_int( argc, argv, "-t", 0 );
 int aGrading = read_int( argc, argv, "-a", -100 );
 int mGrading = read_int( argc, argv, "-m", -100 );
 int BLOCKSIZE = read_int( argc, argv, "-b", 10 );
 int maxQueueSize = read_int( argc, argv, "-q", n_proc );
 int AHEAD_LIMIT = read_int( argc, argv, "-l", n_proc/2 );
 int MAX_ROWS_PER_CYCLE_INPUT = read_int( argc, argv, "-r", 10 );
 int MAX_ROWS_PER_CYCLE = MAX_ROWS_PER_CYCLE_INPUT;

 if( aGrading == -100 || mGrading == -100 ) {
   printf("Must specify in which grading to calculate the boundary image and kernel dimensions with -a and -m\n");
   exit(0);
 }
 if( !saveDir ) {
   printf("Use -i to specify where the output of gradings was written\n");
   exit(0);
 }

 char knotFile[30];
 sprintf(knotFile, "%s/knot", saveDir);
 FILE *f = fopen(knotFile, "r");
 if( !f ) {
   printf("Error opening file %s\n", knotFile);
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

 Factorial = (long long *) malloc( gridsize * sizeof( long long ) );
 Factorial[0] = 1;
 for( int i = 1; i < gridsize; i++ )
   Factorial[i] = i*Factorial[i-1];
   

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
   if( !justTime )
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
     if( !justTime )
       printf("warning, no such file %s\n", rowFile);
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
 
 int imageDimension, kernelDimension;

 if( cols.size() == 0 || rows.size() == 0 ) {
   imageDimension = 0;
   kernelDimension = cols.size();
 } else {
   // Decide who has ownership of what subset of the generators (without communication).
   // The division will be block cyclic, with size set by BLOCKSIZE
   int num_generators = cols.size();
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
   vector<Generator> GraphIn( rows.size() );
   vector<Generator> GraphOut( numCols );
   int globalIndex = BLOCKSIZE*rank-BLOCKSIZE*(n_proc-1)-1;
   for(int index=0; index < numCols; index++) {
     // update globalIndex
     if( index % BLOCKSIZE == 0 )
       globalIndex += BLOCKSIZE*(n_proc-1);
     globalIndex++;
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
#else
	   GraphOut[index].ones.push_back( indexgij );
	   GraphIn[indexgij].ones.push_back(index);
#endif
	 }
       }
     }
   }
   if( !justTime )
     printf("(%d) Matrix filled\n", rank);

   //Unlike previous versions, send a block at a time.  The format of a message is: <origProc> <col1Size> <col1Entry1> ... <col1lastEntry> <col2size> ...
#define TAG_COLUMN 11
#define TAG_FINISHED 12

#ifdef PROFILE
   int colsAhead = 0;
   int numBlocksSkipped = 0;
   double totalDeleteWait = 0;
   double waitTime = 0;
#endif
   int nextProc = rank + 1;
   if( nextProc == n_proc )
     nextProc = 0;
   int prevProc = rank - 1;
   if( prevProc < 0 )
     prevProc = n_proc - 1;
   
   // These variable store the state of the reduction
   int turn = 0; // the id of the matrix whose turn it is
   int globalTurn = 0; // This just keeps increasing.  At any time, turn=globalTurn%n_proc
   std::deque<int*> colsToReduceBy;
   int finishedSignals[2*n_proc];
   int *finishedInSignal = finishedSignals+2;
   int indexInCurrentCols=-1;
   int colInCurrentCols;
   int *currentCols;
   list<int>::iterator indexInCol;
   list<int>::iterator indexInColEndValue = indexInCol;
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
   bool allOthersFinished = false;
   
   while( numProcReportingFinished < n_proc ) { // main loop for the reduction.  The key principle is that we only do one thing per loop iteration
#ifdef PROFILE
     //       totalCycles++;
#endif
     // make sure turn is someone who is actually still working
     while( procsReportingFinished[turn] != -1 && procsReportingFinished[turn] <= globalTurn ) {
       increment(turn, 0, n_proc);
       globalTurn++;
     }
     // this line will produce a lot of output
     //printf("(%d) in reduction, block = %d/%d, turn = %d(%d), numfinished = %d, selfFinished = %d(%d), colsToReduceBy.size = %lu, outRequests.size() = %lu\n", rank, block, numBlocks, turn, globalTurn, numProcReportingFinished, selfFinished, procsReportingFinished[rank], colsToReduceBy.size(), outRequests.size());

     if( turn == rank && // it is our turn to send 
	 colsToReduceBy.size() < maxQueueSize // but don't let the queue get too big
	 ){
       // resolve a block of columns with eachother, (fully? there is no need to do it fully, unless we are going to have multiple threads work on different ones at once)
       /* For full reduction, use these lines
       //for( int source = block*BLOCKSIZE; source < block*BLOCKSIZE+currentBlockSize; source++ )
       //for( int target = block*BLOCKSIZE; target < block*BLOCKSIZE+currentBlockSize; target++ )
       //if( source != target )
       */
       for( int source = block*BLOCKSIZE; source < block*BLOCKSIZE+currentBlockSize-1; source++ )
	 for( int target = source+1; target < block*BLOCKSIZE+currentBlockSize; target++ )
	   resolveColsInternal(GraphOut, GraphIn, source, target);

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
	     for( list<int>::iterator it = GraphOut[source].ones.begin(); it != GraphOut[source].ones.end(); it++, e++ )
	       columns[e] = *it;
	   }
	 }
	 
	 colsToReduceBy.push_back(columns);
       } else { // just send out the empty message, but don't put in on the queue
       }
       increment(turn, 0, n_proc);
       globalTurn++;

       block++;
       if( block == numFullBlocks )
	 currentBlockSize = tailSize;
       // try skipping over some zeros
       bool allZero = true;
       while( allZero && block < numBlocks ) {
	 for( int source = block*BLOCKSIZE; source < block*BLOCKSIZE+currentBlockSize; source++ )
	   if( GraphOut[source].ones.size() ) {
	     allZero = false;
	     break;
	   }
	 if( allZero ) {
#ifdef PROFILE
	   //printf("(%d) skipping over a block that is all zero %d/%d (%d)\n", rank, block, numBlocks, numFullBlocks);
	   numBlocksSkipped++;
#endif
	   block++;
	   if( block == numFullBlocks )
	     currentBlockSize = tailSize;
	 }
       
       
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
       }
       if( block >= numBlocks ) { // we are done, let the other processes know
	   selfFinished = true;
	   procsReportingFinished[rank] = globalTurn;
	   numProcReportingFinished++;
       }
       continue;
     }

     // expedited reduction if our queue is too big
     MAX_ROWS_PER_CYCLE = MAX_ROWS_PER_CYCLE_INPUT;
     if( colsToReduceBy.size() >= maxQueueSize ) {
       MAX_ROWS_PER_CYCLE = 1000000; // really I just want infininty here
     }     

     // do some work on the columns we are working on
     if( colsToReduceBy.size() ) {
       if( indexInCurrentCols == -1 ) {
	 currentCols = colsToReduceBy.front();
	 colInCurrentCols = 0;
	 indexInCurrentCols = 2;
       }
       
       // this is the code to initialize a new column, still in currentCols
       if( indexInCol == indexInColEndValue ) {
	 indexInCol = GraphIn[currentCols[indexInCurrentCols+1]].ones.begin();
	 indexInColEndValue = GraphIn[currentCols[indexInCurrentCols+1]].ones.end();
       }

       //for( list<int>::iterator it = GraphIn[currentCols[indexInCurrentCols+1]].ones.begin(); it != GraphIn[currentCols[indexInCurrentCols+1]].ones.end(); it++ ) {
       for( int count = 0; count < MAX_ROWS_PER_CYCLE && indexInCol != indexInColEndValue; count++, indexInCol++ ) {
	 // *indexInCol is what was called j in the serial code
	 if( *indexInCol < block*BLOCKSIZE+currentBlockSize ) {
	   continue; // these columns have already been dealt with
	 }
	 for( int k = indexInCurrentCols+1; k < indexInCurrentCols+1+currentCols[indexInCurrentCols]; k++ ) {
	   list<int>::iterator search = find( GraphOut[*indexInCol].ones.begin(), GraphOut[*indexInCol].ones.end(), currentCols[k] );
	   if( search != GraphOut[*indexInCol].ones.end() ) {
	     GraphOut[*indexInCol].ones.erase(search);
	     if( k != indexInCurrentCols+1 )
	       GraphIn[currentCols[k]].ones.remove(*indexInCol);
	   } else {
	     GraphOut[*indexInCol].ones.push_back(currentCols[k]);
	     if( k != indexInCurrentCols+1 )
	       GraphIn[currentCols[k]].ones.push_back(*indexInCol);	       
	   }
	 }
       }
       if( indexInCol == indexInColEndValue ) { // we have finished with a single column
	 GraphIn[currentCols[indexInCurrentCols+1]].ones.clear();
	 colInCurrentCols++;
	 indexInCurrentCols += currentCols[indexInCurrentCols]+1;
       }
       if( colInCurrentCols >= currentCols[1] ) { // we have finished with an entire block of columns
	 // make sure we have sent on this column before deleting it
#ifdef PROFILE
	 double deleteWait = read_timer();
#endif
#ifdef PROFILE
	 totalDeleteWait += read_timer()-deleteWait;
#endif
	 
	 free(currentCols);
	 indexInCurrentCols = -1;
	 colsToReduceBy.pop_front();
       }
       //printf("(%d) work cycle finished\n", rank);
       continue;
     }

#ifdef PROFILE
     double waitStart = read_timer();
#endif
#ifdef PROFILE
     waitTime += read_timer()-waitStart;
#endif

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
	 for( list<int>::iterator it = GraphOut[col].ones.begin(); it != GraphOut[col].ones.end(); it++ )
	 printf("%d ", *it);
	 printf("\n");
       */
     }
   unsigned int localID = numCols - localKD;
   
   kernelDimension = localKD;
   imageDimension = localID;
   if( rank == 0 )
     printf("(%d) TAG gradings kdim=%d rank=%d\n", rank, kernelDimension, imageDimension );
   
#ifdef PROFILE
   printf("(%d) number of columns done out of order %d, number of empty blocks skipped %d, time spent waiting to delete %f, time spent waiting for messages %f\n", rank, colsAhead, numBlocksSkipped, totalDeleteWait, waitTime);
   /*
   printf("(%d) profile results, spent %d of %d cycles waiting", rank, totalWaitingCycles, totalCycles);
   printf("(%d) significant ones are:\n", rank);
   for( int i = 0; i < numBlocks+1; i++ )
     if( waitingCycles[i] > totalWaitingCycles/100 )
       printf("(%d) %d : %d\n", rank, i, waitingCycles[i]);
   */
#endif
   if( rank == 0 ) {
     printf("Time to compute and reduce the matrix %f\n", read_timer()-starttime); 
   }
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
