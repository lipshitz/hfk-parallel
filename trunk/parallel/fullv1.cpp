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
 
 char *knotFile = read_string( argc, argv, "-k", NULL );
 bool printMatrices = read_int( argc, argv, "-p", 0 );
 bool justTime = read_int( argc, argv, "-t", 0 );
 if( knotFile ) {
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
 } else {
   white = default_white;
   black = default_black;
 }

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
 if( rank == 0 )
   double starttime = read_timer(); // Used to record how long this takes

 // Record winding numbers around grid points for Alexander grading computations
 // Also record the Alexander grading shift
 // Want this to be non-positive; if it isn't, exchange the white and black dots,
 // multiply all winding numbers by -1 and try again.
 int WN[gridsize][gridsize];
 for(int x=0; x<gridsize; x++) {
  for(int y=0; y<gridsize; y++) {
   WN[x][y] = WindingNumber(x,y);
  }
 }
 // Record the lowest winding number per row
 int lowestWN[gridsize];
 for( int x = 0; x < gridsize; x++ ) {
   int lowest = 100;
   for( int y = 0; y < gridsize; y++ )
     lowest = min(WN[x][y],lowest);
   lowestWN[x] = lowest;
 }

 int temp=0;
 for(int i=0; i<gridsize; i++) {
  temp += WN[i][black[i]];
  temp += WN[i][(black[i]+1) % gridsize];
  temp += WN[(i+1) % gridsize][black[i]];
  temp += WN[(i+1) % gridsize][(black[i]+1) % gridsize];
 } 
 for(int i=0; i<gridsize; i++) {
  temp += WN[i][white[i]];
  temp += WN[i][(white[i]+1) % gridsize];
  temp += WN[(i+1) % gridsize][white[i]];
  temp += WN[(i+1) % gridsize][(white[i]+1) % gridsize];
 }
 
 const int AShift = (temp - 4 * gridsize + 4)/8;
 
 if( rank == 0 && !justTime ) {
   printf("Alexander Grading Shift: %d\n", AShift);
   printf("Matrix of winding numbers and Black/White grid:\n");
   for(int y=gridsize-1; y>=0; y--) {
     for(int x=0; x<gridsize; x++) {
       printf("%2d", WN[x][y]);
     }
     printf("   ");
     for(int x=0; x<gridsize; x++) {
       printf(" ");
       if(black[x]==y) printf("X");
       if(white[x]==y) printf("O");
       if(black[x] != y && white[x] != y) printf(" ");
     }
     printf("\n");
   }
 } 

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

 // Iterate through the generators in lexicographic
 // order and calculate their boundaries
 // Identify each permutation with the integer given by the number of permutations preceding
 // it in lexicographic order.
 // Populate Graph[count].out with a list of integers corresponding to the permutations that
 // are boundaries of the permutation corresponding to count.
 // Populate Graph[count].

 
 // This will hold the generators, sorted by grading, the first index is Agrading+30, the second is Maslov-grading+30
 std::vector<long long> *generators[60][60];
 unsigned int imageDimensions[60][60];
 unsigned int kernelDimensions[60][60];
 int homologyDimensions[60][60];
 for( int i = 0; i < 60; i++ )
   for( int j = 0; j < 60; j++ )
     generators[i][j] = new std::vector<long long>();
 
 if( rank == 0 && !justTime )
   printf("Searching through %lld generators to compute Alexander gradings...\n", Factorial[gridsize]);
 double agStartTime = read_timer();
 double agEndTime;

 if( n_proc == 1 ) {
 int g[gridsize];
 int taken[gridsize];
 int istack[gridsize];
 for( int i = 0; i < gridsize; i++ )
   istack[i] = 0;
 for( int i = 0; i < gridsize; i++ )
   taken[i] = 0;
 int depth = 0; // this is the index of g we are working on
 int AGrading = AShift;
 while( true ) {
   if ( depth == gridsize ) { // we are at the end of the recursion, use the permutation
     if (AGrading >= amin && AGrading <= amax) {
       int MGrading = MaslovGrading(g);
       generators[AGrading+30][MGrading+30]->push_back(getIndex(g));
       //NumGenByAGrading[AGrading+30]++;
     }
     depth--;
     if( depth < 0 )
       break;
     taken[g[depth]] = 0;
     AGrading += WN[depth][g[depth]];
     continue;
   }

   // If there is no hope of getting a non-negative Alexander Grading from here on, decrease the depth
   int maxAGrading = AGrading;
   for( int i = depth; i < gridsize; i++ ) {
     maxAGrading -= lowestWN[i];
   }
   if( maxAGrading < amin ) {
     depth--;
     if( depth < 0 )
       break;
     taken[g[depth]] = 0;
     AGrading += WN[depth][g[depth]];
     continue;
   }

   bool callBack = true;
   for( int i = istack[depth]; i < gridsize; i++ ) { // consider using value i for position depth
     if( !taken[i] ) {
       g[depth] = i;
       taken[i] = 1;
       // push onto the stack
       callBack = false;
       istack[depth] = i+1;
       if( depth < gridsize-1 )
	 istack[depth+1] = 0;
       AGrading -= WN[depth][g[depth]];
       depth++;
       break;
     }
   }
   if( callBack ) {
     depth--;
     if( depth < 0 )
       break;
     taken[g[depth]] = 0;
     AGrading += WN[depth][g[depth]];
   }
 }


 } else {

 const int min_depth = 2; // This is the depth from which the master hands out packets
 int num_generators[60*60];
 if( rank == 0 ) {
   // Assign the partitions and collect the results
   long long prefix_max = 1;
   for( int i = gridsize; i > gridsize-min_depth; i-- )
     prefix_max *= i;
   for( long long i = n_proc-1; i < prefix_max+n_proc-1; i++ ) { // the processors will figure out their initial assignments on their own
     long long permnum = i * Factorial[gridsize-min_depth];
     if( i >= prefix_max )
       permnum = -1;
     MPI_Status status;
     // Wait for a thread to ask for more
     MPI_Recv( num_generators, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
     /* This deals with collecting the data, but let's put that off until after the calculation
     MPI_Recv( num_generators, 60*60, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
     for( int i = 0; i < 60; i++ )
       for( int j = 0; j < 60; j++ )
	 if( num_generators[i*60+j] ) {
	   long long buffer[num_generators[i*60+j]];
	   MPI_Recv( buffer, num_generators[i*60+j], MPI_LONG_LONG_INT, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	   generators[i][j]->reserve(generators[i][j]->size()+num_generators[i*60+j]);
	   for( int k = 0; k < num_generators[i*60+j]; k++ ) {
	     generators[i][j]->push_back(buffer[k]);	 
	   }
	 }
     */
     MPI_Send( &permnum, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD );
   }
 } else {
   int g[gridsize];
   long long startPoint = (rank-1)*Factorial[gridsize-min_depth];
   while( startPoint != -1 ) {
     getPerm(startPoint,g);
     int taken[gridsize];
     int istack[gridsize];
     for( int i = 0; i < gridsize; i++ )
       istack[i] = 0;
     for( int i = 0; i < gridsize; i++ )
       taken[i] = 0;
     for( int i = 0; i < min_depth; i++ )
       taken[g[i]] = 1;
     int depth = min_depth; // this is the index of g we are working on
     int AGrading = AShift;
     for( int i = 0; i < min_depth; i++ )
	 AGrading -= WN[i][g[i]];
     
     while( true ) {
       if ( depth == gridsize ) { // we are at the end of the recursion, use the permutation
	 if (AGrading >= amin && AGrading <= amax) {
	   int MGrading = MaslovGrading(g);
	   generators[AGrading+30][MGrading+30]->push_back(getIndex(g));
	 }
	 depth--;
	 if( depth < min_depth )
	   break;
	 taken[g[depth]] = 0;
	 AGrading += WN[depth][g[depth]];
	 continue;
       }
       
       // If there is no hope of getting a non-negative Alexander Grading from here on, decrease the depth
       int maxAGrading = AGrading;
       for( int i = depth; i < gridsize; i++ ) {
	 maxAGrading -= lowestWN[i];
       }
       if( maxAGrading < amin ) {
	 depth--;
	 if( depth < min_depth )
	   break;
	 taken[g[depth]] = 0;
	 AGrading += WN[depth][g[depth]];
	 continue;
       }
       
       bool callBack = true;
       for( int i = istack[depth]; i < gridsize; i++ ) { // consider using value i for position depth
	 if( !taken[i] ) {
	   g[depth] = i;
	   taken[i] = 1;
	   // push onto the stack
	   callBack = false;
	   istack[depth] = i+1;
	   if( depth < gridsize-1 )
	     istack[depth+1] = 0;
	   AGrading -= WN[depth][g[depth]];
	   depth++;
	   break;
	 }
       }
       if( callBack ) {
	 depth--;
	 if( depth < min_depth )
	   break;
	 taken[g[depth]] = 0;
	 AGrading += WN[depth][g[depth]];
       }
       
     }
     /*
     // Request more to do, send back results, then clear generators
     for( int i = 0; i < 60; i++ )
       for( int j = 0; j < 60; j++ )
	 num_generators[i*60+j] = generators[i][j]->size();
     MPI_Send( num_generators, 60*60, MPI_INT, 0, 0, MPI_COMM_WORLD );
     for( int i = 0; i < 60; i++ )
       for( int j = 0; j < 60; j++ )
	 if( num_generators[i*60+j] ) {
	   MPI_Send( &(generators[i][j]->front()), num_generators[i*60+j], MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD );
	   generators[i][j]->clear();
	 }
     */
     // Just request and receive more
     MPI_Send( num_generators, 0, MPI_INT, 0, 0, MPI_COMM_WORLD );
     MPI_Recv( &startPoint, 1, MPI_LONG_LONG_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
   }

 }
 }
 if( rank == 0 ) {
   agEndTime = read_timer();
   printf("Time to compute all gradings %f\n", agEndTime-agStartTime); 
 }
 if( !justTime ) {
   for( int i = 0; i < 60; i++ )
     for( int j = 0; j < 60; j++ )
       if (generators[i][j]->size()) {
	 printf("Alexander grading %d, Maslov grading %d, rank %d, num generators %lu\n", i-30, j-30, rank, generators[i][j]->size());
	 for( int k = 0; k < generators[i][j]->size(); k++ )
	   if( (*(generators[i][j]))[k] == 265375 or (*(generators[i][j]))[k] == 270415 )
	     printf("(rank %d) %d %d %lld\n", rank, i, j, (*(generators[i][j]))[k]);

       }
 }

 vector<long long> *rows, *cols;
 rows = new vector<long long>;
 cols = new vector<long long>;
 int g[gridsize];
 for( int i = 0; i < 60; i++ )
   for( int j = 0; j < 60; j++ ) {
#define NGEN 5
#define GEN 6
     // Eventually we should do these in parallel, but for now in series
     // in which case these are the correct values
     int firstProcThisMatrix = 0;
     int numProcThisMatrix = n_proc;
     MPI_Comm CommunicatorIJ;
     MPI_Comm_dup(MPI_COMM_WORLD, &CommunicatorIJ);
     //MPI_Comm CommunicatorIJ = MPI_COMM_WORLD;

     int num_generators = generators[i][j]->size();
     delete rows;
     rows = cols;
     cols = new vector<long long>();
     // Gather all and sort.  This could be done more efficiently, but for now just have everyone send all the generators to firstProcThisMatrix.
     if( rank > firstProcThisMatrix && rank < firstProcThisMatrix+numProcThisMatrix ) {
       MPI_Send( &num_generators, 1, MPI_INT, firstProcThisMatrix, NGEN, CommunicatorIJ );
       if( num_generators ) {
	 MPI_Send( &(generators[i][j]->front()), num_generators, MPI_LONG_LONG_INT, firstProcThisMatrix, GEN, CommunicatorIJ );
	 generators[i][j]->clear();
       }
     } else if (rank == firstProcThisMatrix) {
       for( int k = 1; k < numProcThisMatrix; k++ ) {
	 int num_generators;
	 MPI_Status status;
	 //printf("(%d %d) waiting to recieve message %d of %d\n", i, j, k, numProcThisMatrix-1);
	 MPI_Recv( &num_generators, 1, MPI_INT, MPI_ANY_SOURCE, NGEN, CommunicatorIJ, &status );
	 if( num_generators ) {
	   cols->reserve(cols->size()+num_generators);
	   long long buffer[num_generators];
	   MPI_Recv( buffer, num_generators, MPI_LONG_LONG_INT, status.MPI_SOURCE, GEN, CommunicatorIJ, MPI_STATUS_IGNORE );
	   for( int k = 0; k < num_generators; k++ ) {
	     cols->push_back(buffer[k]);	 
	   }
	 }
       }
       if( cols->size() ) {
	 printf("num generators in ranks (%d,%d) is %lu\n", i-30, j-30, cols->size());
       }
       std::sort(cols->begin(), cols->end());
     }

     // Send the sorted list of generators to everyone
     if( rank == firstProcThisMatrix )
       num_generators = cols->size();
     MPI_Bcast( &num_generators, 1, MPI_INT, firstProcThisMatrix, CommunicatorIJ );
     long long *buffer;
     if( rank == firstProcThisMatrix )
       buffer = &(cols->front());
     else
       buffer = (long long*) malloc( num_generators * sizeof(long long) );
     MPI_Bcast( buffer, num_generators, MPI_LONG_LONG_INT, firstProcThisMatrix, CommunicatorIJ );
     if( rank != firstProcThisMatrix )
       cols->assign(buffer, buffer+num_generators);

     // if j is 0, we are done, because there is no calculation with those cols
     if (j == 0)
       continue;

     if( cols->size() == 0 || rows->size() == 0 ) {
       imageDimensions[i][j] = 0;
       kernelDimensions[i][j] = cols->size();
       continue;
     }

     // For Debugging only
     //if( i == 30 && j == 26 )
     //  continue;

     printf("(%d) beginning calculation %d,%d\n", rank, i-30, j-30);
     // Decide who has ownership of what subset of the generators (without communication).
     int colsPerProcMin = num_generators/numProcThisMatrix;
     int extraCols = num_generators - numProcThisMatrix * colsPerProcMin;
     int firstCol, numCols;
     if( rank < firstProcThisMatrix + extraCols ) {
       firstCol = (rank-firstProcThisMatrix)*(colsPerProcMin+1);
       numCols = colsPerProcMin+1;
     } else {
       firstCol = extraCols + (rank-firstProcThisMatrix)*colsPerProcMin;
       numCols = colsPerProcMin;
     }
     printf("(%d) ownership determined, my block is %lu %d starting at %d\n", rank, rows->size(), numCols, firstCol);
     //printf("Here %d %d %d\n", i, j, rank);
     // Each proc calculates the part of the matrix that it owns, storing both rows and columns as in the serial code.  It holds entire columns, but only certain rows
     vector<Generator> GraphIn( rows->size() ); // Will hold boundary data.
     vector<Generator> GraphOut( numCols ); // Will hold boundary data.
     //printf("Filling %d %d\n", I, J);
     for(int index=0; index < numCols; index++) {
       getPerm((*cols)[index+firstCol],g);
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
	     int indexgij = Find(*rows,Indexgij);
	     if(indexgij==-1) {printf("Error with Alexander grading: %lld->%lld\n", (*cols)[index], Indexgij); return 0; }
#ifdef FIELD_Z3
	     printf("Not implemented\n");
	     exit(-1);
	     //GraphOut[index].ones.push_back( indexgij );
	     //GraphIn[indexgij].ones.push_back( index );
	     // This is, of course, wrong
#else
	     GraphOut[index].ones.push_back( indexgij );
	     GraphIn[indexgij].ones.push_back( index );
#endif
	   }
	 }
       }
     }
     if( printMatrices ) {
       char name[15];
       sprintf(name, "mat%d,%dp%d.dat", i, j, rank);
       printMatrix( name, GraphOut );
     }
     printf("(%d) Matrices filled %d %d\n", rank, i-30, j-30);

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

     int turn = firstProcThisMatrix; // the id of the matrix whose turn it is
     int globalTurn = firstProcThisMatrix; // This just keeps increasing
     std::deque<int*> colsToReduceBy;
     std::deque<MPI_Request*> outRequests; // keep these so we can make sure we have sent a given column when we dequeue it befor we delete it
     bool waitingForSize = false;
     bool waitingForCol = false;
     int waitingSource = -1; // the processor from which we are currently waiting
     bool colOfMineOnStack = false; // we don't want to ever have multiple columns of our own on the stack (it means the one we are about to add hasn't been reduced sufficiently yet)
     int nextProc = rank + 1;
     if( nextProc == numProcThisMatrix+firstProcThisMatrix )
       nextProc = firstProcThisMatrix;
     int prevProc = rank - 1;
     if( prevProc < firstProcThisMatrix )
       prevProc = firstProcThisMatrix + numProcThisMatrix - 1;
     MPI_Request *inRequest;
     MPI_Request *finishedInRequest;
     MPI_Request *finishedOutRequest;
     int finishedOutSignal[2];
     int finishedInSignal[2];
     finishedInRequest = new MPI_Request();
     finishedOutRequest = NULL;
     MPI_Irecv( finishedInSignal, 2, MPI_INT, prevProc, TAG_FINISHED, CommunicatorIJ, finishedInRequest );
     int inSize;
     int *inCol=0;
     //bool starting = true;
     list<int>::iterator indexInCol;
     list<int>::iterator indexInColEndValue = indexInCol;
     int *currentCol;
     MPI_Request *currentColOutRequest;
     int c = 0;
     int numProcReportingFinished = 0;
     int procsReportingFinished[numProcThisMatrix]; // this stores the globalTurn at which each proc has reported it is finished, or -1 if it has not reported finished
     for( int p = 0; p < numProcThisMatrix; p++ )
       procsReportingFinished[p] = -1;
     bool selfFinished = false;
     int zero = 0;

     //while( numProcReportingFinished < numProcThisMatrix ) { // main loop for the reduction.  The key principle is that we only do one thing per loop iteration
     while ( true ) {
       bool allDone = true;
       for( int p = 0; p < numProcThisMatrix; p++ )
	 if( procsReportingFinished[p] == -1 || procsReportingFinished[p] > globalTurn ) {
	   allDone = false;
	   break;
	 }
       if( allDone )
	 break;
       // make sure turn is someone who is actually still working
       while( procsReportingFinished[turn-firstProcThisMatrix] != -1 && procsReportingFinished[turn-firstProcThisMatrix] <= globalTurn ) {
	 increment(turn, firstProcThisMatrix, numProcThisMatrix);
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
	   if( numProcReportingFinished < numProcThisMatrix-1 ||
	       numProcReportingFinished == numProcThisMatrix-1 && selfFinished ) {
	     finishedInRequest = new MPI_Request();
	     //printf("(%d) prepared to receive another done signal\n", rank);
	     MPI_Irecv( finishedInSignal, 2, MPI_INT, prevProc, TAG_FINISHED, CommunicatorIJ, finishedInRequest );
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
	   procsReportingFinished[finishedInSignal[0]-firstProcThisMatrix] = finishedInSignal[1];
	   numProcReportingFinished++;
	   if( finishedInSignal[0] != nextProc ) {
	     finishedOutRequest = new MPI_Request();
	     //printf("(%d) sending on finished signal from %d\n", rank, finishedInSignal[0]);
	     MPI_Isend( finishedInSignal, 2, MPI_INT, nextProc, TAG_FINISHED, CommunicatorIJ, finishedOutRequest );
	   } else {
	   // put out a new request if there are any procs left not reporting done
	     if( numProcReportingFinished < numProcThisMatrix-1 ||
		 numProcReportingFinished == numProcThisMatrix-1 && selfFinished ) {
	       finishedInRequest = new MPI_Request();
	       //printf("(%d) prepared to receive another done signal\n", rank);
	       MPI_Irecv( &finishedInSignal, 2, MPI_INT, prevProc, TAG_FINISHED, CommunicatorIJ, finishedInRequest );
	     }	   
	   }
	   continue;
	 }
       }
       // check if we ourselves are done
       if( !selfFinished && c >= numCols ) {
	 selfFinished = true;
	 procsReportingFinished[rank-firstProcThisMatrix] = globalTurn;
	 numProcReportingFinished++;
	 // send a message along the chain
	 MPI_Request req;
	 //printf("(%d) sending out my done signal\n", rank);
	 finishedOutSignal[0] = rank;
	 finishedOutSignal[1] = globalTurn;
	 MPI_Isend( finishedOutSignal, 2, MPI_INT, nextProc, TAG_FINISHED, CommunicatorIJ, &req);
	 continue;
       }

       // check if we should prepare to receive
       if( turn != rank && !waitingForSize && !waitingForCol ) {
	 inRequest = new MPI_Request();
	 MPI_Irecv( &inSize, 1, MPI_INT, prevProc, TAG_COLUMN_SIZE, CommunicatorIJ, inRequest );
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
	     MPI_Irecv( inCol, inSize+2, MPI_INT, prevProc, TAG_COLUMN, CommunicatorIJ, inRequest );
	     waitingForCol = true;
	   } else {
	     waitingSource = -1;
	     globalTurn += rank-turn;
	     if( rank < turn )
	       globalTurn += numProcThisMatrix-firstProcThisMatrix;
	     turn = rank;
	     // a special case: if this process is finished, it needs to pass on a zero size so that the calculation can continue
	     if( selfFinished ) {
	       MPI_Request req;
	       MPI_Isend( &zero, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, CommunicatorIJ, &req);
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
	     if( inCol[e] >= rows->size() ) {
	       printf("(%d) FATAL receive error e=%d inCol[e]=%d numRows=%lu (%d %d, %d)\n", rank, e, inCol[e], rows->size(), i, j, c ); 
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
	     globalTurn += numProcThisMatrix-firstProcThisMatrix;
	   //printf("(%d) received column %d %d %d\n", rank, inCol[0], inCol[1], inCol[2]);
	   turn = inCol[0]; // this should take care of dropped blank columns
	   increment(turn, firstProcThisMatrix, numProcThisMatrix);
	   globalTurn++;
	   if( inCol[1] == 0 ) {
	     printf("(%d) Warning: just received a column of zero size, size should have been %d\n", rank, inSize);
	   }
	   if( inCol[0] != nextProc ) { // if it is, this column has made it all the way around
	     MPI_Request request1, *request2;
	     request2 = new MPI_Request();
	     //printf("(%d) passing along column %d %d %d\n", rank, inCol[0], inCol[1], inCol[2]);
	     MPI_Isend( inCol+1, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, CommunicatorIJ, &request1 ); 
	     MPI_Isend( colsToReduceBy.back(), inCol[1]+2, MPI_INT, nextProc, TAG_COLUMN, CommunicatorIJ, request2 );
	     outRequests.push_back(request2);
	   } else {
	     // need to keep outRequests in sync with colsToReduceBy
	     MPI_Request rnull = MPI_REQUEST_NULL;
	     outRequests.push_back(&rnull);
	   }
	   if( c < numCols ) {
	     resolveCols(GraphOut, GraphIn, inCol+2, inCol[1], c);
	   }
	   continue;
	 }
       }

       if( turn == rank && selfFinished ) {
	 assert(false);
	 /*
	 MPI_Request request1;
	 MPI_Isend( &zero, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, CommunicatorIJ, &request1 );
	 increment(turn, firstProcThisMatrix, numProcThisMatrix);
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
	   if( column[e] >= rows->size() ) {
	     printf("(%d) FATAL error %d %d %lu\n", rank, e, column[e], rows->size() ); 
	   }
	 }
	 //printf("(%d) sending column %d %d %d\n", rank, column[0], column[1], column[2]);
	 increment(turn, firstProcThisMatrix, numProcThisMatrix);
	 globalTurn++;
	 // send this column to the next processor 
	 // first send the size
	 MPI_Request request1, *request2;
	 request2 = new MPI_Request();
	 MPI_Isend( column+1, 1, MPI_INT, nextProc, TAG_COLUMN_SIZE, CommunicatorIJ, &request1 );
	 if( numEntries ) {
	   MPI_Isend( column, numEntries+2, MPI_INT, nextProc, TAG_COLUMN, CommunicatorIJ, request2 );
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
	       resolveCols(GraphOut, GraphIn, columnFromQueue+2, columnFromQueue[1], c);
	     }
	     // also the column we are currently working on
	     if( indexInCol != indexInColEndValue ) {
	       // this should just deal with it, whether or not there is anything to resolve
	       resolveCols(GraphOut, GraphIn, currentCol+2, currentCol[1], c);
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
	   if( *indexInCol <= c ) {
	     continue; // these columns have already been dealt with
	   }
	   for( int k = 2; k < currentCol[1]+2; k++ ) {
	     list<int>::iterator search = find( GraphOut[*indexInCol].ones.begin(), GraphOut[*indexInCol].ones.end(), currentCol[k] );
	     if( search != GraphOut[*indexInCol].ones.end() ) {
	       GraphOut[*indexInCol].ones.erase(search);
	       if( k != 2 )
		 GraphIn[currentCol[k]].ones.remove(*indexInCol);
	     } else {
	       GraphOut[*indexInCol].ones.push_back(currentCol[k]);
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
     printf("(%d) Outside main reduction loop %d %d\n", rank, i-30, j-30);


     printf("(%d) imposing barrier\n", rank);
     // Probably don't need this
     MPI_Barrier(CommunicatorIJ);


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

     // this would probably be more efficient done at the end with all at once
     MPI_Reduce( &localKD, &kernelDimensions[i][j], 1, MPI_UNSIGNED, MPI_SUM, firstProcThisMatrix, CommunicatorIJ );
     MPI_Reduce( &localID, &imageDimensions[i][j], 1, MPI_UNSIGNED, MPI_SUM, firstProcThisMatrix, CommunicatorIJ );
     if( rank == firstProcThisMatrix )
       printf("(%d) TAG gradings %d %d, %d %d\n", rank, i-30, j-30, kernelDimensions[i][j], imageDimensions[i][j] );

     //if( i == 30 && j == 27 ) {
     //  MPI_Finalize();
     //  exit(0);
     //}
     MPI_Comm_free(&CommunicatorIJ);
   }
     
     
 if( rank == 0 ) {
   printf("Time to collect and sort all generators, and compute the matrices %f\n", read_timer()-agEndTime); 
 }

 // Finish off on a single processor

 if( rank == 0 ) {
   for( int i = 0; i < 60; i++ )
     for( int j = 1; j < 60; j++ )
       homologyDimensions[i][j-1] = kernelDimensions[i][j-1] - imageDimensions[i][j];
   printf("Ranks of unshifted homology groups in Alexander grading [%d,%d]:\n", amin, amax);
   for(int a=amax+30; a>=amin+30; a--) {
     for(int m=20; m<40; m++) {
       printf("%4d", homologyDimensions[a][m]);
     }
     printf("\n");
   }

   int HFKRanks [60][60]; // HFKRanks[i][j] will hold rank of HFK^ in Maslov grading=i-30 and Alexander grading=j-30
   for(int a=0; a<60; a++) { for(int m=0; m<60; m++) HFKRanks[m][a]=0; }
 
   // Reproduce HFK^ from HFK^ \otimes K^{gridsize-1} in non-negative Alexander grading
   for(int a=59; a>=0; a--) {
     for(int m=59; m>=0; m--) {
       if( homologyDimensions[a][m] > 0) {
	 HFKRanks[m][a] = homologyDimensions[a][m];
	 for(int i=0; i<=min(gridsize-numcomp,min(a,m)); i++) homologyDimensions[a-i][m-i] -= (HFKRanks[m][a] * Factorial[gridsize-numcomp]) / (Factorial[i] * Factorial[gridsize-numcomp-i]);
       }
     }
   }
   
   // Use symmetry to fill up HFKRanks in negative Alexander gradings
   for(int alex=-1; alex>=-9; alex--){ 
     for(int mas=-20; mas < 12; mas++) {
       HFKRanks[mas+30][alex+30] = HFKRanks[mas-2*alex+30 ][-alex+30];
     }
   }
   if(amin > 0) printf("This Poincare polynomial is only valid in Alexander grading >= %d:\n", amin);
   // Print Results
   bool first=1;
   for(int a=-20; a<19; a++) {
     for(int m=-20; m<19; m++) {
       int rankam = HFKRanks[m+30][a+30];
       if(rankam > 0) {
	 if(!first) printf("+");
	 else first=0;
	 if(rankam > 1 || (rankam==1 && a==0 && m==0) ) printf("%d", rankam);
	 if(m==1) printf("q");
	 if(m != 0 && m != 1) printf("q^{%d}",m);
	 if(a==1) printf("t");
	 if(a != 0 && a != 1) printf("t^{%d}",a);
       }
     }
   }
   printf("\n");
   //time_t endtime = time(NULL);
   //printf("Total time elapsed: %ld seconds.\n", endtime - starttime);
   
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
