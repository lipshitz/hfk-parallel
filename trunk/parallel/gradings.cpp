#include <sys/time.h>
#include <iostream>
#include <stdlib.h>
#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <algorithm>


//#define FIELD_Z3
//#include "matrix-z2z3.h"
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

 
 vector<long long> label; // label[i] will hold the number of perms lexicographically before the i^th generator
 // This will hold the generators, sorted by grading, the first index is Agrading+30, the second is Maslov-grading+30
 std::vector<long long> *generators[60][60];
 unsigned long long imageDimensions[60][60];
 unsigned long long kernelDimensions[60][60];
 int homologyDimensions[60][60];
 for( int i = 0; i < 60; i++ )
   for( int j = 0; j < 60; j++ )
     generators[i][j] = new std::vector<long long>();
 
 if( rank == 0 && !justTime )
   printf("Searching through %lld generators to compute Alexander gradings...\n", Factorial[gridsize]);
 double agStartTime = read_timer();

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

 if( rank == 0 )
   printf("Time to compute all gradings %f\n", read_timer()-agStartTime); 
 if( !justTime ) {
   for( int i = 0; i < 60; i++ )
     for( int j = 0; j < 60; j++ )
       if (generators[i][j]->size()) {
	 printf("Alexander grading %d, Maslov grading %d, rank %d, num generators %lu\n", i-30, j-30, rank, generators[i][j]->size());
       }
 }

 // Just stop here for the moment.

 /*
 // Calculate the homology groups
 for( int I = 0; I < 60; I++ )
   for( int J = 0; J < 59; J++ ) {
     std::vector<long long> &cols = *(generators[I][J+1]);
     std::vector<long long> &rows = *(generators[I][J]);
     if( cols.size() == 0 || rows.size() == 0 ) {
       imageDimensions[I][J] = 0;
       kernelDimensions[I][J] = cols.size();
       continue;
     }
     vector<Generator> GraphIn( rows.size() ); // Will hold boundary data.
     vector<Generator> GraphOut( cols.size() ); // Will hold boundary data.
     long long edges=0;
     printf("Filling %d %d\n", I, J);
     for(int index=0; index < cols.size(); index++) {
       getPerm(cols[index],g);
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
	     if(indexgij==-1) {printf("Error with Alexander grading: %lld\n", Indexgij); return 0; }
#ifdef FIELD_Z3
	     //printf("Not implemented\n");
	     //exit(-1);
	     GraphOut[index].ones.push_back( indexgij );
	     GraphIn[indexgij].ones.push_back( index );
	     // This is, of course, wrong
#else
	     GraphOut[index].ones.push_back( indexgij );
	     GraphIn[indexgij].ones.push_back( index );
#endif
	     edges++;
	   }
	 }
       }
     }
     printf("Size is %lu by %lu\n", rows.size(), cols.size());
     if( printMatrices ) {
       char name[10];
       sprintf(name, "mat%d,%d.dat", I, J);
       printf("Writing to mat%d,%d\n", I, J);
       printMatrix( name, GraphOut );
     }

     printf("Reducing %d %d\n", I, J);
     for( int i = 0; i < GraphOut.size(); i++ ) {
#ifdef FIELD_Z3
       if( (!GraphOut[i].alive) ) continue;
       int target;
       int targetVal; // We plan to use target as a pivot, the value is targetVal
       if( GraphOut[i].ones.size()==0 ) {
	 if( GraphOut[i].twos.size() == 0 )
	   continue;
	 targetVal = 2;
	 target = GraphOut[i].twos.front();
       } else {
	 targetVal = 1;
	 target = GraphOut[i].ones.front();
       }
       GraphOut[i].alive = 0;
       // for all rows with a 1 in column target
       bool add = (targetVal == 2); // whether we add or subtract row i from row j
       for( list<int>::iterator j = GraphIn[target].ones.begin(); j != GraphIn[target].ones.end(); j++ ) {
	 if( !GraphOut[*j].alive ) continue;

	 // for all 1s in row i
	 for( list<int>::iterator k = GraphOut[i].ones.begin(); k != GraphOut[i].ones.end(); k++ ) {
	   // check if that entry in row j is a 1
	   list<int>::iterator search = find( GraphOut[*j].ones.begin(), GraphOut[*j].ones.end(), *k);
	   if( search != GraphOut[*j].ones.end() ) {
	     if( add ) {
	       // Change the 1 to a 2
	       GraphOut[*j].ones.erase(search);
	       GraphOut[*j].twos.push_back(*search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
		 GraphIn[*k].twos.push_back(*j);
	       }
	     } else {
	       // Remove the 1
	       GraphOut[*j].ones.erase(search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
	       }
	     }
	   } else {
	     // check if that entry in row j is a 2
	     search = find( GraphOut[*j].twos.begin(), GraphOut[*j].twos.end(), *k);
	     if( search != GraphOut[*j].twos.end() ) {
	       if( add ) {
		 // Remove the 2
		 GraphOut[*j].twos.erase(search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		 }
	       } else {
		 // Change the 2 to a 1
		 GraphOut[*j].twos.erase(search);
		 GraphOut[*j].ones.push_back(*search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		   GraphIn[*k].ones.push_back(*j);
		 }
	       }
	     }
	     else {
	       // the entry in row j is a 0
	       if( add ) {
		 // Set it to one
		 GraphOut[*j].ones.push_back(*k);
		 if( *k != target )  {
		   GraphIn[*k].ones.push_back(*j);
		 }
	       } else {
		 // Set it to two
		 GraphOut[*j].twos.push_back(*k);
		 if( *k != target )  {
		   GraphIn[*k].twos.push_back(*j);
		 }
	       }
	     }
	   }
	 }
		 
	 // for all 2s in row i
	 for( list<int>::iterator k = GraphOut[i].twos.begin(); k != GraphOut[i].twos.end(); k++ ) {
	   // check if that entry in row j is a 1
	   list<int>::iterator search = find( GraphOut[*j].ones.begin(), GraphOut[*j].ones.end(), *k);
	   if( search != GraphOut[*j].ones.end() ) {
	     if( !add ) {
	       // Change the 1 to a 2
	       GraphOut[*j].ones.erase(search);
	       GraphOut[*j].twos.push_back(*search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
		 GraphIn[*k].twos.push_back(*j);
	       }
	     } else {
	       // Remove the 1
	       GraphOut[*j].ones.erase(search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
	       }
	     }
	   } else {
	     // check if that entry in row j is a 2
	     search = find( GraphOut[*j].twos.begin(), GraphOut[*j].twos.end(), *k);
	     if( search != GraphOut[*j].twos.end() ) {
	       if( !add ) {
		 // Remove the 2
		 GraphOut[*j].twos.erase(search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		 }
	       } else {
		 // Change the 2 to a 1
		 GraphOut[*j].twos.erase(search);
		 GraphOut[*j].ones.push_back(*search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		   GraphIn[*k].ones.push_back(*j);
		 }
	       }
	     }
	     else {
	       // the entry in row j is a 0
	       if( !add ) {
		 // Set it to one
		 GraphOut[*j].ones.push_back(*k);		 
		 if( *k != target )  {
		   GraphIn[*k].ones.push_back(*j);
		 }
	       } else {
		 // Set it to two
		 GraphOut[*j].twos.push_back(*k);
		 if( *k != target )  {
		 GraphIn[*k].twos.push_back(*j);
		 }
	       }
	     }
	   }

	 }
	 
       }
       GraphIn[target].ones.clear();


       // for all rows with a 2 in column target
       add = (targetVal == 1); // whether we add or subtract row i from row j
       for( list<int>::iterator j = GraphIn[target].twos.begin(); j != GraphIn[target].twos.end(); j++ ) {

	 if( !GraphOut[*j].alive ) continue;

	 // for all 1s in row i
	 for( list<int>::iterator k = GraphOut[i].ones.begin(); k != GraphOut[i].ones.end(); k++ ) {
	   // check if that entry in row j is a 1
	   list<int>::iterator search = find( GraphOut[*j].ones.begin(), GraphOut[*j].ones.end(), *k);
	   if( search != GraphOut[*j].ones.end() ) {
	     if( add ) {
	       // Change the 1 to a 2
	       GraphOut[*j].ones.erase(search);
	       GraphOut[*j].twos.push_back(*search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
		 GraphIn[*k].twos.push_back(*j);
	       }
	     } else {
	       // Remove the 1
	       GraphOut[*j].ones.erase(search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
	       }
	     }
	   } else {
	     // check if that entry in row j is a 2
	     search = find( GraphOut[*j].twos.begin(), GraphOut[*j].twos.end(), *k);
	     if( search != GraphOut[*j].twos.end() ) {
	       if( add ) {
		 // Remove the 2
		 GraphOut[*j].twos.erase(search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		 }
	       } else {
		 // Change the 2 to a 1
		 GraphOut[*j].twos.erase(search);
		 GraphOut[*j].ones.push_back(*search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		   GraphIn[*k].ones.push_back(*j);
		 }
	       }
	     }
	     else {
	       // the entry in row j is a 0
	       if( add ) {
		 // Set it to one
		 GraphOut[*j].ones.push_back(*k);
		 if( *k != target )  {
		   GraphIn[*k].ones.push_back(*j);
		 }
	       } else {
		 // Set it to two
		 GraphOut[*j].twos.push_back(*k);
		 if( *k != target )  {
		   GraphIn[*k].twos.push_back(*j);
		 }
	       }
	     }
	   }
	 }
	 
	 // for all 2s in row i
	 for( list<int>::iterator k = GraphOut[i].twos.begin(); k != GraphOut[i].twos.end(); k++ ) {
	   // check if that entry in row j is a 1
	   list<int>::iterator search = find( GraphOut[*j].ones.begin(), GraphOut[*j].ones.end(), *k);
	   if( search != GraphOut[*j].ones.end() ) {
	     if( !add ) {
	       // Change the 1 to a 2
	       GraphOut[*j].ones.erase(search);
	       GraphOut[*j].twos.push_back(*search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
		 GraphIn[*k].twos.push_back(*j);
	       }
	     } else {
	       // Remove the 1
	       GraphOut[*j].ones.erase(search);
	       if( *k != target )  {
		 GraphIn[*k].ones.remove(*j);
	       }
	     }
	   } else {
	     // check if that entry in row j is a 2
	     search = find( GraphOut[*j].twos.begin(), GraphOut[*j].twos.end(), *k);
	     if( search != GraphOut[*j].twos.end() ) {
	       if( !add ) {
		 // Remove the 2
		 GraphOut[*j].twos.erase(search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		 }
	       } else {
		 // Change the 2 to a 1
		 GraphOut[*j].twos.erase(search);
		 GraphOut[*j].ones.push_back(*search);
		 if( *k != target )  {
		   GraphIn[*k].twos.remove(*j);
		   GraphIn[*k].ones.push_back(*j);
		 }
	       }
	     }
	     else {
	       // the entry in row j is a 0
	       if( !add ) {
		 // Set it to one
		 GraphOut[*j].ones.push_back(*k);
		 if( *k != target )  {		 
		   GraphIn[*k].ones.push_back(*j);
		 }
	       } else {
		 // Set it to two
		 GraphOut[*j].twos.push_back(*k);
		 if( *k != target )  {
		   GraphIn[*k].twos.push_back(*j);
		 }
	       }
	     }
	   }

	 }
	 


       }



       GraphIn[target].twos.clear();
#else
       if( (!GraphOut[i].alive) || GraphOut[i].ones.size()==0 ) continue;
       int target = GraphOut[i].ones.front(); // We plan to delete the edge from i to target
       GraphOut[i].alive = 0;
       for( list<int>::iterator j = GraphIn[target].ones.begin(); j != GraphIn[target].ones.end(); j++ ) {
	 if( !GraphOut[*j].alive ) continue;
	 for( list<int>::iterator k = GraphOut[i].ones.begin(); k != GraphOut[i].ones.end(); k++ ) {
	   list<int>::iterator search = find( GraphOut[*j].ones.begin(), GraphOut[*j].ones.end(), *k);
	   if( search != GraphOut[*j].ones.end() ) {
	     GraphOut[*j].ones.erase(search);
	     if( *k != target )  {
	       GraphIn[*k].ones.remove(*j);
	     }
	   } else {
	     GraphOut[*j].ones.push_back(*k);
	     if( *k != target )  {
	       GraphIn[*k].ones.push_back(*j);
	     }
	   }
	 }
       }
       GraphIn[target].ones.clear();
#endif
     }
     kernelDimensions[I][J] = 0;
     for( int i = 0; i < GraphOut.size(); i++ ) {
       if( GraphOut[i].alive )
	 kernelDimensions[I][J]++;
     }
     imageDimensions[I][J] = cols.size() - kernelDimensions[I][J];
     printf("Rank is %llu\n", imageDimensions[I][J]);
     printf("------------\n");
   }
 
 for( int i = 0; i < 60; i++ )
   for( int j = 1; j < 60; j++ ) {
     homologyDimensions[i][j] = kernelDimensions[i][j-1] - imageDimensions[i][j];
     if( homologyDimensions[i][j] )
       printf("Alexander grading %d Maslov grading %d Homology dimension %d\n", i-30, j-30, homologyDimensions[i][j]);
   }
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
 time_t endtime = time(NULL);
 printf("Total time elapsed: %ld seconds.\n", endtime - starttime);

 */
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
  int dwhite[gridsize];
  for (int i =0; i<gridsize; i++){
    dblack[black[i]]=i;
    dwhite[white[i]]=i;
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
