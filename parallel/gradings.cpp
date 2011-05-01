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

int gridsize = 12; // arc-index
int *white;
int *black;
int default_white[12] = {9,5,11,7,8,1,10,4,0,3,2,6};
int default_black[12] = {1,0,4,3,2,6,5,9,8,11,7,10};

// Fill in the big ones later since g++ doesn't seem to like big constants
long long Factorial[16] = {
  1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,
  0,0,0};

void getPerm( long long n, int *P);
void NextPerm(short counter[], int h[]);
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
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

 if( n_proc <= 1 ) {
   printf("Run this with at least 2 threads\n");
   MPI_Finalize();
   exit(0);
 }
 
 char *knotFile = read_string( argc, argv, "-k", NULL );
 bool printMatrices = read_int( argc, argv, "-p", 0 );
 bool justTime = read_int( argc, argv, "-t", 0 );
 char *saveDir = read_string( argc, argv, "-o", NULL );
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

#define NGEN 5
#define GEN 6
 // Gather all the results to the root node and save
 if( saveDir ) {
   for( int i = 0; i < 60; i++ )
     for( int j = 0; j < 60; j++ ) {
       if( rank != 0 ) {
	 int num_generators = generators[i][j]->size();
	 MPI_Send( &num_generators, 1, MPI_INT, 0, NGEN, MPI_COMM_WORLD );
	 if( num_generators ) {
	   MPI_Send( &(generators[i][j]->front()), num_generators, MPI_LONG_LONG_INT, 0, GEN, MPI_COMM_WORLD );
	   generators[i][j]->clear();
	 }
       } else {
	 vector<long long> gens;
	 for( int k = 1; k < n_proc; k++ ) {
	   int num_generators;
	   MPI_Status status;
	   MPI_Recv( &num_generators, 1, MPI_INT, MPI_ANY_SOURCE, NGEN, MPI_COMM_WORLD, &status );
	   if( num_generators ) {
	     gens.reserve(gens.size()+num_generators);
	     //long long buffer[num_generators];
	     long long *buffer = (long long*) malloc( num_generators*sizeof(long long) );
	     MPI_Recv( buffer, num_generators, MPI_LONG_LONG_INT, status.MPI_SOURCE, GEN, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	     for( int k = 0; k < num_generators; k++ ) {
	       gens.push_back(buffer[k]);	 
	     }
	     free(buffer);
	   }
	 }
	 std::sort(gens.begin(), gens.end());
	 if( gens.size() ) {
	   char outFile[30];
	   sprintf(outFile, "%s/gen%d,%d.dat", saveDir, i, j);
	   FILE *f = fopen( outFile, "w" );
	   fprintf(f, "%lu\n", gens.size());
	   for( int k = 0; k < gens.size(); k++ )
	     fprintf(f, "%lld\n", gens[k]);
	   fclose(f);
	 }
       }
       MPI_Barrier(MPI_COMM_WORLD);
     }
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
