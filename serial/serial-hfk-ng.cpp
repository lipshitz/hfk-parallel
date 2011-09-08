#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <algorithm>
#include "command-line-parser.h"
#include "matrix.h"
#include "gradings.h"
using std::list;
using std::vector;

const int NUM_GRADINGS=2;

// Function Prototypes

void NextPerm(short counter[], int h[]);
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
int WindingNumber(int x, int y, int*, int*); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y [], int*, int);
int NumComp(int*, int*, int); //Returns the number of components of the link.
bool ValidGrid(int*, int*, int);

// Main

int main(int argc, char *argv[]){

 int gridsize; //arc index
 int *white;
 int *black;
  
 char *knotFile = read_string( argc, argv, "-k", NULL );
 bool printMatrices = read_int( argc, argv, "-p", 0 );
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
   printf("Please specify a knot file with -k\n");
   exit(-1);
 }

 int  numcomp = NumComp(white, black, gridsize);
 printf("Number of components: %d\n", numcomp);

 if(!ValidGrid(white, black, gridsize)) {printf("Invalid grid!!\n"); return 0;} // Check that the grid is valid
 time_t starttime = time(NULL); // Used to record how long this takes

 // Record winding numbers around grid points for Alexander grading computations
 // Also record the Alexander grading shift
 // Want this to be non-positive; if it isn't, exchange the white and black dots,
 // multiply all winding numbers by -1 and try again.
 int WN[gridsize][gridsize];
 for(int x=0; x<gridsize; x++) {
  for(int y=0; y<gridsize; y++) {
    WN[x][y] = WindingNumber(x,y,white,black);
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
 
 HalfInteger amin = 0;
 double AShiftD = (temp - 4 * gridsize + 4)/8.;
 HalfInteger AShift(AShiftD);
 if( AShiftD != (int)AShiftD ) {
   // Either of the next two lines works; which is more correct?
   amin = -0.5;
   //AShift += 0.5;
   printf("Alexander Grading shift is not integer: %f\n", AShiftD);
 }
 
 //printf("Alexander Grading Shift: %d\n", AShift);
 printf("Alexander Grading Shift: %f\n", AShift.print());
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
 

 // Record for later use whether every possible rectangle has a black or white dot in it
 // This will speed boundary computations.
 printf("Computing which rectangles on the torus have no black or white dots inside.\n");
 initRectangles(gridsize, white, black);

 // Iterate through the generators in lexicographic
 // order and calculate their boundaries
 // Identify each permutation with the integer given by the number of permutations preceding
 // it in lexicographic order.
 

 // The effect of the boundary operator on each grading
 const HalfInteger boundary[2] = {0,-1};
 Gradings generators(boundary);

 printf("Searching through %lld generators to compute Alexander gradings...\n", factorial(gridsize));
 time_t agStartTime = time(NULL);

 int g[gridsize];
 int taken[gridsize];
 int istack[gridsize];
 HalfInteger grad[NUM_GRADINGS];
 for( int i = 0; i < gridsize; i++ )
   istack[i] = 0;
 for( int i = 0; i < gridsize; i++ )
   taken[i] = 0;
 int depth = 0; // this is the index of g we are working on
 HalfInteger AGrading = AShift;
 while( true ) {
   if ( depth == gridsize ) { // we are at the end of the recursion, use the permutation
     if (AGrading >= amin) {
       int MGrading = MaslovGrading(g, white, gridsize);
       grad[0] = AGrading;
       grad[1] = MGrading;
       generators.addGeneratorToGrading( getIndex(g,gridsize), grad );
     }
     depth--;
     if( depth < 0 )
       break;
     taken[g[depth]] = 0;
     AGrading += WN[depth][g[depth]];
     continue;
   }

   // If there is no hope of getting a non-negative Alexander Grading from here on, decrease the depth
   HalfInteger maxAGrading = AGrading;
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
 printf("Time to compute all Alexander gradings %ld\n", time(NULL)-agStartTime); 

 for( gradingIterator it = generators.getGradingIterator(); it != generators.getGradingIteratorEnd(); it++ ) {
   printf("Grading (%2.1f %2.1f): num generators %lu\n", (*it)[0].print(), (*it)[1].print(), generators.getGenerators(*it)->size() );
 }

 for( gradingIterator it = generators.getGradingIterator(); it != generators.getGradingIteratorEnd(); it++ ) {
   vector<generator> &cols = *(generators.getGenerators(*it));
   vector<generator> &rows = *(generators.getBoundaryGenerators(*it));

   int kernelDimension = fillReduceKernel(cols, rows, gridsize);

   int imageDimension = cols.size() - kernelDimension;
   if( kernelDimension != 0 )
     generators.setKernelDimension(*it, kernelDimension);
   if( imageDimension != 0 )
     generators.setImageDimension(*it, imageDimension);
   //printf("Rank is %d\n", imageDimension);
 }
 
 generators.calculateHomology();

 std::map<const HalfInteger*, int, ltgrad> HFKRanks;
 for( gradingIterator it = generators.getGradingIterator(); it != generators.getGradingIteratorEnd(); it++ ) {
   int hd = generators.getHomologyDimension(*it);
   if( hd ) {
     HFKRanks[*it] = hd;
     printf("Alexander grading %2.1f Maslov grading %2.1f homology dimension %d\n", (*it)[0].print(), (*it)[1].print(), hd);
   }
 }
 for( std::map<const HalfInteger*, int, ltgrad>::reverse_iterator rit = HFKRanks.rbegin(); rit != HFKRanks.rend(); rit++ ) {
   for( int i = 1; i <= gridsize-numcomp; i++ ) {
     HalfInteger temp[2];
     temp[0] = (rit->first)[0]-i;
     temp[1] = (rit->first)[1]-i;
     // only try to edit entries that exist
     if( HFKRanks.find((const HalfInteger*)temp) == HFKRanks.end() )
       break;
     HFKRanks[(const HalfInteger*)temp] -= HFKRanks[rit->first]*factorial(gridsize-numcomp) / (factorial(i) * factorial(gridsize-numcomp-i));
   }
 }
 for( std::map<const HalfInteger*, int, ltgrad>::reverse_iterator rit = HFKRanks.rbegin(); rit != HFKRanks.rend(); rit++ ) {
   HalfInteger * temp = new HalfInteger[2];
   if( rit->first[0] < 0 )
     continue;
   temp[0] = amin+amin- (rit->first)[0];
   temp[1] = (rit->first)[1]-rit->first[0]-rit->first[0]+amin+amin;
   HFKRanks[(const HalfInteger*)temp] = HFKRanks[rit->first];
 }
 for( std::map<const HalfInteger*, int, ltgrad>::reverse_iterator rit = HFKRanks.rbegin(); rit != HFKRanks.rend(); rit++ )

 if(amin > 0) printf("This Poincare polynomial is only valid in Alexander grading >= %f:\n", amin.print());
 bool first=1;
 for( std::map<const HalfInteger*, int, ltgrad>::iterator it = HFKRanks.begin(); it != HFKRanks.end(); it++ ) {
   int rankam = it->second;
   HalfInteger m = it->first[1], a = it->first[0];
   if( rankam > 0 ) {
     if(!first) printf("+");
     else first = 0;
     if( rankam > 1 || (rankam==1 && a==0 && m==0) ) printf("%d", rankam);
       if(m==1) printf("q");
       if(m != 0 && m != 1) {
	 double mp = m.print();
	 if (((int)mp)*1. == mp)
	   printf("q^{%d}", (int)mp);
	 else
	   printf("q^{%.1f}",mp);
       }
       if(a==1) printf("t");
       if(a != 0 && a != 1) {
	 double ap = a.print();
	 if(((int)ap)*1. == ap )
	   printf("t^{%d}", (int)ap);
	 else
	   printf("t^{%.1f}",a.print());
       }
   }
 }
 printf("\n");
 time_t endtime = time(NULL);
 printf("Total time elapsed: %ld seconds.\n", endtime - starttime);

 return 0;
}

// Actual Functions

int NumComp(int* white, int* black, int gridsize){
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


int WindingNumber(int x, int y, int* white, int* black){ // Return winding number around (x,y)
 int ret=0;
 for(int i=0; i<x; i++) {
  if ((black[i] >= y) && (white[i] < y)) ret++;
  if ((white[i] >= y) && (black[i] < y)) ret--;
 }
 return ret;
}

int MaslovGrading(int y [], int* white, int gridsize) {

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
 if( P % 4 != 0 )
   printf("Maslov Grading is not integer: %d/4\n", P);
 return (P/4);
}

bool ValidGrid(int* white, int* black, int gridsize) {
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


