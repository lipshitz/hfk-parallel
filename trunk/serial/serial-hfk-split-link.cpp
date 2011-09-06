#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <algorithm>
#include "matrix.h"
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

// Globals

int gridsize = 12; // arc-index
int *white;
int *black;
int default_white[12] = {9,5,11,7,8,1,10,4,0,3,2,6};
int default_black[12] = {1,0,4,3,2,6,5,9,8,11,7,10};


long long Factorial[16] = {
  1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,
  0,0,0};

// Function Prototypes

void getPerm(long long k, int h []); // Fills h with the k^th permutation (in lexico. order)
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
long long getIndex(int y []); // Returns the number of permutations appearing before y lexicographically
long long getIndexSwap(int y [], int , int); // Returns the number of permutations appearing before y lexicographically, if two entries were swapped
int WindingNumber(int x, int y); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y []);
int NumComp(); //Returns the number of components of the link.
bool ValidGrid();
int Find(vector<long long> & V, long long x); // Returns i if V[i]=x or -1 if x isn't V[i] for any i

// Main

int main(int argc, char *argv[]){

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
   white = default_white;
   black = default_black;
 }

 Factorial[13] = 13*Factorial[12];
 Factorial[14] = 14*Factorial[13];
 Factorial[15] = 15*Factorial[14];

 int amin=0;
 int amax=20;

int  numcomp = NumComp();
 printf("Number of components: %d\n", numcomp);

 if(!ValidGrid()) {printf("Invalid grid!!\n"); return 0;} // Check that the grid is valid
 time_t starttime = time(NULL); // Used to record how long this takes

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
 
 //const int AShift = (temp - 4 * gridsize + 4)/8;
 int AShift = (temp - 4 * gridsize + 4)/8;
 double AShiftD = (temp - 4 * gridsize + 4)/8.;
 if( AShiftD != AShift*1. ) {
   printf("Alexander Grading shift is not integer: %f\n", AShiftD);
   //amin -= 1;
   if( AShiftD > 0 )
     AShift += 1;
 }
 
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
 

 // Record for later use whether every possible rectangle has a black or white dot in it
 // This will speed boundary computations.
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

 
 int NumGenByAGrading[60]; // NumGenByAGrading[i] holds number of generators in A Grading i-30
 vector<long long> label; // label[i] will hold the number of perms lexicographically before the i^th generator
 // This will hold the generators, sorted by grading, the first index is Agrading+30, the second is Maslov-grading+30
 std::vector<long long> *generators[60][60];
 unsigned long long imageDimensions[60][60];
 unsigned long long kernelDimensions[60][60];
 int homologyDimensions[60][60];
 for( int i = 0; i < 60; i++ )
   for( int j = 0; j < 60; j++ )
     generators[i][j] = new std::vector<long long>();
 
 for(int i=0; i<60; i++) NumGenByAGrading[i]=0;
 printf("Searching through %lld generators to compute Alexander gradings...\n", Factorial[gridsize]);
 time_t agStartTime = time(NULL);

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
       NumGenByAGrading[AGrading+30]++;
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
 printf("Time to compute all Alexander gradings %ld\n", time(NULL)-agStartTime); 

 for(int i=0;i<60;i++) {
   if(NumGenByAGrading[i]>0) printf("Number of generators in Alexander grading %d: %d\n", (i-30), NumGenByAGrading[i]);
 }
 for( int i = 0; i < 60; i++ )
   for( int j = 0; j < 60; j++ )
     if (generators[i][j]->size())
       printf("Alexander grading %d, Maslov grading %d, num generators %lu\n", i-30, j-30, generators[i][j]->size());
 

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
     vector<GeneratorIn> GraphIn( rows.size() ); // Will hold boundary data.
     vector<GeneratorOut> GraphOut( cols.size() ); // Will hold boundary data.
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
	     GraphOut[index].out.push_back( indexgij );
	     GraphIn[indexgij].in.push_back( index );
	     edges++;
	   }
	 }
       }
     }
     if( printMatrices ) {
       char name[10];
       sprintf(name, "mat%d,%d.dat", I, J);
       printf("Writing to mat%d,%d\n", I, J);
       printf("Size is %lu by %lu\n", rows.size(), cols.size());
       printMatrix( name, GraphOut );
     }

     printf("Reducing %d %d\n", I, J);
     for( int i = 0; i < GraphOut.size(); i++ ) {
       if( (!GraphOut[i].alive) || GraphOut[i].out.size()==0 ) continue;
       int target = GraphOut[i].out.front(); // We plan to delete the edge from i to target
       GraphOut[i].alive = 0;
       for( list<int>::iterator j = GraphIn[target].in.begin(); j != GraphIn[target].in.end(); j++ ) {
	 if( !GraphOut[*j].alive ) continue;
	 for( list<int>::iterator k = GraphOut[i].out.begin(); k != GraphOut[i].out.end(); k++ ) {
	   list<int>::iterator search = find( GraphOut[*j].out.begin(), GraphOut[*j].out.end(), *k);
	   if( search != GraphOut[*j].out.end() ) {
	     GraphOut[*j].out.erase(search);
	     if( *k != target )  {
	       GraphIn[*k].in.remove(*j);
	     }
	   } else {
	     GraphOut[*j].out.push_back(*k);
	     GraphIn[*k].in.push_back(*j);
	   }
	 }
       }
       //for( list<int>::iterator j = GraphOut[i].out.begin(); j != GraphOut[i].out.end(); j++ )
       //	 GraphIn[*j].in.remove(i);
       
       //GraphIn[target].alive = 0;
       GraphIn[target].in.clear();
       //GraphOut[i].out.clear();
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
 if( P % 4 != 0 )
   printf("Maslov Grading is not integer: %d/4\n", P);
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
