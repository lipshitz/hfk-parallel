#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>
#include <string.h>
#include <algorithm>
using std::list;
using std::vector;
using std::cout;

int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

class Generator {
 public:
 list<int> out;
 list<int> in;
 bool alive;
 Generator();
 ~Generator();
};

// Globals

//const int gridsize = 12; // arc-index
//const int gridsize = 12;

int gridsize = 12; // arc-index
int *white;
int *black;
int default_white[12] = {9,5,11,7,8,1,10,4,0,3,2,6};
int default_black[12] = {1,0,4,3,2,6,5,9,8,11,7,10};
//const int gridsize = 12; // arc-index
// Trefoil
//int white[5] = {1, 2, 3, 4, 0};
//int black[5] = {4, 0, 1, 2, 3};

//int white[10] = {8,7,6,5,4,3,2,9,1,0};
//int black[10] = {1,3,9,0,7,5,8,4,6,2};

// Kinoshita-Terasaka KT_{2,1}
//int white[11]={5,10,9,4,8,0,1,6,7,2,3};
//int black[11]={0,6,1,7,10,2,5,9,3,4,8};

//int white[12] = {9,5,11,7,8,1,10,4,0,3,2,6};
//int black[12] = {1,0,4,3,2,6,5,9,8,11,7,10};


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
long long getIndexSwap(int y [], int, int); 
int WindingNumber(int x, int y); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y []);
int NumComp(); //Returns the number of components of the link.
bool ValidGrid();
int Find(vector<long long> & V, long long x); // Returns i if V[i]=x or -1 if x isn't V[i] for any i

// Main

int main(int argc, char *argv[]){

 char *knotFile = read_string( argc, argv, "-k", NULL );
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
 cout<<"Number of components:"<<numcomp<<"\n";

 if(!ValidGrid()) {cout << "Invalid grid!!\n"; return 0;} // Check that the grid is valid
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
 
 const int AShift = (temp - 4 * gridsize + 4)/8;
 
 cout << "Alexander Grading Shift: " << AShift << "\n";
 cout << "Matrix of winding numbers and Black/White grid:\n";
 for(int y=gridsize-1; y>=0; y--) {
  for(int x=0; x<gridsize; x++) {
   if(WN[x][y] >= 0) cout << " ";
   cout << WN[x][y];
  }
  cout << "   ";
  for(int x=0; x<gridsize; x++) {
   cout << " ";
   if(black[x]==y) cout << "X";
   if(white[x]==y) cout << "O";
   if(black[x] != y && white[x] != y) cout << " ";
  }
  cout << "\n";
 }
 

 // Record for later use whether every possible rectangle has a black or white dot in it
 // This will speed boundary computations.
 cout << "Computing which rectangles on the torus have no black or white dots inside.\n";
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
 for(int i=0; i<60; i++) NumGenByAGrading[i]=0;
 cout << "Searching through " << Factorial[gridsize] << " generators to compute Alexander gradings...\n";
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
       label.push_back(getIndex(g));
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
 cout << "Time to compute all Alexander gradings " << time(NULL)-agStartTime << "\n"; 

 for(int i=0;i<60;i++) {
  if(NumGenByAGrading[i]>0) cout << "Number of generators in Alexander grading " << (i-30) << ": "  << NumGenByAGrading[i] << "\n";
 }
 cout << "Total generators: " << label.size() << "\n";
 
 vector<Generator> Graph( label.size() ); // Will hold boundary data.
 cout << "Populating the Graph...\n";
 long long edges=0;

 int gij [gridsize];
 for(int index=0; index < label.size(); index++) {
  getPerm(label[index],g);
  bool firstrect;
  bool secondrect;
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

	  int indexgij = Find(label,Indexgij);
	  if(indexgij==-1) {cout << "Error with Alexander grading: " << Indexgij << "\n"; return 0; }
	  Graph[index].out.push_back( indexgij );
	  Graph[indexgij].in.push_back( index );     
	  edges++;
     }
    }
   }
 }

 cout << "Done computing the graph.  Total edges (boundaries): " << edges << ".\n";
 //PrintGraph(Graph);
 // Kill all the edges in the graph.
 // No live generator should ever have a dead generator on its "out" list
 cout << "Killing edges in the graph...\n";
 for(int i=0; i<Graph.size(); i++) {
  if(i % 1000000 == 0 && i > 0) cout << "Finished " << i << " generators.\n";
  if( (!Graph[i].alive) || Graph[i].out.size()==0) continue;
  int target = Graph[i].out.front(); // We plan to delete the edge from i to target...
 
  // For every m with i in dm, remove i from dm
  for(list<int>::iterator j=Graph[i].in.begin(); j!=Graph[i].in.end(); j++) {
   if(Graph[*j].alive) Graph[*j].out.remove(i);
  }
  Graph[i].alive = 0;
     
  // For every m with target in dm adjust dm appropriately
  for(list<int>::iterator j= Graph[ target ].in.begin(); j != Graph[ target ].in.end(); j++) {
   if( !Graph[*j].alive ) continue;
   for(list<int>::iterator k = Graph[i].out.begin(); k != Graph[i].out.end(); k++) {
    // Search for *k in the boundary of *j
    // If found, remove it; if not, add it to the boundary of *j
    list<int>::iterator search = find( Graph[*j].out.begin(), Graph[*j].out.end(), *k);
    if( search != Graph[*j].out.end() ) {
     Graph[ *j ].out.erase( search );
     if( *k != target) Graph[ *k ].in.remove( *j );
    }
    else {
     Graph[*j].out.push_back(*k);
     Graph[*k].in.push_back(*j);
    } 
   }
  }

  // For each a in di, remove i from the in list of a
  for(list<int>::iterator j=Graph[i].out.begin(); j != Graph[i].out.end(); j++) Graph[*j].in.remove(i);

  Graph[target].alive = 0;
  Graph[target].out.clear();
  Graph[target].in.clear();
  Graph[i].out.clear();
  Graph[i].in.clear();
 }

  
 int HomologyRanks [60][60]; // HomologyRanks[i][j] will hold rank of homology Maslov grading=i-30 and Alexander grading j-30
 for(int a=0; a<60; a++) { for(int m=0; m<60; m++) HomologyRanks[m][a]=0; }
 for(int i=0; i< Graph.size(); i++) {
  if(Graph[i].alive) {
   getPerm(label[i],g);
   int AGrading = AShift;
   for(int j=0; j<gridsize; j++) AGrading -= WN[j][g[j]];
   HomologyRanks[MaslovGrading(g)+30][AGrading+30]++;
  }
 }
 cout << "Ranks of unshifted homology groups in Alexander grading [" << amin << "," << amax << "]:\n";
 for(int a=amax+30; a>=amin+30; a--) {
  for(int m=20; m<40; m++) {
   if(HomologyRanks[m][a] < 10) cout << "   ";
   if(HomologyRanks[m][a] >= 10 && HomologyRanks[m][a] < 100) cout << "  ";
   if(HomologyRanks[m][a] >= 100 && HomologyRanks[m][a] < 1000) cout << " ";
   cout << HomologyRanks[m][a];
  }
  cout << "\n";
 }
 

 int HFKRanks [60][60]; // HFKRanks[i][j] will hold rank of HFK^ in Maslov grading=i-30 and Alexander grading=j-30
 for(int a=0; a<60; a++) { for(int m=0; m<60; m++) HFKRanks[m][a]=0; }
 
 // Reproduce HFK^ from HFK^ \otimes K^{gridsize-1} in non-negative Alexander grading
 for(int a=59; a>=0; a--) {
  for(int m=59; m>=0; m--) {
   if( HomologyRanks[m][a] > 0) {
    HFKRanks[m][a] = HomologyRanks[m][a];
    for(int i=0; i<=min(gridsize-numcomp,min(a,m)); i++) HomologyRanks[m-i][a-i] -= (HFKRanks[m][a] * Factorial[gridsize-numcomp]) / (Factorial[i] * Factorial[gridsize-numcomp-i]);
   }
  }
 }
 // Use symmetry to fill up HFKRanks in negative Alexander gradings
 for(int alex=-1; alex>=-9; alex--){ 
  for(int mas=-20; mas < 12; mas++) {
   HFKRanks[mas+30][alex+30] = HFKRanks[mas-2*alex+30 ][-alex+30];
  }
 }
 if(amin > 0) cout << "This Poincare polynomial is only valid in Alexander grading >= " << amin << ":\n";
 // Print Results
 bool first=1;
 for(int a=-20; a<19; a++) {
  for(int m=-20; m<19; m++) {
   int rankam = HFKRanks[m+30][a+30];
   if(rankam > 0) {
    if(!first) cout << "+";
    else first=0;
    if(rankam > 1 || (rankam==1 && a==0 && m==0) ) cout << rankam;
    if(m==1) cout << "q";
    if(m != 0 && m != 1) cout << "q^{" << m << "}";
    if(a==1) cout << "t";
    if(a != 0 && a != 1) cout << "t^{" << a << "}";
   }
  }
 }
 cout << "\n";
 time_t endtime = time(NULL);
 cout << "Total time elapsed: " << (endtime-starttime) << " seconds.\n";

 return 0;

}

// Class Functions

Generator::Generator(){alive=1;};
Generator::~Generator(){};

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

// Code below added by MC

// Maps a permutation of size n to an integer < n!
// See: Knuth, Volume 2, Section 3.3.2, Algorithm P

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

long long getIndexSwap( int *P, int I, int J ) {
  long long index = 0;
  for( int i = gridsize-2; i >= 0; i-- ) {
    int r = P[i];
    if( i == I )
      r = P[J];
    if( i == J )
      r = P[I];
    int m = 0;
    for( int j = 0; j < i; j++ ) {
      int l = P[j];
      if( j == I )
	l = P[J];
      if( j == J )
	l = P[I];
      if( l < r )
	m++;
    }
    index += Factorial[gridsize-1-i]*(r-m);
  }
  return index;
}

