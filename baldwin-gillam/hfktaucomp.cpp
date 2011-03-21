
// Computes Heegaard Floer knot homology.
// Compile: g++ -o hfk hfk21.cpp
// TO do other knots, edit the global values before compiling

#include <time.h>
#include <iostream>
#include <vector>
#include <list>
using std::list;
using std::vector;
using std::cout;
using std::cin;


class Generator {
 public:
 list<int> out;
 list<int> in;
 long long label;
 bool alive;
 int AGrading;
 Generator();
 ~Generator();
};

// Globals

const int gridsize = 11; // arc-index

//11n31
int white[11] = {9,6,4,5,7,8,3,2,1,10,0};
int black[11] = {2,3,1,8,10,6,7,0,9,4,5};


// Function Prototypes

void GetPerm(long long k, int h []); // Fills h with the k^th permutation (in lexico. order)
bool RectDotFree(int xll, int yll, int xur, int yur, int which); 
bool RectWhiteFree(int xll, int yll, int xur, int yur, int which);
long long Index(int y []); // Returns the number of permutations appearing before y lexicographically
int WindingNumber(int x, int y); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y []);
long long Factorial(int n); // Returns n!
bool ValidGrid();
int Find(vector<long long> & V, long long x); // Returns i if V[i]=x or -1 if x isn't V[i] for any i

// Main

int main(int argc, char *argv[]){

 int MG=0;
 int AG=1;

 cout << "Reading in the graph...\n";
 vector<Generator> Graph;
 long long l;
 cin >> l;
 while(1) {
  Generator Gen;
  Gen.label = l;
  cin >> Gen.AGrading;
  int number;
  cin >> number;
  while(number != -1) {
   Gen.out.push_back(number);
   cin >> number;
  }
  cin >> number;
  while(number != -999 && number != -9999) {
   Gen.in.push_back(number);
   cin >> number;
  }
  Graph.push_back(Gen);
  if(number == -9999) break;
  cin >> l;
 }

 cout << "Number of vertices in the graph: " << Graph.size() << "\n";

 cout << "Adding some E1 boundaries in Maslov grading " << MG << " Alexander grading " << AG << "...\n";
 int previousdim = 0;
 int numberadded=0;
 int mgrading;
 int agrading;
 cin >> mgrading;
 while( mgrading != -9999) {
  Generator Gen;
  bool good = 1;
  if (mgrading != MG) good=0;
  cin >> agrading;
  if(agrading == AG-1) previousdim++;
  if (agrading != AG) good=0;
  Gen.AGrading = agrading;
  int o;
  cin >> o;
  while( o != -999) {
   Gen.out.push_back(o);
   cin >> o;
  }
  if(good) {
   Graph.push_back(Gen);
   numberadded++;
   for(list<int>::iterator j = Gen.out.begin(); j != Gen.out.end(); j++) {
    if( (*j) >= Graph.size() ) cout << "Error.  Trying to look at generator " << (*j) << "\n";
    Graph[(*j)].in.push_back( Graph.size() -1 );
   }
  }

  cin >> mgrading;
   
 }

 cout << "Added " << numberadded << " generators.\n";
 cout << "Before adding d^1 differentials the homology in Alexander grading " << (AG-1) << " has dimension ";
 cout << previousdim << ".\n";
 
 // Kill all the edges in the graph.
 // No live generator should ever have a dead generator on its "out" list
 // Keep track of the isomorphism
 cout << "Killing edges in the graph...\n";
 for(int i=0; i<Graph.size(); i++) {
  if(i % 1000000 == 0 && i > 0) cout << "Finished " << (i/1000000) << " million generators.\n";
  if(!Graph[i].alive || Graph[i].out.size()==0) continue;
  int target = Graph[i].out.front(); // We plan to delete the edge from i to target...
  //cout << "Deleting edge from " << i << " to " << target << "\n";

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
    // cout << "Searching for " << *k << " in out list of " << *j << "\n";
    list<int>::iterator search = find( Graph[*j].out.begin(), Graph[*j].out.end(), *k);
    if( search != Graph[*j].out.end() ) {
     //cout << "Found it.  Erasing " << *k << " from out list of " << *j << " and " << *j << " from the in list of " << *k << ".\n";
     Graph[ *j ].out.erase( search );
     if( *k != target) Graph[ *k ].in.remove( *j );
    }
    else {
     //cout << "Didn't find it... adding " << *k << " to out list of " << *j << " and " << *j << " to the in list of " << *k << ".\n";
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

 int newdimension = 0;
 for(int i=0; i < Graph.size(); i++) {
  if( Graph[i].alive && Graph[i].AGrading == AG-1 ) newdimension++;
 }
 cout << "Rank of homology in Alexander grading " << (AG-1) << " after accounting for d^1 differentials is ";
 cout << newdimension << ".\n";
 cout << "Rank of d^1:E_1(" << MG << "," << AG << ") -> E_1(" << (MG-1) << "," << (AG-1) << ") is ";
 cout << (previousdim-newdimension) << ".\n";


 time_t endtime = time(NULL);
 //cout << "Total time elapsed: " << (endtime-starttime) << " seconds.\n";



 return 0;


}





// Class Functions

Generator::Generator(){alive=1;};
Generator::~Generator(){};

// Actual Functions

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

long long Index(int y []) {
 long long index=0;
 for(int i=0; i<gridsize; i++) {
  int shiftati = y[i];
  for(int j=0; j<i; j++) {
   if(y[i] > y[j]) shiftati--;
  }
  index += shiftati * Factorial(gridsize-i-1);
 }
 return index;
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

long long Factorial(int n) {
 switch(n) {
  case 0: return 1;
  case 1: return 1;
  case 2: return 2;
  case 3: return 6;
  case 4: return 24;
  case 5: return 120;
  case 6: return 720;
  case 7: return 5040;
  case 8: return 40320;
  case 9: return 362880;
  case 10: return 3628800;
  case 11: return 39916800;
  default: return (n*Factorial(n-1));
 }
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

bool RectWhiteFree(int xll, int yll, int xur, int yur, int which) {
 bool dotfree = 1;
 switch (which) {
  case 0: 
   for(int x=xll; x<xur && dotfree; x++) {
    if (white[x] >= yll && white[x] < yur) dotfree = 0;
   }
   return dotfree;
  case 1:
   for(int x=0; x<xll && dotfree; x++) {
    if (white[x] < yll || white[x] >= yur) dotfree = 0;
   }
   for(int x=xur; x<gridsize && dotfree; x++) {
    if (white[x] < yll || white[x] >= yur) dotfree = 0;
   }
   return dotfree;
  case 2:
   for(int x=xll; x<xur && dotfree; x++) {
    if (white[x] < yll || white[x] >= yur) dotfree = 0;
   }
   return dotfree;
  case 3:
   for(int x=0; x<xll && dotfree; x++) {
    if (white[x] >= yll && white[x] < yur) dotfree = 0;
   }
   for(int x=xur; x<gridsize && dotfree; x++) {
    if (white[x] >= yll && white[x] < yur) dotfree = 0;
   }
   return dotfree;
 }
 return 0; //Error!
}

void PrintGraph(vector<Generator> & G) { // This is diagnostic only
 for(int i=0; i<G.size(); i++) {
  if(!G[i].alive) continue;
  cout << i << " in=";
  for(list<int>::iterator j=G[i].in.begin(); j != G[i].in.end(); j++) {
   cout << *j << " ";
  }
  cout << "out=";
  for(list<int>::iterator j=G[i].out.begin(); j != G[i].out.end(); j++) {
   cout << *j << " ";
  }
  cout << "\n";
 }
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


void GetPerm(long long k, int h []) {
 int g[gridsize-1];
 for(int i=0; i<gridsize; i++) h[i]=i;
 for(int i=1; i<gridsize; i++) {
  g[gridsize-1-i] = k/(Factorial(gridsize-i));
  k-=k/(Factorial(gridsize-i))*Factorial(gridsize-i);  
 }
 for (int i=1; i<gridsize; i++){
  int t[gridsize];
  for (int l=0;l<gridsize;l++) t[l] = h[l];
   if (g[gridsize-i-1] != 0) {
    h[i-1]=h[i-1+g[gridsize-i-1]];
    for(int p = 0; p<g[gridsize-i-1];p++) h[i+p] = t[i+p-1];
   }
  }
 return;
}

