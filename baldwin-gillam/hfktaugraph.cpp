// Computes Heegaard Floer knot homology.
// Some speed improvements due to Marc Culler
// Compile: g++ -o hfk -O3 hfk-mc.cpp
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
 list<int> nonreduced;
 int AGrading;
 bool alive;
 Generator();
 ~Generator();
};

// Globals

const int gridsize = 11; // arc-index
// Don't waste time computing factorials.  Look them up.
// Fill in the big ones later since g++ doesn't seem to like big constants
long long Factorial[16] = {
  1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,
  0,0,0};
  
//11n31
int white[11] = {9,6,4,5,7,8,3,2,1,10,0};
int black[11] = {2,3,1,8,10,6,7,0,9,4,5};

//11n_17
//int white[13] = {12,2,8,5,9,6,3,4,11,10,1,7,0};
//int black[13] = {8,10,0,1,6,4,5,2,3,7,9,12,11};

// Function Prototypes

void GetPerm(long long k, int h []); // Fills h with the k^th permutation (in lexico. order)
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
bool RectWhiteFree(int xll, int yll, int xur, int yur, int which);
long long Index(int y []); // Returns the number of permutations appearing before y lexicographically
int WindingNumber(int x, int y); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y []);
bool ValidGrid();
int Find(vector<long long> & V, long long x); // Returns i if V[i]=x or -1 if x isn't V[i] for any i

// Main

int main(int argc, char *argv[]){

 Factorial[13] = 13*Factorial[12];
 Factorial[14] = 14*Factorial[13];
 Factorial[15] = 15*Factorial[14];

 int amin=-1;
 int amax=20;


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
  

 // Record for later use whether every possible rectangle has a black or white dot in it
 // This will speed boundary computations.
 
 bool Rectangles[gridsize][gridsize][gridsize][gridsize][4];
 for(int xll=0; xll < gridsize; xll++) {
  for(int xur=xll+1; xur < gridsize; xur++) {
   for(int yll=0; yll < gridsize; yll++) {
    for(int yur=yll+1; yur < gridsize; yur++) {
     Rectangles[xll][yll][xur][yur][0] = RectWhiteFree(xll,yll,xur,yur,0);
     Rectangles[xll][yll][xur][yur][1] = RectWhiteFree(xll,yll,xur,yur,1);
     Rectangles[xll][yll][xur][yur][2] = RectWhiteFree(xll,yll,xur,yur,2);
     Rectangles[xll][yll][xur][yur][3] = RectWhiteFree(xll,yll,xur,yur,3);     
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
 int g[gridsize]; 
 //cout << "Iterating through " << Factorial[gridsize] << " generators to compute Alexander gradings...\n";
 
 long long count=0;
 // Loop through generators... change the loop for different gridsize
 
 // This array is a factorial base counter used for counting through
 // permutations. - MC
 short counter[gridsize-1];
 for(int i=0; i<gridsize-1; i++) counter[i] = 0;

 for(count = 0; count < Factorial[gridsize]; count++) {
  NextPerm(counter,g);
  int AGrading = AShift;
  for(int i=0; i<gridsize; i++) AGrading -= WN[i][g[i]];
  if (AGrading >= amin && AGrading <= amax) {
   label.push_back(count); 
   NumGenByAGrading[AGrading+30]++;
  }
 }
 
/*
 for(int i=0;i<60;i++) {
  if(NumGenByAGrading[i]>0) cout << "Number of generators in Alexander grading " << (i-30) << ": "  << NumGenByAGrading[i] << "\n";
 }
 cout << "Total generators: " << label.size() << "\n";
 */
 
 vector<Generator> Graph( label.size() ); // Will hold boundary data.
 for(int i=0; i<Graph.size(); i++) {Graph[i].nonreduced.push_back(i);}
 //cout << "Populating the Graph...\n";
 long long edges=0;

 for(int index=0; index < label.size(); index++) {
  GetPerm(label[index],g);
  int AlexanderGrading = AShift;
  for(int i=0; i<gridsize; i++) AlexanderGrading -= WN[i][g[i]];
  Graph[index].AGrading = AlexanderGrading;
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
     if(firstrect != secondrect && -WN[i][g[i]]-WN[j][g[j]]+WN[i][g[j]]+WN[j][g[i]]==0) { 
    // Exactly one rectangle is a boundary and Alex(g)=Alex(gij)) 
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
	  long long Indexgij = Index(gij);
	  int indexgij = Find(label,Indexgij);
	  if(indexgij==-1) {cout << "Error with Alexander grading:\n"; return 0; }
	  Graph[index].out.push_back( indexgij );
	  Graph[indexgij].in.push_back( index );     
	  edges++;
     }
    }
   }
 }

 //cout << "Done computing the graph.  Total edges (boundaries): " << edges << ".\n";
 // Output the whole graph
 
 for(int i=0; i < Graph.size(); i++) {
  cout << label[i] << "\n";
  cout << Graph[i].AGrading << "\n";
  for(list<int>::iterator j = Graph[i].out.begin(); j != Graph[i].out.end(); j++) cout << (*j) << "\n";
  cout << "-1\n";
  for(list<int>::iterator j = Graph[i].in.begin(); j != Graph[i].in.end(); j++) cout << (*j) << "\n";
  if( i==Graph.size()-1) cout << "-9999\n";
  else cout << "-999\n";
 }
 
 //cout << "Killing edges in the graph...\n";
 for(int i=0; i<Graph.size(); i++) {
  //if(i % 1000000 == 0 && i > 0) cout << "Finished " << i << " generators.\n";
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
   
    // Replace Graph[*j].nonreduced with the symmetric difference of
   //  Graph[*j].nonreduced and Graph[i].nonreduced
   for(list<int>::iterator l = Graph[i].nonreduced.begin(); l != Graph[i].nonreduced.end(); l++) {
    list<int>::iterator look = find( Graph[*j].nonreduced.begin(), Graph[*j].nonreduced.end(), *l);
    if(look != Graph[*j].nonreduced.end()) Graph[*j].nonreduced.erase(look);
    else { 
     Graph[*j].nonreduced.push_back(*l);
     //if( Graph[*j].AGrading != Graph[*l].AGrading ) {cout << "Error!!!\n"; return 0;}
    }    
   }
   
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
  Graph[target].nonreduced.clear();
  Graph[i].nonreduced.clear();
  Graph[i].out.clear();
  Graph[i].in.clear();
 }

 
 vector<Generator> E1Homology;
 for(int i=0; i<Graph.size(); i++) {
  if(Graph[i].alive && Graph[i].AGrading > amin) { 
   Graph[i].out.clear();
   E1Homology.push_back( Graph[i] );
  }
 }

  
 // Calculate the boundaries of the generators in one lower Alexander grading
 //cout << "Calculating d^1 boundaries of generators in the E_1 term:\n";
 for(int index=0; index < E1Homology.size(); index++) {
  for(list<int>::iterator actual = E1Homology[index].nonreduced.begin(); actual != E1Homology[index].nonreduced.end(); actual++) {
   GetPerm(label[ *actual ],g);
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
     if(firstrect!=secondrect) { 
     // Exactly one rectangle is a boundary.
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
      int indexgij = Find(label,Index(gij));
      if(indexgij != -1 && (Graph[indexgij].AGrading == E1Homology[index].AGrading || Graph[indexgij].AGrading+1 == E1Homology[index].AGrading)) {
       //if(MaslovGrading(g) != MaslovGrading(gij)+1) {cout << "Error with Maslov Grading!\n"; return 0;}
       list<int>::iterator search = find( E1Homology[index].out.begin(), E1Homology[index].out.end(), indexgij);
       if (search != E1Homology[index].out.end()) E1Homology[index].out.erase(search);
       else {
        E1Homology[index].out.push_back( indexgij );
        //cout << "Added " << indexgij << " to the boundary of " << i << "\n";
       }
      }     
     }
    }
   }
  }
 }


 // Check that each entry of E1Homology is actually closed in F_i/F_{i-1}
 //cout << "Checking that each E1 generator is actually closed in F_i/F_{i-1}...\n";
 for(int i=0; i < E1Homology.size(); i++) {  
  for(list<int>::iterator j = E1Homology[i].out.begin(); j != E1Homology[i].out.end(); j++) {
   if( E1Homology[i].AGrading != Graph[ *j ].AGrading+1) {
    cout << "Error: the " << i << "th generator (Alexander grading " << E1Homology[i].AGrading << ") seems to have boundary in ";
    cout << "Alexander grading " << Graph[ *j ].AGrading << "\n";
    return 0;
   }
  }  
 }
 
 
 //cout << "Outputting E1 boundary data....\n";
 
 for(int i=0; i<E1Homology.size(); i++) {
  GetPerm(label[ E1Homology[i].nonreduced.front() ],g);
  cout << MaslovGrading(g) << "\n";
  cout << E1Homology[i].AGrading << "\n"; 
  for(list<int>::iterator k = E1Homology[i].out.begin(); k != E1Homology[i].out.end(); k++) cout << (*k) << "\n";
  cout << "-999\n";
 }
 cout << "-9999\n"; 
 
  
 return 1;

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
long long Index(int P []) {
 long long index=0;
 int temp, m, r = gridsize;
 while (r > 0){
   for (m=0; m < r; m++){
     if (r - P[m] == 1)
       break;
   }
   index = index*(r) + m;
   r -= 1;
   temp = P[r];
   P[r] = P[m];
   P[m] = temp;
 }
 return index;
}

// Inverse mapping, from integers < n! to permutations of size n
// Writes the permutation corresponding to N into the array P.
void GetPerm(long long N, int P []) {
  int r, m, temp;
  for(int i=0; i<gridsize; i++) P[i]=i;
  r = 1;
  while (r < gridsize)  {
    m = N%(r+1);
    N = N/(r+1);
    temp = P[r];
    P[r] = P[m];
    P[m] = temp;
    r += 1;
  }
  return;
}

// Generator for permutations.  Inputs a factorial based counter and
// an array.  Writes the permutation indexed by the counter into the
// array, and then increments the counter.
void NextPerm(short counter[], int P[]) {
  int r, m, temp, i;
  for(i=0; i<gridsize; i++) P[i]=i;
  r = 1;
  while (r < gridsize) {
    m = counter[r-1];
    temp = P[r];
    P[r] = P[m];
    P[m] = temp;
    r += 1;
  }
  for (i=0; i<gridsize-1; i++) {
    counter[i] += 1;
    if (counter[i] == i+2)
      counter[i] = 0;
    else
      break;
  }
  return;
}
