
// Computes Heegaard Floer knot homology of the Kinoshita-Terasaka knot KT_{2,1}.
// To calculate for different knots, modify the gridsize and the white and black arrays below
// and change the loops at lines 340 and 464.  Won't work for gridsize > 11 unless you have
// about 60GB of RAM, though, with more modification it could be made to work.

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
 bool alive;
 Generator();
 ~Generator();
};

// Globals

const int gridsize = 11; // arc-index

// Kinoshita-Terasaka KT_{2,1}
int white[11]={5,10,9,4,8,0,1,6,7,2,3};
int black[11]={0,6,1,7,10,2,5,9,3,4,8};

// Conway Mutant C_{2,1}
// int white[11]={10,9,3,4,5,8,6,7,1,2,0};
// int black[11]={6,1,7,0,3,10,9,2,4,8,5};

// 8_19
// int white[10] = {2,5,3,4,8,6,7,9,1,0};
// int black[10] = {0,9,1,2,3,4,5,6,7,8};

// 8_20
// int white[9] = {4,5,6,7,3,2,1,0,8};
// int black[9] = {0,1,4,5,8,6,3,7,2};

// 8_21
// int white [10] = {9,6,1,5,2,3,0,4,7,8};
// int black [10] = {1,2,3,7,4,6,5,8,9,0};

// 9_42
// int white[10] = {3,4,5,8,2,1,0,6,7,9};
// int black[10] = {0,1,3,4,9,5,7,8,2,6};

// 9_43
// int black[10] = {2,5,6,1,7,9,0,8,4,3};
// int white[10] = {7,8,5,4,3,6,2,1,0,9};

// 9_44
// int black [9] = {0,6,5,3,1,2,4,8,7};
// int white [9] = {5,4,8,6,7,0,1,3,2};

// 9_45
// int white[10] = {2,4,9,8,5,0,3,6,7,1};
// int black[10] = {0,1,2,3,9,4,7,8,5,6};

// 9_46
// int white [9] = {2,1,5,7,6,0,8,3,4};
// int black [9] = {0,6,8,2,3,5,4,7,1};

// 9_47
// int white[9]= {8,6,3,2,1,4,5,7,0};
// int black[9] = {2,1,0,7,5,6,8,3,4};

// 9_48
// int white[10]= {5,6,2,1,8,7,0,4,3,9};
// int black[10]= {7,1,0,5,3,9,8,2,6,4};

// 9_49
// int white[9]={3,2,8,1,7,0,6,5,4};
// int black[9]={0,6,3,5,2,4,1,8,7};

// 10_124
// int white [9] = {7,6,5,3,4,2,1,0,8};
// int black [9] = {3,1,8,0,7,6,5,2,4};

// 10_125
// int white[10] = {8,6,4,3,7,9,2,1,0,5};
// int black[10] = {4,3,1,0,5,6,8,7,2,9};

// 10_126
// int white[10] = {8,0,5,2,1,3,4,6,7,9};
// int black[10] = {1,7,9,8,5,0,2,3,4,6};

// 10_127
// int white[10] = {7,3,4,5,9,0,8,1,6,2};
// int black[10] = {0,1,2,3,4,6,5,7,9,8};

// 10_128
//int white[10] = {9,3,2,8,6,5,7,4,1,0};
//int black[10] = {2,1,4,5,3,0,9,8,7,6};

// 10_129
// int white[10] = {9,3,2,6,8,7,0,4,1,5};
// int black[10] = {2,1,4,3,5,9,6,8,7,0};

// 10_130
// int white[10] = {7,6,8,5,3,9,2,4,0,1};
// int black[10] = {3,0,4,7,6,5,8,1,2,9};

// 10_131
// int white[10] = {8,1,3,5,6,4,2,7,9,0};
// int black[10] = {2,4,6,7,9,8,0,1,5,3};

// 10_132
// int white[9] = {7,3,4,0,1,8,2,6,5};
// int black[9] = {4,6,1,2,5,3,7,0,8};

// 10_133
// int white[10] = {0,9,1,3,6,5,4,7,8,2};
// int black[10] = {5,4,8,0,9,1,2,3,6,7};

// 10_134
// int white[10] = {6,5,4,9,3,8,7,2,1,0};
// int black[10] = {1,8,7,6,5,4,0,9,3,2};

// 10_135
// int white[10] = {8,0,6,7,9,3,4,2,1,5};
// int black[10] = {2,4,1,5,6,0,8,7,3,9};

// 10_136
// int white[10] = {8,6,3,9,2,0,1,7,5,4};
// int black[10] = {0,4,6,2,1,5,9,3,7,8};

// 10_137
// int white[10] = {5,7,9,8,3,2,4,0,1,6};
// int black[10] = {2,1,6,9,3,4,8,5,7,0};

// 10_138
// int white[10] = {9,5,8,2,6,3,7,1,4,0};
// int black[10] = {2,1,3,5,0,9,4,6,8,7};

// 10_139
// int white[10] = {9,8,6,5,3,7,4,2,1,0};
// int black[10] = {1,5,3,2,0,9,8,7,6,4};

// 10_140
// int white[10] = {5,6,9,4,1,0,7,8,2,3};
// int black[10] = {0,2,5,8,6,4,9,3,7,1};

// 10_141
// int white[10] = {9,7,4,8,0,2,1,6,3,5};
// int black[10] = {2,1,9,6,3,7,5,4,8,0};

// 10_142
// int white[10] = {6,5,4,8,3,2,7,1,0,9};
// int black[10] = {2,1,7,6,5,0,9,8,4,3};

// 10_143
// int white[10] = {8,9,1,4,3,5,2,6,7,0};
// int black[10] = {2,7,8,0,9,1,4,3,5,6};

// 10_144
// int white[10] = {7,0,1,5,4,8,9,2,6,3};
// int black[10] = {1,2,4,3,6,5,7,8,0,9};

// 10_145
// int white[10] = {6,8,3,2,4,5,7,9,0,1};
// int black[10] = {0,4,6,7,1,8,3,2,5,9};

// 10_146
// int white[10] = {7,6,1,2,9,0,8,3,4,5};
// int black[10] = {2,0,3,8,7,5,4,6,1,9};

// 10_147
// int white[10] = {6,0,4,7,5,2,9,8,3,1};
// int black[10] = {2,5,6,3,8,7,4,1,0,9};

// 10_148
// int white[10] = {8,7,5,6,1,9,2,3,0,4};
// int black[10] = {1,0,8,2,3,5,4,7,6,9};

// 10_149
// int white[10] = {4,8,5,9,2,0,1,6,7,3};
// int black[10] = {2,3,1,4,6,5,7,8,9,0};

// 10_150
// int white[10] = {6,5,2,9,3,8,7,4,1,0};
// int black[10] = {3,8,7,6,1,2,0,9,5,4};

// 10_151
// int white[10] = {4,3,5,7,0,2,9,1,8,6};
// int black[10] = {0,9,1,4,5,8,7,6,3,2};

// 10_152
// int white[10] = {5,6,7,8,0,9,1,2,3,4};
// int black[10] = {0,2,1,5,3,7,8,4,6,9};

// 10_153
// int white[10] = {9,1,2,3,5,0,7,6,8,4};
// int black[10] = {3,5,8,7,9,6,4,2,1,0};

// 10_154
// int white[10] = {8,7,6,5,4,3,2,9,1,0};
// int black[10] = {1,3,9,0,7,5,8,4,6,2};

// 10_155
// int white[10] = {9,7,5,8,1,4,2,0,3,6};
// int black[10] = {5,4,1,2,3,0,9,6,7,8};

// 10_156
// int white[10] = {9,0,6,7,3,5,8,1,4,2};
// int black[10] = {4,7,8,1,9,0,2,3,6,5};

// 10_157
// int white[10] = {9,6,4,8,3,7,5,2,0,1};
// int black[10] = {4,3,0,2,1,9,8,6,5,7};

// 10_158
// int white[10] = {7,2,1,8,9,6,0,5,3,4};
// int black[10] = {1,0,5,3,7,2,4,8,6,9};

// 10_159
// int white[10] = {4,1,3,0,8,2,7,9,6,5};
// int black[10] = {2,8,7,6,4,9,1,5,3,0};

// 10_160
// int white[10] = {6,4,2,0,3,8,9,1,7,5};
// int black[10] = {3,9,8,7,1,6,2,5,4,0};

// 10_161
// int white[10] = {5,6,8,7,1,2,9,3,0,4};
// int black[10] = {0,1,5,2,3,4,6,8,7,9};

// 10_162
// int white[10] = {2,7,1,8,9,3,0,4,6,5};
// int black[10] = {0,3,4,2,7,6,5,8,1,9};

// 10_163
// int white[10] = {5,4,2,0,3,7,9,1,8,6};
// int black[10] = {3,9,7,5,8,1,2,6,4,0};

// 10_164
// int white[10] = {7,8,6,5,4,0,2,9,1,3};
// int black[10] = {0,5,2,9,1,3,7,4,6,8};

// 10_165
// int white[10] = {3,2,7,5,6,4,1,8,0,9};
// int black[10] = {0,6,4,9,3,8,7,2,5,1};

// Function Prototypes

bool RectDotFree(int xll, int yll, int xur, int yur, int which); /* Is the rectangle with lowerleft
 corner xll,yll and upper right corner xur,ur (on the torus their are 4 such, indexed by "which" as below) 
 free of black and white dots:

    |      |
 1  |  2   |  1
    |      |
----+------+-------
    |      |
 3  |  0   |  3
    |      |
----+------+-------
    |      |
 1  |  2   |  1
    |      |
 */
int Index(int y []); // Returns the number of permutations appearing before y lexicographically
int WindingNumber(int x, int y); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y []);
int Factorial(int n); // Returns n!
bool ValidGrid();

// Main

int main(int argc, char *argv[]){

 if(!ValidGrid()) return 0; // Check that the grid is valid
 //cout << "This is a valid " << gridsize << " by " << gridsize << " grid diagram.\n";
 int numgen = Factorial(gridsize);
 time_t starttime = time(NULL); // Used to record how long this takes

 // Record winding numbers around grid points for Alexander grading computations
 // Also record the Alexander grading shift
 int temp=0;
 for(int i=0; i<gridsize; i++) {
  temp += WindingNumber(i,black[i]);
  temp += WindingNumber(i,(black[i]+1) % gridsize);
  temp += WindingNumber((i+1) % gridsize, black[i]);
  temp += WindingNumber((i+1) % gridsize, (black[i]+1) % gridsize);
 } 
 for(int i=0; i<gridsize; i++) {
  temp += WindingNumber(i,white[i]);
  temp += WindingNumber(i,(white[i]+1) % gridsize);
  temp += WindingNumber((i+1) % gridsize, white[i]);
  temp += WindingNumber((i+1) % gridsize, (white[i]+1) % gridsize);
 } 
 const int AShift = (temp - 4 * gridsize + 4)/8;
 //cout << "Alexander Grading Shift: " << AShift << "\n";
 //cout << "Calculating winding number around all grid points.\n";
 int WN[gridsize][gridsize];
 for(int x=0; x<gridsize; x++) {
  for(int y=0; y<gridsize; y++) {
   WN[x][y] = WindingNumber(x,y);
  }
 }

 // Record for later use whether every possible rectangle has a black or white dot in it
 // This will speed boundary computations.
 //cout << "Computing which rectangles on the torus have no black or white dots inside.\n";
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

 int g[gridsize]; 
 vector<Generator> Graph(numgen); // Will hold the boundary data.  This is VERY RAM intensive.  Could be improved.
 //cout << "Allocating " << ((sizeof(Generator)*numgen)/1028) << "KB of memory for the boundary graph (more will be needed).\n";
 //cout << "Iterating through " << numgen << " generators to compute boundaries...\n";

 int edges=0;
 int count=0; // Loop through generators... change the loop for different gridsize
 for(g[0]=0; g[0]<gridsize; g[0]++) {
 for(g[1]=0; g[1]<gridsize; g[1]++) { if (g[1]==g[0]) continue;
 for(g[2]=0; g[2]<gridsize; g[2]++) { if (g[2]==g[1] || g[2]==g[0]) continue;
 for(g[3]=0; g[3]<gridsize; g[3]++) { if (g[3]==g[2] || g[3]==g[1] || g[3]==g[0]) continue; 
 for(g[4]=0; g[4]<gridsize; g[4]++) { if (g[4]==g[3] || g[4]==g[2] || g[4]==g[1] || g[4]==g[0]) continue;
 for(g[5]=0; g[5]<gridsize; g[5]++) { if (g[5]==g[4] || g[5]==g[3] || g[5]==g[2] || g[5]==g[1] || g[5]==g[0]) continue;
 for(g[6]=0; g[6]<gridsize; g[6]++) { if (g[6]==g[5] || g[6]==g[4] || g[6]==g[3] || g[6]==g[2] || g[6]==g[1] || g[6]==g[0])   continue;  
 for(g[7]=0; g[7]<gridsize; g[7]++) { if (g[7]==g[6] || g[7]==g[5] || g[7]==g[4] || g[7]==g[3] || g[7]==g[2] || g[7]==g[1]   || g[7]==g[0]) continue;
 for(g[8]=0; g[8]<gridsize; g[8]++) { if (g[8]==g[7] || g[8]==g[6] || g[8]==g[5] || g[8]==g[4] || g[8]==g[3] || g[8]==g[2] || g[8]==g[1] || g[8]==g[0]) continue; 
 for(g[9]=0; g[9]<gridsize; g[9]++) { if (g[9]==g[8] || g[9]==g[7] || g[9]==g[6] || g[9]==g[5] || g[9]==g[4] || g[9]==g[3] || g[9]==g[2] || g[9]==g[1] || g[9]==g[0]) continue; 
 for(g[10]=0; g[10]<gridsize; g[10]++) { if (g[10]==g[9] || g[10]==g[8] || g[10]==g[7] || g[10]==g[6] || g[10]==g[5] || g[10]==g[4] || g[10]==g[3] || g[10]==g[2] || g[10]==g[1] || g[10]==g[0]) continue; 

  int AGrading = AShift;
  for(int i=0; i<gridsize; i++) AGrading -= WN[i][g[i]];
  if (AGrading < 0) { Graph[count].alive=0; count++; continue; }
  
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
     int indexgij = Index(gij);
     Graph[count].out.push_back(indexgij);
     Graph[indexgij].in.push_back(Index(g));     
     edges++;
    }
   }
  }
  count++;
  // if (count % 500000 == 0) cout << "Finished " << count << " generators.\n";
 }}}}}}}}}}}

 //cout << "Done computing the graph.  Total edges (boundaries): " << edges << ".\n";

 // Kill all the edges in the graph.
 // No live generator should ever have a dead generator on its "out" list
 //cout << "Killing edges in the graph...\n";
 for(int i=0; i<numgen; i++) {
  //if(i % 500000 == 1 && i > 1) cout << "Finished " << i-1 << " generators.\n";
  if( (!Graph[i].alive) || Graph[i].out.size()==0) continue;
  int target = Graph[i].out.front(); // We plan to delete the edge from i to target...
  //cout << "Plan to delete edge from " << i << " to " << target << "\n";

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

  //PrintGraph( Graph);
 }
  
 int HomologyRanks [60][60]; // HomologyRanks[i][j] will hold rank of homology Maslov grading=i-30 and Alexander grading j-30
 for(int a=0; a<60; a++) { for(int m=0; m<60; m++) HomologyRanks[a][m]=0; }

 count=0; // Loop through generators...
 for(g[0]=0; g[0]<gridsize; g[0]++) {
 for(g[1]=0; g[1]<gridsize; g[1]++) { if (g[1]==g[0]) continue;
 for(g[2]=0; g[2]<gridsize; g[2]++) { if (g[2]==g[1] || g[2]==g[0]) continue;
 for(g[3]=0; g[3]<gridsize; g[3]++) { if (g[3]==g[2] || g[3]==g[1] || g[3]==g[0]) continue; 
 for(g[4]=0; g[4]<gridsize; g[4]++) { if (g[4]==g[3] || g[4]==g[2] || g[4]==g[1] || g[4]==g[0]) continue;
 for(g[5]=0; g[5]<gridsize; g[5]++) { if (g[5]==g[4] || g[5]==g[3] || g[5]==g[2] || g[5]==g[1] || g[5]==g[0]) continue;
 for(g[6]=0; g[6]<gridsize; g[6]++) { if (g[6]==g[5] || g[6]==g[4] || g[6]==g[3] || g[6]==g[2] || g[6]==g[1] || g[6]==g[0])   continue;  
 for(g[7]=0; g[7]<gridsize; g[7]++) { if (g[7]==g[6] || g[7]==g[5] || g[7]==g[4] || g[7]==g[3] || g[7]==g[2] || g[7]==g[1]   || g[7]==g[0]) continue;
 for(g[8]=0; g[8]<gridsize; g[8]++) { if (g[8]==g[7] || g[8]==g[6] || g[8]==g[5] || g[8]==g[4] || g[8]==g[3] || g[8]==g[2] || g[8]==g[1] || g[8]==g[0]) continue; 
 for(g[9]=0; g[9]<gridsize; g[9]++) { if (g[9]==g[8] || g[9]==g[7] || g[9]==g[6] || g[9]==g[5] || g[9]==g[4] || g[9]==g[3] || g[9]==g[2] || g[9]==g[1] || g[9]==g[0]) continue; 
 for(g[10]=0; g[10]<gridsize; g[10]++) { if (g[10]==g[9] || g[10]==g[8] || g[10]==g[7] || g[10]==g[6] || g[10]==g[5] || g[10]==g[4] || g[10]==g[3] || g[10]==g[2] || g[10]==g[1] || g[10]==g[0]) continue; 

 if(Graph[count].alive) { // This way we only compute Maslov for the surviving generators
  int AGrading = AShift;
  for(int x=0; x<gridsize; x++) AGrading -= WN[x][g[x]];
  HomologyRanks[MaslovGrading(g)+30][AGrading-30]++;  
 }
 count++;
 }}}}}}}}}}}


 int HFKRanks [60][60]; // HFKRanks[i][j] will hold rank of HFK^ in Maslov grading=i-30 and Alexander grading=j-30
 for(int a=0; a<60; a++) { for(int m=0; m<60; m++) HFKRanks[a][m]=0; }
 
 // Reproduce HFK^ from HFK^ \otimes K^{gridsize-1} in non-negative Alexander grading
 for(int a=59; a>=0; a--) {
  for(int m=59; m>=0; m--) {
   if( HomologyRanks[m][a] > 0) {
    HFKRanks[m][a] = HomologyRanks[m][a];
    for(int i=0; i<=gridsize-1; i++) HomologyRanks[m-i][a-i] -= (HFKRanks[m][a] * Factorial(gridsize-1)) / (Factorial(i) * Factorial(gridsize-1-i));
   }
  }
 }
 // Use symmetry to fill up HFKRanks in negative Alexander gradings
 for(int alex=-1; alex>=-9; alex--){ 
  for(int mas=-20; mas < 12; mas++) {
   HFKRanks[mas+30][alex+30] = HFKRanks[mas-2*alex+30 ][-alex+30];
  }
 }
 // Print Results
 cout << "Poincare polynomial of HFK^ is:\n";
 bool first=1;
 for(int a=-20; a<19; a++) {
  for(int m=-20; m<19; m++) {
   int rankam = HFKRanks[m+30][a+30];
   if(rankam > 1 && (!first) && m != 0 && a != 0) cout << " + " << rankam << "q^{" << m << "}t^{" << a << "}";
   if(rankam == 1 && (!first) && m != 0 && a != 0) cout << " + q^{" << m << "}t^{" << a << "}";
   if(rankam > 1 && (!first) && m == 0 && a != 0) cout << " + " << rankam << "t^{" << a << "}";
   if(rankam == 1 && (!first) && m == 0 && a != 0) cout << " + t^{" << a << "}";
   if(rankam > 1 && (!first) && m != 0 && a == 0) cout << " + " << rankam << "q^{" << m << "}";
   if(rankam == 1 && (!first) && m != 0 && a == 0) cout << " + q^{" << m << "}";
   if(rankam > 1 && (!first) && m == 0 && a == 0) cout << " + " << rankam;
   if(rankam == 1 && (!first) && m == 0 && a == 0) cout << " + 1";
   if(rankam > 1 && (first) && m != 0 && a != 0) {cout << rankam << "q^{" << m << "}t^{" << a << "}"; first=0;}
   if(rankam == 1 && (first) && m != 0 && a != 0) {cout << "q^{" << m << "}t^{" << a << "}"; first=0;}
   if(rankam > 1 && (first) && m == 0 && a != 0) {cout << rankam << "t^{" << a << "}"; first=0;}
   if(rankam == 1 && (first) && m == 0 && a != 0) {cout << "t^{" << a << "}"; first=0;}
   if(rankam > 1 && (first) && m != 0 && a == 0) {cout << rankam << "q^{" << m << "}"; first=0;}
   if(rankam == 1 && (first) && m != 0 && a == 0) {cout << "q^{" << m << "}"; first=0;}
   if(rankam > 1 && (first) && m == 0 && a == 0) {cout << rankam; first=0;}
   if(rankam == 1 && (first) && m == 0 && a == 0) {cout << "1"; first=0;}
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

int Index(int y []) {
 int index=0;
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

int Factorial(int n) {
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
