// Computes Heegaard Floer knot homology of all 11 crossing non-alternating knots.
// Output is in latex format.
// Speed improvements due to Marc Culler
// Grid diagrams courtesy of Marc Culler's "gridlink" program
// Compile: g++ -O3 -o hfkbatch11 hfkbatch11.cpp
// Run: ./hfkbatch11

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

int gridsize = 0; // arc-index to be determined...

// Don't waste time computing factorials.  Look them up.
// Fill in the big ones later since g++ doesn't seem to like big constants
long long Factorial[16] = {
  1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,
  0,0,0};

int white[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
int black[12] = {0,0,0,0,0,0,0,0,0,0,0,0};

// Master list of grid diagrams for 11 crossing non-alternating knots
// All have arc-index <= 11.

int Diagrams[185][2][11]={{{1,5,6,2,4,7,9,8,0,3,10},{7,10,0,8,1,2,6,3,5,9,4}},
{{1,3,10,4,9,7,8,5,6,0,2},{6,9,8,0,1,10,3,2,4,5,7}},
{{2,4,8,3,9,6,7,5,0,10,1},{0,9,10,1,5,2,4,8,7,3,6}},
{{0,7,4,1,2,6,10,3,9,5,8},{6,5,9,8,0,1,4,7,2,10,3}},
{{8,1,6,2,3,7,10,4,5,9,0},{5,10,0,9,1,2,6,8,3,4,7}},
{{9,1,8,0,10,3,5,4,2,6,7},{5,7,4,6,2,9,0,1,8,10,3}},
{{2,6,9,3,7,5,0,8,1,4,10},{8,0,4,1,10,2,9,3,6,7,5}},
{{9,10,3,0,4,1,2,7,8,5,6},{4,5,9,7,10,6,8,1,3,0,2}},
{{7,3,4,5,6,9,8,10,0,1,2},{10,6,0,8,1,7,3,5,2,4,9}},   // Pointed out by Lenny Ng
{{2,4,8,5,6,7,9,0,10,3,1},{10,1,3,0,4,5,2,8,6,9,7}},
{{7,9,6,8,10,5,0,1,3,2,4},{5,2,3,4,7,9,8,10,0,6,1}},
{{9,2,3,6,5,0,4,8,7,1,-1},{5,8,1,2,7,6,9,3,0,4,-1}},
{{10,0,2,9,7,1,3,4,5,6,8},{6,8,7,5,10,9,0,1,2,3,4}},
{{4,6,9,10,0,2,8,1,3,5,7},{1,3,5,6,9,10,4,7,0,8,2}},
{{9,7,3,5,4,6,8,10,0,2,1},{4,0,9,2,1,3,5,6,8,10,7}},
{{3,5,6,8,7,10,1,9,2,0,4},{6,1,4,2,3,5,7,0,10,8,9}},
{{10,7,9,3,1,4,0,2,6,5,8},{4,2,6,8,5,9,7,10,1,3,0}},
{{9,3,0,8,1,5,4,6,10,2,7},{1,10,6,2,7,9,8,3,5,0,4}},
{{2,3,0,1,8,4,5,6,7,-1,-1},{5,8,6,7,2,0,1,3,4,-1,-1}},
{{4,5,6,0,1,2,8,9,7,3,-1},{2,9,1,3,7,5,4,6,0,8,-1}},
{{6,10,0,2,3,1,9,7,4,8,5},{4,5,8,10,7,6,3,0,9,2,1}},
{{1,6,4,7,5,9,3,0,8,2,10},{8,10,0,2,1,4,6,5,3,9,7}},
{{9,5,2,4,6,3,7,8,0,1,-1},{3,0,8,9,1,5,2,4,6,7,-1}},
{{9,1,5,0,6,7,3,4,2,8,-1},{4,7,8,2,1,5,9,0,6,3,-1}},
{{8,7,4,5,1,9,10,0,2,3,6},{5,0,8,3,6,4,7,9,10,1,2}},
{{1,8,7,9,2,0,3,6,5,10,4},{5,3,10,4,8,7,9,1,0,2,6}},
{{6,9,0,7,3,8,1,2,4,5,-1},{3,4,5,2,9,6,7,8,0,1,-1}},
{{4,7,5,2,0,8,10,9,1,3,6},{9,2,1,6,5,4,3,7,8,0,10}},
{{6,7,3,10,5,8,0,2,4,9,1},{9,4,0,2,1,3,6,8,10,5,7}},
{{1,4,8,2,7,9,5,10,6,3,0},{9,10,0,6,1,3,2,7,4,8,5}},
{{0,5,2,4,1,7,6,8,10,9,3},{7,1,10,9,3,2,0,4,5,6,8}},
{{4,10,0,2,8,3,9,6,5,7,1},{8,7,9,5,10,1,4,2,0,3,6}},
{{2,6,8,9,1,5,7,0,3,4,10},{8,0,3,4,6,10,2,5,9,1,7}},
{{2,1,0,3,5,9,8,6,7,4,10},{7,8,4,10,1,3,0,2,5,9,6}},
{{0,4,1,6,2,5,7,9,8,10,3},{7,10,8,3,0,1,4,5,2,6,9}},
{{9,7,4,6,8,0,5,1,2,3,10},{3,1,10,2,4,7,9,8,0,6,5}},
{{6,1,2,5,8,7,9,4,3,0,-1},{4,9,0,1,3,2,6,8,7,5,-1}},
{{0,1,8,4,2,3,6,7,5,-1,-1},{6,7,3,0,5,1,2,4,8,-1,-1}},
{{6,1,4,3,8,7,0,2,5,9,10},{2,9,10,5,4,1,6,8,0,3,7}},
{{1,7,8,2,0,9,10,3,5,4,6},{5,3,6,7,4,1,8,9,0,10,2}},
{{7,5,3,4,9,1,10,2,6,0,8},{0,10,8,6,5,7,4,9,1,3,2}},
{{4,0,9,5,6,10,2,8,1,3,7},{8,3,1,0,4,7,9,5,6,10,2}},
{{0,7,9,10,1,2,5,3,4,6,8},{5,4,6,8,9,0,1,7,2,10,3}},
{{1,6,5,2,0,3,10,8,4,7,9},{8,10,9,6,7,1,5,2,0,3,4}},
{{0,10,3,1,2,9,5,6,4,7,8},{7,5,8,4,10,3,1,0,9,2,6}},
{{5,1,4,6,10,8,7,9,2,0,3},{2,7,0,3,5,1,10,4,6,8,9}},
{{9,1,2,5,3,6,4,0,7,10,8},{7,6,8,10,1,2,9,5,3,4,0}},
{{9,0,3,5,4,2,1,6,8,7,-1},{6,8,9,1,0,7,3,2,5,4,-1}},
{{2,5,6,4,3,8,7,1,0,9,-1},{8,9,1,0,7,6,2,5,3,4,-1}},
{{3,1,4,2,6,9,7,5,0,8,-1},{7,5,0,8,1,3,2,9,6,4,-1}},
{{9,8,2,7,4,6,5,3,10,0,1},{3,5,0,1,8,2,9,6,4,7,10}},
{{6,7,9,4,0,10,3,1,2,5,8},{10,5,6,2,8,7,9,4,0,1,3}},
{{6,3,4,0,1,5,9,10,8,2,7},{8,7,1,3,6,2,4,5,0,9,10}},
{{7,2,6,9,8,1,10,3,5,4,0},{3,10,0,2,4,7,6,8,1,9,5}},
{{2,1,6,7,9,8,10,4,5,3,0},{10,3,0,4,6,2,7,9,1,8,5}},
{{1,3,0,2,8,9,4,6,10,7,5},{6,7,5,10,0,3,1,2,8,4,9}},
{{4,6,5,7,8,0,2,9,3,1,-1},{9,2,3,4,5,6,7,1,0,8,-1}},
{{7,4,3,0,2,1,5,6,10,9,8},{10,8,6,4,9,7,0,2,3,5,1}},
{{0,2,4,3,5,8,6,7,9,10,1},{3,8,1,10,2,4,0,5,6,7,9}},
{{1,6,2,3,10,9,4,5,8,7,0},{5,8,7,0,6,1,10,3,2,9,4}},
{{0,1,2,8,3,7,4,6,5,9,-1},{4,6,9,5,0,1,8,2,7,3,-1}},
{{6,4,8,9,10,1,5,3,2,7,0},{1,0,3,7,2,8,10,9,6,4,5}},
{{4,1,0,2,5,6,10,3,8,7,9},{8,3,7,10,9,1,4,0,5,2,6}},
{{7,10,0,1,4,2,9,8,5,3,6},{2,3,7,8,10,5,6,4,9,0,1}},
{{6,0,8,2,7,4,1,3,5,9,-1},{1,5,3,6,9,8,7,0,2,4,-1}},
{{9,2,3,7,4,5,6,0,1,10,8},{6,8,1,10,2,9,3,5,7,4,0}},
{{5,2,7,9,0,8,3,10,1,4,6},{10,8,4,5,6,1,9,2,7,0,3}},
{{2,5,0,4,10,8,1,3,9,7,6},{8,10,7,9,6,3,5,0,2,1,4}},
{{6,8,4,9,5,7,10,0,1,3,2},{3,0,1,7,8,2,4,6,9,10,5}},
{{1,4,8,5,6,7,0,2,3,9,-1},{7,2,3,9,1,4,5,8,0,6,-1}},
{{8,6,5,7,10,2,9,0,3,1,4},{2,1,8,3,6,7,4,5,10,9,0}},
{{9,5,6,10,1,0,2,7,3,4,8},{7,0,4,5,8,6,9,1,10,2,3}},
{{7,8,10,9,0,6,4,2,3,1,5},{4,5,3,1,7,10,8,9,0,6,2}},
{{10,4,5,1,9,2,7,6,0,8,3},{6,0,3,8,4,10,1,9,5,2,7}},
{{7,0,9,8,10,3,1,6,5,2,4},{3,5,2,4,7,8,9,0,10,6,1}},
{{9,8,4,6,2,5,7,3,1,10,0},{2,1,10,9,8,0,4,6,5,3,7}},
{{7,8,5,9,0,10,2,1,6,3,4},{3,6,2,4,5,7,9,8,10,0,1}},
{{10,9,0,8,3,7,6,2,4,5,1},{4,2,5,10,9,1,0,8,7,3,6}},
{{5,8,9,7,1,2,6,0,4,3,-1},{1,2,3,0,8,9,4,5,7,6,-1}},
{{4,0,6,7,8,1,5,9,2,10,3},{9,5,2,3,4,6,10,0,8,7,1}},
{{9,2,3,6,4,5,7,8,0,1,-1},{4,5,0,1,8,9,2,3,6,7,-1}}, // Pointed out by Lenny Ng
{{2,3,1,6,7,0,8,9,4,5,-1},{8,9,4,2,3,6,1,5,7,0,-1}},
{{10,7,5,4,0,2,3,1,8,9,6},{3,1,8,6,5,7,9,10,2,4,0}},
{{0,2,1,5,3,4,6,7,9,8,-1},{4,9,7,0,8,2,3,5,6,1,-1}},
{{2,4,10,0,9,1,8,3,7,5,6},{5,7,3,6,2,10,0,9,1,8,4}},
{{0,5,6,9,7,8,2,4,3,1,-1},{2,1,4,5,3,6,7,0,9,8,-1}},
{{1,5,6,3,2,9,4,7,8,0,10},{8,9,0,10,5,3,1,2,6,7,4}},
{{8,9,7,3,0,1,2,5,4,6,-1},{4,5,2,8,6,7,9,0,1,3,-1}},
{{3,4,0,6,7,1,8,9,10,2,5},{10,1,3,2,5,6,4,7,8,9,0}},
{{0,8,7,9,10,1,2,5,4,6,3},{5,3,2,4,6,7,0,9,8,1,10}},
{{7,1,10,4,3,6,2,5,8,9,0},{3,8,5,7,9,0,10,1,4,6,2}},
{{0,5,6,8,9,7,3,4,1,2,-1},{3,1,2,4,5,0,6,9,7,8,-1}},
{{10,1,9,0,4,2,6,3,7,5,8},{6,7,5,8,9,10,0,1,4,2,3}},
{{3,4,5,2,1,6,0,7,10,8,9},{8,0,10,9,3,2,5,1,6,4,7}},
{{4,7,1,8,0,2,5,6,3,-1,-1},{2,3,5,4,6,7,8,1,0,-1,-1}},
{{8,5,7,1,6,2,3,4,0,9,-1},{3,0,4,5,8,7,9,1,6,2,-1}},
{{3,6,7,8,5,2,10,9,1,0,4},{10,0,1,4,9,6,7,3,8,5,2}},
{{9,3,2,4,6,5,1,10,7,8,0},{1,0,8,10,3,9,6,2,4,5,7}},
{{4,7,5,8,6,0,2,9,1,3,-1},{2,3,9,4,1,5,7,6,8,0,-1}},
{{7,6,10,3,4,2,0,1,8,5,9},{2,8,7,9,1,5,3,10,4,0,6}},
{{1,3,5,6,0,7,2,4,10,8,9},{4,10,8,2,5,1,9,0,6,3,7}},
{{7,8,0,6,9,3,1,5,4,2,-1},{1,5,3,2,4,7,6,0,8,9,-1}},
{{8,9,10,1,2,5,3,6,4,0,7},{6,7,4,8,10,1,0,2,9,5,3}},
{{4,5,6,7,8,9,3,2,0,1,-1},{7,0,9,2,4,5,1,6,3,8,-1}},  // Pointed out by Lenny Ng
{{3,1,5,4,6,8,7,9,0,10,2},{10,7,2,0,1,5,3,6,8,4,9}},
{{3,8,6,7,4,0,5,1,9,2,-1},{0,1,9,2,8,6,3,7,4,5,-1}},
{{6,2,4,5,8,9,3,0,7,1,-1},{9,7,0,1,2,4,6,5,3,8,-1}},
{{10,6,9,0,4,1,5,2,3,7,8},{7,3,5,8,6,10,0,9,1,2,4}},
{{0,6,2,8,3,4,5,7,10,1,9},{8,9,10,1,0,2,3,4,6,7,5}},
{{9,8,2,1,5,3,4,6,0,7,-1},{6,1,0,7,9,8,2,3,5,4,-1}},
{{8,5,9,6,7,1,2,3,0,4,-1},{3,2,4,0,5,6,8,9,7,1,-1}},
{{2,4,5,7,6,10,3,0,1,9,8},{9,1,3,2,4,5,8,7,10,6,0}},
{{7,6,9,8,1,3,10,0,2,5,4},{2,0,5,3,7,9,4,8,6,10,1}},
{{2,5,10,7,9,3,6,8,0,4,1},{7,3,4,1,0,10,2,5,6,9,8}},
{{0,2,5,4,6,8,3,10,7,9,1},{3,8,9,7,1,5,6,4,2,0,10}},
{{9,8,5,1,4,7,6,2,0,3,-1},{4,2,0,6,8,3,9,5,7,1,-1}},
{{0,1,7,8,3,4,6,2,5,9,-1},{6,4,2,5,7,9,1,8,0,3,-1}},
{{6,8,0,1,2,4,3,5,7,-1,-1},{2,3,5,8,7,0,6,1,4,-1,-1}},
{{0,10,7,2,1,6,3,4,9,5,8},{6,5,0,9,8,10,7,2,3,1,4}},
{{6,7,8,5,4,3,9,2,10,1,0},{10,2,6,1,0,7,5,8,4,3,9}},
{{0,3,6,1,4,7,8,5,10,2,9},{8,5,9,10,2,0,3,1,6,7,4}},
{{5,6,8,7,9,1,0,3,2,4,-1},{3,4,5,0,2,8,6,9,7,1,-1}},
{{3,1,10,6,9,0,4,7,2,5,8},{7,4,5,3,2,8,9,1,10,0,6}},
{{9,4,10,0,3,5,6,1,2,8,7},{3,1,2,8,10,9,4,5,7,6,0}},
{{5,9,10,1,0,8,7,2,3,6,4},{3,4,6,9,7,5,10,8,1,2,0}},
{{0,2,1,3,6,4,7,5,8,9,-1},{7,8,5,9,1,0,3,2,4,6,-1}},
{{4,9,7,0,8,10,1,6,3,5,2},{10,2,3,5,4,7,8,0,9,1,6}},
{{10,7,8,9,1,4,5,0,3,2,6},{5,0,6,2,7,10,1,8,9,4,3}},
{{1,3,2,5,6,4,0,7,8,10,9},{8,10,0,1,3,9,5,4,6,7,2}},
{{7,6,0,4,1,2,10,9,3,8,5},{3,4,5,8,7,0,6,1,10,2,9}},
{{1,2,4,5,7,6,0,10,8,9,3},{5,6,10,3,1,9,4,2,0,7,8}},
{{9,7,4,0,2,6,3,1,5,8,-1},{6,1,8,5,9,0,7,4,3,2,-1}},
{{3,4,6,8,7,5,0,9,1,2,-1},{0,1,2,3,9,8,6,4,5,7,-1}},
{{5,7,6,9,2,8,0,4,1,3,-1},{2,4,1,3,7,5,6,9,8,0,-1}},
{{0,2,6,1,4,8,3,5,7,9,-1},{4,8,9,5,7,6,0,2,1,3,-1}},  // Lenny Ng
{{4,5,6,7,8,9,10,2,0,1,3},{2,0,3,5,1,6,8,9,4,7,10}},
{{10,6,2,5,3,4,7,8,0,9,1},{5,0,9,1,10,2,3,4,7,6,8}},
{{0,8,7,2,6,4,1,3,9,5,-1},{4,3,9,5,8,7,6,0,2,1,-1}},
{{9,4,7,2,0,6,8,3,1,5,-1},{2,1,3,5,8,9,4,6,7,0,-1}},
{{3,6,4,9,1,7,10,8,2,0,5},{1,2,8,0,5,3,4,6,9,7,10}},
{{1,8,0,6,7,5,2,3,4,9,-1},{7,2,5,1,4,8,9,6,0,3,-1}},
{{1,2,0,7,9,6,8,3,4,5,-1},{6,9,3,1,5,2,4,7,0,8,-1}},
{{7,4,6,1,5,8,0,9,2,3,-1},{2,0,9,7,3,4,6,5,8,1,-1}},
{{5,8,9,0,1,2,4,3,6,7,10},{1,3,7,5,10,8,0,9,2,4,6}},
{{1,9,3,8,4,6,7,0,5,2,-1},{7,5,0,1,9,3,2,4,8,6,-1}},
{{8,0,9,4,5,6,10,7,2,1,3},{2,7,1,8,0,3,5,4,9,6,10}},
{{4,2,0,6,3,5,7,8,9,1,-1},{8,5,4,1,9,0,2,3,6,7,-1}},
{{8,2,0,9,5,4,1,3,6,7,10},{6,7,3,2,10,8,5,9,0,1,4}},
{{1,5,6,8,7,4,3,9,2,10,0},{9,10,1,2,5,0,6,4,7,3,8}},
{{2,10,5,8,6,7,0,9,4,3,1},{7,8,9,1,10,3,5,2,0,6,4}},
{{4,9,0,2,8,5,6,10,1,3,7},{1,6,7,9,3,0,4,5,8,10,2}},
{{9,10,5,2,0,1,3,8,4,7,6},{4,8,9,6,7,10,0,5,2,1,3}},
{{0,7,10,8,9,3,4,2,1,5,6},{5,2,6,1,7,8,0,10,9,3,4}},
{{9,5,2,3,7,4,6,8,10,1,0},{4,10,9,1,2,0,3,5,7,8,6}},
{{0,6,10,7,8,9,3,5,4,2,1},{3,2,4,1,6,5,7,0,9,10,8}},
{{4,9,0,3,5,10,6,8,2,7,1},{8,2,7,9,1,4,0,5,10,3,6}},
{{3,6,8,2,0,9,10,5,4,1,7},{10,1,4,7,3,6,8,2,9,5,0}},
{{3,9,8,4,7,6,1,5,0,2,10},{1,0,2,10,3,9,8,7,4,6,5}},
{{4,6,10,0,9,7,8,2,5,1,3},{2,3,1,5,4,10,6,7,9,8,0}},
{{3,0,8,9,2,10,1,4,6,5,7},{10,6,3,7,8,5,9,0,2,1,4}},
{{6,10,7,4,0,5,1,3,2,8,9},{3,5,2,8,6,9,7,10,4,0,1}},
{{5,3,6,2,4,7,10,8,0,9,1},{10,9,1,8,0,3,2,6,5,4,7}},
{{5,3,6,4,1,7,0,8,9,10,2},{10,0,2,9,5,3,6,1,7,4,8}},
{{3,5,7,0,8,1,9,2,4,6,-1},{8,9,4,5,2,3,6,7,0,1,-1}}, // Lenny Ng
{{4,6,10,7,3,9,8,1,5,0,2},{8,1,5,0,6,4,10,9,2,3,7}},
{{2,3,4,8,0,9,10,7,6,5,1},{10,6,1,3,5,4,8,2,0,9,7}},
{{10,0,2,7,3,5,4,9,6,8,1},{6,9,10,1,0,2,8,7,3,5,4}},
{{4,8,10,5,6,0,7,3,1,2,9},{0,3,1,9,4,5,2,10,6,8,7}},
{{8,9,1,10,0,4,2,5,6,3,7},{2,6,7,3,8,10,9,0,4,1,5}},
{{0,4,2,8,6,7,5,9,1,10,3},{8,10,7,3,1,4,0,2,5,6,9}},
{{3,5,1,4,2,6,8,7,0,9,10},{7,9,8,10,0,3,5,1,4,2,6}},
{{9,0,7,3,8,6,4,5,10,2,1},{2,8,9,10,1,0,7,3,4,6,5}},
{{8,7,3,5,4,6,0,2,1,9,-1},{2,0,9,1,8,3,5,7,6,4,-1}},
{{0,1,10,3,5,2,4,6,7,9,8},{4,9,7,8,0,10,1,3,5,6,2}},
{{0,1,7,4,8,9,2,5,6,10,3},{5,9,3,10,0,6,7,1,4,2,8}},
{{9,7,8,1,4,2,3,10,5,0,6},{5,10,6,7,9,8,0,2,1,4,3}},
{{2,10,0,3,5,7,1,9,8,6,4},{7,6,8,10,2,4,5,3,1,0,9}},
{{1,9,0,2,8,3,5,4,6,7,10},{6,4,7,10,1,9,2,8,0,5,3}},
{{2,6,3,4,10,7,9,5,1,8,0},{9,1,8,2,6,3,0,10,7,4,5}},
{{6,3,4,5,7,10,8,9,0,2,1},{0,9,2,3,4,5,1,6,7,8,10}},
{{2,4,10,1,5,7,6,8,0,9,3},{7,1,2,9,0,3,10,4,6,5,8}},
{{8,10,0,4,2,3,1,9,7,5,6},{5,7,9,10,8,0,6,4,2,1,3}},
{{1,3,2,4,6,5,7,9,8,0,-1},{5,8,9,1,3,0,2,6,4,7,-1}},
{{0,2,6,4,1,3,5,9,7,10,8},{9,7,10,8,6,0,2,4,1,3,5}},
{{5,1,3,4,6,9,7,10,8,2,0},{10,8,0,2,3,5,1,4,6,7,9}}};


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
long long Index(int y []); // Returns the number of permutations appearing before y lexicographically
int WindingNumber(int x, int y); // Return winding number of the knot projection around (x,y)
int MaslovGrading(int y []);
int Find(vector<long long> & V, long long x); // Returns i if V[i]=x or -1 if x isn't V[i] for any i

// Main

int main(int argc, char *argv[]){

 Factorial[13] = 13*Factorial[12];
 Factorial[14] = 14*Factorial[13];
 Factorial[15] = 15*Factorial[14];

 // Select the proper gridsize
 for(int knot=1; knot <= 185; knot++) {

 if(Diagrams[knot-1][0][9] == -1) { 
  gridsize = 9;
 }
 if(Diagrams[knot-1][0][9] != -1 && Diagrams[knot-1][0][10] == -1) { 
  gridsize = 10;
 }
 if(Diagrams[knot-1][0][10] != -1 ) { 
  gridsize = 11;
 }
 if( knot == 9 ) gridsize = 12;

 int amin=0;  // To compute knot Floer homology it suffices to compute in non-negative Alexander grading
 int amax=20;

 // Fill out the white and black arrays from the master list
 for(int i=0; i < gridsize; i++) {white[i]=Diagrams[knot-1][0][i]; black[i]=Diagrams[knot-1][1][i];}

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
 //cout << "Populating the Graph...\n";
 long long edges=0;

 for(int index=0; index < label.size(); index++) {
  GetPerm(label[index],g);
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
 
 // Kill all the edges in the graph.
 // No live generator should ever have a dead generator on its "out" list
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
   GetPerm(label[i],g);
   int AGrading = AShift;
   for(int j=0; j<gridsize; j++) AGrading -= WN[j][g[j]];
   HomologyRanks[MaslovGrading(g)+30][AGrading+30]++;
  }
 }
 //cout << "Ranks of unshifted homology groups in Alexander grading [" << amin << "," << amax << "]:\n";
 /*
 for(int a=amax+30; a>=amin+30; a--) {
  for(int m=20; m<40; m++) {
   if(HomologyRanks[m][a] < 10) cout << "   ";
   if(HomologyRanks[m][a] >= 10 && HomologyRanks[m][a] < 100) cout << "  ";
   if(HomologyRanks[m][a] >= 100 && HomologyRanks[m][a] < 1000) cout << " ";
   cout << HomologyRanks[m][a];
  }
  cout << "\n";
 }
 */

 int HFKRanks [60][60]; // HFKRanks[i][j] will hold rank of HFK^ in Maslov grading=i-30 and Alexander grading=j-30
 for(int a=0; a<60; a++) { for(int m=0; m<60; m++) HFKRanks[m][a]=0; }
 
 // Reproduce HFK^ from HFK^ \otimes K^{gridsize-1} in non-negative Alexander grading
 for(int a=59; a>=0; a--) {
  for(int m=59; m>=0; m--) {
   if( HomologyRanks[m][a] > 0) {
    HFKRanks[m][a] = HomologyRanks[m][a];
    for(int i=0; i<=gridsize-1; i++) HomologyRanks[m-i][a-i] -= (HFKRanks[m][a] * Factorial[gridsize-1]) / (Factorial[i] * Factorial[gridsize-1-i]);
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
 cout << "11n_{" << knot << "} & ";
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
 cout << " \\\\ \n";
 time_t endtime = time(NULL);
 //cout << "Total time elapsed: " << (endtime-starttime) << " seconds.\n";

}
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
