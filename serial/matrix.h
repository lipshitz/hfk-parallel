#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include <algorithm>
#include "gradings.h"

class GeneratorIn {
 public:
  std::list<int> in;
  bool alive;
  GeneratorIn();
  ~GeneratorIn();
};

class GeneratorOut {
 public:
  std::list<int> out;
  bool alive;
  GeneratorOut();
  ~GeneratorOut();
};


void getPerm(generator k, int h [], int gridsize); // Fills h with the k^th permutation (in lexico. order)
int rectIndex(int i, int j, int k, int l, int m, int gridsize);
generator getIndex( int *P, int gridsize );
int Find(std::vector<generator> & V, generator x);
generator factorial( int );
void printMatrix( char*, std::vector<GeneratorOut> );
int fillReduceKernel( std::vector<generator>&, std::vector<generator>&, int);
// Decides whether one of the four rectangles on the torus with corners at(xll,yll) and (xur,yur) contains no white or black dots
/*
 1 | 2 | 1
---+---+---
 3 | 0 | 3
---+---+---
 1 | 2 | 1
*/
bool RectDotFree(int xll, int yll, int xur, int yur, int which, int*, int*, int); 
void initRectangles(int gridsize, int *white, int *black);

#endif
