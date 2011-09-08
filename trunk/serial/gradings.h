#ifndef GRADINGS_H
#define GRADINGS_H

#include "half-integer.h"
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <set>

/*
  This class holds information about the gradings, and the list of generators in each.  Gradings arrays of HalfInteger's, all of which better be the same length, set by NUM_GRADINGS.
 */

// this needs to be set by whatever program uses this to the number of gradings that we are keeping track of
extern const int NUM_GRADINGS;

// change this if you want to represent generators by something other than long long's
typedef long long generator;

struct ltgrad {
  bool operator()(const HalfInteger* g1, const HalfInteger* g2) {
    for( int i = 0; i < NUM_GRADINGS; i++ ) {
      if( g1[i] < g2[i] )
	return 1;
      if( g1[i] > g2[i] )
	return 0;
    }
    return 0;
  }
};

typedef std::set<const HalfInteger*, ltgrad>::iterator gradingIterator;

class Gradings {
 public:
  Gradings(const HalfInteger *bm);

  void addGeneratorToGrading( const generator gen, const HalfInteger *grad );

  std::vector<generator>* getGenerators( const HalfInteger* grad );

  std::vector<generator>* getBoundaryGenerators( const HalfInteger* grad );

  gradingIterator getGradingIterator();
  gradingIterator getGradingIteratorEnd();

  void setKernelDimension(const HalfInteger* grad, int dim);

  void setImageDimension(const HalfInteger* grad, int dim);

  int getHomologyDimension( const HalfInteger* grad );

  void calculateHomology();

  ~Gradings();

 private:
  std::map<const HalfInteger*, std::vector<generator>*, ltgrad> generatorsOfGrading;
  std::set<const HalfInteger*, ltgrad> gradings; // this will hold the gradings that have been defined, so we know which to save, and what to clean up
  std::map<const HalfInteger*, int, ltgrad> kernelDimensions;
  std::map<const HalfInteger*, int, ltgrad> imageDimensions;
  std::map<const HalfInteger*, int, ltgrad> homologyDimensions;
  const HalfInteger *boundaryMap;
  std::vector<generator> emptyGenerators; // this vector should always be empty.

  int getKernelDimension( const HalfInteger* grad );

  int getImageDimension( const HalfInteger* grad );

  // be sure to delete[] the output of this when done with it.
  HalfInteger* getBoundaryGrading(const HalfInteger* grad ) const;
  // be sure to delete[] the output of this when done with it.
  HalfInteger* getInverseBoundaryGrading(const HalfInteger* grad ) const;
};

#endif
