#include "gradings.h"
#include <stdio.h>

const int NUM_GRADINGS = 2;

int main( int argc, char **argv ) {
  HalfInteger a[] = {-1,0};
  Gradings gradings( a );
  HalfInteger b[] = {2.5,3};
  gradings.addGeneratorToGrading(1,b);
  gradings.addGeneratorToGrading(2,b);
  printf("%lu\n", gradings.getGenerators(b)->size());
  b[0] = 0;
  gradings.addGeneratorToGrading(3,b);
  
  printf("A %lu\n", gradings.getGenerators(b)->size());
  HalfInteger c[] = {2.5,3};
  printf("C %lu\n", gradings.getGenerators(c)->size());
  return 0;
}
