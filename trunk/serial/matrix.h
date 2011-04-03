#include <stdio.h>
#include <list>
#include <vector>

using std::list;
using std::vector;

class GeneratorIn {
 public:
 list<int> in;
 bool alive;
 GeneratorIn();
 ~GeneratorIn();
};

class GeneratorOut {
 public:
 list<int> out;
 bool alive;
 GeneratorOut();
 ~GeneratorOut();
};

// Class Functions

GeneratorIn::GeneratorIn(){alive=1;};
GeneratorIn::~GeneratorIn(){};
GeneratorOut::GeneratorOut(){alive=1;};
GeneratorOut::~GeneratorOut(){};

void printMatrix( char* fileName, vector<GeneratorOut> mat ) {
  FILE *f = fopen(fileName, "w");
  for( int i = 0; i < mat.size(); i++ )
    for( list<int>::iterator j = mat[i].out.begin(); j != mat[i].out.end(); j++ )
      fprintf(f, "%d %d %e\n", i+1, *j+1, 1.);
  fclose(f);
}
