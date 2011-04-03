#include <stdio.h>
#include <list>
#include <vector>

using std::list;
using std::vector;

class Generator {
 public:

 list<int> ones;
#ifdef FIELD_Z3
 list<int> twos;
#endif
 bool alive;
 Generator();
 ~Generator();
};

// Class Functions

Generator::Generator(){alive=1;};
Generator::~Generator(){};

void printMatrix( char* fileName, vector<Generator> mat ) {
  FILE *f = fopen(fileName, "w");
  for( int i = 0; i < mat.size(); i++ ) {
    for( list<int>::iterator j = mat[i].ones.begin(); j != mat[i].ones.end(); j++ )
      fprintf(f, "%d %d %e\n", i+1, *j+1, 1.);
#ifdef FIELD_Z3
    for( list<int>::iterator j = mat[i].twos.begin(); j != mat[i].twos.end(); j++ )
      fprintf(f, "%d %d %e\n", i+1, *j+1, -1.);
#endif
  }
  fclose(f);
}
