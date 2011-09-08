#include "matrix.h"

bool *Rectangles;

GeneratorIn::GeneratorIn(){alive=1;};
GeneratorIn::~GeneratorIn(){};
GeneratorOut::GeneratorOut(){alive=1;};
GeneratorOut::~GeneratorOut(){};

void printMatrix( char* fileName, std::vector<GeneratorOut> mat ) {
  FILE *f = fopen(fileName, "w");
  for( int i = 0; i < mat.size(); i++ )
    for( std::list<int>::iterator j = mat[i].out.begin(); j != mat[i].out.end(); j++ )
      fprintf(f, "%d %d %e\n", i+1, *j+1, 1.);
  fclose(f);
}

int fillReduceKernel( std::vector<generator> &cols, std::vector<generator> &rows, int gridsize ) {
  if( cols.size() == 0 || rows.size() == 0 )
    return cols.size();
  int g[gridsize];
  std::vector<GeneratorIn> GraphIn( rows.size() ); // Will hold boundary data.
  std::vector<GeneratorOut> GraphOut( cols.size() ); // Will hold boundary data.
  generator edges=0;
  for(int index=0; index < cols.size(); index++) {
    getPerm(cols[index],g,gridsize);
    bool firstrect, secondrect;
    for(int i=0; i<gridsize; i++) {
      for(int j=i+1; j<gridsize; j++) {
	if(g[i]<g[j]) {
	  firstrect = Rectangles[rectIndex(i,g[i],j,g[j],0,gridsize)];
	  for(int k=i+1; k<j && firstrect; k++) {
	    if(g[i] < g[k] && g[k] < g[j]) firstrect=0;
	  }
	  secondrect = Rectangles[rectIndex(i,g[i],j,g[j],1,gridsize)];
	  for(int k=0; k<i && secondrect; k++) {
	    if(g[k]<g[i] || g[k] > g[j]) secondrect=0;
	  }
	  for(int k=j+1; k<gridsize && secondrect; k++) {
	    if(g[k]<g[i] || g[k] > g[j]) secondrect=0;
	  }
	}
	if(g[j]<g[i]) {
	  firstrect = Rectangles[rectIndex(i,g[j],j,g[i],2,gridsize)];
	  for(int k=i+1; k<j && firstrect; k++) {
	    if(g[k]<g[j] || g[k] > g[i]) firstrect=0;
	  }
	  secondrect = Rectangles[rectIndex(i,g[j],j,g[i],3,gridsize)];
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
	  generator Indexgij = getIndex(gij, gridsize);
	  int indexgij = Find(rows,Indexgij);
	  if(indexgij==-1) {printf("Error with Alexander grading: %lld\n", Indexgij); 
	    exit(-1); }
	  GraphOut[index].out.push_back( indexgij );
	  GraphIn[indexgij].in.push_back( index );
	  edges++;
	}
      }
    }
  }
  
  for( int i = 0; i < GraphOut.size(); i++ ) {
    if( (!GraphOut[i].alive) || GraphOut[i].out.size()==0 ) continue;
    int target = GraphOut[i].out.front(); // We plan to delete the edge from i to target
    GraphOut[i].alive = 0;
    for( std::list<int>::iterator j = GraphIn[target].in.begin(); j != GraphIn[target].in.end(); j++ ) {
      if( !GraphOut[*j].alive ) continue;
      for( std::list<int>::iterator k = GraphOut[i].out.begin(); k != GraphOut[i].out.end(); k++ ) {
	std::list<int>::iterator search = find( GraphOut[*j].out.begin(), GraphOut[*j].out.end(), *k);
	if( search != GraphOut[*j].out.end() ) {
	  GraphOut[*j].out.erase(search);
	  if( *k != target )  {
	    GraphIn[*k].in.remove(*j);
	  }
	} else {
	  GraphOut[*j].out.push_back(*k);
	  GraphIn[*k].in.push_back(*j);
	}
      }
    }
    GraphIn[target].in.clear();
  }
  int kernelDimension = 0;
  for( int i = 0; i < GraphOut.size(); i++ ) {
    if( GraphOut[i].alive )
      kernelDimension++;
  }
  return kernelDimension;
}


// Inverse mapping, from integers < gridsize! to permutations of size n
// Writes the permutation corresponding to n into the array P.
void getPerm( generator n, int *P, int gridsize) {
  int taken[gridsize];
  int offset;
  for( int i = 0; i < gridsize; i++ )
    taken[i] = 0;
  for( int i = 0; i < gridsize; i++ ) {
    offset = n / factorial(gridsize-1-i);
    n -= offset*factorial(gridsize-1-i);
    for( int j = 0; j <= offset; j++ )
      if( taken[j] )
	offset++;
    P[i] = offset;
    taken[P[i]] = 1;
  }
}


int rectIndex(int i, int j, int k, int l, int m, int gridsize) {
  return (m+4*(l+gridsize*(k+gridsize*(j+gridsize*(i)))));
}

void initRectangles(int gridsize, int *white, int *black) {
  Rectangles = (bool*) malloc(gridsize*gridsize*gridsize*gridsize*4*sizeof(bool));
  for(int xll=0; xll < gridsize; xll++) {
    for(int xur=xll+1; xur < gridsize; xur++) {
      for(int yll=0; yll < gridsize; yll++) {
	for(int yur=yll+1; yur < gridsize; yur++) {
	  Rectangles[rectIndex(xll,yll,xur,yur,0,gridsize)] = RectDotFree(xll,yll,xur,yur,0,white,black,gridsize);
	  Rectangles[rectIndex(xll,yll,xur,yur,1,gridsize)] = RectDotFree(xll,yll,xur,yur,1,white,black,gridsize);
	  Rectangles[rectIndex(xll,yll,xur,yur,2,gridsize)] = RectDotFree(xll,yll,xur,yur,2,white,black,gridsize);
	  Rectangles[rectIndex(xll,yll,xur,yur,3,gridsize)] = RectDotFree(xll,yll,xur,yur,3,white,black,gridsize);
	}
      }
    }
  }
  
}

generator getIndex( int *P, int gridsize ) {
  generator index = 0;
  for( int i = gridsize-2; i >= 0; i-- ) {
    int r = P[i];
    int m = 0;
    for( int j = 0; j < i; j++ )
      if( P[j] < r )
	m++;
    index += factorial(gridsize-1-i)*(r-m);
  }
  return index;
}

int Find(std::vector<generator> & V, generator x) {
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

generator factorial( int i ) {
  static std::vector<generator> results;
  if( results.size() == 0 )
    results.push_back(1);
  for( int j = results.size(); j <= i; j++ )
    results.push_back(j*results[j-1]);
  return results[i];
}

bool RectDotFree(int xll, int yll, int xur, int yur, int which, int *white, int *black, int gridsize) {
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

