
#ifndef __MATRIX_REP_H
#define __MATRIX_REP_H

#include <linbox/linbox-config.h>
#include <linbox/blackbox/blackbox-interface.h>
//#include <linbox/matrix/sparse.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/vector/vector-traits.h>
#include <linbox/vector/stream.h>
#include <linbox/util/field-axpy.h>
#include <linbox/field/hom.h>
#include <linbox/field/rebind.h>

#include <stdio.h>

extern int gridsize;
extern int Find(std::vector<long long> & V, long long x);

namespace LinBox
{

template <class _Field>
class MatrixRep : public BlackboxInterface
{
public:

	typedef _Field Field;
	typedef typename Field::Element Element;
        typedef MatrixRep<_Field> Self_t;
	typedef MatrixCategories::BlackboxTag MatrixCategory;
    

#ifdef __LINBOX_PARALLEL
	BB_list_list sub_list;
#endif

	//	FileFormatTag Format;

	// Constructor 
	MatrixRep (Field F, std::vector<long long> *cols, std::vector<long long> *rows, bool *rect):_F(F), _VD(F), _MD(F) {
	  M = rows->size();
	  N = cols->size();
	  colIndices = cols;
	  rowIndices = rows;
	  Rectangles = rect;

	}

	~MatrixRep () {
#ifdef __LINBOX_PARALLEL

		BB_list_list::iterator p;

		BB_list::iterator e_p;

		for (p = sub_list. begin(); p != sub_list. end(); ++ p)
			for (e_p = p -> second. begin(); 
			     e_p != p -> second. end(); ++ e_p) {

				Thread::terminate_thread (*e_p);

				delete (*e_p);
			}
#endif
	}

	/** Matrix-vector product
	 * y = A x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	template <class OutVector, class InVector>
	  OutVector &apply (OutVector &y, const InVector &x) const {
	  //#ifdef __LINBOX_PARALLEL
	  //	return BlackboxParallel (y, *this, x, BBBase::Apply);
	  //#else
	  //y.clear();
	  typename InVector::const_iterator I;
	  //printf("Input x: ");
	  //for( I = x.begin(); I != x.end(); ++I )
	  //  printf("%d ", *I);
	  //printf("\n");
	  Element one;
	  _F.init(one, 1);
	  for( I = y.begin(); I != y.end(); ++I )
	    _F.init(*(y.begin()+(I-y.begin())), 0);
	  int g[gridsize];
	  for( I = x.begin(); I != x.end(); ++I ) 
	    if( *I ){
	      //printf("Here %d %d\n", *I, int(I-x.begin()));
	      OutVector w;
	      getPerm((*colIndices)[int(I-x.begin())],g);
	      bool firstrect;
	      bool secondrect;
	      for(int i=0; i<gridsize; i++) {
		for(int j=i+1; j<gridsize; j++) {
		  if(g[i]<g[j]) {
		    firstrect = Rectangles[((((i*gridsize+g[i])*gridsize)+j)*gridsize+g[j])*4+0];
		    for(int k=i+1; k<j && firstrect; k++) {
		      if(g[i] < g[k] && g[k] < g[j]) firstrect=0;
		    }
		    secondrect = Rectangles[((((i*gridsize+g[i])*gridsize)+j)*gridsize+g[j])*4+1];
		    for(int k=0; k<i && secondrect; k++) {
		      if(g[k]<g[i] || g[k] > g[j]) secondrect=0;
		    }
		    for(int k=j+1; k<gridsize && secondrect; k++) {
		      if(g[k]<g[i] || g[k] > g[j]) secondrect=0;
		    }
		  }
		  else {
		    firstrect = Rectangles[((((i*gridsize+g[j])*gridsize)+j)*gridsize+g[i])*4+2];
		    for(int k=i+1; k<j && firstrect; k++) {
		      if(g[k]<g[j] || g[k] > g[i]) firstrect=0;
		    }
		    secondrect = Rectangles[((((i*gridsize+g[j])*gridsize)+j)*gridsize+g[i])*4+3];
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
		    long long Indexgij = getIndex(gij);
		    int indexgij = Find(*rowIndices,Indexgij);
		    
		    if(indexgij==-1) {printf("Error with Alexander grading: %lld -> %lld by switching %d,%d and %d,%d\n", (*colIndices)[*I], Indexgij, i, g[i], j, g[j]); exit(-1); }
		    //			w.push_back(std::pair<size_t, Element> (indexgij, 1));
		    //_F.assign(*(y.begin()+indexgij), 1^*(y.begin()+indexgij));
		    _F.add(*(y.begin()+indexgij), *(y.begin()+indexgij), one);
		    //printf("Flipping %d\n", indexgij);
		    //_F.assign(*(w.begin()+indexgij), 1);
		  }
		}
	      }
	      // We need to somehow add vectors here
	      //y += w;
	    }
	  //printf("output y nonzero at: ");
	  //for( I = y.begin(); I != y.end(); ++I )
	  //  if( *I )
	  //    printf("%d ", int(I-y.begin()));
	  //printf("\n");
	  return y;
	  //#endif
	}
	
	/** Transpose matrix-vector product
	 * y = A^T x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	template <class OutVector, class InVector>
	  OutVector &applyTranspose (OutVector& y, const InVector &x) const { 
	  //#ifdef __LINBOX_PARALLEL
	  //	return BlackboxParallel (y, *this, x, BBBase::ApplyTranspose);
	  //#else
	  typename InVector::const_iterator I;
	  //printf("(transp) Input x: ");
	  //for( I = x.begin(); I != x.end(); ++I )
	  //  printf("%d ", *I);
	  //printf("\n");
	  Element one;
	  _F.init(one, 1);
	  for( I = y.begin(); I != y.end(); ++I )
	    _F.init(*(y.begin()+(I-y.begin())), 0);
	  int g[gridsize];
	  for( I = x.begin(); I != x.end(); ++I ) 
	    if( *I ){
	      //printf("Here %d %d\n", *I, int(I-x.begin()));
	      OutVector w;
	      getPerm((*rowIndices)[int(I-x.begin())],g);
	      bool firstrect;
	      bool secondrect;
	      for(int i=0; i<gridsize; i++) {
		for(int j=i+1; j<gridsize; j++) {
		  if(g[i]<g[j]) {
		    firstrect = Rectangles[((((i*gridsize+g[i])*gridsize)+j)*gridsize+g[j])*4+2];
		    for(int k=i+1; k<j && firstrect; k++) {
		      if(g[i] > g[k] || g[k] > g[j]) firstrect=0;
		    }
		    secondrect = Rectangles[((((i*gridsize+g[i])*gridsize)+j)*gridsize+g[j])*4+3];
		    for(int k=0; k<i && secondrect; k++) {
		      if(g[k] > g[i] && g[k] < g[j]) secondrect=0;
		    }
		    for(int k=j+1; k<gridsize && secondrect; k++) {
		      if(g[k] > g[i] && g[k] < g[j]) secondrect=0;
		    }
		  }
		  else {
		    firstrect = Rectangles[((((i*gridsize+g[j])*gridsize)+j)*gridsize+g[i])*4+0];
		    for(int k=i+1; k<j && firstrect; k++) {
		      if(g[k] > g[j] && g[k] < g[i]) firstrect=0;
		    }
		    secondrect = Rectangles[((((i*gridsize+g[j])*gridsize)+j)*gridsize+g[i])*4+1];
		    for(int k=0; k<i && secondrect; k++) {
		      if(g[k] < g[j] || g[k] > g[i]) secondrect=0;
		    }
		    for(int k=j+1; k<gridsize && secondrect; k++) {
		      if(g[k] < g[j] || g[k] > g[i]) secondrect=0;
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
		    long long Indexgij = getIndex(gij);
		    int indexgij = Find(*colIndices,Indexgij);
		    
		    if(indexgij==-1) {printf("(transp) Error with Alexander grading: %lld -> %lld by switching %d,%d and %d,%d\n", (*rowIndices)[*I], Indexgij, i, g[i], j, g[j]); exit(-1); }
		    //			w.push_back(std::pair<size_t, Element> (indexgij, 1));
		    //_F.assign(*(y.begin()+indexgij), 1^*(y.begin()+indexgij));
		    _F.add(*(y.begin()+indexgij), *(y.begin()+indexgij), one);
		    //printf("Flipping %d\n", indexgij);
		    //printf("Flipping %d\n", indexgij);
		    //_F.assign(*(w.begin()+indexgij), 1);
		  }
		}
	      }
	      // We need to somehow add vectors here
	      //y += w;
	    }
	  //printf("output y nonzero at: ");
	  //for( I = y.begin(); I != y.end(); ++I )
	  //  if( *I )
	  //    printf("%d ", int(I-y.begin()));
	  //printf("\n");
	  //exit(-1);
	  return y;
	  //#endif
	}


	/** Retreive row dimensions of Sparsemat matrix.
	 * @return integer number of rows of SparseMatrix0Base matrix.
	 */
	size_t rowdim () const { return M; }

	/** Retreive column dimensions of Sparsemat matrix.
	 * @return integer number of columns of SparseMatrix0Base matrix.
	 */
	size_t coldim () const { return N; }

	const Field& field () const { return _F;}

    protected:
	size_t M; // Sizes of the matrix
	size_t N;
	std::vector<long long> *colIndices;
	std::vector<long long> *rowIndices;
	bool *Rectangles;
	Field                             _F;      // Field used for all arithmetic
	VectorDomain<Field>                     _VD;     // Vector domain for matrix operations
	MatrixDomain<Field>                     _MD;     // Matrix domain for matrix operations

	//TransposeMatrix<MatrixRep<_Field> > _AT;

};


  /*struct MatrixTraits< MatrixRep<Field> >
{
  typedef MatrixRep<Field> MatrixType;
  typedef MatrixCategories::BlackboxTag MatrixCategory;
  };*/

#endif // __MATRIX_REP_H
}
