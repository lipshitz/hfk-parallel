#include "gradings.h"

Gradings::Gradings(const HalfInteger *bm) { boundaryMap = bm; }

void Gradings::addGeneratorToGrading( const generator gen, const HalfInteger *grad ) {
  const HalfInteger *gradInternal;
  gradingIterator it = gradings.find(grad);
  if( it != gradings.end() )
    gradInternal = *it;
  else {
    HalfInteger *gtemp = new HalfInteger[NUM_GRADINGS];
    for( int i = 0; i < NUM_GRADINGS; i++ )
      gtemp[i] = grad[i];
    gradInternal = gtemp;
    gradings.insert(gradInternal);
    generatorsOfGrading[gradInternal] = new std::vector<generator>();
  }
  generatorsOfGrading[gradInternal]->push_back(gen);
}

std::vector<generator>* Gradings::getGenerators( const HalfInteger* grad ) {
  if( generatorsOfGrading.find(grad) == generatorsOfGrading.end() )
    return &emptyGenerators;
  return generatorsOfGrading[grad];
}

std::vector<generator>* Gradings::getBoundaryGenerators( const HalfInteger* grad ) {
  HalfInteger *bGrad = getBoundaryGrading(grad);
  std::vector<generator> *ret = getGenerators( bGrad );
  delete [] bGrad;
  return ret;
}

gradingIterator Gradings::getGradingIterator() {
  return gradings.begin();
}

gradingIterator Gradings::getGradingIteratorEnd() {
  return gradings.end();
}

void Gradings::setKernelDimension(const HalfInteger* grad, int dim) {
  assert( gradings.find(grad) != gradings.end() );
  kernelDimensions[grad] = dim;
}

void Gradings::setImageDimension(const HalfInteger* grad, int dim) {
  assert( gradings.find(grad) != gradings.end() );
  imageDimensions[grad] = dim;
}

int Gradings::getHomologyDimension( const HalfInteger* grad ) {
  if( homologyDimensions.find(grad) == homologyDimensions.end() )
    return 0;
  return homologyDimensions[grad];
}

void Gradings::calculateHomology() {
  for( gradingIterator it = getGradingIterator(); it != getGradingIteratorEnd(); it++ ) {
    int kd = getKernelDimension(*it);
    HalfInteger *bGrad = getInverseBoundaryGrading(*it);
    int id = getImageDimension(bGrad);
    delete [] bGrad;
    homologyDimensions[*it] = kd - id;
  }
}

Gradings::~Gradings() {
  for( gradingIterator it = gradings.begin(); it != gradings.end(); it++ ) {
    delete generatorsOfGrading[*it];
    delete [] (*it);
  }
}

int Gradings::getKernelDimension( const HalfInteger* grad ) {
  if( kernelDimensions.find(grad) == kernelDimensions.end() )
    return 0;
  return kernelDimensions[grad];
}

int Gradings::getImageDimension( const HalfInteger* grad ) {
  if( imageDimensions.find(grad) == imageDimensions.end() )
    return 0;
  return imageDimensions[grad];
}

HalfInteger* Gradings::getBoundaryGrading(const HalfInteger* grad ) const {
  HalfInteger *bGrad = new HalfInteger[NUM_GRADINGS];
  for( int i = 0; i < NUM_GRADINGS; i++ )
    bGrad[i] = grad[i]+boundaryMap[i];
  return bGrad;
}

HalfInteger* Gradings::getInverseBoundaryGrading(const HalfInteger* grad ) const {
  HalfInteger *bGrad = new HalfInteger[NUM_GRADINGS];
  for( int i = 0; i < NUM_GRADINGS; i++ )
    bGrad[i] = grad[i]-boundaryMap[i];
  return bGrad;
}

