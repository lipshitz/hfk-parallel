all:	gradings fullv1 boundaries boundaries-bc boundaries-bc-set

gradings:	gradings.cpp
	mpicxx -o gradings gradings.cpp

fullv1:	fullv1.cpp matrix-z2z3.h reduction-functions.h
	mpicxx -o fullv1 fullv1.cpp

boundaries:	boundaries.cpp matrix-z2z3.h reduction-functions.h
	mpicxx -o boundaries boundaries.cpp

boundaries-bc:	boundaries-bc.cpp matrix-z2z3.h reduction-functions.h
	mpicxx -o boundaries-bc boundaries-bc.cpp

boundaries-openmp:	boundaries-openmp.cpp matrix-z2z3.h reduction-functions.h
	mpicxx -o boundaries-openmp boundaries-openmp.cpp -fopenmp

boundaries-bc-set:	boundaries-bc-set.cpp matrix-z2z3-set.h reduction-functions-set.h
	mpicxx -o boundaries-bc-set boundaries-bc-set.cpp