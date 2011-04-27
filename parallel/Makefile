all:	gradings fullv1

gradings:	gradings.cpp
	mpicxx -o gradings gradings.cpp

fullv1:	fullv1.cpp matrix-z2z3.h reduction-functions.h
	mpicxx -o fullv1 fullv1.cpp