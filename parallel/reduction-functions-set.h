using namespace std;

void increment( int &turn, int fP, int nP ) {
  turn++;
  if( turn == fP+nP )
    turn = fP;
}

// At the moment, everything below this line is specialized to Z2


// add column to c in the graph if the first one of column also appears in c
// it would be very inefficient to use this for full matrix reduction, but it should be fine if used only on a single column
void resolveCols(vector<Generator> &GraphOut, vector<Generator> &GraphIn, int* column, int size, int c) {
  int pivot = column[0];
  set<int> &col2 = GraphOut[c].ones;
  set<int>::iterator location = col2.find(pivot);
  if( location == col2.end() )
    return;
  for( int k = 0; k < size; k++ ) {
    pair<set<int>::iterator,set<int>::iterator> search = GraphOut[c].ones.equal_range(column[k]);
    if( search.first != search.second ) {
      GraphOut[c].ones.erase(search.first);
      if( column[k] != pivot ) 
	GraphIn[column[k]].ones.erase(c);
    } else {
      GraphOut[c].ones.insert(search.first, column[k]);
      if( column[k] != pivot )
	GraphIn[column[k]].ones.insert(c);
    }
  }
}

// for resolving with a column in GraphOut
void resolveColsInternal(vector<Generator> &GraphOut, vector<Generator> &GraphIn, int column, int c) {
  set<int>::iterator k = GraphOut[column].ones.begin();
  if( k == GraphOut[column].ones.end() )
    return; // column is empty
  int pivot = *k;
  set<int> &col2 = GraphOut[c].ones;
  set<int>::iterator location = col2.find(pivot);
  if( location == col2.end() )
    return;
  for( ; k != GraphOut[column].ones.end(); k++ ) {
    pair<set<int>::iterator,set<int>::iterator> search = GraphOut[c].ones.equal_range(*k);
    if( search.first != search.second ) {
      GraphOut[c].ones.erase(search.first);
      //if( *k != pivot )  // I don't think we want these conditions
	GraphIn[*k].ones.erase(c);
    } else {
      GraphOut[c].ones.insert(search.first, *k);
      //if( *k != pivot )
	GraphIn[*k].ones.insert(c);
    }
  }
}

// now that we are using sets, we should have some multiple-column operations

void displayMatrix(vector<Generator> &GraphOut, vector<Generator> &GraphIn, int rows, int cols, int rank) {
  for( int r = 0; r < rows; r++ ) {
    printf("(%d): ", rank);
    for( int c = 0; c < cols; c++ ) {
      bool foundOut = false;
      bool foundIn = false;
      for( set<int>::iterator it = GraphOut[c].ones.begin(); it != GraphOut[c].ones.end(); it++ )
	if( *it == r )
	  foundOut = true;
      for( set<int>::iterator it = GraphIn[r].ones.begin(); it != GraphIn[r].ones.end(); it++ ) {
	if( *it == c )
	  foundIn = true;
      }
      if( foundOut && foundIn ) 
	printf("1 ");
      else if( !foundOut && !foundIn )
	printf("0 ");
      else if( foundOut )
	printf("+ ");
      else
	printf("- ");
    }
    printf("\n");
  }
}
