using namespace std;

void increment( int &turn, int fP, int nP ) {
  turn++;
  if( turn == fP+nP )
    turn = fP;
}

// At the moment, everything below this line is specialized to Z2


// add column to c in the graph if the first one of column also appears in c
// it would be very inefficient to use this for full matrix reduction, but it should be fine if used only on a single column
void resolveCols(vector<Generator> &GraphOut, vector<Generator> &GraphIn, int* column, int size, int c, int firstCol) {
  int pivot = column[0];
  list<int> &col2 = GraphOut[c].ones;
  list<int>::iterator location = find(col2.begin(), col2.end(), pivot);
  if( location == col2.end() )
    return;
  for( int k = 0; k < size; k++ ) {
    list<int>::iterator search = find( GraphOut[c].ones.begin(), GraphOut[c].ones.end(), column[k] );
    if( search != GraphOut[c].ones.end() ) {
      GraphOut[c].ones.erase(search);
      if( column[k] != pivot ) 
	GraphIn[column[k]].ones.remove(c+firstCol);
    } else {
      GraphOut[c].ones.push_back(column[k]);
      if( column[k] != pivot )
	GraphIn[column[k]].ones.push_back(c+firstCol);
    }
  }
}

void displayMatrix(vector<Generator> &GraphOut, vector<Generator> &GraphIn, int rows, int cols, int rank, int firstCol) {
  for( int r = 0; r < rows; r++ ) {
    printf("(%d): ", rank);
    for( int c = 0; c < cols; c++ ) {
      bool foundOut = false;
      bool foundIn = false;
      for( list<int>::iterator it = GraphOut[c].ones.begin(); it != GraphOut[c].ones.end(); it++ )
	if( *it == r )
	  foundOut = true;
      for( list<int>::iterator it = GraphIn[r].ones.begin(); it != GraphIn[r].ones.end(); it++ ) {
	if( *it == c+firstCol )
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
