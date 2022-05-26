#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include <list>

using std::copy;
using std::cout;
using std::endl;
using std::vector;

using namespace std;

vector< vector<int> > comb;

void combinations_r_recursive(vector<int> elems, int req_len,
			    vector<int> pos, int depth,
			    int margin)
{
	vector<int> pair;
	// Have we selected the number of required elements?
	if (depth >= req_len) {
		for (int ii = 0; ii < pos.size(); ++ii)
			// cout << elems[pos[ii]];
			pair.push_back(elems[pos[ii]]);
			comb.push_back(pair);
		// cout << endl;
		return;
	}

	// Try to select new elements to the right of the last selected one.
	for (int ii = margin; ii < elems.size(); ++ii) {
		pos[depth] = ii;
		// pair.push_back(pos[depth]);
		combinations_r_recursive(elems, req_len, pos, depth + 1, ii);
		
	}
	return;

	
}

void combinations_r(vector<int> elems, int req_len)
{
	// assert(req_len > 0 && req_len <= elems.size());
	vector<int> positions(req_len, 0);
	combinations_r_recursive(elems, req_len, positions, 0, 0);
}


void printMultiVec(const vector< vector<int> > vec) {
	cout << "vec.size()" << vec.size() << endl;
	
    for (int i = 0; i < vec.size(); i++) {
		// cout << "vec[i].size()  " << vec[i].size() << endl;
        for (int j = 0; j < vec[i].size(); j++)
            cout << vec[i][j] << " ";
        cout << endl;
    }
}


int main()
{	
	vector <int> swag = {1,2,3,4};
	int num_elements = swag.size();

	vector<int> elements(num_elements);
	list<int> SwagList;
	copy( swag.begin(), swag.end(), std::back_inserter( SwagList) );

	int elements_str[num_elements + 1] = {1,2,3,4};

	copy(swag.begin(), swag.begin() + num_elements, elements.begin());

	combinations_r(elements, 2);

	printMultiVec(comb);


	return 0;
}
