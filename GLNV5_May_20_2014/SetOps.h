//
//  SetOps.h
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//

#ifndef __gln__SetOps__
#define __gln__SetOps__

#include <iostream>
#include <vector>

using namespace std;

vector<int> unionSet(const vector<int> & firstSet, const vector<int> & secondSet);
vector<int> complementSet(const vector<int> & currentSet, const vector<int> & completeSet);
vector<vector<int> > powerSet(const vector<int> & parents);
vector<int> intersectionSet(const vector<int> & firstSet, const vector<int> & secondSet);

#endif /* defined(__gln__SetOps__) */
