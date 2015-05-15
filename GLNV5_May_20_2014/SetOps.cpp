//
//  SetOps.cpp
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//
#include <algorithm>
#include <cassert>
#include <cmath>
#include "SetOps.h"

//added by Yang Zhang 2.16.2010
vector<int> unionSet(const vector<int> & firstSet, const vector<int> & secondSet)
{
	vector<int> unionset(firstSet);  // Joe Song 10/8/2011. Use copy constructor instead of copy
    
	//Joe Song 10/8/2011 copy( firstSet.begin(),firstSet.end(),back_insert_iterator<vector<int> >( unionset ) );
    
	for(size_t i=0; i<secondSet.size(); i++)
	{
		if(find(firstSet.begin(),firstSet.end(),secondSet[i])
           == firstSet.end())
		{
			unionset.push_back(secondSet[i]);
		}else{
			continue;
		}
	}
	//copy( secondSet.begin(),secondSet.end(),back_insert_iterator<vector<int> >( firstSet ) );
	sort(unionset.begin(), unionset.end());
	return unionset;
}



vector<int> complementSet(const vector<int> & currentSet, const vector<int> & completeSet)
{
	assert(completeSet.size() >= currentSet.size());
	vector<int> complementset;
	for(size_t i=0; i<completeSet.size(); i++)
	{
		if(find(currentSet.begin(),currentSet.end(),completeSet[i])
           == currentSet.end())
		{
			complementset.push_back(completeSet[i]);
		}else{
			continue;
		}
	}
	return complementset;
}

vector<vector<int> > powerSet(const vector<int> & parents)
{
	vector<vector<int> > powerset((size_t) pow(2.0, (int) parents.size()));
    
	for(size_t i=1; i<powerset.size(); i++) {
		int mask = 1;
		for(size_t j=0; j<parents.size(); j++)
		{
			if((i & mask) != 0)
			{
				powerset[i].push_back(parents[j]);
			}
			mask <<= 1;
		}
	}
    
	return powerset;
}

vector<int> intersectionSet(const vector<int> & firstSet, const vector<int> & secondSet)
{
	vector<int> intersectionset;
	for(size_t i=0; i<firstSet.size(); i++)
	{
		if(find(secondSet.begin(),secondSet.end(),firstSet[i])
           == secondSet.end())
		{
			continue;
		}else{
			intersectionset.push_back(firstSet[i]);
		}
	}
	return intersectionset;
}

void removeDuplicates(std::vector<int>& vec)
{
    
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}