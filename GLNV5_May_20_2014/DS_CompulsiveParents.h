#pragma once

#include <vector>
using std::vector;

#include "DS_AdjacencyList.h"

class DS_CompulsiveParents :
	public DS_AdjacencyList
{
public:
	DS_CompulsiveParents(void);
	~DS_CompulsiveParents(void);

public:
	//vector<int> m_polarityRelations;

    //Oct 28, 2012
    //In order to support non-integer relationship between parent and child
    vector<double> m_polarityRelations;
};

