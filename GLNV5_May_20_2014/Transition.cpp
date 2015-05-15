// Transition.cpp
//
// Joe Song
// Created: November 30, 2008

#include <cstdlib>
#include "Transition.h"

void Transition::initialize(const int child, const int base,
							const vector< int > & parents, const vector<int> & parentBases, 
							const vector< int > & delays,  const string & nameChild, 
                            const vector<string> & nameParents)
{
	m_child = child;
    m_nameChild = nameChild;
    
	m_base = base;

    /*
	m_parents = parents;
    m_nameParents = nameParents;
    */
    
    setParents(parents, nameParents);
    
	m_parentBases = parentBases;

	m_delays = delays;

	allocate();
}

size_t Transition::calculateSize() const
{
	size_t nRows = 0;

	if(m_parentBases.size() > 0) {
		nRows=1;
		for(size_t i=0; i<m_parentBases.size(); i++) {
			nRows *= m_parentBases[i];
		}
	}

	return nRows;
}

void Transition::setParents(const vector<int> & parents, const vector<string> & names)
{ 
    m_parents = parents; 
    
    if(parents.size() > 0 && names == vector<string>(0) ) {
        m_nameParents = vector<string>(parents.size(), "");
    } else {
        m_nameParents = names;
    }
}
