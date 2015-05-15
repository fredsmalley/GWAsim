// JobEnvironment.h
//
// Curtis Luce, Joe Song
// Created: July 2007.
// Last modified: Nov 18, 2011. Added virtual function getMaxNumDiffParents()

#pragma once

#include <vector>
using std::vector;

#include "EnumMethod.h"
#include "GLNConfig.h"

class JobEnvironment : public EnumMethod {

public:
	
    virtual	void finalize() = 0;
    
	virtual size_t getGLNSize() const = 0;
	virtual char getNodeType(int i) const = 0;
    virtual size_t getMaxNumDiffParents() const = 0;
    
    virtual bool withinRange(const vector<vector<int> > & parentSets) const 
    {
        for ( size_t k=0; k < parentSets.size(); k++ ) {
            if((int) parentSets[k].size() < m_gc.m_min_parents_allowed) {
                return false;
            } else if((int) parentSets[k].size() > m_gc.m_max_parents_allowed) {
                return false;
            }
        }
        return true;
    }

public:

	GLNConfig	m_gc; 
	
	//added by YangZhang 2/23/2009 to record all results with p-value smaller than threshold
	vector< vector<string> >   m_recordresult;
	//added by CMHJ 4/18/2009 to record the counts found in contingency tables of the m_recordresults
	vector< vector<string> >   m_recordcounts;
};

