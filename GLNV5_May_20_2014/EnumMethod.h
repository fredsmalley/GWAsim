// EnumMethod.h
//
// Joe Song
// Created: November 30, 2008
// Modified: October 1, 2011. Joe Song.  Added this header.
// Last modified: Novermber 18, 2011. Added enumerateParentsSubsets()

#pragma once

#include <vector>
using std::vector;

#include "GLNConfig.h"

class EnumMethod
{
public:

    static int next_N_Elements_To_K_Sets(vector<unsigned long> & a, unsigned K);
    static int nextAssignmentWithRepeat(vector<unsigned> & a, unsigned K);
    static int nextAssignmentWithRepeat(vector<int> & a, const vector<int> & K);
    static int nextAssignmentNoRepeat(vector<int> & a, int K);
    
	// Add parameter parentList by Yang Zhang on 7.25.2010 in order to only enumerate parents on candidate file list
	bool enumerateParents(const int child, const vector<int> & parentsBegin, const vector<int> & parentsEnd,
                          const int nNodes, const GLNConfig & gc, const vector<int>  & parentList);

	// Add parameter parentList by Yang Zhang on 7.25.2010 in order to only enumerate parents on candidate file list
	bool enumerateParents2(const int child, const vector<int> & parentsBegin, const vector<int> & parentsEnd, 
                           const int nNodes, const GLNConfig & gc, const vector<int> & parentList);

    bool enumerateParentSubsets(const int child, const vector<int> & parentsBegin, const vector<int> & parentsEnd,	
                                const int nNodes, const GLNConfig & gc, const vector<int> & candidateParents, const vector<int> & mustIncludeParents = vector<int>(0) );

	bool enumerateDelays(const int child, const vector<int> & parents, const GLNConfig & gc);

	bool enumerateDelays2(const int child, const vector<vector<int> > & parents, const GLNConfig & gc);

	//added by Yang Zhang 2.24.2011
	bool generateDelays (vector<int> & delays, const int min_Markov_order, const int max_Markov_order, bool first); // Joe Song 08/07/2008

	virtual bool processOneEnumeration(const int child, const vector< int > & parents, 
		const vector< int > & delays) = 0;

	virtual bool processOneEnumeration(const int child, const vector<vector< int > > & parents, 
		const vector<vector< int > > & delays) = 0;
    
    virtual bool withinRange(const vector<vector<int> > & parentSets) const = 0;

};


class Parent
{
public:
    Parent(): m_id(0), m_delay(0) { }
    Parent(int id, int delay): m_id(id), m_delay(delay) { } 
    
    int getId() const { return m_id; }
    int getDelay() const { return m_delay; }
    
    //operator< defined in "EnvGLNCmp.cpp";
    friend bool operator< (const Parent & p1, const Parent & p2);
    
private:
    int m_id;
    int m_delay;
};

void findCommonParentsWithSameDelays(const vector<vector< int > > & parents, 
                                     const vector<vector< int > > & delays,
                                     vector<int> & commonParents,
                                     vector<int> & commonDelays);

void findCommonParents(const vector<vector< int > > & parents, 
                       const vector<vector< int > > & delays,
                       vector<int> & commonParents,
                       vector<int> & commonDelays);

void mergeParentsWithSameDelays(const vector<vector< int > > & parents,
                                const vector<vector< int > > & delays,
                                vector<int> & parentSuperset,
                                vector<int> & delaySuperset);