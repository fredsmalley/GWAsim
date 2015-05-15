// EnvGLNRec.h -- The job environment for transition and enumeration based GLN reconstruction
//
// Joe Song
// Created: November 28, 2008
// Modified:
//   July 13, 2014. MS. Added saveTopology() member function

#pragma once

#include "JobEnvironment.h"
#include "EnumMethod.h"
#include "TrajectoryCollection.h"
#include "TransitionTable.h"
#include "GLNConfig.h"
#include "GLN.h"

class EnvGLNRec :	public JobEnvironment
{
public:

	void initialize(const GLNConfig & gc);
    
	virtual void finalize();
	
    /* Core enumeration and processing functions */
	virtual bool processOneEnumeration(const int child, const vector< int > & parents, const vector< int > & offset);

	virtual bool processOneEnumeration(const int child, const vector<vector< int > > & parents, 
                                       const vector<vector< int > > & offset);

    bool compareParentSets(const TransitionTable & tt1, const TransitionTable & tt2);
    
	void calculateOverallPval(int pValMode);

	/* get/set functions */
    
    size_t getMaxNumDiffParents() const { 
        // To make sure that each node can have any combination of m_max_parents_allowed 
        return m_gc.m_max_parents_allowed; }

	const GeneralizedLogicalNetwork & getGLN() const { return m_gln; }
	const TransitionTable & getTransTable(int i) const { return m_transTables[i]; }
	void setTransTable(int i, const TransitionTable & tt) { m_transTables[i] = tt; }

	//added by Yang Zhang 06/12/2009
	const vector<TransitionTable> getTransTable() const {return m_transTables;}

	virtual size_t getGLNSize() const { return m_gln.size(); }

	virtual char getNodeType(int i) const { return m_gln.getNodeType(i); }

	/* I/O */
	void print(void) const;
    void recordResults(int child, const TransitionTable & tt);
    void saveTopology(const string & fileTopology) const;

    //Added by Haizhou Wang, Jun 21, 2012
    TrajectoryCollection getTraCol() const { return m_trajCol; }
    //End of Haizhou's Code

    void calculatePathwayChisqStats() const; // MS Feb 7, 2014. Added.

private:

	GeneralizedLogicalNetwork	m_gln; 
	TrajectoryCollection		m_trajCol;

	vector<TransitionTable>		m_transTables;

	double m_p_val; 		//the result after calling calculate_pvalue
	double m_chisq;		// Joe Song 08/03/06.  The chisq value after calling calculate_pvalue

};
