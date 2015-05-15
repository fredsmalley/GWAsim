// EnvGLNCmp.h -- The job environment for transition and enumeration based GLN comparison
//                Multiple trajectories can be compared
//
// Joe Song
// Created: November 30, 2008
// Updated: October 1, 2011. Joe Song.  Changed "parentSelect()" to "compareParentSets"
//   Added "cmpParentSetsbyType"
// Last modified: October 31, 2011. Joe Song. Added several inline functions to simplify the 
//   Class implementation file.

#pragma once

#include <algorithm>
#include <set>
#include <cassert>

using namespace std;

#include "JobEnvironment.h"
#include "EnumMethod.h"
#include "TrajectoryCollection.h"
#include "TransitionTable.h"
#include "GLNConfig.h"
#include "GLN.h"

#include "StatDistributions.h"
#include "PathwayStats.h"

// 0: null interation, 1: relatively differential, 2: absolutely differential, 3: homogeneous
enum CmpIntxType {NULL_INTX=0, REL_DIFF=1, ABS_DIFF=2, CONSERVED=3};

enum WorkingZone {UNCHANGED=0, CHANGED=1};

class EnvGLNCmp :	public JobEnvironment
{
public:
    
	void initialize(const GLNConfig & gc);
	virtual void finalize();
    
	/* Core enumeration and processing functions */
	virtual bool processOneEnumeration(const int child, const vector< int > & parents, 
                                       const vector< int > & offset);
    
	virtual bool processOneEnumeration(const int child, const vector<vector< int > > & parents, 
                                       const vector<vector< int > > & delays);
    
	/* get/set functions */
    virtual size_t getMaxNumDiffParents() const;
    
	virtual size_t getGLNSize() const { return m_glnDiff.size(); }
	virtual char getNodeType(int i) const { return m_glnDiff.getNodeType(i); }
	virtual const string & getNodeName(int i) const { return m_glnDiff.getNodeName(i); }
    
	// const GeneralizedLogicalNetwork & getGLN() const { return m_gln; }
    
	const TransitionTable & getPooledTransTable(int i) const { return m_pooledTransTables[i]; }
	void setPooledTransTable(int i, const TransitionTable & tt) { m_pooledTransTables[i] = tt; }
    
	const TransitionTable & getDiffTransTable(int i) const { return m_diffTransTables[i]; }
	void setDiffTransTable(int i, const TransitionTable & tt) { m_diffTransTables[i] = tt; }
    
    const vector<TransitionTable> & getTransTables(int i) const { return m_TransTables[i]; }
	void setTransTables(int i, const vector<TransitionTable> & tts) { m_TransTables[i] = tts; }
    
	double getChisq_t(int i) const{ return  m_chisqs_t[i]; }
	void setChisq_t(int i, double Chi_t){ m_chisqs_t[i] = Chi_t; }
    
	const int getDf_t(int i) const{ return  m_dfs_t[i]; }
	void setDf_t(int i, int Df_t){  m_dfs_t[i] = Df_t; }
    
	const double getpChisq_t(int i) const{ return  m_pchisqs_t[i]; }
	void setpChisq_t(int i, double pChi_t){ m_pchisqs_t[i] = pChi_t; }
    
	const int getType(int i) const{ return  m_types[i]; }
	void setType(int i, int type){  m_types[i] = type; }
    
	double getChisq_z(int i) const{ return  m_chisqs_z[i]; }
	void setChisq_z(int i, double Chi_z){ m_chisqs_z[i] = Chi_z; }
    
	const int getDf_z(int i) const{ return  m_dfs_z[i]; }
	void setDf_z(int i, int Df_z){  m_dfs_z[i] = Df_z; }
    
	const double getpChisq_z(int i) const{ return  m_pchisqs_z[i]; }
	void setpChisq_z(int i, double pChi_z){ m_pchisqs_z[i] = pChi_z; }
    
	double getChildChisq_z(int i) const{ return  m_childchisqs_z[i]; }
	void setChildChisq_z(int i, double Chi_z){ m_childchisqs_z[i] = Chi_z; }
    
	const int getChildDf_z(int i) const{ return  m_childdfs_z[i]; }
	void setChildDf_z(int i, int Df_z){  m_childdfs_z[i] = Df_z; }
    
	const double getChildpChisq_z(int i) const{ return  m_childpchisqs_z[i]; }
	void setChildpChisq_z(int i, double pChi_z){ m_childpchisqs_z[i] = pChi_z; }
       
	const GeneralizedLogicalNetwork getDiffGLN() const {return m_glnDiff;}
	const GeneralizedLogicalNetwork getPooledGLN() const {return m_glnPooled;}
    
    void setMarkovOrder (int min, int max) {
        
        m_gc.m_min_Markov_order = min;
        m_gc.m_max_Markov_order = max;
        
        m_glnPooled.setMaxMarkovOrder(max);  // m_gc.m_max_Markov_order
        m_glnDiff.setMaxMarkovOrder(max);
    }
    
    const EnvGLNCmp & modify(int min_Markov_order, int max_Markov_order,
                             const vector<TrajectoryCollection> & trajCols);
    
	// added by YangZhang 2/7/2009 to denote whether or not the networks are different
	int checkdiff(const vector<TransitionTable> & diffTransTables, 
                  const vector<TransitionTable> & pooledTransTables, double alpha) const;
    
	// added by YangZhang 5/9/2009
	// Returns a Boolean value indicating whether this parent & delay enumeration 
    //   has improved and if so update the parent and delay sets
	bool decideInteractionType(int child, const vector< int > & parents, 
                               const vector< int > & delays);
	
	bool decideInteractionType(int child, const vector<vector< int > > & parents, 
                               const vector<vector< int > > & delays);

    void collectStats(int child, const vector<vector< int > > & parents,
                      const vector<vector< int > > & delays,
                      vector<TransitionTable> & tts, TransitionTable & ttTot,
                      TransitionTable & ttHom, TransitionTable & ttHet) const;
    
    void collectStatsPooled(int child, const vector<vector< int > > & parents,
                      const vector<vector< int > > & delays, 
                      vector<TransitionTable> & tts, TransitionTable & ttTot, 
                      TransitionTable & ttHom, TransitionTable & ttHet) const;

	void collectStatsCPFunX2(int child, const vector<vector< int > > & parents,
		const vector<vector< int > > & delays,
		vector<TransitionTable> & tts, TransitionTable & ttTot,
		TransitionTable & ttHom, TransitionTable & ttHet) const;

    void collectStatsIndep(int child, const vector<vector< int > > & parents,
                      const vector<vector< int > > & delays,
                      vector<TransitionTable> & tts, TransitionTable & ttTot,
                      TransitionTable & ttHom, TransitionTable & ttHet) const;

    void collectDiffCorPooled(int child, const vector<vector< int > > & parents, const vector<vector< int > > & delays, vector<TransitionTable> & tts, TransitionTable & ttTot, TransitionTable & ttHom, TransitionTable & ttHet) const;
    
    bool select(int child, const vector<TransitionTable> & tts, const TransitionTable & ttTot, const TransitionTable & ttHom, const TransitionTable & ttHet);
    
    bool analyze(int child, const vector<vector< int > > & parents, 
                 const vector<vector< int > > & delays);
    
	bool compareParentSets(int child, int type, const vector <TransitionTable> & tts,
                           const TransitionTable & ttTot,
                           const TransitionTable & ttHom, 
                           const TransitionTable & ttHet) const;
    
    bool cmpParents_BY_HOM(int child, const vector <TransitionTable> & tts,
                           const TransitionTable & ttTot,
                           const TransitionTable & ttHom, 
                           const TransitionTable & ttHet) const;
    
    bool cmpParents_BY_EACH_COND(int child, const vector <TransitionTable> & tts,
                                 const TransitionTable & ttTot,
                                 const TransitionTable & ttHom, 
                                 const TransitionTable & ttHet) const;
    
    bool cmpParents_BY_T_IMP_PER_DF(int child, const vector <TransitionTable> & tts,
                                    const TransitionTable & ttTot,
                                    const TransitionTable & ttHom, 
                                    const TransitionTable & ttHet) const;
    
    bool cmpParents_BY_T_C_IMP(int child, const vector <TransitionTable> & tts,
                               const TransitionTable & ttTot,
                               const TransitionTable & ttHom, 
                               const TransitionTable & ttHet) const;
    
    bool cmpParents_BY_INTX_TYPE(int child, int type, 
                                 const vector <TransitionTable> & tts,
                                 const TransitionTable & ttTot,
                                 const TransitionTable & ttHom, 
                                 const TransitionTable & ttHet) const;
	
    bool analyze_BY_BEST_HOM(int child, const vector<vector< int > > & parents, 
                             const vector<vector< int > > & delays);
    
    bool select_BY_BEST_HOM(int child, 
                            const vector<TransitionTable> & tts,
                            const TransitionTable & ttTot, 
                            const TransitionTable & ttHom);
    
	void updateEnvironment(int child,int type, const vector <TransitionTable> & tts,
                           const TransitionTable & ttTot,
                           const TransitionTable & ttHom, 
                           const TransitionTable & ttHet, 
                           double p_z, int df_z, double chi_z);
    
	/* Child working zone statistics */
	void compareChildWorkingZones();
    
	/* Pathway statistics */
	void comparePathway();
    void gatherPathwayStats();
    void savePathwayStats() const;
    
	/* I/O */
	void print(const GeneralizedLogicalNetwork & gln, const vector<TransitionTable> & tts) const;
	void print(int i) const;
	void print(void) const;
	void printk2(void) const;
    
	void saveDotInteraction(const string & filename);
	void saveDotFile(const string & filename);
	void saveDotFilek2(const string & filename); //added by Yang Zhang 10.25.2012, support different parents
    void saveDOTNodes(ofstream & out, const vector<bool> & used,
                                    const vector<bool> & isParent) const;

    void saveDOTCondition(ofstream & out, size_t cond, const string colors[],
                          bool showAllNodes, bool mergeNodes,
                          vector<bool> & used, vector<bool> & isParent) const; // May 6, 2014. Added by MS
    void saveDOTFilebyCondition(const string & filename) const;  // May 6, 2014. Added by MS

	void save(const string & file) const;

	void saveTopology(const string & fileToplogy, int typeIntx, int typeParents, int typeChild) const;
	void saveTopology(const string & fileToplogy) const;
    
	void saveAllInteractionToDot();
    
    // added by CMHJ (4/17/09): Records results to Header file:	
    void recordResults(int child, const TransitionTable &tt, string description);  
   
	// added by Yang Zhang (4/15/14): Add homogeneity and total strength
	void recordResults(int child, const TransitionTable &ttHet, const TransitionTable &ttHom, const TransitionTable &ttTot, string description);

    const vector<TrajectoryCollection> & getTrajCols() const { return m_trajCols; }
    void setTrajCols(const vector<TrajectoryCollection> & trajCols)
    { m_trajCols = trajCols; }

    const MultiplePathwayStats & getPathwayStats() const { return m_pathwayStats; }
    void setPathwayStats(const MultiplePathwayStats & stats) { m_pathwayStats=stats; }
    
private:
    
	GeneralizedLogicalNetwork	m_glnDiff, m_glnPooled;
    
	vector<TrajectoryCollection> m_trajCols;
    
	vector<TransitionTable>		m_pooledTransTables;
	vector<TransitionTable>		m_diffTransTables;
    
	//following members added for different parent selecting strategy
	vector< vector<TransitionTable> > m_TransTables;
    
	vector<double> m_chisqs_t;
	vector<int> m_dfs_t;
	vector<double> m_pchisqs_t;
    
	//these correspond to parent working zone
	vector<double> m_chisqs_z;
	vector<int> m_dfs_z;
	vector<double> m_pchisqs_z;
    
	//these correspond to child working zone
	vector<double> m_childchisqs_z;
	vector<int> m_childdfs_z;
	vector<double> m_childpchisqs_z;
    
	vector<int> m_types;
    
	vector<enum WorkingZone> m_workingZones;
	vector<enum WorkingZone> m_childWorkingZones;
    
	// double m_p_val; 		//the result after calling calculate_pvalue
	// double m_chisq;		// Joe Song 08/03/06.  The chisq value after calling calculate_pvalue
    
    MultiplePathwayStats  m_pathwayStats;
};

inline bool isSelfCycleZeroDelay(int child, const vector<vector< int > > & parents, 
                                 const vector<vector< int > > & delays)
{
    size_t m;
    
	bool flag = false; //Check if child itself at delay=0 (current time) has been picked as its own parent.
    
	for(size_t i=0; i<parents.size(); i++)
	{
		for(m = 0; m<parents[i].size(); m++) {
			// Check if child itself at delay=0 (current time) has been picked as its own parent.
			if(child == parents[i][m]-1 && delays[i][m] == 0) {
                break;
			}
		}
		if(m < parents[i].size())
		{
			flag = true;
			break;
		}
	}
    return flag;
}

inline CmpIntxType labelIntxType(double p_t, double p_c, double p_d, double alpha)
{
    CmpIntxType type = NULL_INTX;  	//int type = 0;
    
    if(p_d <= alpha)
    {
        if(p_t > alpha) {
            type = REL_DIFF; // 1; //relative differential interaction
        } else {
            type = ABS_DIFF; //2; //absolute differential interaction
        }
        
    } else if(p_c <= alpha) {
        
        type = CONSERVED; // 3;
        
    } else {
        type = NULL_INTX; // 0;
    }
    
    return type;
}
