// TransitionTable.h
//
// Joe Song
// Created: September 27, 2008

#pragma once

#include <vector>
using std::vector;

#include <iostream>
using std::ofstream;

#include "Transition.h"
#include "TrajectoryCollection.h"
#include "PermStatTable.h"

// History:
//   Added by Haizhou Wang, Feb 13, 2013
//   #define INITIAL_P_VALUE 1.0 // MS 2/23/2013 Changed from 1.01 to 1.0
//   MS 6/1/2013 Changed to the following 
const double INITIAL_P_VALUE = 1.01;

class TransitionTable : public Transition
{
	// friend bool compare(const TransitionTable & tt1, const TransitionTable & tt2);
        
	friend TransitionTable operator - (const TransitionTable & tt1, const TransitionTable & tt2);
	friend TransitionTable & abs(TransitionTable & tt);
    
public:
    
    TransitionTable(void) : m_chisq(0), m_df(0), m_pValue(INITIAL_P_VALUE), m_adjustedpValue(INITIAL_P_VALUE) { }
    
	TransitionTable(const int child, const vector< int > & parents,
                    const vector< int > & delays, const TrajectoryCollection & trajCol);
    
    TransitionTable(const int child, const int base, const vector< int > & parents,
                    const vector< int > & parentBases, const vector< int > & delays,
                    const string & nameChild="",
                    const vector<string> & nameParents=vector<string>(0));
    
	void initialize(const int child, const vector< int > & parents,
                    const vector< int > & delays, const TrajectoryCollection & trajCol);
    
	void allocate(void);
    
	vector<int> getRowSums() const;

	//added by Yang Zhang 2013.1.2
	vector<int> getColSums() const;
	int getTotalSum() const;
    
	const TransitionTable & operator += (const TransitionTable & tt);
    
	int get(int r, int c) const { return m_transitionTable[r][c]; }
	void set(int r, int c, int val) { m_transitionTable[r][c] = val; }
    
	const vector< vector<int> > & getTransitionTable() const { return m_transitionTable; }
    
	void setTransitionTable(vector< vector<int> > transtable) { m_transitionTable = transtable; }
    
	void accumulate(int value, const vector<int> & parentValues);
    
	void setpValue(const double pvalue) { m_pValue = pvalue; }
	void setAdjustedpValue(const double pvalue) { m_adjustedpValue = pvalue; }
    
	void setChisq(const double chisq) { m_chisq = chisq; }
    
	void setDf(const int df) { m_df = df; }
    
	//added by yangzhang 11.26.2008
    
	double getpValue() const { return m_pValue; }
	double getAdjustedpValue() const { return m_adjustedpValue; }
    
	double getChisq() const { return m_chisq; }
    
	int getDf() const { return m_df; }
    
	void fill(const TrajectoryCollection & trajCol);
    
    //added by Tyler Hunt 07.31.2012 //
    // Fill the table based on condition conditionNodeIndex when it is at value level
    void condFill(const TrajectoryCollection & trajCol, int conditionNodeIndex, int conditionDelay, int level);
    
    // added by Tyler Hunt 08/01.2012 //
    // reset the table to contain only zeros
    void reset();
    
	void applyChisquareTest(int pValueMode,
                            const vector< vector<double> > & null_prob
                            = vector< vector<double> >(0));
	
	void applyChisqPermTest(int pValueMode, const CPermStatTable & permtable,int maxDelay, int minDelay);
    
	// Communication related functions:
	size_t numBytes() const;
	size_t pack(unsigned char buffer[], size_t size) const;
	size_t unpack(unsigned char buffer[], size_t size);
    
	void print(void) const;
	void heading(void) const;
    
	string countsToString() const;
    
    bool scan4by2(const string & file);
    void saveTable(const string & file) const;
    
    bool scan(ifstream & ifs);
    void save(ofstream & ofs) const;
    void saveDOT(ofstream & ofs) const;
    
    //Added by Haizhou Wang, Aug 2, 2012
    bool pure() const; //This function check to see if each row in m_transitionTable has at most one non-zero cell.
    //End of Haizhou's code
    
    TransitionTable subtable(const vector <int> & parentSubset) const;
    
    size_t nrow () const { return m_transitionTable.size(); }
    size_t ncol () const {
        if(m_transitionTable.size()>0) {
            return m_base;
        } else {
            return 0;
        }
    }
    
    double cor() const;
    
protected:
    /*
     Comment added by Haizhou Wang, July 23, 2012
     This m_transitionTable actually is a contingency table.
     It records how many times each child value show up under different parent combination
     */
	vector< vector<int> > m_transitionTable;
    
	double m_chisq;		// chisquare value
	int m_df;			// degrees of freedom
	double m_pValue;	// p-value
    
	double m_adjustedpValue; // p-value adjusted for multiple comparison
    
	//int m_pValueMode; // chisquare test method
    
	// vector<int> m_truthtable;
    
};

