// PermStatTable.h
//
// Joe Song
// Created: Septebmer 07, 2008

#pragma once

#include<vector>
using std::vector;

class CPermStatTable {

public: // public member functions

	CPermStatTable(void);
	CPermStatTable(const char * filePermStatTable);

	void setNumNodes(int numNodes) 
	{
		m_numNodes = numNodes;
	}

	int getNumNodes() const
	{
		return m_numNodes;
	}

	void setNumPermutations(int numPermutations) 
	{
		m_numPermutations = numPermutations;
	}

	int getNumPermutations() const
	{
		return m_numPermutations;
	}

	void setMaxNumParents(int maxNumParents)
	{
		m_maxNumParents = maxNumParents;
	}

	int getMaxNumParents() const
	{
		return m_maxNumParents;
	}

	void setNumMarkovianOrderRanges(int numMarkovianOrderRanges)
	{
		m_numMarkovianOrderRanges = numMarkovianOrderRanges;
	}

	int getNumMarkovianOrderRanges() const
	{
		return m_numMarkovianOrderRanges;
	}

	int scan(const char * filePermStatTable);

	void setNullStat(int node, int nParents, int highestMarkovianOrder, int lowestMarkovianOrder, int iPerm, double nullStat)
	{
		int whichRange = findMarkovianOrderRange(highestMarkovianOrder, lowestMarkovianOrder);
		m_permStatTable[node][nParents-1][whichRange][iPerm] = nullStat;
	}

	void setNullStat(int node, int nParents, int whichRange, int iPerm, double nullStat)
	{
		m_permStatTable[node][nParents-1][whichRange][iPerm] = nullStat;
	}

	const vector<double> & getNullStats(int node, int nParents, int highestMarkovianOrder, int lowestMarkovianOrder) const
	{
		int whichRange = findMarkovianOrderRange(highestMarkovianOrder, lowestMarkovianOrder);
		return m_permStatTable[node][nParents-1][whichRange];
	}

	double compute_pValue(int node, int nParents, int highestMarkovianOrder, int lowestMarkovianOrder, double stat) const;

private: // private data members

	vector< vector< vector< vector<double> > > > m_permStatTable;
	// m_permStatTable[node][numParents][theMOrange][iPerm]

	int m_numNodes; // Number of nodes in the table

	int m_maxNumParents; // Maximum number of parents in the table
	
	int m_numMarkovianOrderRanges; // Number of different Markovian ranges in the table

	vector<int> m_highestMarkovianOrders;  // The highest Markovian order in each range
	vector<int> m_lowestMarkovianOrders; // The lowest Markovian order in each range

	int m_numPermutations;  // The number of permutations for each parameter set

private: // private functions

	void allocate();

	int findMarkovianOrderRange(int highestMarkovianOrder, int lowestMarkovianOrder) const
	{
		int whichRange=0;

		for(int i=0; i<getNumMarkovianOrderRanges(); i++) {
			if(m_highestMarkovianOrders[i] == highestMarkovianOrder 
				&& m_lowestMarkovianOrders[i] == lowestMarkovianOrder) {
					whichRange = i;
					break;
			}
		}
		return whichRange;
	}

};
