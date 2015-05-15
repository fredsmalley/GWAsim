// PermStatTable.cpp
//
// Joe Song
// Created: September 7, 2008

#include <iostream>
using std::ios;
using std::cerr;
using std::endl;

#include <fstream>
using std::ifstream;

#include <vector>
using std::vector;

#include <algorithm>

#include <cstring>

#include "PermStatTable.h"

CPermStatTable::CPermStatTable(void)
{
	setNumNodes(0);
	setNumPermutations(0);
	setMaxNumParents(0);
	setNumMarkovianOrderRanges(0);
}

CPermStatTable::CPermStatTable(const char * filePermStatTable)
{
	scan(filePermStatTable);
}

void CPermStatTable::allocate() 
{
	m_permStatTable.resize(getNumNodes());

	for(int n=0; n<getNumNodes(); n++) {
		m_permStatTable[n].resize(getMaxNumParents());
		for(int i=0; i < getMaxNumParents(); i ++) {
			m_permStatTable[n][i].resize(getNumMarkovianOrderRanges());
			for(int j=0; j < getNumMarkovianOrderRanges(); j++) {
				m_permStatTable[n][i][j].resize(getNumPermutations());
			}
		}
	}

	m_highestMarkovianOrders.resize(getNumMarkovianOrderRanges());
	m_lowestMarkovianOrders.resize(getNumMarkovianOrderRanges());
}

int CPermStatTable::scan(const char * filePermStatTable)
{
	if(strlen(filePermStatTable) > 0) {

		// Open file for reading
		ifstream ifs(filePermStatTable, ios::in);

		if(!ifs.is_open()) {
			cerr << "ERROR: openning file \"" << filePermStatTable << "\"" << endl;
			exit(EXIT_FAILURE);
		}

		// File format validation by tag

		// Load header
		int numNodes, numPermutations, maxNumParents, numMarkovianOrderRanges;
		ifs >> numNodes >> numPermutations >> maxNumParents >> numMarkovianOrderRanges;

		// Set the parameters
		setNumNodes(numNodes);
		setNumPermutations(numPermutations);
		setMaxNumParents(maxNumParents);
		setNumMarkovianOrderRanges(numMarkovianOrderRanges);

		// Allocate memory
		allocate();

		for(int iMORange=0; iMORange < numMarkovianOrderRanges; iMORange ++) {
			ifs >> m_highestMarkovianOrders[iMORange];
			ifs >> m_lowestMarkovianOrders[iMORange];
		}

		// Read the table
		for(int n=0; n < numNodes; n++) {
			for(int iPerm=0; iPerm < numPermutations; iPerm++) {
				for(int theMORange=0; theMORange < numMarkovianOrderRanges; theMORange ++) {
					for(int nPars=1; nPars <= maxNumParents; nPars ++) {
						double nullStat;
						if(ifs.good()) {
							ifs >> nullStat;
							setNullStat(n, nPars, theMORange, iPerm, nullStat);
						}
					}
				}
			}
		}

		// Close file
		ifs.close();
	}
	return 0;
}

double CPermStatTable::compute_pValue(int node, int nParents, int highestMarkovianOrder, int lowestMarkovianOrder, double stat) const
{
	double pValue;
	int rank;

	vector<double> vecNullStats = getNullStats(node, nParents, highestMarkovianOrder, lowestMarkovianOrder);

	// Perform binary search 
	rank = (int) ( lower_bound(vecNullStats.begin(), vecNullStats.end(), stat) - vecNullStats.begin());

	pValue = 1.0 - rank / (double) getNumPermutations();

	// Return the p-value
	return pValue;
}
