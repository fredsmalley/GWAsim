// TransitionTable.cpp
//
// Joe Song
// Created: September 27, 2008
// Modified:
//   September 12, 2011. Added different options for parent comparison
//   January 24, 2013.  Extracted friend functions for comparative chisq
//     analysis to ComparativeChisq.cpp

#include <iostream>
#include <fstream>
#include <sstream>
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;
//using std::exit;

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "TransitionTable.h"
#include "ChisqTests.h"
#include "GLN.h"

//inline
int digitsToNumber(const vector<int> & digits, const vector<int> & bases);

TransitionTable operator - (const TransitionTable & tt1, const TransitionTable & tt2)
{
	TransitionTable tt = tt1;
    
	for(size_t r=0; r < tt1.m_transitionTable.size(); r++) {
		for(size_t c=0; c < tt1.m_transitionTable[r].size(); c++) {
			tt.m_transitionTable[r][c] = (tt1.m_transitionTable[r][c] - tt2.m_transitionTable[r][c]); // absolute value?
		}
	}
	
	return tt;
}

TransitionTable & abs(TransitionTable & tt)
{
	for(size_t r=0; r < tt.m_transitionTable.size(); r++) {
		for(size_t c=0; c < tt.m_transitionTable[r].size(); c++) {
			tt.m_transitionTable[r][c] = abs(tt.m_transitionTable[r][c]);
		}
	}
	return tt;
}

TransitionTable::TransitionTable(const int child, const vector< int > & parents,
                                 const vector< int > & delays,
								 const TrajectoryCollection & trajCol)
: m_chisq(0), m_df(0), m_pValue(INITIAL_P_VALUE), m_adjustedpValue(INITIAL_P_VALUE)
{
	initialize(child, parents, delays, trajCol);
    
	/*
     m_child = child;
     m_base = trajCol.base( child );
     
     m_parents = parents;
     
     m_parentBases.resize(parents.size());
     for(size_t m=0; m < m_parentBases.size(); m++){
     m_parentBases[m] = trajCol.base( parents[m]-1 );
     }
     
     m_delays = delays;
     allocate();
     */
}

TransitionTable::TransitionTable
(const int child, const int base, const vector< int > & parents,
 const vector< int > & parentBases, const vector< int > & delays,
 const string & nameChild, const vector<string> & nameParents)
: m_chisq(0), m_df(0), m_pValue(INITIAL_P_VALUE), m_adjustedpValue(INITIAL_P_VALUE)
{
    Transition::initialize(child, base, parents, parentBases,
                           delays, nameChild, nameParents);
}

/*There is a potential off-by-one error when using this function.
It assumes IDs stored in the incoming vector 'parents' start from 1
*/
void TransitionTable::initialize(const int child, const vector< int > & parents, const vector< int > & delays, 
								 const TrajectoryCollection & trajCol) 

{
	m_child = child;
    
    m_nameChild = trajCol.getNodeNames(vector<int>(1, child+1))[0];
    
	m_base = trajCol.base( child );
    
	m_parents = parents;
    
    m_nameParents = trajCol.getNodeNames(parents);
    
	m_parentBases.resize(parents.size());
	for(size_t m=0; m < m_parentBases.size(); m++){
		m_parentBases[m] = trajCol.base( parents[m]-1 );
	}
    
	m_delays = delays;
	allocate();
}

const TransitionTable &
TransitionTable::operator += (const TransitionTable & tt)
{
	// check compatibility of the two tables
	if( m_transitionTable.size() > 0 ) {

        if(m_transitionTable.size() != tt.m_transitionTable.size() ||
           m_transitionTable[0].size() != tt.m_transitionTable[0].size()) {
            
            cerr << "ERRRO in TransitionTable::operator +=: incompatible tables!"
            << endl;
            
            exit(EXIT_FAILURE);
        }
        
		for(size_t r = 0; r < nrow(); r++) {
			for(size_t c = 0; c < m_transitionTable[r].size(); c++) {
				m_transitionTable[r][c] += tt.m_transitionTable[r][c];
			}
		}
        
	} else {
        
        *this = tt;
        
    }
    
	return *this;
}


void TransitionTable::allocate()
{
	size_t nRows = calculateSize();
    
	m_transitionTable.resize(nRows);
    
	if(nRows > 0) {
		for(size_t i=0; i<nRows; i++) {
			m_transitionTable[i].resize(m_base);
			for(int j=0; j<m_base; j++) {
				m_transitionTable[i][j] = 0;
			}
		}
	}
}

vector<int> TransitionTable::getRowSums() const
{
    return ::getRowSums(m_transitionTable);
}

//added by Yang Zhang 1.2.2013
vector<int> TransitionTable::getColSums() const
{
    return ::getColSums(m_transitionTable);
}

int TransitionTable::getTotalSum() const
{
    return ::getTotalSum(m_transitionTable);
}

void TransitionTable::accumulate(int value, const vector<int> & parentValues)
{
	// Calculate the transition table index (share code)
	int row = digitsToNumber(parentValues, m_parentBases);
    
	// Fill in the right entry in the transition table
	m_transitionTable[row][value]++;
    
}

void TransitionTable::applyChisquareTest(int pValueMode,
                                         const vector< vector< double > > & null_prob)
{
    size_t df;
	m_pValue = ChisqTest(m_transitionTable, m_base, pValueMode,
                         m_chisq, df, null_prob);
    m_df = (size_t) df;
}

void TransitionTable::applyChisqPermTest(int pValueMode,
                                         const CPermStatTable & permtable,
                                         int maxDelay, int minDelay)
{
    size_t df;
	ChisqTest(m_transitionTable, m_base, pValueMode, m_chisq, df);
    m_df = (size_t) df;
    
	m_pValue = permtable.compute_pValue(m_child, (int) m_parents.size(),
                                        maxDelay, minDelay, m_chisq);
}

void TransitionTable::print(void) const
{
	cout << "Child=" << m_nameChild << '[' << m_child + 1 << ']'
    << ", with " << m_parents.size() << " parents = ";
    
    cout << '{';
	for(size_t i=0; i < m_parents.size(); i++) {
		cout << m_nameParents[i] << '[' << m_parents[i] << "]" << m_delays[i];
		if(i < m_parents.size() - 1) {
			cout << ';';
		}
	}
	cout << '}' << " (name[ID]delay)" << endl;
    
	for(size_t j=0; j < m_transitionTable.size(); j++){
		for(size_t k=0; k < m_transitionTable[j].size(); k++){
			cout << m_transitionTable[j][k] << "\t";
		}//end for
		cout << endl;
	}//end for
	cout << "\n";
}

void TransitionTable::heading(void) const
{
	cout << m_parents.size() << "\t";
    
	for(size_t j=0; j < getParents().size(); j++) {
		cout << m_parents[j] << ",";
	}
	cout << "\t";
    
	for (size_t j = 0; j < m_delays.size(); j++) {
        cout << m_delays[j] << ",";
	}//end if
	cout << "\t";
    
	cout << getpValue() << "\t";
	cout << getChisq() << "\t";
	cout << getDf() << endl;
}

string TransitionTable::countsToString() const
{
    stringstream ttStream;
    for(size_t j=0; j < m_transitionTable.size(); j++){
        for(size_t k=0; k < m_transitionTable[j].size(); k++){
            ttStream << m_transitionTable[j][k] << " ";
        }//end for
        ttStream << endl;
    }//end for
    ttStream << "\n";
    
    return ttStream.str();
}

void TransitionTable::fill(const TrajectoryCollection & trajCol)
{
	size_t k = m_parents.size();

    
    
    /*If not return when there is no parent,
     int tBegin = - * ( min_element(m_delays.begin(), m_delays.end()) );
     will crash.
     Since the dereference of a vector::end() is invalid.
     */
    if( k == 0 )
        return;


	vector<int> parentValues(k);
    
	// Compute the starting time of the trajectory
    
	//int maxMO = 0;
	//for(size_t i=0; i<k; i++) {
	//	if(m_delays[i] < maxMO) {
	//		maxMO = m_delays[i];
	//	}
	//}
	// int tBegin = - maxMO;
    
    // Simplified by using min_element()  MS Nov 21, 2011
    int tBegin = - * ( min_element(m_delays.begin(), m_delays.end()) );
    
	for(size_t j = 0; j < trajCol.getNumTrajectories(); j ++) {
        
        //if the child node is excluded from this trajectory,
        //then we simply pass this one by letting j increment 1
		if(!trajCol.getexcludeTrajectory().empty())
		{
            if(trajCol.getexcludeTrajectory()[j].end()
               != find(trajCol.getexcludeTrajectory()[j].begin(),
                       trajCol.getexcludeTrajectory()[j].end(),
                       m_child)
               ) { // MS added Nov 19, 2011
                continue;
            }
            
            /* MS commented out on Nov 19, 2011.  The following code
             seems to have a logical error with i++
             
             bool flag = false;
             for(size_t i=0; i < trajCol.getexcludeTrajectory()[j].size(); i++)
             {
             if((trajCol.getexcludeTrajectory())[j][i]==m_child)
             {
             flag = true;
             }
             i++;
             }
             
             if(flag)
             {
             continue;
             }
             */
		}
		for(size_t t = tBegin; t < (size_t) trajCol.getTrajectoryLength((int) j); t++){
			
			bool incomplete = false;
			
            if( trajCol.hasMissingValues() ) { // MS added Nov 20, 2011
                
                if(trajCol.missing( (int) j, (int) t, m_child ))
                {
                    incomplete = true;
                }
                else
                {
                    for(size_t i=0; i < k; i ++){
                        if(trajCol.missing( (int) j, (int) (t + m_delays[i]), m_parents[i]-1))
                        {
                            incomplete = true;
                            break;
                        }
                        parentValues[i] = trajCol.get( (int)j, (int) (t + m_delays[i]), m_parents[i]-1 );
                    }
                }
                
            } else {
                
                for(size_t i=0; i < k; i ++){
                    parentValues[i] = trajCol.get( (int)j, (int) (t + m_delays[i]), m_parents[i]-1 );
                }
                
            }
			
			if(!incomplete)
			{
				accumulate( trajCol.get((int) j, (int) t, m_child), parentValues );
			}
            
		}
	}
    
}

void TransitionTable::condFill(const TrajectoryCollection & trajCol, int condition, int conditionDelay, int level)
{
	size_t k = m_parents.size();
    
    if (k == 0)
        return;
    
	vector<int> parentValues(k);
    
	// Compute the starting time of the trajectory
    
    int tBegin = - * ( min_element(m_delays.begin(), m_delays.end()) );
    
	for(size_t j = 0; j < trajCol.getNumTrajectories(); j ++) {
        
        //if the child node is excluded from this trajectory,
        //then we simply pass this one by letting j increment 1
		if(!trajCol.getexcludeTrajectory().empty()){
            if(trajCol.getexcludeTrajectory()[j].end()
               != find(trajCol.getexcludeTrajectory()[j].begin(),
                       trajCol.getexcludeTrajectory()[j].end(),
                       m_child)
               )
                continue;
		}
		for(size_t t = tBegin; t < (size_t) trajCol.getTrajectoryLength((int) j); t++){
            
            // THIS FUNCTION IS IDENTICAL TO fill(...) except for the following lines
            if(trajCol.missing((int) j, (int) (t + conditionDelay), condition)) continue;
            
            if(trajCol.get((int) j, (int) (t + conditionDelay), condition) != level ) continue;
			
			bool incomplete = false;
			
            if( trajCol.hasMissingValues() ) {
                
                
                if(trajCol.missing( (int) j, (int) t, m_child )) incomplete = true;
                
                else
                    for(size_t i=0; i < k; i ++){
                        if(trajCol.missing( (int) j, (int) (t + m_delays[i]), m_parents[i]-1)) {
                            incomplete = true;
                            break;
                        }
                        parentValues[i] = trajCol.get( (int)j, (int) (t + m_delays[i]), m_parents[i]-1 );
                    }
                
            }
            else
                for(size_t i=0; i < k; i ++)
                    parentValues[i] = trajCol.get( (int) j,(int) (t + m_delays[i]), m_parents[i]-1 );
			
			if(!incomplete) accumulate( trajCol.get((int) j, (int) t, m_child), parentValues );
            
		}
	}
    
}

void TransitionTable::reset()
{
    for(size_t i = 0; i < m_transitionTable.size(); i++)
        for(size_t j = 0; j < m_transitionTable[i].size(); j++)
            m_transitionTable[i][j] = 0;
}

size_t TransitionTable::numBytes() const
{
	//child + base + bestparents +  parentbases + bestoffsets + numparents + row + transtable + df
	size_t intelementsSize = 1 + 1 + m_parents.size()
    + m_parentBases.size() + m_delays.size()
    + 1 + 1 + (m_base *  m_transitionTable.size()) + 1;
	
    // child name + parent names
    size_t sizeParentNames = 0;
    
    for(size_t i=0; i<m_nameParents.size(); i++) {
        sizeParentNames += m_nameParents[i].size() + 1;
    }
    
    size_t charelementSize = m_nameChild.size() + 1 + sizeParentNames;
    
	size_t doubleelementsSize = 1 + 1; //chisq + pvalue
    
	size_t totalSize = sizeof(int) * intelementsSize + sizeof(char) * charelementSize + sizeof(double) * doubleelementsSize;
    
	return totalSize;
}

size_t TransitionTable::pack(unsigned char buffer[], size_t size) const
{
	size_t nBytes = numBytes();
    
	assert(buffer);
    
	if( (! buffer) || size < nBytes ) {
		cerr << "ERROR: the buffer is either empty or insufficient for packing!" << endl;
		exit(EXIT_FAILURE);
	}
    
	unsigned char * p = buffer;
    
	*(int*) p = m_child;
	p += sizeof(int);
    
    for(size_t i=0; i<m_nameChild.size(); i++) {
        * p ++ = m_nameChild[i];
    }
    *p ++ = '\0';
    
	*(int*) p = m_base;
	p += sizeof(int);
    
    
	*(int*) p = (int) m_parents.size();
	p += sizeof(int);
    
	*(int*) p = (int) m_transitionTable.size();
    
	for (size_t i = 0; i < m_transitionTable.size(); i++){
		for (size_t j = 0; j < m_transitionTable[i].size(); j++){
			p += sizeof(int);
			*(int*) p = m_transitionTable[i][j];
		}
	}
    
	for (size_t i = 0; i < m_parents.size(); i++){
		p += sizeof(int);
		*(int*) p = m_parents[i];
	}
    p += sizeof(int);
    
	for (size_t i = 0; i < m_nameParents.size(); i++){
        for (size_t j=0; j < m_nameParents[i].size(); j++) {
            * p ++ = m_nameParents[i][j];
        }
        * p ++ = '\0';
	}
    
	for (size_t i = 0; i < m_parentBases.size(); i++){
		*(int*) p = m_parentBases[i];
		p += sizeof(int);
	}
    
	for (size_t i = 0; i < m_delays.size(); i++){
		*(int*) p = m_delays[i];
		p += sizeof(int);
	}
    
	// p += sizeof(int);
    
	*(double*) p = m_chisq;
	p += sizeof(double);
    
	*(int *) p = m_df;
	p += sizeof(int);
    
	*(double*) p = m_pValue;
    
	return nBytes;
}

size_t TransitionTable::unpack(unsigned char buffer[], size_t size)
{
	unsigned char * p = buffer;
    
	m_child = *(int*) p;
	p += sizeof(int);
    
    m_nameChild = (char *) p;
    p += m_nameChild.size() + 1;
    
	m_base = *(int*) p;
	p += sizeof(int);
    
	int num_parents = *(int*) p;
	p += sizeof(int);
    
	int ttrows = *(int*) p;
	
	// trans_table = createTable(size_transtable, penv->m_trajCol.base(nodeId-1) );
	m_transitionTable.resize(ttrows);
	for(int i=0; i<ttrows; i++) {
		m_transitionTable[i].resize( m_base );
	}
    
	for (int i = 0; i < ttrows; i++){
		for (int j = 0; j < m_base; j++){
			p += sizeof(int);
			m_transitionTable[i][j] =*(int*) p;
		}
	}
    
	m_parents.resize(num_parents);
	
	for (int i = 0; i < num_parents; i++){
		p += sizeof(int);
		m_parents[i] = *(int*) p;
	}
    
    p += sizeof(int);
    
	m_nameParents.resize(num_parents);
	for (int i = 0; i < num_parents; i++){
		m_nameParents[i] = (char *) p;
        p += m_nameParents[i].size() + 1;
	}
    
	m_parentBases.resize(num_parents);
	
	for (int i = 0; i < num_parents; i++){
		m_parentBases[i] = *(int*) p;
		p += sizeof(int);
	}
    
	m_delays.resize(num_parents);
    
	for (int i = 0; i < num_parents; i++){
		m_delays[i] = *(int*) p;
		p += sizeof(int);
	}
	
    // p += sizeof(int);
    
	m_chisq = *(double*) p;
	p += sizeof(double);
    
	m_df = *(int *) p;
	p += sizeof(int);
    
	m_pValue = *(double*) p;
	p += sizeof(double);
    
	return (size_t) (p - buffer);
}



/*Added by Haizhou Wang, Aug 2, 2012
 This function check to see if each row in m_transitionTable 
 has at most one non-zero cell, i.e. given a combination of
 parents, the can only be a unique child value.
 If more than one cells is non-zero, return false;
 */
bool TransitionTable::pure() const
{
    for(size_t i=0; i<m_transitionTable.size(); ++i)
    {
        bool nonZero_exist = false;
        for(size_t j=0; j<m_transitionTable[i].size(); ++j){
            if( m_transitionTable[i][j] != 0 ) {
                if(nonZero_exist == true) {
                    return false;
                } else {
                    nonZero_exist = true;
                }
            }
        }
    }
    return true;
}

TransitionTable
TransitionTable::subtable(const vector <int> & parentSubset) const
// Extract a subtable from the transition table that only covers
//   the given parent subset.
// ASSUMPTIONS:
//   1. Parents in the original table are expected to be unique
{
    if (parentSubset.empty()) {
        return TransitionTable();
    }
    
    // Identify relative parent index in the original parent list
    vector<int> subsetIndices(parentSubset.size());
    for(size_t i = 0; i<parentSubset.size(); i++) {
        subsetIndices[i] = find(m_parents.cbegin(), m_parents.cend(),
                                parentSubset[i]) - m_parents.cbegin();
        if(subsetIndices[i] == (int) m_parents.size()) {
            cerr << "ERROR in TransitionTable.sub(): invalid parent subset!"
            << endl;
            exit(EXIT_FAILURE);
        }
    }

    vector<int> parentSubsetBases(parentSubset.size());
    vector<int> parentSubsetDelays(parentSubset.size());
    vector<string> parentSubsetNames(parentSubset.size());

    for(size_t i = 0; i<parentSubset.size(); i++) {
        parentSubsetBases[i] = m_parentBases[subsetIndices[i]];
        parentSubsetDelays[i] = m_delays[subsetIndices[i]];
        parentSubsetNames[i] = m_nameParents[subsetIndices[i]];
    }
    
    // Initialize a new transition table with id, base, and delay of
    //   the parent subset
    TransitionTable sub(m_child, m_base, parentSubset, parentSubsetBases,
                        parentSubsetDelays, m_nameChild, parentSubsetNames);

    // Go through each row and map each row in the original table
    //   to the sub-table
    vector<int> parentValues(m_parents.size());
    vector<int> subsetParentValues(parentSubset.size());
    
    for (long long unsigned row=0; row<m_transitionTable.size(); row++) {
        LinearIndexToArrayIndex(m_parentBases, row, parentValues);
        for(size_t i=0; i<parentSubset.size(); i++) {
            subsetParentValues[i] = parentValues[subsetIndices[i]];
        }
        long long unsigned rowSub = ArrayIndexToLinearIndex(parentSubsetBases,
                                                            subsetParentValues);
        
        for (size_t col=0; col<m_transitionTable[row].size(); col++) {
            sub.m_transitionTable[rowSub][col] += m_transitionTable[row][col];
        }
    }
    return sub;
}

double TransitionTable::cor() const
{
    double corCoeff;
    int n = 0;
    long sxy=0, sx=0, sy=0, sxsq=0, sysq=0;
    for (size_t x=0; x<m_transitionTable.size(); x++) {
        for(size_t y=0; y<m_transitionTable[x].size(); y++) {
            sxy += m_transitionTable[x][y] * x * y;
            sx += m_transitionTable[x][y] * x;
            sy += m_transitionTable[x][y] * y;
            sxsq += m_transitionTable[x][y] * x * x;
            sysq += m_transitionTable[x][y] * y * y;
            n += m_transitionTable[x][y];
        }
    }
    
    double ux = sx / (double) n;
    double uy = sy / (double) n;
    double uxy = sxy / (double) n;
    double ux2 = sxsq / (double) n;
    double uy2 = sysq / (double) n;
    
    corCoeff = (uxy - ux * uy) / ( sqrt(ux2-ux*ux) * sqrt(uy2-uy*uy) );
    
    return corCoeff;
}
