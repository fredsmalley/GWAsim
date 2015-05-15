//FisherComparative.cpp
//
//created by Yang Zhang 1.18.2013
//
// Modified:
//   MS 2/3/2013.
//    1. Changed the code to handle more than two conditions
//    2. Improved the memory efficiency by removing unnecessary intermediate
//       data structures.
//    3. Reuse the ArrayIndexToLinearIndex() and LinearIndexToArrayIndex(),
//       two inline functions already written in GLN.h.

#include <cmath>
#include <cassert>

using namespace std;

#include "SetOps.h"
#include "TransitionTable.h"
#include "GLNGlobals.h"
#include "ChisqTests.h"
#include "ComparativeChisq.h"

vector<vector<TransitionTable> >
getPooledAndUnique(const vector<TransitionTable> & tts);

vector<vector<vector<double> > >
getExpectedTables(const vector<TransitionTable> & tts);

//compute the fisher's heterogeity exact test p-value
void exact_comparative_chisq_1(const vector<TransitionTable> & tts,
                               double &pd, double &pc, double &pt,
                               int pvaluemode)
{
	double qd = 1.0;
	double qt = 1.0;
    
	vector<vector<TransitionTable> > pooledandunique = getPooledAndUnique(tts);
	TransitionTable pooledTable = pooledandunique[0][0];
    
	vector<vector<vector<double> > > expectedTable = getExpectedTables(tts);
    
	for(size_t i=0; i<tts.size(); i++)
	{
		vector<vector<int> > originaltt = tts[i].getTransitionTable();
		
		if(originaltt.size() > 0)
		{
			vector<vector<double> > expectedtt = getExpectedTable(originaltt);
            
			double pdi = fisher_exact_test(originaltt, expectedTable[i]);
            
			double pti = 1.0;
            
			if(pvaluemode == 8)
			{
				pti = fisher_exact_test(originaltt);
			}
			else if(pvaluemode == 9)
			{
				pti = fisher_exact_test(originaltt, expectedtt);
			}
			qd *= abs(1 - pdi);
			qt *= abs(1 - pti);
		}
		
	}
	pd = 1 - qd;
	pt = 1 - qt;
    
	if(pvaluemode == 8)
	{
		pc = fisher_exact_test(pooledTable.getTransitionTable());
	}
	else if(pvaluemode == 9)
	{
		vector<vector<int> > pooled = pooledTable.getTransitionTable();
		if(pooled.size() > 0)
		{
			vector<vector<double> > expectedpc = getExpectedTable(pooled);
			pc = fisher_exact_test(pooled, expectedpc);
		}else
		{
			pc = 1;
		}
	}
    
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// TO BE DELETED FROM BELOW

//convert row index in contingency table into parent values, index starts from 0
vector<int> indexToValues(int index, const vector<int> & bases)
{
	vector<int> bitsigs; //denote significance bits for each parent position
    
    bitsigs.resize(bases.size(), 1);
    
	vector<int> parentvalues;
    
	parentvalues.resize(bases.size());
    
	for(size_t i=0; i<bitsigs.size(); i++)
	{
		for(size_t j=i+1; j<bases.size(); j++)
		{
			bitsigs[i] = bitsigs[i]*bases[j];
		}
        
		parentvalues[i] = index/bitsigs[i];
		index -= bitsigs[i]*parentvalues[i];
	}
    
	return(parentvalues);
}

//convert parent values into row index in contingency table, index starts from 0.
int valueToIndex(const vector<int> & parentvalues, const vector<int> & bases)
{
	assert(parentvalues.size()==bases.size());
    
	vector<int> bitsigs; //denote significance bits for each parent position
    
    bitsigs.resize(parentvalues.size(), 1);
    
	int rowindex = 0;
    
	for(size_t i=0; i<bitsigs.size(); i++)
	{
		for(size_t j=i+1; j<bases.size(); j++)
		{
			bitsigs[i] = bitsigs[i]*bases[j];
		}
        
		rowindex += parentvalues[i]*bitsigs[i];
	}
    
	return(rowindex);
    
}

//compute pooled transition table and transition table with unique parents for each condition
//the first dimension is pooled, second dimension is for unique
vector<vector<TransitionTable> >
getPooledAndUniqueOld(const vector<TransitionTable> & tts,
                      const vector<TrajectoryCollection> & m_trajCols)
{
    vector<vector<TransitionTable> > ttspoolandhetero;
    ttspoolandhetero.resize(2);
    ttspoolandhetero[0].resize(1);
    ttspoolandhetero[1].resize(2);
    
    vector<vector<int> > parentsvec;
    parentsvec.resize(tts.size());
    
    for(size_t i=0; i<tts.size(); i++)
    {
        parentsvec[i] = tts[i].getParents();
    }
    vector<int> intersecParents = intersectionSet(parentsvec[0],parentsvec[1]);
    vector<vector<int> > uniqueParents;
    uniqueParents.resize(tts.size());
    uniqueParents[0] = complementSet(intersecParents, parentsvec[0]);
    uniqueParents[1] = complementSet(intersecParents, parentsvec[1]);
    vector<vector<int> > uniqueDelays;
    uniqueDelays.resize(tts.size());
    vector<vector<int> > pooledDelays;
    pooledDelays.resize(tts.size());
    vector<vector<int> > uniqueBases;
    uniqueBases.resize(tts.size());
    vector<vector<int> > pooledBases;
    pooledBases.resize(tts.size());
    
    for(size_t i=0; i<tts.size(); i++)
    {
        vector<int> tempDelays = tts[i].getDelays();
        vector<int> tempParents = tts[i].getParents();
        vector<int> tempParentBases = tts[i].getParentBases();
        for(size_t parentsindex = 0; parentsindex<tempParents.size(); parentsindex++)
        {
            bool matchflag = false;
            for(size_t pooledindex = 0; pooledindex<intersecParents.size(); pooledindex++)
            {
                if(intersecParents[pooledindex] == tempParents[parentsindex])
                {
                    matchflag = true;
                    pooledDelays[i].push_back(tempDelays[parentsindex]);
                    pooledBases[i].push_back(tempParentBases[parentsindex]);
                }
            }
            if(!matchflag)
            {
                uniqueDelays[i].push_back(tempDelays[parentsindex]);
                uniqueBases[i].push_back(tempParentBases[parentsindex]);
            }
        }
        
    }
    
    //generate transition tables with pooled parents for each condition
    TransitionTable ttpooled = TransitionTable(tts[0].getChild(),intersecParents, pooledDelays[0], m_trajCols[0]);
    ttpooled.fill(m_trajCols[0]);
    for(size_t i=1; i<tts.size(); i++)
    {
        ttpooled += TransitionTable(tts[i].getChild(),intersecParents, pooledDelays[i], m_trajCols[i]);
        ttpooled.fill(m_trajCols[i]);
    }
    ttspoolandhetero[0][0] = ttpooled;
    
    //generate transition tables with unique parents for each condition
    for(size_t i=0; i<tts.size(); i++)
    {
        ttspoolandhetero[1][i] = TransitionTable(tts[i].getChild(),uniqueParents[i], uniqueDelays[i], m_trajCols[i]);
        ttspoolandhetero[1][i].fill(m_trajCols[i]);
    }
    return(ttspoolandhetero);
}

//compute expect contingency table for each condition
//currently restricted to 2 conditions
vector<vector<vector<double> > >
getExpectedTableOld(const vector<TransitionTable> & tts,
                    const vector<TrajectoryCollection> & m_trajCols)
{
    vector<vector<vector<double> > > ttsexpected; //only contains the contingency table
    ttsexpected.resize(tts.size());
    
    int counts = tts[0].getTotalSum(); //total number of observations
    for(size_t i=1; i<tts.size(); i++)
    {
        counts += tts[i].getTotalSum();
    }
    
    vector<vector<TransitionTable> > pooledandunique = getPooledAndUniqueOld(tts, m_trajCols);
    TransitionTable ttpooled = pooledandunique[0][0];
    vector<vector<int> > ctpooled = ttpooled.getTransitionTable(); //pooled contingency table
    
    vector<TransitionTable> ttshetero = pooledandunique[1];
    
    //compute overall expected table
    for(size_t i=0; i<tts.size(); i++)
    {
        int tempcounts = tts[i].getTotalSum();
        vector<int> parentbases = tts[i].getParentBases();
        vector<vector<int> > ttint = tts[i].getTransitionTable();
        vector< vector<double> > temptt;
        temptt.resize(ttint.size());
        if(ttint.size()>0)
        {
            for(size_t rowindex=0; rowindex<temptt.size();rowindex++)
            {
                temptt[rowindex].resize(ttint[rowindex].size());
                for(size_t colindex=0;colindex<temptt[rowindex].size();colindex++)
                {
                    temptt[rowindex][colindex] = ttint[rowindex][colindex];
                }
                
            }
            
            vector<int> rowsums = tts[i].getRowSums();
            vector<int> colsums = tts[i].getColSums();
            vector<int> heterorowsums = ttshetero[i].getRowSums();
            vector<int> heterocolsums = ttshetero[i].getColSums();
            
            vector<int> pooledparentvalues;
            vector<int> pooledparents = ttpooled.getParents();
            vector<int> pooledBases = ttpooled.getParentBases();
            pooledparentvalues.resize(pooledparents.size());
            
            //compute unique parent values
            vector<int> uniqueParents = ttshetero[i].getParents();
            vector<int> uniqueBases = ttshetero[i].getParentBases();
            vector<int> uniqueparentvalues;
            uniqueparentvalues.resize(ttshetero[i].getParents().size());
            
            for(size_t rowindex=0; rowindex<temptt.size();rowindex++)
            {
                vector<int> allparentvalues = indexToValues(rowindex, parentbases);  //bit-wise, not decimal
                
                //compute pooled parent values for each condition
                vector<int> tempParents = tts[i].getParents();
                for(size_t parentsindex = 0; parentsindex<tempParents.size(); parentsindex++)
                {
                    for(size_t pooledindex = 0; pooledindex<pooledparentvalues.size(); pooledindex++)
                    {
                        if(pooledparents[pooledindex] == tempParents[parentsindex])
                        {
                            //uniqueparentindices.push_back(parentsindex);
                            pooledparentvalues[pooledindex] = allparentvalues[parentsindex];
                        }
                    }
                }
                
                //convert parent values into indices in pooled matrices
                int pooledrowindex = valueToIndex(pooledparentvalues, pooledBases);
                
                for(size_t parentsindex = 0; parentsindex<tempParents.size(); parentsindex++)
                {
                    for(size_t uniqueindex = 0; uniqueindex<uniqueParents.size(); uniqueindex++)
                    {
                        if(uniqueParents[uniqueindex] == tempParents[parentsindex])
                        {
                            //uniqueparentindices.push_back(parentsindex);
                            uniqueparentvalues[uniqueindex] = allparentvalues[parentsindex];
                        }
                    }
                }
                
                //convert parent values into indices in hetero matrices
                int uniquerowindex = valueToIndex(uniqueparentvalues, uniqueBases);
                
                for(size_t colindex=0; colindex<temptt[0].size();colindex++)
                {
                    double p_pooled = (ctpooled.size()>0)?(double)ctpooled[pooledrowindex][colindex]/counts:(double)colsums[colindex]/tempcounts;
                    double p_hetero = (heterorowsums.size()>0)?(double)heterorowsums[uniquerowindex]/tempcounts:(1.0);
                    temptt[rowindex][colindex] = tempcounts*p_pooled*p_hetero;
                }
            }
            ttsexpected[i].resize(temptt.size());
            ttsexpected[i] = temptt;
        }
    }
    return(ttsexpected);
    
}

