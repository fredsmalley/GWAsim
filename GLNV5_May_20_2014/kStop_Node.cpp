#include "kStop_Node.h"

#include <deque>
using std::deque;

#include <functional>

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <numeric>
using std::accumulate;

#include <vector>
using std::vector;


#include "Transition.h"   //In order to use inline function digitToNumber()



//This function converts a number into bit wise, with each bit's radix given in bitRadixs
const vector<int> getBits(int number, const vector<int> & bitRadixs)
{
    if(number < 0)
    {
        cerr<<"getBits(): Can not convert number less than 0"<<endl;
        exit(-1);
    }

    int sum = accumulate( bitRadixs.begin(), bitRadixs.end(), 1, std::multiplies<int>() );

    if( sum < number )
    {
        cerr<<"getBits(): Can not convert number "<<number<<" into bit form. It exceeds the maximum of the bit form you want to convert to."<<endl;
        exit(-1);
    }

    vector<int> result;
    for(size_t t=0; t<bitRadixs.size(); ++t)
    {
        sum = sum/bitRadixs[t];
        int quotient = number / sum;
        result.push_back(  quotient  );
        number -= quotient * sum;
    }

    return result;
}



const vector<int> nextVertexContainsPattern(const vector<int> & pattern, const vector<int> & unfixedRadix, int k)
{
    vector<int> bitsK( getBits(k, unfixedRadix) );

	vector<int> result;
	for(size_t r=0,t=0; t<pattern.size(); ++t)
	{
		(pattern[t] == -1)?
			result.push_back( bitsK[r++] ):
			result.push_back( pattern[t] );
	}
	
    return result;
}


const deque<int> getVertexContainPattern(const vector<int> & pattern, const vector<int> & radixs)
{
     //Will enumerate bits with undefined value
    vector<int> unfixedRadix;
	int numberOfEnum = 1;
    for(size_t i=0; i<pattern.size(); ++i) 
	{
		if( pattern[i] == -1 )
		{
			numberOfEnum *= radixs[i];
            unfixedRadix.push_back( radixs[i] );
		}
	}


    deque<int> vertices;
    if( numberOfEnum == 1 )
        vertices.push_back( digitsToNumber(pattern, radixs) );
    else
    {
        for(int k=0; k<numberOfEnum; ++k)
            vertices.push_back( digitsToNumber( nextVertexContainsPattern(pattern, unfixedRadix, k), radixs)  ); 
    }

	return vertices;
}


template<typename T>
const bool sameSign( const vector<T> & v)
{

    if( v.empty() ) return true;

    auto sign = v.front();
    for( int i=1; i<(int)v.size(); ++i)
    {
        if( v[i] == 0 ) continue; //Numeric '0' has no sign, or always have the same sign with others.

        sign *= v[i];
        if( sign < 0 )
            return false;
    }
    return true;
}



kStop_Node::kStop_Node(const kStop_Node & rhs):m_ID(-1),m_color(WHITE),m_mark(-1),m_ParentID(-1),m_adjList()
{
    m_ID = rhs.m_ID;
    m_color = rhs.m_color;
    m_mark = rhs.m_mark;
    m_ParentID = rhs.m_ParentID;
    m_adjList = rhs.m_adjList;
}


kStop_Node & kStop_Node::operator=(const kStop_Node & rhs)
{
    if( m_ID == rhs.m_ID )
        return *this;

    m_ID = rhs.m_ID;
    m_color = rhs.m_color;
    m_mark = rhs.m_mark;
    m_ParentID = rhs.m_ParentID;
    m_adjList = rhs.m_adjList;
    
    return *this;
}


void kStop_Node::buildAdjList( const DS_LocalConstraint & lc )
{
    /*
    1. Convert vertex ID into bit form
	2. For each bit, based on local constraints, compute its next value. If no local constraint for this bit, its next value will be unknown. All bits consist the pattern for next vertices
	3. Call getVertexContainPattern() to get a list of vertices containing the pattern.
	4. Get the ID of those vertices from step3, store them in vertex ID's adjacent list.
	*/

    const vector<int> bitV( getBits(m_ID, lc.m_radixs) );
    vector<int> nextPattern(lc.m_radixs.size(), -1); 
    
    //For each bit (representing a node in the original topological graph), compute its next value
	for(int j=0; j<(int)bitV.size(); ++j)
	{
        if(lc.m_compulsiveParents[j].m_nodeList.empty() ) //No parent, this bit remains -1
            continue;

        if( sameSign( lc.m_compulsiveParents[j].m_polarityRelations ) )
        {
            int parentSum = 0;

		    //Traverse parents of node j. Compute their effects on node j
		    for(int i=0; i<(int)lc.m_compulsiveParents[j].m_nodeList.size(); ++i)
		    {
                int parentID = lc.m_compulsiveParents[j].m_nodeList[i]-1;   //In m_nodeList, parents ID start from 1. Need to subtract 1 here.
                parentSum += bitV[parentID];
            
		    }
		    
            if( parentSum == (int)lc.m_compulsiveParents[j].m_nodeList.size() )
                lc.m_compulsiveParents[j].m_polarityRelations.front() > 0 ?
                    nextPattern[j] = 1 :
                    nextPattern[j] = 0;
            else if( parentSum == 0 )
                lc.m_compulsiveParents[j].m_polarityRelations.front() > 0 ?
                    nextPattern[j] = 0 :
                    nextPattern[j] = 1;
            else
                nextPattern[j] = -1;

        }
        else
        {
            int positiveParentSum = 0;
            int negativeParentSum = 0;

            int numOfPositiveParent = 0, numOfNegativeParent = 0;

		    //Traverse parents of node j. Compute their effects on node j
		    for(int i=0; i<(int)lc.m_compulsiveParents[j].m_nodeList.size(); ++i)
		    {
                int parentID = lc.m_compulsiveParents[j].m_nodeList[i]-1;   //In m_nodeList, parents ID start from 1. Need to subtract 1 here.

                if( lc.m_compulsiveParents[j].m_polarityRelations[i]>0 )
                {
                    ++numOfPositiveParent;
                    positiveParentSum += bitV[parentID];
                }
                else
                {
                    ++numOfNegativeParent;
                    negativeParentSum += bitV[parentID];
                }
		    }
            

            if( positiveParentSum == numOfPositiveParent && negativeParentSum == 0)
                nextPattern[j] = 1;
            else if( positiveParentSum == 0 && negativeParentSum == numOfNegativeParent )
                nextPattern[j] = 0;
            else
                nextPattern[j] = -1;
        }
        
	}

    //Based on the next pattern, enumerate all possible next vertices
	//const deque<int> nextVertices( getVertexContainPattern(nextPattern, lc.m_radixs) );

    //Store this node's the adjacent list
    //m_adjList.insert( nextVertices.begin(), nextVertices.end() );
    m_adjList = getVertexContainPattern(nextPattern, lc.m_radixs);
}




void kStop_Node::showAdj() const
{
    cout<<"Node ID = "<<m_ID<<endl<<"List:"<<endl;
    for(auto t=m_adjList.begin(); t!=m_adjList.end(); ++t)
    {
		  cout<<(*t)<<" ";
    }
    cout<<endl<<endl;
}
