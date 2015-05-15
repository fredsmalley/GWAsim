// GeneralizedTruthTable.cpp -- Generalized truth table class 
//
// Joe Song
// Created: November 30, 2008.
// Last modified: November 12, 2009.  Added the probability transition table

#include "GeneralizedTruthTable.h"


#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>    //In order to use ptrdiff_t
#include <ctime>      //In order to use rand()
#include <numeric>    //In order to use accumulate()
#include <functional> //In order to use multiplies<>()

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <set>
using std::set;

#include <string>
using std::string;

#include <vector>
using std::vector;



//Random number generator to be used in random_shuffle()
inline ptrdiff_t RNG( ptrdiff_t i){return rand()%i;}



GeneralizedTruthTable::GeneralizedTruthTable(void)
{
}

GeneralizedTruthTable::GeneralizedTruthTable(const int child, const int base, const vector< int > & parents, 
		const vector<int> & parentBases, const vector< int > & delays, const vector<int> & gtt)
{
	initialize(child, base, parents, parentBases, delays);
	m_truthValues = gtt;

	m_transProbs.resize(gtt.size());

	for(size_t i = 0; i < gtt.size(); i++) {

		m_transProbs[i].resize(base);

		for(int j=0; j < base; j ++) {
			m_transProbs[i][j] = 0;
		}

		m_transProbs[i][ gtt[i] ] = 1.0;

	}
}

GeneralizedTruthTable::~GeneralizedTruthTable(void)
{
}

void GeneralizedTruthTable::allocate() 
{
	m_truthValues.resize( size() );	
}

int GeneralizedTruthTable::lookup(const vector<int> & parentValues, int previousValue) const
{
	int value;

	if (parentValues.size() == 0) { // repeat prevoius value if no truth table is available

        value = previousValue;
	
	} else {

		// Calculate the transition table index (share code)	
		int row = digitsToNumber(parentValues, m_parentBases);
	
		value = m_truthValues[row];
	}

	return value;
}

void GeneralizedTruthTable::convert(const TransitionTable & transTable, double alpha)
{
	// Obtain the probabiliy transition table.  Added Nov 11, 2009.  Joe Song
	* (Transition *) this = (Transition &) transTable;

	m_transProbs.resize( transTable.nrow() );

	for(int j=0; j < (int) transTable.nrow(); j++) {

		int row_sum=0;

		for(int p=0; p < m_base; p++) {
			
			row_sum += transTable.get(j, p);

		}

		m_transProbs[j].resize(m_base, 0);

		if(row_sum > 0) {
			for(int p=0; p < m_base; p++) {
			
				m_transProbs[j][p] = transTable.get(j, p) / (double) row_sum;

			}
		}

	}

	if(transTable.getpValue() <= alpha) {

		// * (Transition *) this = (Transition &) transTable;

		m_truthValues.resize( transTable.nrow() );

		for(int j=0; j < (int) transTable.nrow(); j++) {

			// Go down the truth table and compare the columns
			int highest = transTable.get(j, 0);	//Make the first position the highest by default

			int indexat = 0;

			// Go through the columns
			// The following algorithm to derive a truth table needs to be changed
			for(int p=1; p < m_base; p++){
				if(transTable.get(j, p) > highest){
					highest = transTable.get(j, p);
					indexat = p;
				}//end if
			}//end for

			m_truthValues[j] = indexat;

		}

	} else {

		m_child = transTable.getChild();
		m_base = transTable.getBase();
		
		m_parents.clear();
		m_parentBases.clear();
		m_delays.clear();

		m_truthValues.clear();

	}

}

vector< vector<int> > convertTruthTables(const vector<TransitionTable> & transTables)
{ 
	int highest;
	int indexat;			//index of where highest is

	int size = (int) transTables.size();

	vector< vector<int> > converted_truth_tables(size);

	//Allocate the 2nd dimension
	for(int i=0; i < size; i++){

		converted_truth_tables[i].resize(transTables[i].nrow());

	}//end for

	/** for testing only, print out table of truth table sizes**/
	
	// cout << "printing out the table of truth table sizes" << endl;
	cout << endl;
	for(int i=0; i < size; i++)
		cout << transTables[i].nrow() << " ";
	cout << endl << endl;
	
	for(int i=0; i < size; i++){
		for(int j=0; j < (int) transTables[i].nrow(); j++){
			//Go down the truth table and compare the columns
			highest= transTables[i].get(j, 0);	//Make the first position the highest by default
			indexat=0;
			//Go through the columns
			// The following algorithm to derive a truth table needs to be changed
			for(int p=1; p < transTables[i].getBase(); p++){
				if(transTables[i].get(j, p) > highest){
					highest = transTables[i].get(j, p);
					indexat = p;
				}//end if
				converted_truth_tables[i][j] = indexat;
			}//end for

		}//end for

	}//end outside for

	return converted_truth_tables;

}//end convertTruthTable

bool sameChildParentsDelays(const GeneralizedTruthTable & g1, const GeneralizedTruthTable & g2)
{
	bool identical = true;

	if( g1.getChild() != g2.getChild() ) {

		identical = false; 

	} else if (g1.getParents().size() != g2.getParents().size()) {

		identical = false;

	} else {

		for(size_t j=0; j<g1.getParents().size(); j++) {

			if ( ( g1.getParents()[j] != g2.getParents()[j] ) 
				|| ( g1.getDelays()[j] != g2.getDelays()[j] ) ) {
					identical = false;
					break;
			}
		}
	}

	return identical;
}

double computeKLDivergence(const vector< double > & px, const vector< double > & qx)
{
	double divergence = 0;

	for(size_t i=0; i< px.size() && i < qx.size(); i++) {

		if( px[i] != 0 && qx[i] != 0 ) {
			divergence += px[i] * log( px[i] / qx[i] );
		}
	}

	return divergence;
}

double computeMeanSquareDistance(const vector< double > & px, const vector< double > & qx)
{
	double divergence = 0;

	size_t length = px.size() > qx.size() ? qx.size() : px.size();

	for(size_t i=0; i < length; i++) {

		divergence += ( px[i] - qx[i] ) * ( px[i] - qx[i] );

	}

	divergence /= length;

	return divergence;
}

double distanceTransProbTables(const GeneralizedTruthTable & g1, const GeneralizedTruthTable & g2)
{
	double distance = 0;

	if( sameChildParentsDelays(g1, g2) ) {

		const vector< vector<double> > & tpt1 = g1.getTransProbTable();
		const vector< vector<double> > & tpt2 = g2.getTransProbTable();
		
		for(size_t r=0; r < tpt1.size(); r++) {
			// distance += (computeKLDivergence(tpt1[r], tpt2[r]) + computeKLDivergence(tpt2[r], tpt1[r])) / 2;
			distance += computeMeanSquareDistance(tpt1[r], tpt2[r]);
		}

		distance /= g1.size();

	} else {

		distance = 1.0; // maximum possible

	}

	return distance;
}

bool compareGenTruthTables(const GeneralizedTruthTable & g1, const GeneralizedTruthTable & g2)
{
	bool identical = true;

	if( sameChildParentsDelays(g1, g2) ) {

		for(size_t k=0; k < g1.getTruthValues().size(); k++)
		{
			if( (g1.getTruthValues())[k] != (g2.getTruthValues())[k] )
			{
				identical = false;
				break;
			}
		}

	} else {

		identical = false;

	}

	return identical;
}





//Added by Haizhou Wang, March 22, 2012
//Return 'True' if parent at position 'index' is fictitious
const bool GeneralizedTruthTable::isFictitious( const int index ) const
{
	vector<bool> checked(m_truthValues.size(), false);
	int rowsToCompare; //indicates number of rows with the same pattern
	
	//Compute the interval between two rows with the same pattern
    int interval = accumulate( m_parentBases.cbegin()+index+1, m_parentBases.cend(), 1, std::multiplies<int>() );

	//Traverse 'm_truthValues' vector, compare rows with the same pattern
	for(int i=0; i<(int)m_truthValues.size(); ++i)
	{
		if( !checked[i] )
		{
			checked[i] = true;
			rowsToCompare = m_parentBases[index] - 1;
			int location = i;
			do
			{	
				if( (location + interval)>(int)m_truthValues.size() ){
					cerr<<"Error!\nGeneralizedTruthTable::isFictitious(): Out of m_truthValues's bound!";
                    exit(-1);
				}

                if( m_truthValues[location] != -1 && 
                    m_truthValues[location+interval] != -1 &&                        
                    m_truthValues[location] != m_truthValues[location+interval]
                  )
                    return false; //This parent is NOT fictitious
				checked[location + interval] = true;
				--rowsToCompare;
				location = location + interval;
			}while(rowsToCompare>0);
		}
	}
	return true;
}


//Retrun parent index in 'm_nameParents'
//Return '-1' if the given name is not found in 'm_nameParents'
//NOTE: the return value is INDEX, not the actual ID of the parent. 
//(Acutal ID is: m_nameParents[RETURN])
int GeneralizedTruthTable::loopUpParentIndexByName(const string & name) const
{
    for(int t=0; t<(int)m_nameParents.size(); ++t)
    {
        if( m_nameParents[t] == name )
        {
            if( find( m_nameParents.begin()+t+1, m_nameParents.end(), name ) != m_nameParents.end() )
            {
                cerr<<"Error!\nGeneralizedTruthTable::loopUpParentIndexByName(): When eliminating fictitious parents for node "<<m_nameChild<<", there are more than one node has name \""<<name<<"\"\n"<<endl;
                exit(-1);
            }
            return t;
        }
    }    

    return -1;
}



//m_parentBases stores base for each parent
//m_parentBases[0] for parent1, m_parentBases[1] for parent2 and so on.
/*
Oct 30, 2012
Now this function supports eliminating fictitious parent(s) from incomplete truth tables.
NOTE that unknown values are represented by -1.
This is to say, a node's "normal" value can NOT be -1.

A parent node is assumed fictitious by default.
For known entries, if no child change is observed when parent changed, the parent will be eliminated.
*/
void GeneralizedTruthTable::getRidOfFictitious( const vector<int> & mustIncludeParentList )
{
    if( m_truthValues.size() == 0 || m_parentBases.size() == 0 || m_parents.size() == 0)
		return;

	//Check if 'm_truthValues' is a full truth table
    size_t sum = accumulate( m_parentBases.cbegin(), m_parentBases.cend(), 1, std::multiplies<int>() );
	
	if( m_truthValues.size() != sum )
	{
		cerr<<"Error: truth values vector of node "<<m_child<<" is NOT a full list according to the information given in m_parentBases (Lengths don't match)\n";
		exit( EXIT_FAILURE );
	}

	//Make histogram of 'm_truthValues', check to see if all parents are fictitious
    vector<unsigned> truthValuesHistogram( m_base, 0 );
    for( const auto & i: m_truthValues )
    {
        if( i >= 0 ) //Child value may be '-1', if incomplete 'm_truthValues' is given.
        {
            ++truthValuesHistogram[i];
        }
    }
		

    //If any entry in 'truthValuesHistogram' equals to the size of 'm_truthValues',
    //then all parents are fictitious.
	//'m_truthValues' will be resized to 1.
	//And this child will have its own base in m_parentBases (But m_delays will be empty)
	for(size_t i=0; i<truthValuesHistogram.size(); ++i)
	{	
        if( truthValuesHistogram[i] == m_truthValues.size() )
		{
			int finalValue = m_truthValues[0];
            assert( finalValue >= 0 );

			m_truthValues.resize(1);
			m_truthValues[0] = finalValue;

			m_transProbs.resize(1);
			m_transProbs[0] = vector<double>(m_base, 0);
			m_transProbs[0][ finalValue ] = 1.0;

            m_parents.clear();
            m_nameParents.clear();
            m_delays.clear();

			m_parentBases.resize(1);
			m_parentBases.push_back(m_base);
			
			return;
		}
     }



    
    //Check fictitious-ness of each parent (All parents are assumed fictitious by default.)
    //By shffuling the order of parents to be check
    //different result may vary from run to run
    //For example, assume parent A and C are both ficitous.
    //However, eliminating either of them will make the other 'un'-fictious. 
    //(All parents are assumed ficitous by default)
    //In this case, two different result will be produced.
    vector<string> parentsToBeChecked( m_nameParents );
    random_shuffle( parentsToBeChecked.begin(), parentsToBeChecked.end(), RNG);  //RNG: Random number generater. Defined in the beginning of this file.

    for(const auto & name: parentsToBeChecked)
	{
        int i = loopUpParentIndexByName( name ); //It's parent's index, actuall parent ID is m_parent[i]
        if( i == -1 )
        {
            cerr<<"GeneralizedTruthTable::getRidOfFictitious(): Error!\nWhen eliminating node "<<m_nameChild<<"'s fictitious parents, can not find parent "<<name<<" in child's parent list\n"<<endl;
            exit(-1);
        }


        //If a node is in the must-include parent list, keep it anyway.
        if( find( mustIncludeParentList.cbegin(), mustIncludeParentList.cend(), m_parents[i] ) != mustIncludeParentList.cend() )
            continue;

		if( isFictitious(i) )
		{	
            //Parent i is fictitious, delete this parent and shrink 'm_truthValues'
			vector<int> oldTruthValue (m_truthValues);
			vector<bool> checked(m_truthValues.size(), false);
			int rowsToCompare;
			
			//Computer the distance between two rows with the same pattern
            int interval = accumulate( m_parentBases.cbegin()+i+1, m_parentBases.cend(), 1, std::multiplies<int>() );


			m_truthValues.resize( oldTruthValue.size()/m_parentBases[i] );
			m_transProbs.resize (  m_transProbs.size()/m_parentBases[i] );
			
			//Erase the fictitious parent
			for(int j=0, k=0; k<(int)oldTruthValue.size(); ++k)
			{
				if( !checked[k])
				{
					checked[k] = true;
					rowsToCompare = m_parentBases[i] - 1;
					int location = k;
                    int valueToBeAssigned = oldTruthValue[k];
					do
					{
						if( (location + interval)>(int)oldTruthValue.size() )
						{
							cerr<<"GeneralizedTruthTable::getRidOfFictitious(): Error!\nWhen eliminating fictitious parent, index out of bound\n";
							exit( EXIT_FAILURE );
						}

                        if( oldTruthValue[location] != -1 && oldTruthValue[location+interval] != -1 )
                            assert( oldTruthValue[location] == oldTruthValue[location+interval] );

                        if( oldTruthValue[location+interval] != -1 )
                            valueToBeAssigned = oldTruthValue[location+interval];
						
						checked[location + interval] = true;
						--rowsToCompare;
						location = location + interval;
					}while( rowsToCompare > 0 );

                    m_truthValues[j++] = valueToBeAssigned;
				}
			}
			
			m_parentBases.erase( m_parentBases.begin() + i  );
			m_parents.erase    ( m_parents.begin()     + i  );
			m_nameParents.erase( m_nameParents.begin() + i  );
			m_delays.erase     ( m_delays.begin()      + i  );
		}//Processed to next parent
	}//Finished checking all parents


    //After eliminating fictitious parents, if still unknown (-1) cells exist, replace them with 0
    replace( m_truthValues.begin(), m_truthValues.end(), -1, 0);

    //And set the corresponding possibility transition table
    for(size_t x = 0; x < m_transProbs.size(); ++x) 
	{
		for(int y = 0; y < m_base; ++y) 
		{
			m_transProbs[x][y] = 0;
		}

		m_transProbs[x][ m_truthValues[x] ] = 1.0;
	}

}



//Added by Haizhou Wang, Oct 25, 2012
bool GeneralizedTruthTable::operator==( const GeneralizedTruthTable & gtt ) const
{
    if( m_transProbs != gtt.m_transProbs ||
        m_truthValues != gtt.m_truthValues ||
        m_child != gtt.m_child ||
        m_nameChild != gtt.m_nameChild ||
        m_base != gtt.m_base ||
        m_parents.size() != gtt.m_parents.size() ||
        m_nameParents.size() != gtt.m_nameParents.size() ||
        m_parentBases.size() != gtt.m_parentBases.size() ||
        m_delays.size() != gtt.m_delays.size()
        )
        return false;

    for(unsigned i=0; i<m_parents.size(); ++i)
    {
        if( find( m_parents.cbegin(), m_parents.cend(), gtt.m_parents[i] ) == m_parents.cend() || find( gtt.m_parents.cbegin(), gtt.m_parents.cend(), m_parents[i] ) == gtt.m_parents.cend()) 
            return false;
    }


    for(unsigned i=0; i<m_nameParents.size(); ++i)
    {
        if( find( m_nameParents.cbegin(), m_nameParents.cend(), gtt.m_nameParents[i] ) == m_nameParents.cend() || find( gtt.m_nameParents.cbegin(), gtt.m_nameParents.cend(), m_nameParents[i] ) == gtt.m_nameParents.cend()) 
            return false;
    }


    for(unsigned i=0; i<m_parentBases.size(); ++i)
    {
        if( find( m_parentBases.cbegin(), m_parentBases.cend(), gtt.m_parentBases[i] ) == m_parentBases.cend() || find( gtt.m_parentBases.cbegin(), gtt.m_parentBases.cend(), m_parentBases[i] ) == gtt.m_parentBases.cend()) 
            return false;
    }


    for(unsigned i=0; i<m_delays.size(); ++i)
    {
        if( find( m_delays.cbegin(), m_delays.cend(), gtt.m_delays[i] ) == m_delays.cend() || find( gtt.m_delays.cbegin(), gtt.m_delays.cend(), m_delays[i] ) == gtt.m_delays.cend()) 
            return false;
    }


    return true;

}
//End of Haizhou's code