// GeneralizedTruthTable.h -- Generalized truth table class 
//
// Joe Song
// Created: November 30, 2008
// Last modified: November 11, 2009.  Added the ProbGenTruthTable class

#pragma once

#include "Transition.h"
#include "TransitionTable.h"

class GeneralizedTruthTable :
	public Transition
{
public:

	friend double distanceTransProbTables(const GeneralizedTruthTable & g1, const GeneralizedTruthTable & g2);
	friend bool compareGenTruthTables(const GeneralizedTruthTable & g1, const GeneralizedTruthTable & g2);
	friend bool sameChildParentsDelays(const GeneralizedTruthTable & g1, const GeneralizedTruthTable & g2);

	GeneralizedTruthTable(void);

	GeneralizedTruthTable(const int child, const int base, const vector< int > & parents, 
		const vector<int> & parentBases, const vector< int > & delays, const vector<int> & gtt);

	GeneralizedTruthTable(const int child, const int base, const vector< int > & parents, 
		const vector<int> & parentBases, const vector< int > & delays, const vector<int> & gtt,
		const vector< vector<double> > & pgtt)
	{
		initialize(child, base, parents, parentBases, delays);
		m_truthValues = gtt;
		setTransProbTable(pgtt);
	}

	~GeneralizedTruthTable(void);

	void allocate();

	int lookup(const vector<int> & parentValues, int previousValue) const;
	
	void convert(const TransitionTable & transTable, double alpha);

	size_t size() const { return m_truthValues.size(); }

	const vector< int > & getTruthValues() const { return m_truthValues; }
	
	void setTruthValues( const vector< int > & truthvalues ) { m_truthValues = truthvalues ;}

	const vector< vector<double> > & getTransProbTable () const { return m_transProbs; }
	void setTransProbTable ( const vector< vector<double> > & transProbs) { m_transProbs = transProbs; }


	//Added by Haizhou Wang, March 22, 2012
	void getRidOfFictitious(const vector<int> & mustIncludeParentList = vector<int>() );
	const bool isFictitious( const int index) const;
	//End of Haizhou's code


    //Added by Haizhou Wang, Oct 25, 2012
    bool operator==( const GeneralizedTruthTable & gtt ) const;
    //End of Haizhou's code


private:
    int loopUpParentIndexByName(const string & name) const;

protected:
	vector< int > m_truthValues; // This is a vector holding the truth table values for each combination of parent values.

	vector< vector<double> > m_transProbs; // This is a matrix holding the transition probabilities of the child given parent values
};

/*
class ProbGenTruthTable : 
	public GeneralizedTruthTable
{
public:

	ProbGenTruthTable(const int child, const int base, const vector< int > & parents, 
		const vector<int> & parentBases, const vector< int > & delays, const vector<int> & gtt)
	{
		GeneralizedTruthTable(child, base, parents, parentBases, delays, gtt);
	}

	ProbGenTruthTable(const int child, const int base, const vector< int > & parents, 
		const vector<int> & parentBases, const vector< int > & delays, const vector<int> & gtt,
		const vector< vector<double> > & pgtt)
	{
		GeneralizedTruthTable(child, base, parents, parentBases, delays, gtt);
		setTransProbTable(pgtt);
	}

	void setTransProb( const vector<int> parentValues, int childValue, double prob ) 
	{	
		int row = digitsToNumber(parentValues, m_parentBases);
		m_transProbs[row][childValue] = prob; 
	}

	double getTransProb( const vector<int> parentValues, int childValue ) const
	{ 
		int row = digitsToNumber(parentValues, m_parentBases);
		return m_transProbs[row][childValue]; 
	}

	const vector< vector<double> > & getTransProbTable () const { return m_transProbs; }
	void setTransProbTable ( const vector< vector<double> > & transProbs) { m_transProbs = transProbs; }

protected:
	vector< vector<double> > m_transProbs;

};
*/
