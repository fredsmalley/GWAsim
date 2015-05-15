// Transition.h
//
// Joe Song
// Created: November 30, 2008

#pragma once

#include <vector>
using std::vector;
#include <string>
using std::string;

class Transition
{
public:

    Transition(void) : m_child(0), m_base(0) { }
  
	void initialize(const int child, const int base, const vector< int > & parents, 
                    const vector<int> & parentBases, const vector< int > & delays, 
                    const string & nameChild="", 
                    const vector<string> & nameParents=vector<string>(0));

	virtual void allocate() = 0;
    
	void setChild(int child, int base, const string & name="")
    { m_child = child; m_base = base; m_nameChild=name; }

    void setChildName(const string & name) { m_nameChild=name; }

    void setParents(const vector<int> & parents, const vector<string> & names = vector<string>(0));
    
    void setParentNames(const vector<string> & names) { m_nameParents = names; }

	void setParentBases(const vector<int> & parentBases) { m_parentBases = parentBases; } 

	void setDelays(const vector<int> & delays) { m_delays = delays; }

	int getChild() const { return m_child; }

    //Added by Haizhou Wang, Oct 11, 2012
    const string getChildName() const {return m_nameChild;}
    //End of Haizhou's Code

	int getBase() const { return m_base; }

	const vector<int> & getParents() const { return m_parents;}

	const vector<int> & getParentBases() const { return m_parentBases; }

	const vector<int> & getDelays() const { return m_delays; }

    
protected:
	
    size_t calculateSize() const;


protected:

	int m_child; // child ID, zero based
    string m_nameChild; // child name
    
	int m_base;

	vector<int> m_parents; // parent IDs, zero based
    vector<string> m_nameParents; // parent names
    
	vector<int> m_parentBases;

	vector<int> m_delays; // delays of the effect of each parent on the child
    
};


inline long long unsigned ArrayIndexToLinearIndex(const vector<int> & bases,
                                                  const vector<int> & ai)
{
    // Converts a multi-dimensional array index to a linear index.
    // ai[0] is the LEAST significant dimension
    
	long long unsigned li=0; // linear index
	for(int j=bases.size()-1; j > 0; --j){
		li = (li + ai[j]) * bases[j-1];
	} //end for
	li += ai[0]; // add the last digit in the array
	return li;
}

inline void LinearIndexToArrayIndex(const vector<int> & bases,
                                    const long long unsigned li,
                                    vector<int> & ai)
{
    // Converts a linear index to a multi-dimensional array index.
    // ai[0] is the LEAST significant dimension
    
    long long unsigned current = li;
    for(size_t j = 0; j < ai.size(); ++j){
        ai[j] = current % bases[j];
        current = current / bases[j];
    }//end for
    
}

inline long long unsigned ArrayIndexToLinearIndexMSB0(const vector<int> & bases,
                                                      const vector<int> & ai)
{
    // Converts a multi-dimensional array index to a linear index.
    // ai[0] is the MOST significant dimension
    
	long long unsigned li=0; // linear index
	for(int j=0; j < (int) bases.size()-1; j++){
		li = (li + ai[j]) * bases[j+1];
	} //end for
	li += ai.back(); // add the last digit in the array
	return li;
}

inline void LinearIndexToArrayIndexMSB0(const vector<int> & bases,
                                        const long long unsigned li,
                                        vector<int> & ai)
{
    // Converts a linear index to a multi-dimensional array index.
    // ai[0] is the MOST significant dimension
    
    long long unsigned current = li;
    for(int j = ai.size() - 1; j >=0; j--){
        ai[j] = current % bases[j];
        current = current / bases[j];
    }//end for
    
}

inline int digitsToNumber(const vector<int> & digits, const vector<int> & bases)
{ //digits[0] is the LEAST significant dimension
    
	int numDigits = (int) digits.size();
    
	int x = digits[numDigits-1];
	int B = bases[numDigits-1];
    
	for(int i=numDigits-2; i>=0; i--) {
		x += digits[i] * B ;
		B *= bases[i];
	}
	return x;
}

inline void numberToDigits(const int x, const vector<int> & bases,
                           vector<int> & digits)
{ //digits[0] is the LEAST significant dimension
    int current = x;
    for(size_t j = 0; j < digits.size(); j++) {
        digits[j] = current % bases[j];
        current = current / bases[j];
    }//end for
}

