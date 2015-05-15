#pragma once

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <valarray>
using std::valarray;

//#include "DS_GLN.h"
#include "DS_ConstraintFile.h"
#include "DS_CompulsiveParents.h"

//#include "DS_RecognizeState.h"
//#include "DS_Feedback.h"

class DS_LocalConstraint : public DS_ConstraintFile
{
public:
	DS_LocalConstraint(void);
	~DS_LocalConstraint(void);
public:
	// initial value for simulation
	vector<int> m_initValue;

	// Local constraint checking
	vector<DS_CompulsiveParents> m_compulsiveParents;
        vector<string> m_ExcludedChild;  //Name of nodes who can NOT HAVE parents.
        vector<string> m_ExcludedParent; //Name of nodes who can NOT BE parent.

public:
	// read input file
	DS_LocalConstraint(const string fileName);
        void scan(const string & fileName);

	static valarray<int> convert(const vector<int> & vec);
   	static float correlation(const vector<int> & x, const vector<int> & y);
	static float correlation(const valarray<int> & x, const valarray<int> & y);
	static float signCorrelation(const vector<int> & x, const vector<int> & y);
	static float signCorrelation(const valarray<int> & x, const valarray<int> & y);
	static float expectedValue(const valarray<int> & v);

    //Added by Haizhou Wang, Oct 15, 2012
    //NOTE: The output local constraint file format can ONLY be used by kstop
    void saveLocalConstraint(const string output_file_name) const;
    //End of Haizhou's code


    void saveConstraints( const string outputFileName ) const;
    void saveLocalConstraintPair( const string output) const;


    //Added by Haizhou Wang, Jan 29, 2013
    //Output local constraints as topology file which can be used by GLN Reconstrucion.
    void saveAsTopologyFile(const string & output_file_name) const;
    //End of Haizhou's code


public:
	static void generateLocalConstraintTestFile();
	static void generateLocalConstraintTestGLNFile();


    //Added by Haizhou Wang, Jan 23, 2013
    //Return the number of parents of the node which has the most.
    int getMaxNumParents() const;
    bool empty() const { return m_compulsiveParents.empty(); }
};


