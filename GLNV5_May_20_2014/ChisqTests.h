//
//  ChisqTests.h
//  gln
//
//  Created by Joe Song on 1/24/13.
//
//

#ifndef gln_ChisqTests_h
#define gln_ChisqTests_h

#include <vector>
using namespace std;

class Chisq {
public:
    
    friend bool operator != (const Chisq & c1, const Chisq & c2);
    friend bool operator == (const Chisq & c1, const Chisq & c2);
    
    Chisq(): m_x(0), m_df(0), m_pval(1) {}
    
    Chisq(double x, size_t df, double pval): m_x(x), m_df(df), m_pval(pval) {};
    
    double m_x;
    size_t m_df;
    double m_pval;
};

double chisq(const vector<vector<int> > & O, const vector<vector<double> > & E);
double chisq(const vector<vector<int> > & O);

double ChisqTest(const vector< vector<int> > & truth_table, int base, int mode, double & chisq, size_t & df, const vector< vector<double> > & null_prob = vector< vector<double> >(0));

double ChisquareTest(const vector< vector<int> > & table_obs,
                     double & chisq, size_t & df,
                     const vector< vector<double> > & null_prob
                     = vector< vector<double> >(0) );

double ChisqDirTest(const vector< vector<int> > & table_obs,
                    double & chisquare, size_t & df, const string & method="");

void ChisquareTest1DNoPValue(const vector<int> & x_obs,
                             const vector<double> & p_null,
                             int K, double & chisq,
                             size_t & df);//Added by Hua Apr 14 2014

/*
double fisher_exact_test(const vector<vector<int > > & O,
                         const vector<vector<double > > & E = vector<vector<double > >(0));
*/

double exact_functional_test(const vector<vector<int > > & C,
                             const string & discrepancy_measure);//Add by Hua Apr 14 2014

//double multinomial_exact_test(const vector< vector<int> > & C,
//                              const string & discrepancy_measure);//Add by Hua Apr 18 2014

double fisher_exact_test(const vector<vector<int > > & C,
                         const string & discrepancy_measure);

double FisherProb(const vector<vector<int> > & A);

vector<vector<double> > getExpectedTable(const vector<vector<int> > & contingencytable);

vector<vector<double> > getExpectedTable(const vector<vector<int> > & contingencytable);

int getTotalSum(const vector<vector<int> > & tt);
vector<int> getColSums(const vector<vector<int> > & tt);
vector<int> getRowSums(const vector<vector<int> > & tt);

double chisq_upper_bound(const vector<int> & a, const vector<int> & b,
                         const string & method="");

void applyFunctionalChisqTest(const string & file, int pValMode);//Hua Apr 14 2014

#endif
