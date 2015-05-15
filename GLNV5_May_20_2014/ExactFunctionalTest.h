//
//  EMTFisherTests.h
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//

#ifndef __gln__EMTFisherTests__
#define __gln__EMTFisherTests__

//#include <iostream>
#include <algorithm>
#include "ExactMultiTableTest.h"

double marginalDistribution(const vector<TransitionTable> & Cs,  int index,  vector<vector<int> > & marginal, string mode);//index:0->row, 1->column//Hua added, May 19 2014

double FisherProb(const vector<vector<int> > & A);

//Add by Hua Apr 14 2014
class EMTFunctionalChisq : public EMTEvaluator
{
public:
    
    void initialize_customized_row_col_sum(const vector<TransitionTable> & Cs);
    
    virtual void initialize(const vector<TransitionTable> & Cs);
    
    virtual void processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs);
    
    virtual bool isMoreExtreme() const;
    
    double evaluate(const EMTEnumerator & e, const vector<TransitionTable> & Cs);
    
	vector<TransitionTable> generateTables(const vector<TransitionTable> & Cs) const;


    //Added by Hua Mar 20 2014
    virtual BOUND_CHECK bound(size_t k, size_t i, size_t j, const EMTEnumerator & e, const vector<TransitionTable> & Cs);
    //////
    //Hua added, Apr 21 2014
    void setRowAndColSum(vector<vector<int> > &row, vector<vector<int> > &col, const vector<TransitionTable> & Cs);
    
    //bool getExtremeness() const {return m_moreOrLessExtreme;};
    ////
    
    virtual double add(size_t k, size_t i, size_t j,
               const EMTEnumerator & e,
               const vector<TransitionTable> & Cs);
    
protected:
    vector<double> m_observedChisq; // observed statistics
    double m_totalObservedChisq;
    vector<double> m_colSumChisq; //chisq of columTotles of the tables
    vector<double> m_nullChisq; // null statistics
    vector<double> m_boundChisqs; // a bound on chisq under each condition
    //    double m_boundHeteroChisq;
    //    int m_row;
    bool m_skip;
//bool m_moreOrLessExtreme;//Hua added, May 17 2014. Test observed table's FunChisq, if <=0.5, false; else true. If two observed tables, 0.1 and 0.7, compare 0.1 and 0.3, so choose the 0.1, choose false; if 0.3 and 0.9, compare 0.3 and 0.1, choose 0.9, choose true.
};


////
//Add by Hua Apr 14 2014
double multi_table_functional_test(const vector<TransitionTable> & Cs,
                                   const string & discrepancy_measure);

double exact_functional_test(const TransitionTable & C,
                             const string & discrepancy_measure);

double exact_functional_test(const vector<vector<int > > & C,
                             const string & discrepancy_measure);
////
////

#endif /* defined(__gln__EMTFisherTests__) */
