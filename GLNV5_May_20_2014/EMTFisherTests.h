//
//  EMTFisherTests.h
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//

#ifndef __gln__EMTFisherTests__
#define __gln__EMTFisherTests__

#include <iostream>
#include <algorithm>
#include "ExactMultiTableTest.h"

double marginalDistribution(const vector<TransitionTable> & Cs,  int index,  vector<vector<int>> & marginal, string mode);//index:0->row, 1->column//Hua added, May 19 2014

class EMTFisherProb : public EMTEvaluator
{
public:
    virtual void initialize(const vector<TransitionTable> & Cs);
    virtual vector<TransitionTable>generateTables(const vector<TransitionTable> & Cs) const;
    
    virtual double evaluate(const EMTEnumerator & e,
                            const vector<TransitionTable> & Cs);
    
    virtual void processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs);
    
    virtual bool isMoreExtreme() const;
    
protected:
    
    vector<double> m_observedFisherProb; // observed statistics for each table
    vector<double> m_nullFisherProb; // null statistics for each table
};

class EMTFisherChisq : public EMTFisherProb
{
public:
    virtual void initialize(const vector<TransitionTable> & Cs)
    {
        EMTFisherProb::initialize(Cs);
        m_observedChisq.resize(Cs.size());
        for (size_t k=0; k<Cs.size(); k++) {
            m_observedChisq[k] = chisq(Cs[k].getTransitionTable());
        }
        m_nullChisq.resize(Cs.size());
    }
    
    virtual void processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs)
    {
        EMTFisherProb::processTable(k, e, Cs);
        m_nullChisq[k] = chisq(e.As[k].getTransitionTable());
    }
    
    virtual bool isMoreExtreme() const
    {
        double nullTotalChisq = 0, observedTotalChisq = 0;

        for(size_t k=0; k < m_nullFisherProb.size(); k++) {
            nullTotalChisq += m_nullChisq[k];
            observedTotalChisq += m_observedChisq[k];
        }
        return nullTotalChisq >= observedTotalChisq;

        /*
        for(size_t k=0; k < m_nullFisherProb.size(); k++) {
            if( m_nullChisq[k] >= m_observedChisq[k] ) {
                return true;
            }
        }
        return false;
        */
    }
    
protected:
    
    vector<double> m_observedChisq; // observed statistics
    vector<double> m_nullChisq; // null statistics
};


double multi_table_fisher_test(const vector<TransitionTable> & Cs,
                               const string & discrepancy_measure);

double fisher_exact_test(const TransitionTable & C,
                         const string & discrepancy_measure);

double fisher_exact_test(const vector<vector<int > > & C,
                         const string & discrepancy_measure);


#endif /* defined(__gln__EMTFisherTests__) */
