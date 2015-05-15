//
//  ExactCPChisq3.h
//  gln
//
//  Created by Joe Song on 2/11/13.
//
//

#ifndef __gln__ExactComparativeChisq__
#define __gln__ExactComparativeChisq__

#include <iostream>
#include "ExactMultiTableTest.h"
#include "ComparativeChisq.h"

class ProbHetero: public EMTEvaluator
{
public:
    ProbHetero(const TransitionTable & stackedProjectedAs,
               const TransitionTable & Cpooled,
               const vector<vector<int> > & uniqueParents,
               const vector<vector<int> > & uniqueParentBases,
               const vector<vector<int> > & uniqueParentDelays,
               double observedHeteroProb,
               double nullPooledProb)
    : m_stackedCommonAs(stackedProjectedAs),
    m_Cpooled(Cpooled),
    m_uniqueParents(uniqueParents),
    m_uniqueParentBases(uniqueParentBases),
    m_uniqueParentDelays(uniqueParentDelays),
    m_observedHeteroProb(observedHeteroProb),
    m_nullPooledProb(nullPooledProb)
    {  m_nullFisherProb.resize(uniqueParents.size()); }
    
    virtual void initialize(const vector< TransitionTable > & Cs);
    virtual void processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs);
    virtual double evaluate(const EMTEnumerator & e,
                            const vector< TransitionTable > & Cs);
    virtual bool isMoreExtreme() const;
    virtual double add(size_t k, size_t i, size_t j, const EMTEnumerator & e,
                       const vector<TransitionTable> & Cs);

    vector< TransitionTable > generateTables(const vector< TransitionTable > & Cs) const;
    
    double computeNullHeteroProb() const;
    
protected:
    
    const TransitionTable & m_stackedCommonAs; // A($\Pi0$): given projections to common parents A
        
    const TransitionTable & m_Cpooled;
    
    const vector<vector<int > > & m_uniqueParents; // for each condition
    const vector<vector<int > > & m_uniqueParentBases; // for each condition
    const vector<vector<int > > & m_uniqueParentDelays; // for each condition

    double m_observedHeteroProb;
    double m_nullPooledProb;

    double m_nullHeteroProb;
    vector<double> m_nullFisherProb;
};

class ProbHeteroChisq: public ProbHetero
{
public:

    ProbHeteroChisq(const TransitionTable & stackedProjectedAs,
                    double observedHeteroChisq, double nullHomoChisq,
                    const TransitionTable & Cpooled,
                    const vector<vector<int> > & uniqueParents,
                    const vector<vector<int> > & uniqueParentBases,
                    const vector<vector<int> > & uniqueParentDelays,
                    double observedHeteroProb,
                    double nullFisherProbPooled,
                    CMPX_MARGINAL marginal)
    :
    ProbHetero(stackedProjectedAs, Cpooled,
               uniqueParents, uniqueParentBases, uniqueParentDelays,
               observedHeteroProb, nullFisherProbPooled),
    m_observedHeteroChisq(observedHeteroChisq), m_nullHomoChisq(nullHomoChisq)
    {
        size_t K = stackedProjectedAs.nrow();
        m_chisqs.resize(K);
        m_boundChisqs.resize(K);
        m_cmpxMarginal = marginal;
    }

    virtual void initialize(const vector< TransitionTable > & Cs);
    
    virtual BOUND_CHECK bound(size_t k, size_t i, size_t j, const EMTEnumerator & e,
                         const vector<TransitionTable> & Cs);
    
    virtual void processTable(size_t k, const EMTEnumerator & e, const vector<TransitionTable> & Cs);
    bool isMoreExtreme() const;

    TransitionTable backProject(size_t k, const EMTEnumerator & e,
                                const vector< TransitionTable > & Cs);
        
protected:

    const double m_observedHeteroChisq; // $\chi^2d$
    const double m_nullHomoChisq; // $\chi^2c from As$

    vector<double> m_chisqs;  // save chisquares of individual tables
    vector<double> m_boundChisqs; // a bound on chisq under each condition
    
    CMPX_MARGINAL m_cmpxMarginal;
    
    bool m_skip;
};

////////////////////////////////////////////////////////////////////////////////

class ProbHeteroGivenPooled: public EMTEvaluator
{
public:
    
    virtual void initialize(const vector< TransitionTable > & Cs);
    
    vector< TransitionTable >
    generateTables(const vector< TransitionTable > & Cs) const;
    
    double evaluate(const EMTEnumerator & e, const vector< TransitionTable > & Cs);
    
    bool isMoreExtreme() const { return true; }
    
protected:

    double m_observedHeteroProb; // observed $\chi^2d$
    
    TransitionTable m_Cpooled;
    
    vector<vector<int > > m_uniqueParents; // for each condition
    vector<vector<int > > m_uniqueParentBases; // for each condition
    vector<vector<int > > m_uniqueParentDelays; // for each condition
    
    TransitionTable m_stackedCommonCs; // C($\Pi0$)
};

class ProbHeteroChisqGivenPooled: public ProbHeteroGivenPooled
{
public:
    
    virtual void initialize(const vector< TransitionTable > & Cs);
    
    void setCmpxMarginal(const CMPX_MARGINAL & marginal)
    { m_cmpxMarginal = marginal; }
    
    virtual double evaluate(const EMTEnumerator & e, const vector< TransitionTable > & Cs)
    {
        // ProbHetero v(A($\Pi0$), observed $\chi^2d$, common parents, unique parents);
        ProbHeteroChisq vSub(e.As[0], m_observedHeteroChisq, m_observedHomoChisq,
                             m_Cpooled, m_uniqueParents, m_uniqueParentBases,
                             m_uniqueParentDelays, m_observedHeteroProb,
                             FisherProb(e.As[0].getTransitionTable()), m_cmpxMarginal);
        EMTEnumerator eSub;
        
        return exact_multi_table_test(Cs, vSub, eSub);
    }
    
protected:
    
    double m_observedHeteroChisq; // observed $\chi^2d$
    double m_observedHomoChisq; // observed $\chi^2c$
    
    CMPX_MARGINAL m_cmpxMarginal;
};

#endif /* defined(__gln__ExactComparativeChisq__) */
