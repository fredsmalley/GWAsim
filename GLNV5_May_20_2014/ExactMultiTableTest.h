//
//  ExactMultiTableTest.h
//  gln
//
//  Created by Joe Song on 2/10/13.
//
//

#ifndef __gln__ExactMultiTableTest__
#define __gln__ExactMultiTableTest__

#include <vector>
#include <iostream>

using namespace std;

#include "TransitionTable.h"
#include "ChisqTests.h"

//Hua added, Mar 9, 2015
enum BOUND_CHECK {
    NOT_TO_SKIP,
    TO_SKIP_ENTIRE_BRANCH,
    TO_KEEP_ENTIRE_BRANCH
};
////

class EMTEvaluator;

class EMTEnumerator
{
public:
    virtual void initialize(const EMTEvaluator & v,
                            const vector<TransitionTable> & Cs);

    virtual double traverse(size_t k, size_t i, size_t j,
                    EMTEvaluator & v, const vector< TransitionTable > & Cs);
    virtual void next(size_t & k, size_t & i, size_t & j);
    virtual void limits(size_t k, size_t i, size_t j, const EMTEvaluator & v,
                        int & lij, int & uij);

    virtual void update(size_t k_next, size_t i_next, size_t j_next,
                        size_t k, size_t i, size_t j,
                        EMTEvaluator & v,
                        const vector<TransitionTable> & Cs);
    
public:
    vector<TransitionTable> As;
    vector<vector<int> > ARowsums;
    vector<vector<int> > AColsums;
};

class EMTEvaluator
{
    friend class EMTEnumerator;
    
public:
    virtual void initialize(const vector<TransitionTable> & Cs)=0;
    virtual vector<TransitionTable>generateTables(const vector<TransitionTable> & Cs) const = 0;

    virtual double evaluate(const EMTEnumerator & e,
                            const vector<TransitionTable> & Cs)=0;

    virtual BOUND_CHECK bound(size_t k, size_t i, size_t j,
                         const EMTEnumerator & e,
                         const vector<TransitionTable> & Cs);
    virtual double add(size_t k, size_t i, size_t j,
                       const EMTEnumerator & e,
                       const vector<TransitionTable> & Cs);

    virtual void processTable(size_t k, const EMTEnumerator & e,
                         const vector<TransitionTable> & Cs);
    
    virtual bool isMoreExtreme() const = 0;
    
protected:
    // required row sums
    vector<vector<int> > m_requiredRowSums;
    
    // required col sums
    vector<vector<int> > m_requiredColSums;
    
};

double exact_multi_table_test(const vector< TransitionTable > & Cs,
                              EMTEvaluator & v, EMTEnumerator & e);

#endif /* defined(__gln__ExactMultiTableTest__) */
