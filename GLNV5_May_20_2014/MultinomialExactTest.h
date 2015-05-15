//
//  MultinomialExactTest.h
//  GLNV5_Apr_14_2014
//
//  Created by Hua Zhong on 4/17/14.
//  Copyright (c) 2014 ZH. All rights reserved.
//

#ifndef __GLNV5_Apr_14_2014__MultinomialExactTest__
#define __GLNV5_Apr_14_2014__MultinomialExactTest__

#include <iostream>
#include "ExactFunctionalTest.h"
#include "TransitionTable.h"
#include "ExactMultiTableTest.h"
#include <math.h>
class MultinomialExactTest;

class oneWayEnumerator : public EMTEnumerator
{
public:
    
    virtual void initialize(const EMTEvaluator & v, const vector<TransitionTable> & Cs){return;};
    
    void initialize(const MultinomialExactTest & v, const vector<TransitionTable> & Cs);
    
    double traverse(size_t k, size_t i, size_t j, MultinomialExactTest & v, const vector< TransitionTable > & Cs, const oneWayEnumerator & oe);
    
    virtual void next(size_t & k, size_t & i, size_t & j);
    
    virtual void limits(size_t k, size_t i, size_t j, const EMTEvaluator & v, int & lij, int & uij){return;};
    
    void limits(size_t k, size_t i, size_t j, const MultinomialExactTest & v, int & lij, int & uij);
    
    virtual void update(size_t k_next, size_t i_next, size_t j_next, size_t k, size_t i, size_t j, EMTEvaluator & v, const vector<TransitionTable> & Cs){return;};
    
    void setIndex(int index){m_index=index;};
    
public:
    //vector<TransitionTable> As;
    //vector<vector<int> > ARowsums;
    int m_index; //0 for doing row enumeration, 1 for doing column enumeration.
};

class MultinomialExactTest : public EMTEvaluator
{
    friend class oneWayEnumerator;
    
public:
    void initializeRequirement(const oneWayEnumerator e, const vector<TransitionTable> & Cs);
    
    virtual void initialize(const vector<TransitionTable> & Cs);
    
    virtual vector<TransitionTable> generateTables(const vector<TransitionTable> & Cs) const{
        vector<TransitionTable> As(Cs);
        
        for (size_t k=0; k<Cs.size(); k++) {
            As[k].reset();
        }
        
        return As;
    };
    
//    double evaluate(const oneWayEnumerator & e, const vector<TransitionTable> & Cs);
    
    double evaluate(const oneWayEnumerator & e, const vector<TransitionTable> & Cs, const oneWayEnumerator & oe);
    
    virtual double evaluate(const EMTEnumerator & e, const vector<TransitionTable> & Cs){return 0;};
    
//    virtual string bound(size_t k, size_t i, size_t j,
//                         const EMTEnumerator & e,
//                         const vector<TransitionTable> & Cs){return "not-to-skip";};
//    
//    virtual double add(size_t k, size_t i, size_t j,
//                       const EMTEnumerator & e,
//                       const vector<TransitionTable> & Cs){return 1;};
    
    virtual void processTable(size_t k, const EMTEnumerator & e, const vector<TransitionTable> & Cs){return;};
    
    virtual bool isMoreExtreme() const{return true;};
    
    vector<int> getNum() const {return m_num;};
    
protected:
    vector<int> m_totalSampleSize;

    vector<int> m_num;//for enumeration

};


//////

double multi_table_multinomial_exact_test(const vector<TransitionTable> & Cs,
                                          const string & discrepancy_measure);
double multinomial_exact_test(const TransitionTable & C,
                              const string & discrepancy_measure);
double multinomial_exact_test(const vector< vector<int> > & C,
                              const string & discrepancy_measure);

#endif /* defined(__GLNV5_Apr_14_2014__MultinomialExactTest__) */
