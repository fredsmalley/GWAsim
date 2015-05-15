//
//  ExactMultiTableTest.cpp
//  gln
//
//  Created by Joe Song on 2/10/13.
//
//

#include "ExactMultiTableTest.h"
#include "ChisqTests.h"

void EMTEnumerator::limits(size_t k, size_t i, size_t j, const EMTEvaluator & v,
                          int & lij, int & uij)
{
    // lower and upper bounds for Ak[i,j]
    int bri = v.m_requiredRowSums[k][i] - ARowsums[k][i];
    int bcj = v.m_requiredColSums[k][j] - AColsums[k][j];
    lij = 0;
    if(i == As[k].nrow() - 1) { // last row
        lij = bcj;
    }
    if(j == As[k].ncol() - 1) { // last column
        lij = max(lij, bri);
    }
    uij = min(bcj, bri);
}

double
EMTEnumerator::traverse(size_t k, size_t i, size_t j,
                       EMTEvaluator & v,
                       const vector< TransitionTable > & Cs)
{
    double P;

    // Base case:
    if(k == As.size()) {
        // k being out of bound indciates (k,i,j) is beyond a leaf
        P = v.evaluate(*this, Cs);
        return P;
    } else {
        P = 0;
    }
    
    // branch-and-bound: decide whether to skip the entire branch starting
    //   at and including Ak[i,j]
    BOUND_CHECK check = v.bound(k, i, j, *this, Cs);
    if( check == TO_SKIP_ENTIRE_BRANCH ) {
        return 0;
    } else if( check == TO_KEEP_ENTIRE_BRANCH ) {
        P = v.add(k, i, j, *this, Cs);
        return P;
    }
    
    size_t k_next=k, i_next=i, j_next=j;
    next(k_next, i_next, j_next);

    if (0) {
        cout << "k=" << k << ", i=" << i << ", j=" << j
        << "; next k=" << k_next << ", i=" << i_next << ", j=" << j_next
        << endl;
    }

    if(As[k].nrow() == 0 || As[k].ncol() == 0) {

        // As[k] is an empty matrix.  Move on to next one
        update(k_next, i_next, j_next, k, i, j, v, Cs);
        P += traverse(k_next, i_next, j_next, v, Cs);
        
    } else {
        // lower and uppler limits for Ak[i,j]
        int lij, uij;
        limits(k, i, j, v, lij, uij);
        
        for( int a=lij; a <= uij; a++) { // will enter loop if and only if lij<=uij
            
            As[k].set(i, j, a); // As[k][i][j] = a;
            
            ARowsums[k][i] += a;
            AColsums[k][j] += a;
            
            //size_t k_next=k, i_next=i, j_next=j;
            //next(k_next, i_next, j_next);
            
            update(k_next, i_next, j_next, k, i, j, v, Cs);
            
            // cout << "k=" << k << ", i=" << i << ", j=" << j << ", a=" << a << endl;
            
            P += traverse(k_next, i_next, j_next, v, Cs);
            
            ARowsums[k][i] -= a;
            AColsums[k][j] -= a;
            
            As[k].set(i, j, 0); // As[k][i][j] = 0; // clear up the entry to return
            
        }
    }
    
    return P;
}

void EMTEnumerator::next(size_t & k, size_t & i, size_t & j)
{
    // k: contingency table index from 0 to K-1 (#matrices)
    // i: row index from 0 to I-1 (I: #rows)
    // j: column index from 0 to J-1 (J: #columns)
    
    if(As[k].nrow() == 0 || As[k].ncol() == 0) {
        k = k+1;
        i = 0;
        j = 0;
    } else if( j < As[k].ncol()-1 ) {  // # columns in Ak
        j = j+1;
    } else if ( i < As[k].nrow()-1 ) { // #rows in Ak
        j = 0;
        i = i+1;
    } else {
        k = k+1;
        i = 0;
        j = 0;
    }
}

void EMTEnumerator::update(size_t k_next, size_t i_next, size_t j_next,
                          size_t k, size_t i, size_t j,
                          EMTEvaluator & v,
                          const vector<TransitionTable> & Cs)
{
    if(k_next != k) {
        // gather statistics on table k
        v.processTable(k, *this, Cs);
    }
}

void EMTEnumerator::initialize(const EMTEvaluator & v,
                              const vector<TransitionTable> & Cs)
{
    As = v.generateTables(Cs);

    ARowsums.resize(As.size());
    AColsums.resize(As.size());

    for(size_t k=0; k<As.size(); k++) {
        ARowsums[k].resize(As[k].nrow());
        AColsums[k].resize(As[k].ncol());
    }
}

////////////////////////////////////////////////////////////////////////////////

BOUND_CHECK EMTEvaluator::bound(size_t k, size_t i, size_t j, const EMTEnumerator & e,
                          const vector<TransitionTable> & Cs)
// (virtual default)
{
    return NOT_TO_SKIP;
}

double EMTEvaluator::add(size_t k, size_t i, size_t j, const EMTEnumerator & e,
                        const vector<TransitionTable> & Cs)
// (virtual default)
{
    return 0;
}

void EMTEvaluator::processTable(size_t k, const EMTEnumerator & e,
                                 const vector<TransitionTable> & Cs)
{
    return;
}

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
double exact_multi_table_test(const vector< TransitionTable > & Cs,
                              EMTEvaluator & v, EMTEnumerator & e)
{
    v.initialize(Cs);
    
    e.initialize(v, Cs);
    
    return e.traverse(0, 0, 0, v, Cs);
}

