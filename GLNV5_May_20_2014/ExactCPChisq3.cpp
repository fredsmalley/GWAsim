//
//  ExactCPChisq3.cpp -- Exact comparative chisq analysis based on Fisher's
//    probability
//
//  gln
//
//  Created by Joe Song on 2/11/13.
//
//
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;

#include "SetOps.h"
#include "GLN.h"
#include "ChisqTests.h"
#include "ComparativeChisq.h"
#include "ExactCPChisq3.h"
#include "ExactMultiTableTest.h"
#include "EMTFisherTests.h"

vector<TransitionTable> project_to_unique_parents(const vector<TransitionTable> & tts,
                                                  const TransitionTable & pooled)
// Compute transition table projected to the space of unique parents
//   for each condition
// OUTPUT:
//   a list of transition tables projected from original tables to unique parents
// ASSUMPTIONS:
//   1. The same parent ID will have the same delay
//   2. Common parents show up in the same order in parent list for each
//      condition
{
    vector<TransitionTable> projected(tts.size());
    
    for(size_t i=0; i<tts.size(); i++) {
        projected[i] = tts[i].subtable(
                                       complementSet(pooled.getParents(),
                                                  tts[i].getParents()
                                                  )
                                       );
    }
    return projected;
}

vector<TransitionTable> project_to_unique_parents(const vector<TransitionTable> & tts)
{
    TransitionTable pooled = pool(tts);
    return project_to_unique_parents(tts, pooled);
}

void findParentIndices(const vector<int> & uniqueParents,
                       const vector<int> & pooledParents,
                       const vector<int> & parents,
                       vector<int> & uniqueParentIndices,
                       vector<int> & pooledParentIndices)
{
    for(size_t parentindex = 0; parentindex<parents.size(); parentindex++) {
        
        size_t pooledindex = find(pooledParents.cbegin(), pooledParents.cend(),
                                  parents[parentindex]) - pooledParents.cbegin();
        
        if(pooledindex < pooledParents.size()) {
            pooledParentIndices[pooledindex] = parentindex;
            continue;
        }
        
        size_t uniqueindex = find(uniqueParents.cbegin(), uniqueParents.cend(),
                                  parents[parentindex]) - uniqueParents.cbegin();
        
        if(uniqueindex < uniqueParents.size()) {
            uniqueParentIndices[uniqueindex] = parentindex;
            continue;
        }
    }
}

void assigneParentValues(const vector< int > & uniqueParentValues,
                         const vector< int > & commonParentValues,
                         const vector< int > & uniqueParentIndices,
                         const vector< int > & commonParentIndices,
                         vector< int > & parentValues)
{
    for(size_t i=0; i<uniqueParentValues.size(); i++) {
        parentValues[uniqueParentIndices[i]] = uniqueParentValues[i];
    }
    for(size_t i=0; i<commonParentValues.size(); i++) {
        parentValues[commonParentIndices[i]] = commonParentValues[i];
    }
}

TransitionTable forward_transform(TransitionTable s,
                                  const vector<int> & rowParents,
                                  const vector<int> & colParents)
// Transform table s to t as follows:
// s:
//   (rowParents, colParents) x child
// t:
//    rowParents x (colParents , child)
{
    if( colParents.size() == 0 )
        return s;
    
    vector<int> rowParentIndices(rowParents.size());
    vector<int> rowParentBases(rowParents.size());
    
    vector<int> colParentIndices(colParents.size());
    vector<int> colParentBases(colParents.size());
    
    findParentIndices(rowParents,colParents, s.getParents(),
                      rowParentIndices, colParentIndices);
    
    for(size_t i=0; i<rowParents.size(); i++) {
        rowParentBases[i] = s.getParentBases()[rowParentIndices[i]];
    }
    
    size_t prodColParentBases = 1;
    for(size_t i=0; i<colParents.size(); i++) {
        colParentBases[i] = s.getParentBases()[colParentIndices[i]];
        prodColParentBases *= colParentBases[i];
    }
    
    size_t t_ncol = prodColParentBases * s.getBase();
    
    TransitionTable t(0, t_ncol, rowParents, rowParentBases,
                      vector<int>(rowParents.size(), 0));
    
    vector<int> parentValues(s.getParents().size());
    
    for(size_t i=0; i < t.nrow(); i++) {
        
        vector<int> rowParentValues(rowParents.size());
        
        LinearIndexToArrayIndex(rowParentBases, i, rowParentValues);
        
        for(size_t j=0; j < t.ncol(); j++) {
            
            size_t liColParents = j / prodColParentBases; // linear index of column parents
            size_t c = j % prodColParentBases;
            
            vector<int> colParentValues(colParents.size());
            LinearIndexToArrayIndex(colParentBases, liColParents,
                                    colParentValues);
            
            // arrange unique and common parents in the order of parents in A
            assigneParentValues(rowParentValues, colParentValues,
                                rowParentIndices, colParentIndices,
                                parentValues);
            
            size_t r = ArrayIndexToLinearIndex(s.getParentBases(), parentValues);
            t.set(i, j, s.get(r, c));
        }
    }
    
    return t;
}

void backward_transform(const TransitionTable & projected,
                        TransitionTable & A,
                        const TransitionTable & pooled,
                        const vector< int > & uniqueParents,
                        const vector< int > & uniqueParentBases)
// Transform table s to t as follows:
// Input:
//   projected (rowParents, colParents) x child
// Output:
//   A: rowParents x (colParents , child)
{
    // Find for each row and column in projecte
    //   the corresponding values of parents and child in
    //   the original table
    vector<int> parentValues(A.getParents().size());
    
    vector<int> uniqueParentIndices(uniqueParents.size());
    vector<int> commonParentIndices(pooled.getParents().size());
    
    findParentIndices(A.getParents(), uniqueParents, pooled.getParents(),
                      uniqueParentIndices, commonParentIndices);

    vector<int> uniqueParentValues(uniqueParents.size());
    vector<int> commonParentValues(pooled.getParents().size());
    
    for(size_t i=0; i < projected.nrow(); i++) {
        
        LinearIndexToArrayIndex(uniqueParentBases, i, uniqueParentValues);
        
        for(size_t j=0; j < projected.ncol(); j++) {
            
            size_t liCommonParents = j / pooled.nrow(); // linear index of common parents
            size_t c = j % pooled.nrow();
            
            LinearIndexToArrayIndex(pooled.getParentBases(), liCommonParents,
                                    commonParentValues);
            
            // arrange unique and common parents in the order of parents in A
            assigneParentValues(uniqueParentValues, commonParentValues,
                                uniqueParentIndices, commonParentIndices,
                                parentValues);
            
            size_t r = ArrayIndexToLinearIndex(A.getParentBases(), parentValues);
            A.set(r, c, projected.get(i,j));
        }
    }
}

void backproject(const TransitionTable & stacked, size_t k, TransitionTable & A)
{
    // Find for each row and column in projecte
    //   the corresponding values of parents and child in
    //   the original table
    for(size_t j=0; j < stacked.ncol(); j++) {
        
        size_t liCommonParents = j / A.nrow(); // linear index of common parents
        size_t c = j % A.nrow();
        A.set(liCommonParents, c, stacked.get(k, j));
    }
}


TransitionTable projectAndStack(const vector< TransitionTable > & Cs,
                                const vector< int > & commonParents,
                                int ncols)
{
    // Create matrix C($\Pi0$) of size K $\times$ Q($\Pi0$)Q(X)
    TransitionTable target = TransitionTable(0, ncols,
                                             vector<int>(1, Cs.size()+1),
                                             vector<int> (1, Cs.size()),
                                             vector<int> (1, 0));
    
    for( size_t k=0; k<Cs.size(); k++) {
        // by projecting Ck in the space of common parents $\Pi0$
        // to row k of C($\Pi0$) for k=1...K
        TransitionTable CkProjected(Cs[k].subtable(commonParents));
        for(size_t r=0; r<CkProjected.nrow(); r++) {
            for(size_t c=0; c<CkProjected.ncol(); c++) {
                size_t l = r * CkProjected.ncol() + c;
                target.set(k, l, CkProjected.get(r,c));
            }
        }
    }
    return target;
}

////////////////////////////////////////////////////////////////////////////////

TransitionTable ProbHeteroChisq::backProject(size_t k, const EMTEnumerator & e,
                                             const vector< TransitionTable > & Cs)
{
    TransitionTable Dk(Cs[k]);
    
    if(m_uniqueParents[k].size() > 0 && m_Cpooled.nrow() > 0) {
        // having both unique and common parents
        // convert e.As[k].getTransitionTable() to the size of Ck
        backward_transform(e.As[k], Dk, m_Cpooled,
                           m_uniqueParents[k], m_uniqueParentBases[k]);
    } else if (m_Cpooled.nrow() > 0) {
        // having only common parents
        backproject(m_stackedCommonAs, k, Dk);
    } else if (m_uniqueParents[k].size() > 0 && m_Cpooled.nrow() == 0) {
        // having only unique parents
        Dk = e.As[k];
    } else {
        Dk.reset();
    }
    return Dk;
}

void ProbHeteroChisq::processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs)
{
    ProbHetero::processTable(k, e, Cs);
    TransitionTable tt = backProject(k, e, Cs);
    
    if (m_uniqueParents[k].size() == 0 && m_Cpooled.nrow() > 0) {
        
        switch (m_cmpxMarginal) {
            case INDEP_SUB:
                m_chisqs[k] = chisq( tt.getTransitionTable() );
                break;
            case INDEP_SUP:
                cerr << "ERROR: INDEP_SUP is not supported!" << endl;
                exit(EXIT_FAILURE);
                break;
            case POOLED:
            case POOLED_NO_ROWCOL:
            {
                // calculate expected count using the pooled row and col sum
                vector< vector<double> > E = getExpectedTable(m_Cpooled.getTransitionTable());
                
                size_t nk = tt.getTotalSum(), nPool = m_Cpooled.getTotalSum();
                
                for (size_t r=0; r<E.size(); r++) {
                    for (size_t c=0; c<E[r].size(); c++) {
                        E[r][c] = E[r][c] * nk / nPool;
                    }
                }
                
                m_chisqs[k] = chisq( tt.getTransitionTable(), E );
                break;
            }
                
            default:
                break;
        }
        
    } else {
        
        m_chisqs[k] = chisq( tt.getTransitionTable() );
        
    }
    // assert(m_chisqs[k] <= m_boundChisqs[k]);
}

bool ProbHeteroChisq::isMoreExtreme() const
{
    // Compute discrepancy measure
    double nullHeteroChisq = 0;
    
    for (size_t k=0; k<m_chisqs.size(); k++) {
        
        // nullHeteroChisq += chisq(k, e, Cs);
        nullHeteroChisq += m_chisqs[k];
    }
    
    nullHeteroChisq -= m_nullHomoChisq;
    
    assert(!(m_skip==true && abs(nullHeteroChisq) >= abs(m_observedHeteroChisq)));
    
    return abs(nullHeteroChisq) >= abs(m_observedHeteroChisq);
}

void ProbHeteroChisq::initialize(const vector< TransitionTable > & Cs)
{
    ProbHetero::initialize(Cs);
    for (size_t k=0; k<Cs.size(); k++) {
        
        if(m_uniqueParents[k].size() > 0 && m_Cpooled.nrow() > 0) {
            // having unique and common parents
            // m_boundChisqs[k] = Cs[k].getTotalSum() * (min(Cs[k].nrow(), Cs[k].ncol()) - 1);
            
            // INCORRECT:
            m_boundChisqs[k] = chisq_upper_bound(m_requiredColSums[k], m_requiredRowSums[k]);
            
        } else if ( m_Cpooled.nrow() > 0 ) {
            // having only common parents
            TransitionTable Dk(Cs[k]);
            backproject(m_stackedCommonAs, k, Dk);
            m_boundChisqs[k] = chisq(Dk.getTransitionTable());
            
        } else if ( m_uniqueParents[k].size() > 0 && m_Cpooled.nrow() == 0 ) {
            // having only unique parents
            m_boundChisqs[k] = chisq_upper_bound(Cs[k].getColSums(), Cs[k].getRowSums());
        } else {
            m_boundChisqs[k] = 0.0;
        }
    }
}

BOUND_CHECK ProbHeteroChisq::
bound(size_t k, size_t i, size_t j, const EMTEnumerator & e,
      const vector<TransitionTable> & Cs)
{
    BOUND_CHECK result = NOT_TO_SKIP; // "to keep entire branch"

    m_skip = false;

    return result;
    
    if( i == 0 && j == 0) {
        double boundHeteroChisq = 0;
        
        for(size_t m=0; m < k; m++) {
            boundHeteroChisq += m_chisqs[m]; 
        }
        
        if(boundHeteroChisq - m_nullHomoChisq >= abs(m_observedHeteroChisq) ) {
            result = TO_SKIP_ENTIRE_BRANCH;
            m_skip = true;
            
        } else {
        
            for(size_t m=k; m<m_boundChisqs.size(); m++) {
                boundHeteroChisq += m_boundChisqs[m];
            }
        
            boundHeteroChisq -= m_nullHomoChisq;
        
            result = (abs(boundHeteroChisq) < abs(m_observedHeteroChisq)) ?
                TO_SKIP_ENTIRE_BRANCH : NOT_TO_SKIP;
        }
    
        if(result == TO_SKIP_ENTIRE_BRANCH) {
            // cout << "Suppose to skip";
            m_skip = true;
            // result = "not-to-skip";
            
            /*
            cout << "skip: actual chisq="
                << m_boundChisqs[0] + m_boundChisqs[1] - m_nullHomoChisq
                << " bound=" << boundHeteroChisq
                << " observed=" << m_observedHeteroChisq
                << endl;
            */
        }
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////

void ProbHetero::initialize(const vector< TransitionTable > & Cs)
{
    m_requiredRowSums.resize(Cs.size());
    m_requiredColSums.resize(Cs.size());
    
    if(m_stackedCommonAs.ncol() == 0) {
        for (size_t k=0; k<Cs.size(); k++) {
            m_requiredRowSums[k] = Cs[k].getRowSums();
            m_requiredColSums[k] = Cs[k].getColSums();
        }
    } else {
        for (size_t k=0; k<Cs.size(); k++) {
            // Compute required row sums for Bk from Ck(${\Pi}$k-$\Pi0$)
            m_requiredRowSums[k] = Cs[k].subtable(m_uniqueParents[k]).getRowSums();
            // Get required column sums of Bk from row k of A($\Pi0$)
            m_requiredColSums[k] = m_stackedCommonAs.getTransitionTable()[k];
        }
    }
}

vector< TransitionTable >
ProbHetero::generateTables(const vector< TransitionTable > & Cs) const
{
    if(m_stackedCommonAs.ncol() == 0 ) {
        return Cs;
    } else {
        vector<TransitionTable> Bs(Cs.size());

        for (size_t k=0; k<Cs.size(); k++) {
            // Bk size: Q(${\Pi}$k-$\Pi0$) $\times$ Q($\Pi0$)Q(X)
            Bs[k]
            = TransitionTable(0, m_stackedCommonAs.ncol(),
                              m_uniqueParents[k],
                              m_uniqueParentBases[k],
                              m_uniqueParentDelays[k]);
        }
        return Bs;
    }
}

void ProbHetero::processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs)
{
    m_nullFisherProb[k] = FisherProb(e.As[k].getTransitionTable());
}

double ProbHetero::computeNullHeteroProb() const
{
    double P = m_nullPooledProb;
    for (size_t k=0; k<m_nullFisherProb.size(); k++) {
        P *= m_nullFisherProb[k];
    }
    return P;
}

bool ProbHetero::isMoreExtreme() const
{
    return m_nullHeteroProb <= m_observedHeteroProb;
}

double ProbHetero::evaluate(const EMTEnumerator & e,
                            const vector< TransitionTable > & Cs)
{
    double P = 0;

    m_nullHeteroProb = computeNullHeteroProb();
    
    if ( isMoreExtreme() ) {
        P = m_nullHeteroProb;
    }
    return P;
}

double ProbHetero::
add(size_t k, size_t i, size_t j, const EMTEnumerator & e,
    const vector<TransitionTable> & Cs)
{
    double P=0.0;
    if(i == 0 && j==0) { // This is only valid when i=j=0
        P = m_nullPooledProb;
        for (size_t m=0; m<k; m++) {
            P *= m_nullFisherProb[m];
        }
    }
    assert(P>0);
    return P;
}

////////////////////////////////////////////////////////////////////////////////

void ProbHeteroGivenPooled::initialize(const vector< TransitionTable > & Cs)
{
    // Identify common parents $\Pi0$ from C1...CK
    m_Cpooled = pool(Cs); // pooledandunique[0][0];

    m_uniqueParents.resize(Cs.size());
    m_uniqueParentBases.resize(Cs.size());
    m_uniqueParentDelays.resize(Cs.size());
        
    vector<TransitionTable> projected = project_to_unique_parents(Cs, m_Cpooled);
    
    for(size_t k=0; k<Cs.size(); k++) {
        if(projected[k].nrow() > 0) {
            // Identify unique parents for each condition
            m_uniqueParents[k] = projected[k].getParents();
            m_uniqueParentBases[k] = projected[k].getParentBases();
            m_uniqueParentDelays[k] = projected[k].getDelays();
        }
    }
    
    // Computed C($\Pi0$) by projecting C1...CK to common parents
    m_stackedCommonCs = projectAndStack(Cs, m_Cpooled.getParents(),
                                        m_Cpooled.nrow()*m_Cpooled.ncol());

    m_observedHeteroProb = FisherProb(m_stackedCommonCs.getTransitionTable());
    for (size_t k=0; k<Cs.size(); k++) {
        m_observedHeteroProb
            *= FisherProb(forward_transform(Cs[k], m_uniqueParents[k],
                                            m_Cpooled.getParents()).getTransitionTable());
    }

    // Compute required row sums
    m_requiredRowSums = vector<vector<int > >(1, m_stackedCommonCs.getRowSums());
    
    // Compute required column sums
    m_requiredColSums = vector<vector<int > >(1, m_stackedCommonCs.getColSums());
}

vector< TransitionTable >
ProbHeteroGivenPooled::generateTables(const vector< TransitionTable > & Cs) const
{
    vector< TransitionTable > As(1, m_stackedCommonCs); // C($\Pi0$)
    As[0].reset(); // reset A to zero
    return As;
}

double ProbHeteroGivenPooled::evaluate(const EMTEnumerator & e,
                                       const vector< TransitionTable > & Cs)
{
    // ProbHetero v(A($\Pi0$), observed $\chi^2d$, common parents, unique parents);
    ProbHetero vSub(e.As[0], m_Cpooled, m_uniqueParents, m_uniqueParentBases,
                    m_uniqueParentDelays, m_observedHeteroProb,
                    FisherProb(e.As[0].getTransitionTable()));
    
    EMTEnumerator eSub;
        
    return exact_multi_table_test(Cs, vSub, eSub);
}

////////////////////////////////////////////////////////////////////////////////
void ProbHeteroChisqGivenPooled::
initialize(const vector< TransitionTable > & Cs)
{
    ProbHeteroGivenPooled::initialize(Cs);

    if (m_Cpooled.nrow() > 0 && m_Cpooled.ncol() > 0) {
        m_Cpooled.applyChisquareTest(3);
    }

    vector< TransitionTable > tts = Cs;
    
    switch (m_cmpxMarginal) {
        case INDEP_SUB:
            computeIndividualChisq(tts, 3);
            break;
        case INDEP_SUP:
            cerr << "ERROR: INDEP_SUP is not supported!" << endl;
            exit(EXIT_FAILURE);
            break;
        case POOLED:
        case POOLED_NO_ROWCOL:
            m_observedHomoChisq = m_Cpooled.getChisq();
            computeIndividualChisq(tts, 3, m_Cpooled);
            break;
            
        default:
            break;
    }
    
    TransitionTable ttTot = computeTotalChisq(tts, 3);
    
    TransitionTable ttHet = computeHetChisq(tts, ttTot, m_Cpooled, 3, m_cmpxMarginal);
    
    m_observedHeteroChisq = ttHet.getChisq();
    
    /*
    double chisq_t = 0;

    for( size_t k=0; k<Cs.size(); k++) {

     chisq_t += chisq(Cs[k].getTransitionTable());
     
    }
    
    m_observedHomoChisq = chisq(m_Cpooled.getTransitionTable());
    
    // Compute observed $\chi^2d$
    m_observedHeteroChisq = chisq_t - m_observedHomoChisq;
    
    // cout << "chi2d=" << m_observedHeteroChisq << ", chi2c=" << m_observedHomoChisq << endl;
    */
}

///////////////////////////////////////////////////////////////////////////////

static string mapPvalMode(int pvaluemode)
{
    string discrepancy_measure;
    switch (pvaluemode) {
        case 8:
            discrepancy_measure = "fisher";
            break;
            
        case 9:
            discrepancy_measure = "chisq";
            break;
            
        default:
            break;
    }
    return discrepancy_measure;
}

double exact_heterogeneity_test(const vector< TransitionTable > & Cs,
                                const string & discrepancy_measure,
                                const CMPX_MARGINAL & marginal)
{
    EMTEnumerator e;

    if( discrepancy_measure == "fisher") {
        ProbHeteroGivenPooled v;
        return exact_multi_table_test(Cs, v, e);
    } else if ( discrepancy_measure == "chisq") {
        ProbHeteroChisqGivenPooled v;
        v.setCmpxMarginal(marginal);
        
        return exact_multi_table_test(Cs, v, e);
    } else {
        return 1;
    }
}

double exact_heterogeneity_test(const vector<TransitionTable> & tts,
                                int pvaluemode, const CMPX_MARGINAL & marginal)
{
    return exact_heterogeneity_test(tts, mapPvalMode(pvaluemode), marginal);
}

double exact_homogeneity_test(const vector<TransitionTable> & tts,
                               const string & discrepancy_measure)
{ // Exact homogeneity test
    
    double pc;
        
    TransitionTable pooled = pool(tts);
    
    pc = fisher_exact_test(pooled, discrepancy_measure);
    
    return pc;
}

double exact_homogeneity_test(const vector<TransitionTable> & tts,
                              int pvaluemode)
{
    return exact_homogeneity_test(tts, mapPvalMode(pvaluemode));
}

double exact_total_strength_test(const vector<TransitionTable> & tts,
                                 const string & discrepancy_measure)
{ // Exact total strength test

    return multi_table_fisher_test(tts, discrepancy_measure);

    /*
    double qt = 1.0;
    
    
    for(size_t i=0; i<tts.size(); i++) {
        qt *= 1 - fisher_exact_test(tts[i], discrepancy_measure);
    }
    
    return 1-qt;
    */
}

double exact_total_strength_test(const vector<TransitionTable> & tts,
                                 int pvaluemode)
{
    // Exact total strength test
    return exact_total_strength_test(tts, mapPvalMode(pvaluemode));
}

//compute the fisher's heterogeity exact test p-value
void exact_comparative_chisq_3(const vector<TransitionTable> & tts,
                             double &pd, double &pc, double &pt,
                             int pvaluemode, const CMPX_MARGINAL & marginal)
{
    string discrepancy_measure = mapPvalMode(pvaluemode);
    
    // Exact total strength test
    pt = exact_total_strength_test(tts, discrepancy_measure);
    
    // Exact homogeneity test
    pc = exact_homogeneity_test(tts, discrepancy_measure);
    
    // Exact heterogeneity test
    pd = exact_heterogeneity_test(tts, discrepancy_measure, marginal);
}

void test_exact_hetero_3 ()
{
    cout << endl << ">>>>> Testing exact heterogeneity test 3 ..." << endl;
    vector<TransitionTable> Cs(2);
    
    const size_t m1=2, m2=4, n=2;
    
    int mat1[m1][n] = {
        {8,0},
        {0,8}
    };
    
    int mat2[m2][n] = {
        {8,0},
        {0,8},
        {8,0},
        {8,0}
    };
    
    vector<int> parents1(1);
    parents1[0] = 3;
    
    Cs[0]=TransitionTable(0, 2, parents1, vector<int>(1,2), vector<int>(1,0));
    for(size_t i=0; i<m1; i++)
        for(size_t j=0; j<n; j++)
            Cs[0].set(i, j, mat1[i][j]);
    
    vector<int> parents2(2);
    parents2[0] = 2, parents2[1]=3;
    
    Cs[1]=TransitionTable(0, 2, parents2, vector<int>(2,2), vector<int>(2,0));
    for(size_t i=0; i<m2; i++)
        for(size_t j=0; j<n; j++)
            Cs[1].set(i, j, mat2[i][j]);
    
    double pd = exact_heterogeneity_test(Cs, "chisq");
    cout << "Example 1: pd=" << pd << endl;
    
    TransitionTable applyHeteroChisqTest(const TransitionTable & tt1,
                                         const TransitionTable & tt2);
    TransitionTable ttChisq = applyHeteroChisqTest(Cs[0], Cs[1]);
    cout << "chisq approximation: p.val=2.497998e-5, chisq_d=16+32-24=24, df=1x1+3x1-1x1=3"
    // << ttChisq.getpValue() << ", " << "chisq_d=" << ttChisq.getChisq()
    << endl << endl;
    
    int matAtoC[m1][n] = {
        {10,0},
        {5,5}
    };
    
    int matABtoC[m2][n] = {
        {5,0},
        {5,0},
        {5,0},
        {0,5}
    };
    
    vector<int> parentsAtoC(1);
    parentsAtoC[0] = 1;
    
    Cs[0]=TransitionTable(0, 2, parentsAtoC, vector<int>(1,2), vector<int>(1,0));
    for(size_t i=0; i<m1; i++)
        for(size_t j=0; j<n; j++)
            Cs[0].set(i, j, matAtoC[i][j]);
    
    vector<int> parentsABtoC(2);
    parentsABtoC[0] = 1, parentsABtoC[1]=2;
    
    Cs[1]=TransitionTable(0, 2, parentsABtoC, vector<int>(2,2), vector<int>(2,0));
    for(size_t i=0; i<m2; i++)
        for(size_t j=0; j<n; j++)
            Cs[1].set(i, j, matABtoC[i][j]);
    
    pd = exact_heterogeneity_test(Cs, "chisq");
    cout << "Example 2: pd=" << pd << endl;
    cout << "chisq approximation: p.val=0.00397, chisq_d=13.3333, df=1x1+3x1-1x1=3"
    << endl << endl;
    
    vector<TransitionTable> CsRev(2);
    CsRev[0] = Cs[1];
    CsRev[1] = Cs[0];
    pd = exact_heterogeneity_test(CsRev, "chisq");
    cout << "Example 3: pd=" << pd << endl;
    cout << "   (reverse of the two tables in Ex 2. should be same.)";
    cout << endl << endl;
    
}
