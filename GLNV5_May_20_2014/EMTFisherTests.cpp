//
//  EMTFisherTests.cpp
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//
#include <cmath>
#include <ctime>
#include <numeric>

using namespace std;

#include "EMTFisherTests.h"

double factorial (double n);
double Fexact(const vector< vector<int> > & X);
double factorial_vector(vector<int>);//Hua Added, May 2 2014

//----------------------------------------------------------------------------
// Fisher exact statsitics related functions
double FisherProbConstant(const vector<vector<int> > & A)
{
    // Input:
    //   A -- contingency table
    
	vector<int> rowsums = getRowSums(A);
	double rowprod = 1.0;
	for(size_t i=0; i<rowsums.size(); i++)
	{
		rowprod *= factorial(rowsums[i]);
	}
    
	vector<int>	colsums = getColSums(A);
	double colprod = 1.0;
	for(size_t i=0; i<colsums.size(); i++)
	{
		colprod *= factorial(colsums[i]);
	}
    
	// int totalcounts = getTotalSum(A);
    int total = accumulate(rowsums.begin(), rowsums.end(), 0);
    return rowprod * colprod / factorial(total);
}

//compute fisher's probability
double FisherProb(const vector<vector<int> > & A)
{
    // Input:
    //   A -- contingency table
    
	double product = 1.0;
    
    //Modified by Hua, May 2 2014
	vector<int> rowsums = getRowSums(A);
	vector<int>	colsums = getColSums(A);
	vector<int> v;
	v.reserve(rowsums.size() + colsums.size());
	v.insert(v.end(), rowsums.begin(), rowsums.end());
	v.insert(v.end(), colsums.begin(), colsums.end());
	int total = accumulate(rowsums.begin(), rowsums.end(), 0);
	v.push_back(-total);
    
	for(size_t i = 0; i<A.size(); i++)
	{
		for(size_t j=0; j<A[i].size(); j++)
		{
			v.push_back(-A[i][j]);//Added by Hua, May 2 2014
		}
	}
    
	return factorial_vector(v);//Added by Hua, May 2 2014
    ////
    
//Hua comment it, May 2 2014
//	for(size_t i = 0; i<A.size(); i++)
//	{
//		for(size_t j=0; j<A[i].size(); j++)
//		{
//			product *= factorial(A[i][j]);
//		}
//	}
//    
//	return FisherProbConstant(A) / product;
////
}

void EMTFisherProb::processTable(size_t k, const EMTEnumerator & e,
                                 const vector<TransitionTable> & Cs)
{
    m_nullFisherProb[k] = FisherProb(e.As[k].getTransitionTable());
}

void EMTFisherProb::initialize(const vector<TransitionTable> & Cs)
{
    // compute required row sums
    m_requiredRowSums.resize(Cs.size());
    for(size_t k=0; k < Cs.size(); k++) {
        m_requiredRowSums[k] = Cs[k].getRowSums();
    }
    
    // compute required column sums
    m_requiredColSums.resize(Cs.size());
    for(size_t k=0; k < Cs.size(); k++) {
        m_requiredColSums[k] = Cs[k].getColSums();
    }
    
    m_observedFisherProb.resize( Cs.size() );
    // compute observed statistics
    for(size_t k=0; k < Cs.size(); k++) {
        m_observedFisherProb[k] = FisherProb(Cs[k].getTransitionTable());
    }
    
    m_nullFisherProb.resize(Cs.size());
}

bool EMTFisherProb::isMoreExtreme() const
{
    /*
    for(size_t k=0; k < m_nullFisherProb.size(); k++) {
        if( m_nullFisherProb[k] <= m_observedFisherProb[k] ) {
            return true;
        }
    }
    return false;
    */
   
    double nullProb = 1, observedProb = 1;
    
    for(size_t k=0; k < m_nullFisherProb.size(); k++) {
        nullProb *= m_nullFisherProb[k];
        observedProb *= m_observedFisherProb[k];
    }
    
    return nullProb <= observedProb;
    
}

double EMTFisherProb::evaluate(const EMTEnumerator & e,
                               const vector<TransitionTable> & Cs)
{
    // Multiple independent Fisher's exact tests
    double P;
    
    if( isMoreExtreme() ) {
        P = 1;
        for(size_t k=0; k < Cs.size(); k++) {
            P *= m_nullFisherProb[k];
        }
    } else {
        P = 0;
    }
    
    return P;
}

vector<TransitionTable>
EMTFisherProb::generateTables(const vector<TransitionTable> & Cs) const
{
    vector<TransitionTable> As(Cs);
    
    for (size_t k=0; k<Cs.size(); k++) {
        As[k].reset();
    }
    
    return As;
}

////////////////////////////////////////////////////////////////////////////////
double multi_table_fisher_test(const vector<TransitionTable> & Cs,
                               const string & discrepancy_measure)
{
    EMTEnumerator e;
    
    if( discrepancy_measure == "fisher" ) {
        EMTFisherProb v;
        return exact_multi_table_test(Cs, v, e);
    } else if ( discrepancy_measure == "chisq" ) {
        EMTFisherChisq v;
        return exact_multi_table_test(Cs, v, e);
    } else {
        cerr << "ERROR: unknown discrepancy measure in multi_table_fisher_test!"
        << endl;
        return 1;
    }
}

double fisher_exact_test(const TransitionTable & C,
                         const string & discrepancy_measure)
{
    vector<TransitionTable> Cs(1, C);
    return multi_table_fisher_test(Cs, discrepancy_measure);
}

double fisher_exact_test(const vector< vector<int> > & C,
                         const string & discrepancy_measure)
{
    double pval;
    
    if( discrepancy_measure == "fisher" ) {
 
        pval = Fexact(C);
        
    } else {

        TransitionTable tt(0, C[0].size(), vector<int>(1,0),
                           vector<int>(1, C.size()), vector<int>(1,0));
        tt.setTransitionTable(C);
        pval = fisher_exact_test(tt, discrepancy_measure);
    }
    
    return pval;
}


////////////////////////////////////////////////////////////////////////////////
// =====================================================================
// Testing function
void test_EMT_fisher()
{
    cout << ">>>>> Testing EMT Fisher's exact test ..." << endl;
    
    /*
     const int m = 3, n = 3;
     int X[3][3] = {
     {10,0,0},
     {0,10,0},
     {0,0,10}
     };
     */
    
    const int m = 3, n = 3;
    int X[m][n] = {
        {3,0,4},
        {1,3,0},
        {1,3,0}
    };
    
    double true_p_vals[2] = {0.0, 0.0273726}; // = {0.02235542, 0.0273726};
    
    
    const char * tests[] = {"Fisher exact test", "Fisher (Chisq) exact test"};
    
    vector<vector<int> > O(m, vector<int>(n));
    
    for(int r=0; r<m; r++) {
        for(int c=0; c<n; c++) {
            O[r][c] = X[r][c];
        }
    }

    true_p_vals[0] = Fexact(O);
    
	time_t seconds1 = time (NULL);
    
    double calculated_p_vals[2];
    
	// for(int s = 0; s<20000; s++)
	{
        calculated_p_vals[0] = fisher_exact_test(O, "fisher");
        calculated_p_vals[1] = fisher_exact_test(O, "chisq");
    }
	
    time_t seconds2 = time (NULL);
	cout << "Seconds taken " << difftime (seconds2, seconds1) << endl;
	
    for(int i=0; i<2; i++) {
        
        cout<< tests[i] << "p.val = " << calculated_p_vals[i]
        << " (truth=" << true_p_vals[i] << ")" << endl;
        
        if(abs(calculated_p_vals[i]-true_p_vals[i]) > 1e-7) {
            cerr << "ERROR: failed!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    cout << "PASSED Fisher exact test testing." << endl << endl;
    
    const int K = 2, M = 9;
    int Y[K][M] = {
        {0,0, 0,0,0,0,6,2,0},
        {0,1,13,0,0,2,0,0,0}
    };
    
    vector<vector<int> > YO(K, vector<int>(M));
    
    for(int r=0; r<K; r++) {
        for(int c=0; c<M; c++) {
            YO[r][c] = Y[r][c];
        }
    }
        
	{
        calculated_p_vals[0] = fisher_exact_test(YO, "fisher");
        calculated_p_vals[1] = fisher_exact_test(YO, "chisq");
    }
    
    for(int i=0; i<2; i++) {
        cout<< tests[i] << "p.val = " << calculated_p_vals[i] << endl;
    }
    
}