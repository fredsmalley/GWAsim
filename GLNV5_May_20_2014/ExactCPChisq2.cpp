//
//  ExactCPChisq2.cpp
//  gln
//
//  Created by Joe Song on 1/25/13.
//
//
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

#include "SetOps.h"
#include "ChisqTests.h"
#include "TransitionTable.h"
#include "GLN.h"

vector<vector<TransitionTable> >
getPooledAndUnique(const vector<TransitionTable> & tts)
// Compute pooled transition table and transition table with unique parents
//   for each condition
// OUTPUT:
//   The first dimension is pooled, second dimension is for unique
// ASSUMPTIONS:
//   1. The same parent ID will have the same delay
//   2. Common parents show up in the same order in parent list for each
//      condition
{
    vector<vector<TransitionTable> > ttspoolandhetero(2);
    ttspoolandhetero[0].resize(1);
    ttspoolandhetero[1].resize(tts.size());
    
    vector<int> intersecParents = tts[0].getParents();
    for(size_t i=1; i<tts.size(); i++) {
        if(intersecParents.empty()) {
            break;
        }
        intersecParents = intersectionSet(intersecParents, tts[i].getParents());
    }
    
    for(size_t i=0; i<tts.size(); i++) {
        
        vector<int> uniqueParents = complementSet(intersecParents, tts[i].getParents());
        ttspoolandhetero[1][i] = tts[i].subtable(uniqueParents); // hetero
        ttspoolandhetero[0][0] += tts[i].subtable(intersecParents); // pooled
    }
    
    return ttspoolandhetero;
}


vector<vector<vector<double> > >
getExpectedTables(const vector<TransitionTable> & tts)
// Compute expect contingency table for each condition
{
    vector<vector<vector<double> > > ttsexpected(tts.size()); //only contains the contingency table
    
    int sum = 0; //total number of observations
    for(size_t i=0; i<tts.size(); i++) {
        sum += tts[i].getTotalSum();
    }
    
    vector<vector<TransitionTable> > pooledandunique = getPooledAndUnique(tts);
    
    const vector<vector<int> > & ctpooled = pooledandunique[0][0].getTransitionTable(); //pooled contingency table
    const vector<int> & pooledparents = pooledandunique[0][0].getParents();
    const vector<int> & pooledBases = pooledandunique[0][0].getParentBases();
    
    vector<int> pooledparentvalues(pooledparents.size());
    
    const vector<TransitionTable> & ttshetero = pooledandunique[1];
    
    for(size_t i=0; i<tts.size(); i++) { //compute expected tables for each condition
        
        const vector<vector<int > > & tt_i = tts[i].getTransitionTable();
        const vector<int> & parentbases = tts[i].getParentBases();
        
        int sum_i = tts[i].getTotalSum();
        
        if( tt_i.size() > 0 ) {
            
            ttsexpected[i].resize(tt_i.size(), vector<double>(tt_i[0].size(), 0) );
            
            if(sum_i == 0) continue;
            
            vector<int> rowsums = tts[i].getRowSums();
            vector<int> colsums = tts[i].getColSums();
            vector<int> heterorowsums = ttshetero[i].getRowSums();
            vector<int> heterocolsums = ttshetero[i].getColSums();
            
            //compute unique parent values
            const vector<int> & uniqueParents = ttshetero[i].getParents();
            const vector<int> & uniqueBases = ttshetero[i].getParentBases();
            vector<int> uniqueparentvalues(uniqueParents.size());
            
            vector<int> allparentvalues(parentbases.size()); //bit-wise, not decimal
            
            for(size_t row=0; row < tt_i.size(); row ++) {
                
                LinearIndexToArrayIndex(parentbases, row, allparentvalues);
                //compute pooled parent values for each condition
                const vector<int> & parents = tts[i].getParents();
                
                for(size_t parentindex = 0; parentindex<parents.size(); parentindex++) {
                    
                    size_t pooledindex = find(pooledparents.cbegin(), pooledparents.cend(),
                                              parents[parentindex]) - pooledparents.cbegin();
                    
                    if(pooledindex < pooledparents.size()) {
                        pooledparentvalues[pooledindex] = allparentvalues[parentindex];
                        continue;
                    }
                    
                    size_t uniqueindex = find(uniqueParents.cbegin(), uniqueParents.cend(),
                                              parents[parentindex]) - uniqueParents.cbegin();
                    
                    if(uniqueindex < uniqueParents.size()) {
                        uniqueparentvalues[uniqueindex] = allparentvalues[parentindex];
                        continue;
                    }
                }
                
                double p_hetero, p_pooled;
                
                if(uniqueparentvalues.size() > 0) {
                    //convert parent values into indices in hetero matrices
                    long long unsigned uniquerowindex
                    = ArrayIndexToLinearIndex(uniqueBases, uniqueparentvalues);
                    p_hetero = heterorowsums[uniquerowindex] / (double) sum_i;
                    
                } else {
                    p_hetero = 1.0;
                }
                
                long long unsigned pooledrowindex;
                if(pooledparentvalues.size() > 0) {
                    //convert parent values into indices in pooled matrices
                    pooledrowindex = ArrayIndexToLinearIndex(pooledBases, pooledparentvalues);
                }
                
                for(size_t col=0; col<tt_i[0].size(); col++) {
                    
                    p_pooled = pooledparentvalues.size() > 0 ?
                    ctpooled[pooledrowindex][col] / (double) sum :
                    colsums[col] / (double) sum_i;
                    
                    ttsexpected[i][row][col] = sum_i * p_pooled * p_hetero;
                }
            }
        }
    }
    return ttsexpected;
}

vector<vector<vector<double > > >
getExpectedTables(const vector<vector<vector<int > > > & Os)
{
    vector<vector<vector<double > > >
    Es(Os.size(),
       vector<vector<double > >(Os[0].size(),
                             vector<double>(Os[0][0].size(), 0)));
    
    vector<vector<int > > pooled(Os[0].size(),
                                 vector<int>(Os[0][0].size(), 0));
    
    vector<int> ns(Os.size(), 0);
    int n=0;
    
    for(size_t k=0; k<Os.size(); k++) {
        for(size_t i=0; i<Os[k].size(); i++) {
            for(size_t j=0; j<Os[k][i].size(); j++) {
                pooled[i][j] += Os[k][i][j];
                ns[k] += Os[k][i][j];
            }
        }
        n += ns[k];
    }

    for(size_t k=0; k<Os.size(); k++) {
        for(size_t i=0; i<Os[k].size(); i++) {
            for(size_t j=0; j<Os[k][i].size(); j++) {
                if(n > 0) {
                    Es[k][i][j] = pooled[i][j] * ns[k] / (double) n;
                }
            }
        }
    }
    
    return Es;
}

double evaluate(const vector<vector<vector<int > > > & As,
                const vector<vector<vector<int > > > & Cs)
{
    double P=0.0;

    // Compute expected contingency tables A0(1), A0(2), ..., A0(K) for A1...AK
    vector<vector<vector<double > > > A0 = getExpectedTables(As);
    vector<vector<vector<double > > > C0 = getExpectedTables(Cs);
    
    bool more_heterogeneous = false;

    for(size_t k=0; k<As.size(); k++) {
        
       if ( chisq( As[k], A0[k] ) >= chisq( Cs[k], C0[k] ) ) {
            more_heterogeneous = true;
            break;
       }
    }
    
    if(more_heterogeneous) {
        P = 1.0;
        
        for(size_t k=0; k<As.size(); k++) {
            // P <- Prob(A1) * Prob(A2) * ... * Prob(AK)
            P *= FisherProb(As[k]);
        }
    }
            
    return P;
}

double evaluate(const vector< TransitionTable > & As,
                const vector< TransitionTable > & Cs)
{
    double P=0.0;
    
    // Compute expected contingency tables A0(1), A0(2), ..., A0(K) for A1...AK
    vector< vector<vector<double> > > A0 = getExpectedTables(As);
    vector< vector<vector<double> > > C0 = getExpectedTables(Cs);
    
    bool more_heterogeneous = false;
    
    for(size_t k=0; k<As.size(); k++) {
        
        if ( chisq( As[k].getTransitionTable(), A0[k] )
            >= chisq( Cs[k].getTransitionTable(), C0[k] ) ) {
            more_heterogeneous = true;
            break;
        }
    }
    
    if(more_heterogeneous) {
        P = 1.0;
        
        for(size_t k=0; k<As.size(); k++) {
            // P <- Prob(A1) * Prob(A2) * ... * Prob(AK)
            P *= FisherProb(As[k].getTransitionTable());
        }
    }
    
    return P;
}

void next(size_t & k, size_t & i, size_t & j, size_t I, size_t J)
{
    // k: contingency table index from 0 to K-1 (#matrices)
    // i: row index from 0 to I-1 (I: #rows)
    // j: column index from 0 to J-1 (J: #columns)

    if( j < J-1 ) {  // J: # columns in Ak
        j = j+1;
    } else if ( i < I-1 ) { // I: #rows in Ak
        j = 0;
        i = i+1;
    } else {
        k = k+1;
        i = 0;
        j = 0;
    }
    // cout << "k=" << k << ", i=" << i << ", j=" << j << endl;
}

double traverse(size_t k, size_t i, size_t j,
                vector<vector<vector<int > > > & As,
                vector<int> & AkRowsums, vector<int> & AkColsums,
                const vector<vector<vector<int > > > & Cs,
                const vector<vector<int> > & CRowsums,
                const vector<vector<int> > & CColsums)
{
    double P;
 
    vector<int> * pAkRowsums = & AkRowsums;
    vector<int> * pAkColsums = & AkColsums;
    // Branch-and-bound
    if(k == As.size()) {
        // k being out of bound indciates (k,i,j) is beyond a leaf
        P = evaluate(As, Cs);
        return P;
    } else {
        P = 0;
    }

    if(i == 0 && j==0) {
        pAkRowsums = new vector<int>(As[k].size());
        pAkColsums = new vector<int>(As[k][0].size());
    }

    // lower and upper bounds for Ak[i,j]
    int bri = CRowsums[k][i] - (*pAkRowsums)[i]; // Ui(Ck)-Ui(Ak);
    int bcj = CColsums[k][j] - (*pAkColsums)[j]; // Vj(Ck)-Vj(Ak);
    int lij = 0;
    if(i == As[k].size()-1) { // last row
        lij = bcj;
    }
    if(j == As[k][i].size() -1) { // last column
        lij = max(lij, bri);
    }
    int uij = min(bcj, bri);

    for( int a=lij; a <= uij; a++) { // will enter loop if and only if lij<=uij

        As[k][i][j] = a;
    
        size_t k_next=k, i_next=i, j_next=j;
        next(k_next, i_next, j_next, As[k].size(), As[k][0].size());
        
        (*pAkRowsums)[i] += a;
        (*pAkColsums)[j] += a;
        
        P += traverse(k_next, i_next, j_next, As, *pAkRowsums, *pAkColsums,
                      Cs, CRowsums, CColsums);

        (*pAkRowsums)[i] -= a; // = Ui(Ak) - Ak[i,j];
        (*pAkColsums)[j] -= a; // = Vj(Ak) - Ak[i,j];

        As[k][i][j] = 0; // clear up the entry to return
               
    }

    if(i == 0 && j==0) {
        delete pAkRowsums;
        delete pAkColsums;
    }
    
    return P;
}

double traverse(size_t k, size_t i, size_t j,
                vector< TransitionTable > & As,
                vector<int> & AkRowsums, vector<int> & AkColsums,
                const vector< TransitionTable > & Cs,
                const vector<vector<int> > & CRowsums,
                const vector<vector<int> > & CColsums)
{
    double P;
    
    vector<int> * pAkRowsums = & AkRowsums;
    vector<int> * pAkColsums = & AkColsums;
    // Branch-and-bound
    if(k == As.size()) {
        // k being out of bound indciates (k,i,j) is beyond a leaf
        P = evaluate(As, Cs);
        return P;
    } else {
        P = 0;
    }
    
    if(i == 0 && j==0) {
        pAkRowsums = new vector<int>(As[k].nrow());
        pAkColsums = new vector<int>(As[k].getBase());
    }
    
    // lower and upper bounds for Ak[i,j]
    int bri = CRowsums[k][i] - (*pAkRowsums)[i]; // Ui(Ck)-Ui(Ak);
    int bcj = CColsums[k][j] - (*pAkColsums)[j]; // Vj(Ck)-Vj(Ak);
    int lij = 0;
    if(i == As[k].nrow()-1) { // last row
        lij = bcj;
    }
    if((int) j == As[k].getBase() -1) { // last column
        lij = max(lij, bri);
    }
    int uij = min(bcj, bri);
    
    for( int a=lij; a <= uij; a++) { // will enter loop if and only if lij<=uij
        
        As[k].set(i, j, a); // As[k][i][j] = a;
        
        size_t k_next=k, i_next=i, j_next=j;
        next(k_next, i_next, j_next, As[k].nrow(), As[k].getBase());
        
        (*pAkRowsums)[i] += a;
        (*pAkColsums)[j] += a;
     
        // cout << "k=" << k << ", i=" << i << ", j=" << j << ", a=" << a << endl;
        
        P += traverse(k_next, i_next, j_next, As, *pAkRowsums, *pAkColsums,
                      Cs, CRowsums, CColsums);
        
        (*pAkRowsums)[i] -= a; // = Ui(Ak) - Ak[i,j];
        (*pAkColsums)[j] -= a; // = Vj(Ak) - Ak[i,j];
        
        As[k].set(i, j, 0); // As[k][i][j] = 0; // clear up the entry to return
        
    }
    
    if(i == 0 && j==0) {
        delete pAkRowsums;
        delete pAkColsums;
    }
    
    return P;
}

double exact_heterogeneity_test_2(const vector<vector<vector<int > > > & Cs)
{
    // Compute pooled contingency table C0 from C1...CK
    // Compute expected C0(k) as described in Section 4.1
    
    vector<vector<int> > CRowsums(Cs.size());
    vector<vector<int> > CColsums(Cs.size());
    
    for(size_t k=0; k<Cs.size(); k++) {
        CRowsums[k] = getRowSums(Cs[k]);
        CColsums[k] = getColSums(Cs[k]);
    }
    
    vector<vector<vector<int > > > As(Cs.size());
    for(size_t k=0; k<As.size(); k++) {
        // A1...AK <- zero matrices of identical sizes with C1...CK
        As[k] = vector<vector<int> >(Cs[k].size(), vector<int>(Cs[k][0].size()));
    }
    
    vector<int> Uk; // row sums: Ui(A1) <- [0...0]
    vector<int> Vk; // col sums: Vj(A1) <- [0...0]
    
    double pd = traverse(0, 0, 0, As, Uk, Vk, Cs, CRowsums, CColsums);
    
    return pd;
}

double exact_heterogeneity_test_2(const vector< TransitionTable > & Cs)
{
    // Compute pooled contingency table C0 from C1...CK
    // Compute expected C0(k) as described in Section 4.1
    
    vector<vector<int> > CRowsums(Cs.size());
    vector<vector<int> > CColsums(Cs.size());
    
    for(size_t k=0; k<Cs.size(); k++) {
        CRowsums[k] = Cs[k].getRowSums();
        CColsums[k] = Cs[k].getColSums();
    }
    
    vector< TransitionTable > As(Cs.size());
    for(size_t k=0; k<As.size(); k++) {
        // A1...AK <- zero matrices of identical sizes with C1...CK
        As[k] = Cs[k];
        As[k].reset();
    }
    
    vector<int> Uk; // row sums: Ui(A1) <- [0...0]
    vector<int> Vk; // col sums: Vj(A1) <- [0...0]
    
    double pd = traverse(0, 0, 0, As, Uk, Vk, Cs, CRowsums, CColsums);
    
    return pd;
}

double exact_total_strength_test2(const vector<TransitionTable> & tts,
                                 int pvaluemode)
{ // Exact total strength test
    double qt = 1.0;
    
	for(size_t i=0; i<tts.size(); i++) {
        
		const vector<vector<int> > & originaltt = tts[i].getTransitionTable();
		
		if(originaltt.size() > 0) {
            
			double pti=1.0;
            
			if(pvaluemode == 8) {
				pti = fisher_exact_test(originaltt, "fisher");
			} else if(pvaluemode == 9) {
                vector<vector<double> > expectedtt = getExpectedTable(originaltt);
                
				pti = fisher_exact_test(originaltt, "chisq"); // expectedtt);
			}
			qt *= abs(1 - pti);
		}
	}
    
	double pt = 1 - qt;
    
    return pt;
}

double exact_homogeneity_test2(const vector<TransitionTable> & tts,
                              int pvaluemode)
{ // Exact homogeneity test
    
    double pc=1.0;
    
    vector<vector<TransitionTable> > pooledandunique = getPooledAndUnique(tts);
    const TransitionTable & pooledTable = pooledandunique[0][0];
    
    const vector<vector<int> > & pooled = pooledTable.getTransitionTable();
    
	if(pvaluemode == 8) {
        
		pc = fisher_exact_test(pooled, "fisher");
        
	} else if(pvaluemode == 9) {
        
		if(pooled.size() > 0) {
			vector<vector<double> > expectedpc = getExpectedTable(pooled);
			pc = fisher_exact_test(pooled, "chisq"); //expectedpc);
		} else {
			pc = 1;
		}
	}
    
    return pc;
}

//compute the fisher's heterogeity exact test p-value
void exact_comparative_chisq_2(const vector<TransitionTable> & tts,
                               double &pd, double &pc, double &pt,
                               int pvaluemode)
{
    // Exact total strength test
    pt = exact_total_strength_test2(tts, pvaluemode);
    
    // Exact homogeneity test
    pc = exact_homogeneity_test2(tts, pvaluemode);
    
    // Exact heterogeneity test
    pd = exact_heterogeneity_test_2(tts);
    
}

vector<vector<int> > convert(int x[3][3], const int m, const int n)
{
    vector<vector<int> > X(m, vector<int>(n));

    for(int r=0; r<m; r++) {
        for(int c=0; c<n; c++) {
            X[r][c] = x[r][c];
        }
    }
    
    return X;
}

void test_exact_hetero()
{
    double pd;
    
    vector<vector<vector<int > > > Cs(2);
    vector<vector<vector<int > > > Ds(2);
    const int m = 3, n = 3;

    int c1[m][n] = {
        {3,0,4},
        {1,3,0},
        {1,3,0}
    };

    int c2[m][n] = {
        {3,0,4},
        {1,3,0},
        {1,3,0}
    };

    Cs[0] = convert(c1, m, n);
    Cs[1] = convert(c2, m, n);
    
    pd = exact_heterogeneity_test_2(Cs);
    cout << "Exact heterogeneity test p-value=" << pd << endl;

    int c3[m][n] = {
        {3,0,0},
        {0,3,0},
        {0,0,3}
    };
    
    int c4[m][n] = {
        {0,0,3},
        {0,3,0},
        {3,0,0}
    };
    
    Ds[0] = convert(c3, m, n);
    Ds[1] = convert(c4, m, n);
    
    pd = exact_heterogeneity_test_2(Ds);
    cout << "Exact heterogeneity test p-value=" << pd << endl;
    
    cout << "evaluate(Cs, Ds)=" << evaluate(Cs, Ds) << endl;  // shall be zero

    cout << "evaluate(Ds, Cs)=" << evaluate(Ds, Cs) << endl;  // shall be positive

}
