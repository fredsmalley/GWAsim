//FisherExactTest.cpp -- Implements Fisher's exact test of contigency tables
//
//created by Yang Zhang 1.2.2013

#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <numeric>

using namespace std;

#include "ChisqTests.h"

double factorial (double n);

double FisherProbConstant(const vector<vector<int> > & A);

double evaluate_fisher(vector<vector<int> > & A, double AStat, double OStat);

double update_fisher(double ARunStat, const vector<vector<int> > & A, int i, int j,
                     const vector<vector<double> > & E);

double evaluate_chisq(vector<vector<int> > & A, double AStat, double OStat);
double update_chisq(double ARunStat, const vector<vector<int> > & A, int i, int j,
                    const vector<vector<double> > & E);

double traverse(vector<vector<int> > & A, int i, int j,
                vector<int> & ARowRunSums,
                vector<int> & AColRunSums,
                double ARunStat,
                const vector<vector<double> > & E,
                const vector<vector<int> > & O,
                const vector<int> & ORowSums,
                const vector<int> & OColSums,
                const double OStat,
                double (* evaluate) (vector<vector<int> > & A, double AStat, double OStat),
                double (* update) (double ARunStat, const vector<vector<int> > & A,
                                   int i, int j,
                                   const vector<vector<double> > & E)
                );

double fisher_exact_test(const vector<vector<int > > & O,
                         const vector<vector<double > > & E)
{
    // Input:
    //   O - observed contigency table
    //   E - expected contigency table
    
    double p_val;
    
    if(O.size() <= 1) return 1;
    
    size_t nrows = O.size(), ncols = O[0].size();
    
    // Auxiliary data structures for observed contingency table O
    vector<int> ORowSums = getRowSums(O);
    vector<int> OColSums = getColSums(O);
    
    // Enumeration matrix A and auxiliary data structures
    vector<vector<int > > A(nrows, vector<int>(ncols, 0));
    
    vector<int> ARowRunSums(nrows, 0);
    vector<int> AColRunSums(ncols, 0);
    
    // Statistics associated with O and A
    double OStat;
    double ARunStat;
    
    if(E.size() == 0) {
        // Using standard Fisher's probability for discrepancy measure
        OStat = FisherProb(O);
        ARunStat = FisherProbConstant(O);
        
        p_val = traverse(A, 0, 0, ARowRunSums, AColRunSums, ARunStat, E,
                         O, ORowSums, OColSums, OStat,
                         evaluate_fisher, update_fisher);
        
    } else {
        
        // Using chi-square discrepancy measure
        OStat = chisq(O, E);
        ARunStat = 0;
        
        p_val = traverse(A, 0, 0, ARowRunSums, AColRunSums, ARunStat, E,
                         O, ORowSums, OColSums, OStat,
                         evaluate_chisq, update_chisq);
    }
    
    return p_val;
}




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
	for(size_t i = 0; i<A.size(); i++)
	{
		for(size_t j=0; j<A[i].size(); j++)
		{
			product *= factorial(A[i][j]);
		}
	}
    
	return FisherProbConstant(A) / product;
}

double evaluate_fisher(vector<vector<int> > & A, double AStat, double OStat)
{
	return AStat <= OStat ? AStat : 0;
}

double update_fisher(double ARunStat, const vector<vector<int> > & A, int i, int j,
              const vector<vector<double> > & E)
{
    return ARunStat / factorial( static_cast<double>(A[i][j]) );
}


//----------------------------------------------------------------------------
// chi-square related functions
double evaluate_chisq(vector<vector<int> > & A, double AStat, double OStat)
{
    return AStat >= OStat ? FisherProb(A) : 0;
}

double update_chisq(double ARunStat, const vector<vector<int> > & A, int i, int j,
              const vector<vector<double> > & E)
{
    if(E[i][j] != 0) {
        double d = A[i][j]-E[i][j];
        ARunStat +=  d * d / E[i][j];
    }

    return ARunStat;
}

//----------------------------------------------------------------------------
// The general branch-and-bound function to enumerate all contingency tables

void next(int & rowindex, int & colindex, int numbofrows, int numbofcols)
{ // return row-wise the next index in the matrix
    
	if(colindex < numbofcols-1)
	{
		colindex += 1;
	}else
	{
		colindex = 0;
		rowindex += 1;
	}
}

double traverse(vector<vector<int> > & A, int i, int j,
                 vector<int> & ARowRunSums,
                 vector<int> & AColRunSums,
                 double ARunStat,
                 const vector<vector<double> > & E,
                 const vector<vector<int> > & O,
                 const vector<int> & ORowSums,
                 const vector<int> & OColSums,
                 const double OStat,
                 double (* evaluate) (vector<vector<int> > & A, double AStat, double OStat),
                 double (* update) (double ARunStat, const vector<vector<int> > & A,
                                    int i, int j,
                                    const vector<vector<double> > & E)
                )
{
	double p;
    
	if(i==(int)O.size()) { // beyond a leaf node. A is complete

        p = evaluate(A, ARunStat, OStat);
        
	} else {

        int nrows = O.size(), ncols = O[0].size();

		int bcj = OColSums[j] - AColRunSums[j];
		int bri = ORowSums[i] - ARowRunSums[i];
        
        int ubound = bcj > bri ? bri : bcj; //upperbound
		int lbound = 0; //lowerbound
        
        if(i == nrows - 1)
		{
			lbound = bcj;
		}
		
		if(j== ncols - 1)
		{
			lbound = lbound > bri ? lbound : bri;
		}
        
        p = 0.0;
		for(int aij = lbound; aij <= ubound; aij++)
		{
			A[i][j] = aij;
			ARowRunSums[i] += aij;
			AColRunSums[j] += aij;
            double ARunStatBackup = ARunStat;
            ARunStat = update(ARunStat, A, i, j, E);

            int i_next = i, j_next = j;
			next(i_next, j_next, nrows, ncols);
			p += traverse(A, i_next, j_next, ARowRunSums, AColRunSums, ARunStat,
                          E, O, ORowSums, OColSums, OStat, evaluate, update);
            
            ARunStat = ARunStatBackup;
			ARowRunSums[i] -= aij;
			AColRunSums[j] -= aij;
		}
        
		A[i][j] = 0;
	}
	return p;
}

// =====================================================================
// Testing function
void test_fisher_exact()
{
    cout << ">>>>> Testing Fisher's exact test ..." << endl;
    
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
    double true_p_vals[] = {0.02235542, 0.0273726};
    const char * tests[] = {"Fisher exact test", "Fisher (Chisq) exact test"};
    
    vector<vector<int> > O(m, vector<int>(n));

    for(int r=0; r<m; r++) {
        for(int c=0; c<n; c++) {
            O[r][c] = X[r][c];
        }
    }
    
	time_t seconds1 = time (NULL);
    
    double calculated_p_vals[2];
    
    vector<vector<double > > E = getExpectedTable(O);
    
	// for(int s = 0; s<20000; s++)
	{
        calculated_p_vals[0] = fisher_exact_test(O);
        calculated_p_vals[1] = fisher_exact_test(O, E);
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
	
}