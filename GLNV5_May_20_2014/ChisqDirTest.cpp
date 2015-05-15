// ChisqDirTest.cpp
//
//  Test whether the column variable is some function
//  of the row variable, in a contingency table.  The row variable can be
//  a combination of multiple parents.  The column variable can be considered
//  as a child.
//  
//  This test is unsymmetric, i.e., it does not give the same test statistics
//  if we rotate the row and column into another contigency matrix.
//
//  If the interaction is presented by x -> y, we decide if y is a function of x, 
//  i.e., if y=f(x) for some f.  
//
//  The chi square statistic is computed as 
//  
//                chisq = sum_x( chisq(y|x) ) - chisq(y)
//
//  with degrees of freedom of 
//				  df = (|x|-1)(|y|-1)
//  where |x| is the radix of x and |y| is the radix of y.
//
// Joe Song
// Created: September 23, 2010
// Modified:
//   October 4, 2011. MS. Name changed from "ChisqFunTest.cpp"
//   April 7, 2014. MS. Added an option parameter "method" to ChisqDirTest()
//     to compute normalized chisq

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

using namespace std;

#include "StatDistributions.h"

void ChisquareTest1DNoPValue(const vector<int> & x_obs,
                             const vector<double> & p_null, int K,
                             double & chisquare, size_t & df);

double ChisqDirTest(const vector< vector<int> > & table_obs, double & chisquare,
                    size_t & df, const string & method)
{
    double pval;
    
	int nrow = (int) table_obs.size();
	int K = (int) table_obs[0].size();

	vector<double> p_null(K, 1.0/K); 
	vector<int> n_y(K);  // the histogram of child y
	
	for(int j = 0; j < K; j++){

		for(int i = 0; i < nrow; i++){

			n_y[j] += table_obs[i][j];

		}//end for

	}//end for

	chisquare = 0;

	for(int i=0; i<nrow; i++) {
		double chisq_row = 0;
		ChisquareTest1DNoPValue(table_obs[i], p_null, K, chisq_row, df);
		chisquare += chisq_row;
	}

	double chisq_y;

	ChisquareTest1DNoPValue(n_y, p_null, K, chisq_y, df);

	chisquare -= chisq_y;

	df = (nrow-1) * (K-1);

    if(method == "normalized") {
        
        chisquare = (chisquare - df) / sqrt(2*df);
        pval = NormalPvalue(chisquare, 0, 1, false);
        
    } else {
        pval = ChisqPvalue(chisquare, (int) df); // qchisq((int) df, chisquare);
    }
    
	return pval;
}

void test_ChisqDirTest()
{
    int x0[3][3] = {
        {4,0,1},
        {0,4,0},
        {4,0,1}
    };
    
    vector< vector<int> > x(3, vector<int>(3));

    for (size_t i=0; i<x.size(); ++i) {
        for (size_t j=0; j<x[i].size(); ++j) {
            x[i][j] = x0[i][j];
        }
    }
    
    double chisq, chisqNormalized;
    double pval, pvalNormalized;
    size_t df;
    
    pval = ChisqDirTest(x, chisq, df, "");

    if (abs(pval - 0.006122) < 0.000001) {
        cout << "ChisqDirTest using chi-square distribution passed." << endl;
    } else {
        cerr << "ChisqDirTest using chi-square distribution failed!" << endl;
    }
    
    pvalNormalized = ChisqDirTest(x, chisqNormalized, df, "normalized");

    if (abs(pvalNormalized - 0.000118) < 0.000001) {
        cout << "ChisqDirTest using normalized chi-square passed." << endl;
    } else {
        cerr << "ChisqDirTest using normalized distribution failed!" << endl;
    }

}
