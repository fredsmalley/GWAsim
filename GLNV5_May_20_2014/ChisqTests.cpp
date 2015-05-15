// ChisqTests.cpp -- Chisquare tests on contingency tables
//
// Joe Song
//
// Created: June 16, 2006
// Updated: September 23, 2006.  Added MultinomialChisquareTestTable()
// Modified: September 5, 2008.  Corrected the degrees of freedom
// Last modified: Name changed from multinom.cpp to ChisqTests.cpp
// Modified: September 11, 2011.  Add a new chisq comparison function to make the 
//           comparison more stringent than penalizing by degrees of freedom.

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <numeric>
#include <vector>
#include <list>
#include <string>
#include <algorithm>

//Added by Hua Apr 29 2014
#if defined _WIN32
#include <time.h>
#endif
////

using namespace std;

#include "StatDistributions.h"
#include "ChisqTests.h"

vector< vector<int> > greedyconst(const vector<int> & a, const vector<int> & b)
{
    vector< vector<int> > x(a.size(), vector<int>(b.size(), 0));
    
    vector<int> a1(a), b1(b);
    
    double T = 0;

    for (size_t m=0; m<a.size(); m++) {
        T += a[m];
    }
    
    list<int> M, N;
    
    for(int m=0; m<(int)a.size(); m++) {
        M.push_back(m);
    }

    for(int n=0; n<(int)b.size(); n++) {
        N.push_back(n);
    }
    
    while (M.size() > 0 && N.size() > 0) {
        
        list<int>::const_iterator iM, iN;
        
        int i, j, u=0;
        double rmax=-1;

        for(iM = M.cbegin(); iM != M.cend(); iM++) {
            for(iN = N.cbegin(); iN != N.cend(); iN++) {
                int uhk = min(a1[*iM], b1[*iN]);
                double e = a[*iM] * b[*iN] / T;
                double r;
                if( e == 0 ) {
                    r = 0;
                } else {
                    r = (uhk - e) * (uhk - e) / e;
                }
                if(r > rmax) {
                    rmax = r;
                    i = *iM;
                    j = *iN;
                    u = uhk;
                }
            }
        }
        x[i][j] = u;
        a1[i] -= u;
        b1[j] -= u;
        if (a1[i] == 0) {
            M.remove(i);
        }
        if (b1[j] == 0) {
            N.remove(j);
        }
    }
    return x;
}

vector< vector<int> > greedyvar(const vector<int> & a, const vector<int> & b)
{
    vector< vector<int> > x(a.size(), vector<int>(b.size(), 0));
    
    vector<int> a1(a), b1(b);
    
    double T = 0.0, T1;
    
    for (size_t m=0; m<a.size(); m++) {
        T += a[m];
    }
    T1 = T;
    
    list<int> M, N;
    
    for(int m=0; m<(int)a.size(); m++) {
        M.push_back(m);
    }
    
    for(int n=0; n<(int)b.size(); n++) {
        N.push_back(n);
    }
    
    while (M.size() > 0 && N.size() > 0) {
        
        list<int>::const_iterator iM, iN;
        
        int i, j, u=0;
        double rmax=-1;
        
        for(iM = M.cbegin(); iM != M.cend(); iM++) {
            for(iN = N.cbegin(); iN != N.cend(); iN++) {
                int uhk = min(a1[*iM], b1[*iN]);
                double e = a1[*iM] * b1[*iN] / T1;
                double r;
                if( e == 0 ) {
                    r = 0;
                } else {
                    r = (uhk - e) * (uhk - e) / e;
                }
                if(r > rmax) {
                    rmax = r;
                    i = *iM;
                    j = *iN;
                    u = uhk;
                }
            }
        }
        x[i][j] = u;
        a1[i] -= u;
        b1[j] -= u;
        T1 -= u;
        if (a1[i] == 0) {
            M.remove(i);
        }
        if (b1[j] == 0) {
            N.remove(j);
        }
    }
    return x;
}

double theoretical(const vector<int> & a, const vector<int> & b)
{
    int M1=0, N1=0;  // non zero rows and columns
    int T = 0;
    for (size_t m=0; m<a.size(); m++) {
        M1 = ( a[m] != 0 ) ? M1+1 : M1;
        T += a[m];
    }
    for (size_t n=0; n<b.size(); n++) {
        N1 = ( b[n] != 0 ) ? N1+1 : N1;
    }

    return ((double)T) * min(M1-1, N1-1);
}

double chisq_upper_bound(const vector<int> & a, const vector<int> & b,
                         const string & method)
{
    double upper;
    
    vector< vector<int> > x;
    
    if (method == "GREEDYCONST") {
        x = greedyconst(a, b);
        upper = chisq(x);
    } else if (method == "GREEDYVAR") {
        x = greedyvar(a, b);
        upper = chisq(x);
    } else {
        x = greedyconst(a, b);
        double upper1 = chisq(x);
        x = greedyvar(a, b);
        upper = chisq(x);
        upper = min(upper, upper1);
        upper1 = theoretical(a, b);
        // assert(upper <= upper1);
        upper = min(upper, upper1);
    }
    
    return upper;
}

void test_chisq_upper_bound()
{
    vector<int> a(3), b(3);
    a[0]=13, a[1]=22, a[2]=57;
    b[0]=28, b[1]=35, b[2]=29;
    
    double upper;
    
    upper = chisq_upper_bound(a, b, "GREEDYCONST");
    
    if(abs(upper - 94.3415435139573) > 1e-6) {
        cout << "Function greedyconst() failed!" << endl;
        exit(EXIT_FAILURE);
    } else {
        cout << "Function greedyconst() passed the test." << endl;
    }

    upper = chisq_upper_bound(a, b, "GREEDYVAR");
    
    if(abs(upper - 84.7711) > 1e-4) {
        cout << "Function greedyvar() failed!" << endl;
        exit(EXIT_FAILURE);
    } else {
        cout << "Function greedyvar() passed the test." << endl;
    }

    upper = chisq_upper_bound(a, b);
    
    if(abs(upper - 84.7711) > 1e-4) {
        cout << "Function chisq_upper_bound() failed!" << endl;
        exit(EXIT_FAILURE);
    } else {
        cout << "Function chisq_upper_bound() passed the test." << endl;
    }

}

bool operator != (const Chisq & c1, const Chisq & c2)
{
    return ! (c1 == c2);
}

bool operator == (const Chisq & c1, const Chisq & c2)
{
    if (c1.m_x != c2.m_x) {
        return false;
    } else if (c1.m_df != c2.m_df) {
        return false;
    } else if (c1.m_pval != c2.m_pval) {
        return false;
    } else {
        return true;
    }
}

vector<int> getRowSums(const vector<vector<int> > & tt)
{
	size_t nRows = tt.size();
	vector<int> rowSums(nRows, 0);
	for(size_t row = 0; row < nRows; row ++) {
		for(size_t col = 0; col < tt[row].size(); col ++) {
			rowSums[row] += tt[row][col];
		}
	}
    
	return rowSums;
}

vector<int> getColSums(const vector<vector<int> > & tt)
{
	size_t nRows = tt.size();
    
    vector<int> colSums;
    
    if(nRows > 0) {
        
        colSums.resize(tt[0].size(), 0);
        
        for(size_t col = 0; col < tt[0].size(); col ++) {
            for(size_t row = 0; row < nRows; row ++) {
                colSums[col] += tt[row][col];
            }
        }
    }
    
	return colSums;
}

int getTotalSum(const vector<vector<int> > & tt)
{
	int totalsum = 0;
	for(size_t row = 0; row < tt.size(); row ++) {
		for(size_t col = 0; col < tt[0].size(); col ++) {
			totalsum += tt[row][col];
		}
	}
    
	return totalsum;
}

double chisq(const vector<vector<int> > & O, const vector<vector<double> > & E)
{
    double chisq = 0.0;
    
    for(size_t i=0; i<O.size(); i++) {
        for(size_t j=0; j<O[i].size(); j++) {
            if(E[i][j] != 0) {
                double d = O[i][j]-E[i][j];
                chisq += d * d / E[i][j];
            }
        }
    }
    
    return chisq;
}

double chisq(const vector<vector<int> > & O)
{
    if (O.size() > 0 ) {
        vector<vector<double> > E = getExpectedTable(O);
        return chisq(O, E);
    } else {
        return 0.0;
    }
}

vector<vector<double> > getExpectedTable(const vector<vector<int> > & contingencytable)
{
    assert(contingencytable.size()>0);
    vector<vector<double> > expectedtable(contingencytable.size(),
                                          vector<double>(contingencytable[0].size(), 0)
                                          );
    
    vector<int> rowsums = getRowSums(contingencytable);
    vector<int> colsums = getColSums(contingencytable);
    // int totalcounts = getTotalSum(contingencytable);
    int total = accumulate(rowsums.begin(), rowsums.end(), 0);
    for(size_t i=0; i<expectedtable.size(); i++)
    {
        for(size_t j=0; j<expectedtable[i].size(); j++)
        {
            expectedtable[i][j] = (double)rowsums[i]*(double)colsums[j]/(double)total;
        }
    }
    return expectedtable;
}

void ChisquareTest1DNoPValue(const vector<int> & x_obs, 
                             const vector<double> & p_null, 
                             int K, double & chisq, 
                             size_t & df)
{
	int N = 0;

	chisq = 0;
    df = K-1;

	for(int k=0; k<K; k++) {
		N += x_obs[k];
	}
    
	if(N <= 0) {
		return;
	}
    
	for(int k=0; k<K; k++) {
		double x_exp = N * p_null[k];
		if(x_exp != 0) {
			chisq += (x_obs[k] - x_exp)*(x_obs[k] - x_exp)/ x_exp;
		} else if(x_obs[k] != 0) {
			cerr << "ERROR: expected is zero, but observed is not. Impossible!" << endl;
			exit(EXIT_FAILURE); 
		}
	}
    
	//legacy code
	//if(N % (int) pow(7.0, 3) == 0) {  // Quick hack to adjust for duplicated sample
	//	chisquare /= pow(7.0, 3);
	//	// cout << ".";
	//}
}

double ChisquareTest(const vector<int> & x_obs, const vector<double> & p_null, int K, double & chisq, size_t & df)
{
    ChisquareTest1DNoPValue(x_obs, p_null, K, chisq, df);
	// return qchisq(K-1, chisq);
    return ChisqPvalue(chisq, K-1);
}

double ChisquareTestPrior2011_10_01(const vector< vector<int> > & table_obs, double & chisquare, size_t & df)
{
	int nrow = (int) table_obs.size();
	int K = (int) table_obs[0].size();
	int N = 0;
	vector<double> p_null(K); 

	for(int j = 0; j < K; j++){

		for(int i = 0; i < nrow; i++){

			p_null[j] += table_obs[i][j];

		}//end for
		N += (int) p_null[j];

	}//end for

	//Divide the total in each column by the total sum of each column for null hypothesis
	for(int j = 0; j < K; j++){
		p_null[j] /= N;
	}//end for
	
	chisquare = 0;
	for(int i=0; i<nrow; i++) {
		double chisq_row = 0;
		ChisquareTest1DNoPValue(table_obs[i], p_null, K, chisq_row, df);
		chisquare += chisq_row;
	}

	df = (nrow-1) * (K-1);

	//return 1 - pchisq((int) df, chisquare);  // Joe Song. September 5, 2008. Changed nrow to (nrow-1)
	// return qchisq((int) df, chisquare);  // Yang Zhang. 2/23/2010
    return ChisqPvalue(chisquare, (int) df);
}

double ChisquareTest(const vector< vector<int> > & table_obs, double & chisq, size_t & df,
                     const vector< vector<double> > & null_prob)
{
	int nrow = (int) table_obs.size();
	int K = (int) table_obs[0].size();
	
    int N = 0;
    
	vector<int> row_sum(nrow), col_sum(K); 

    for(int i = 0; i < nrow; i++) {
    	for(int j = 0; j < K; j++) {
            
			col_sum[j] += table_obs[i][j];
            row_sum[i] += table_obs[i][j]; 
            
		}
		N += row_sum[i];
	}
    	
	chisq = 0;
	
    if(N > 0) {
        for(int i = 0; i < nrow; i++) {
            for(int j = 0; j < K; j++) {
                
                double Eij;
                if (! null_prob.empty()) {
                    Eij = N * null_prob[i][j];
                } else {
                    Eij = row_sum[i] * col_sum[j] / static_cast<double>(N);
                }
                
                if ( Eij > 0 ) {
                    double d = table_obs[i][j] - Eij;
                    
                    chisq += d * d / Eij;
                }                
            }
        }
    }
    
    if( ! null_prob.empty() ) {
        df = nrow * K - 1;  // MS 6/15/2013
    } else {
        df = (nrow-1) * (K-1);
    }
	// return qchisq((int) df, chisq);
    return ChisqPvalue(chisq, (int) df);
}


double GTest(const vector<int> & x_obs, const vector<double> & p_null, int K, double & g, size_t & df)
{
	int N=0;

	g=0;
	df = K-1;
    
	for(int k=0; k<K; k++) {
		N += x_obs[k];
	}

	if(N <= 0) {
		return 1;
	}

	for(int k=0; k<K; k++) {
		double x_exp = N * p_null[k];
		if(x_exp != 0) {
			if(x_obs[k] != 0)
				g += x_obs[k] * log(x_obs[k] / x_exp);
		} else if(x_obs[k] != 0) {
			cerr << "ERROR: expected is zero, but observed is not. Impossible!" << endl;
			exit(EXIT_FAILURE); 
		}
	}
	
	g *= 2;

	// return qchisq(K-1, g);
    return ChisqPvalue(g, K-1);
}

double GTest(const vector< vector<int> > & table_obs, double &g, size_t & df)
{
	int nrow = (int) table_obs.size();
	int K = (int) table_obs[0].size();
	int N = 0;
	vector<double> p_null(K); 

	for(int j = 0; j < K; j++){

		for(int i = 0; i < nrow; i++){

			p_null[j] += table_obs[i][j];

		}//end for
		N += (int) p_null[j];

	}//end for

	if(N > 0) {
        //Divide the total in each column by the total sum of each column for null hypothesis
        for(int j = 0; j < K; j++){
            p_null[j] /= N;
        }//end for
	}
    
	g = 0;
	for(int i=0; i<nrow; i++) {
		double g_row = 0;
		GTest(table_obs[i], p_null, K, g_row, df);
		g += g_row;
	}

	df = (nrow-1) * (K-1);

	//return 1 - pchisq((int) df, chisquare);  // Joe Song. September 5, 2008. Changed nrow to (nrow-1)
	// return qchisq((int) df, g);  // Yang Zhang. 2/23/2010
    return ChisqPvalue(g, (int) df);
}

double hTest(const vector<vector<int> >& table_obs, double &condEn, size_t &df)
{
	int nrow = (int) table_obs.size();
	int K = (int) table_obs[0].size();
	int N = 0;
	vector< vector<double> > table_prob(nrow, vector<double>(K) );
	vector< double > p_x(nrow);

	// change table values to probability
	for(int j = 0; j < K; j++)
		for(int i = 0; i < nrow; i++)
			N += table_obs[i][j];

    if( N > 0 ) {
        for(int j = 0; j < K; j++)
            for(int i = 0; i < nrow; i++)
                table_prob[i][j] = (double)table_obs[i][j] / N;
    }

	// calculate conditional entropy ( H(Y|X) = sum_x,y p(x,y)log(p(x)/p(x,y)) )
	for(int x = 0; x < nrow; x++)
		for(int y = 0; y < K; y++)
			p_x[x] += table_prob[x][y];

	for(int i = 0; i < nrow; i++) {
		for(int j = 0; j < K; j++) {
			if(p_x[i] != 0 || table_prob[i][j] != 0) {
                if( table_prob[i][j] > 0 ) {
                    condEn += table_prob[i][j] * log(p_x[i] / table_prob[i][j]);
                }
            }
        }
    }

    double h = K > 0 ? condEn / log(static_cast<double>(K)) : 0;
    
	return h;
}

void computeColumnProb(const vector< vector<int> > & trans_table, vector<double> & null_hypoth) 
{
	int base = (int)null_hypoth.size();

	int total = 0;

	/****************************************
	*Get the total for each row
	*****************************************/
	for(int j = 0; j < base; j++){

		int total_in_col = 0;

		for(size_t i = 0; i < trans_table.size(); i++){

			total_in_col = total_in_col + trans_table[i][j];	//total is the total sum of the current column
		}//end for

		null_hypoth[j] = total_in_col;				//Store the total into null hypothesis
		total = total + total_in_col;			//total all is for total number in all column
	}//end for

	/*****************************************************************************************
	*Divide the total in each column by the total sum of each column for null hypothesis
	******************************************************************************************/

	/** Test this out!! make sure it is a correct double answer **/
	//Divide the total in each column by the total sum of each column for null hypothesis
	if (total > 0) {
        for(int i = 0; i < base; i++){
            null_hypoth[i] = null_hypoth[i]/total;
        }//end for
    }
}

//modified by Yang Zhang 1.17.2013 to add fisher exact test mode 8 and 9
double ChisqTest(const vector< vector<int> > & trans_table, int base, int mode, double & chisq_max, size_t & df, const vector< vector<double> > & null_prob)
{  // Joe Song

	// int total_in_col;		//The total you get in each column
	// int total_all=0;
	double p_value=1.0;

	vector<double> null_hypoth(base);

	double chisq; // Joe Song
	chisq_max = 0; // Joe Song

    vector<vector<double> > expected = getExpectedTable(trans_table);    
	/*******************************************************
	*This loop will call MultinomialChisquareTest
	********************************************************/

	switch(mode) {
		case 1:  // pick the smallest p-value or the largest chisq as the table p-value
				// p-value = min (p_row)

			computeColumnProb(trans_table, null_hypoth);

			for(size_t i = 0; i < trans_table.size(); i++){
				double p_row = ChisquareTest(trans_table[i], null_hypoth, base, chisq, df);
				// Joe Song if(temp < p_value){ 
				if(p_row < p_value || chisq > chisq_max) { // Joe Song
					p_value = p_row;
					chisq_max = chisq; // Joe Song
				}//end if
			}//end for

			break;

		case 2:  //  p-value = 1 - prod ( 1 - p_row )

			computeColumnProb(trans_table, null_hypoth);

			for(size_t i = 0; i < trans_table.size(); i++){

				double p_row = ChisquareTest(trans_table[i], null_hypoth, base, chisq, df);
				p_value = p_value * (1 - p_row);

			}//end for

			p_value = 1 - p_value;

			break;

		case 3:
		case 4:
            // p_value = MultinomialChisquareTestTable(trans_table, null_hypoth, base, chisq);
			p_value = ChisquareTest(trans_table, chisq, df, null_prob);
			chisq_max = chisq;

			break;

		case 5: // Joe Song. Case added on October 1, 2010
			p_value = ChisqDirTest(trans_table, chisq, df);
			chisq_max = chisq;
			
			break;

        case 10: // Joe Song. Case added on Feb 28, 2015
            p_value = ChisqDirTest(trans_table, chisq, df, "normalized");
            chisq_max = chisq;
            
            break;
            
		case 6: // Tyler Hunt. Case added on September 21, 2011
			p_value = GTest(trans_table, chisq, df);
			chisq_max = chisq;

			break;

		case 7:
			p_value = hTest(trans_table, chisq, df);
			chisq_max = chisq;
			break;

		case 8://fisher exact test using fisher probability
			ChisquareTest(trans_table, chisq, df, null_prob); //still save chisquare stats but use fisher's p-value
			chisq_max = chisq;
			// p_value = fisher_exact_test(trans_table);
            p_value = fisher_exact_test(trans_table, "fisher");
			break;
            
        case 9://fisher exact test using chisq stats
			ChisquareTest(trans_table, chisq, df, null_prob); //still save chisquare stats but use fisher's p-value
			chisq_max = chisq;
            // p_value = fisher_exact_test(trans_table, expected);
            p_value = fisher_exact_test(trans_table, "chisq");
			break;
            
        case 15://Add by Hua Apr 14 201m4
            ChisqDirTest(trans_table, chisq, df);
            chisq_max = chisq;
            p_value = exact_functional_test(trans_table, "chisq");
            break;
            
//        case 11://Add by Hua Apr 18 2014, multinomial distribution of exact functional test
//            ChisqDirTest(trans_table, chisq, df);
//            chisq_max = chisq;
//            p_value = multinomial_exact_test(trans_table, "chisq");
//            break;
            
		default:
			cerr << "ERROR: unrecognized p_value mode!\n" << endl;

	}

	return p_value;
}

////
//Add by Hua Apr 14 2014, for -M test to do simulation study on contingency table directly.
#include "TransitionTable.h"
void applyFunctionalChisqTest(const string & file, int pValMode)
{
    TransitionTable tt;
    
    tt.scan4by2(file);
    double chisq;
    size_t df=tt.getDf();
    double p_value = 1.0;
    int base=tt.getBase();
    double chisq_max=0.0;
    
    //Hua Mar 22 2014
    vector< vector<double> > null(0,vector<double>(0,0.0));
    
    clock_t start;
    double duration;
    
    start = clock();
    //////
    
    switch (pValMode) {
        case 3:
        {
            vector< vector<double> > null(0,vector<double>(0,0.0));
            
            p_value = ChisqTest(tt.getTransitionTable(), base, pValMode, chisq_max, df, null);
        }
            break;
            
        case 5:
            p_value = ChisqDirTest(tt.getTransitionTable(), chisq, df);
			break;
            
        case 8:
            p_value = fisher_exact_test(tt.getTransitionTable(), "fisher");
			break;
        case 15:
        {
			ChisqDirTest(tt.getTransitionTable(), chisq, df);//Still output the FunChisq
            chisq_max = chisq;
            p_value = exact_functional_test(tt.getTransitionTable(), "chisq");
            break;
        }
//        case 11://Add by Hua Apr 18 2014, multinomial distribution of exact functional test
//        {
//			ChisqDirTest(tt.getTransitionTable(), chisq, df);//Still output the FunChisq
//            chisq_max = chisq;
//            p_value = multinomial_exact_test(tt.getTransitionTable(), "chisq");
//            break;
//        }
        default:
            cerr<<"Error:ChisqTest.cpp:723!"<<endl;
            break;
    }
    
    //Hua Mar 22 2014
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    
//cout<<duration<<endl;
    //////
    
    cout << "pd\tchisqd\tdf" << endl;
    cout << p_value << "\t" << tt.getChisq()
    << "\t" << tt.getDf() << endl;
}
/////////

void test_Homo() 
{
	// Test HomoChisquareTest()
}

void test_ChisqPvalue()
{
    double chisqs[] = {9.6, 16, 16, 1000, 511};
    int dfs[] = {511, 1023, 511, 511, 511};
    
    for(int i=0; i<5; i++) {
        double pv = ChisqPvalue(chisqs[i], dfs[i]);
        cout << "chisq=" << chisqs[i] << ", df=" << dfs[i] << ", p-value=" << pv << endl;
    }
}

void test_Multinomial() 
{
    cout << ">>>> Testing ChisqPvalue() and ChisquareTest() ..." << endl;
    
	// int x_obs[][4] = {{1, 3, 2, 4}, {10, 3, 20, 4}};
	// double p_null[] = {0.2, 0.3, 0.4, 0.1};

	vector< vector<int> > x_obs(2, vector<int>(4));

	x_obs[0][0] = 1, x_obs[0][1] = 3, x_obs[0][2] = 2, x_obs[0][3] = 4; 
	x_obs[1][0] = 10, x_obs[1][1] = 3, x_obs[1][2] = 20, x_obs[1][3] = 4; 
		
	vector<double> p_null(4);
	p_null[0] = 0.2, p_null[1] = 0.3, p_null[2] = 0.4, p_null[3] = 0.1;


	int K=4;
	double p_value[] = {0.01476090, 0.03392869};

	double chisq;
	size_t df=0;

	cout << "Chisquare accumulative probability with" << endl
		<< "  4 degrees of freedom from x=3 to infinity: " 
		<< ChisqPvalue(3, 4) << " (0.557825 expected)" << endl;

	cout << "Multinomial chisquare test: " << endl;
	for(int m=0; m<2; m++) {
        
        double p = ChisquareTest(x_obs[m], p_null, K, chisq, df);
		cout << p << " (" << p_value[m] << " expected)" << endl;
        if(abs((p - p_value[m])/p_value[m]) > 1e-6) {
            cerr << "ERROR: The relative difference between the computed p-value and the correct one is greater than 1e-6!" << endl;
            exit(EXIT_FAILURE);
        }
	}

	cout << "Done." << endl << endl;
}

