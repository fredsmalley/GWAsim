//
//  ComparativeChisq.cpp
//  gln
//
//  Created by Joe Song on 1/24/13.
//
//  Modified:
//
//  May 3, 2014. (MS) Amended KSampleChisqTest() to prevent a divide-by-zero error.


#include <iostream>
#include <cmath>
#include <cassert>
#include <set>

using namespace std;

#include "ChisqTests.h"
#include "ComparativeChisq.h"
#include "StatDistributions.h"
#include "SetOps.h"

inline size_t sizeOfSuperset(const vector<vector<int> > & subsets)
{ // MS added Nov 20, 2011
    
    set<int> superset;
    for(size_t i=0; i < subsets.size(); i++) {
        superset.insert( subsets[i].begin(), subsets[i].end() );
    }
    return superset.size();
}

vector<TransitionTable> fillIndividualTables
(int child, const vector< vector<int> > & parents,
 const vector< vector<int> > & delays,
 const vector<TrajectoryCollection> & trajCols)
{
    vector<TransitionTable> tts(trajCols.size());
    
    tts.resize(trajCols.size());
    for(size_t k=0; k<parents.size(); k++) {
        if(parents[k].size()>0) {
            tts[k].initialize(child, parents[k], delays[k], trajCols[k]);
            tts[k].fill(trajCols[k]);
        }
    }
    return tts;
}

const vector<TransitionTable> &
computeIndividualChisq(vector<TransitionTable> & tts, int pValMode,
                       const TransitionTable & ttHom)
/*vector<TransitionTable> computeIndividualChisq
(int child, const vector< vector<int> > & parents,
 const vector< vector<int> > & delays,
 const vector<TrajectoryCollection> & trajCols,
 int pValMode, const TransitionTable & ttHom)*/
{
    // vector<TransitionTable> tts(trajCols.size());
    
    vector< vector<double> > null_prob;
    
    if (ttHom.nrow() > 0 && ttHom.ncol() > 0) {
        
        double total = ttHom.getTotalSum();
        if (total > 0) {
            
            null_prob = getExpectedTable(ttHom.getTransitionTable());
            
            for (size_t r=0; r<null_prob.size(); r++) {
                for (size_t c=0; c<null_prob[r].size(); c++) {
                    null_prob[r][c] /= total;
                }
            }
            
        } else {
            null_prob.clear();
        }
    }

    for(size_t k=0; k<tts.size(); k++) {
        if(tts[k].nrow()>0) {
            tts[k].applyChisquareTest(pValMode, null_prob);
            
            /*
            if( ! null_prob.empty() ) { // MS June 15, 2013. Increase degrees of
                // freedom if null probabilies are not estimated from contingency
                // table tts[k]
                tts[k].setDf( tts[k].getDf() + tts[k].nrow() - 1
                             + tts[k].ncol() - 1 );
            }
            */
        }
    }

    /*
    tts.resize(trajCols.size());
    for(size_t k=0; k<parents.size(); k++) {
        if(parents[k].size()>0) {
            tts[k].initialize(child, parents[k], delays[k], trajCols[k]);
            tts[k].fill(trajCols[k]);
            tts[k].applyChisquareTest(pValMode, null_prob);
        }
    }
    */
    
    return tts;
}

TransitionTable computeTotalChisq(const vector<TransitionTable> & tts, int pValMode)
{    
    TransitionTable ttTot;
    
    double chisq_t=0;
    int df_t=0;
    double p_t;
        
    double q = 1.0; //p-value for fisher test
    for(size_t k=0; k<tts.size(); k++) {
        if(tts[k].nrow() >0) {
			q *= (1-tts[k].getpValue());
            df_t += tts[k].getDf();
            chisq_t += tts[k].getChisq();
        }
    }
    
    switch (pValMode) {
        case 8:
        case 9:
            p_t = exact_total_strength_test(tts, pValMode);
            break;
            
        default:
            p_t = ChisqPvalue(chisq_t, df_t);
            break;
    }
    
    // MS added Nov 20, 2011.  Still under testing
    // df_t = (sizeOfSuperset(parents) - 1) * (tts[0].getBase() - 1);
    ttTot.setChisq( chisq_t );
    ttTot.setDf( df_t );
    ttTot.setpValue( p_t );
    
    return ttTot;
}

//modified by Yang Zhang 2.14.2014
//add parent information for ttHet
//only useful when parents are the same
//modified by Yang Zhang 6.3.2014
//CPFunX2 computes heterogeneity by adding row-wise FunChisqs
//modified by Yang Zhang 7.2.2014
//Allow differential parents by compute CPFunX2 on pooled table
TransitionTable computeHetChisq(const vector<TransitionTable> & tts,
	const TransitionTable & ttTot,
	const TransitionTable & ttHom,
	int pValMode, const CMPX_MARGINAL & marginal)
{ // Heterogeneity test

	TransitionTable ttHet = ttHom;

	int df_d = 1;
	double p_d;	// p-value for chisq_d
	double chisq_d = 0.0;

	double associationChisq = ttTot.getChisq();

	switch (marginal) {
	case POOLED_NO_ROWCOL:
	{
							 size_t df_row, df_col;
							 double p_z;
							 double rowChisq = applyKSampleRowChisqTest(tts, df_row, p_z);
							 double colChisq = applyKSampleColChisqTest(tts, df_col, p_z);
							 associationChisq -= rowChisq + colChisq;
							 df_d = ttTot.getDf() - df_row - df_col - ttHom.getDf();
							 break;
	}
	case CPFUNX2:
	{
					for (size_t i = 0; i < ttHom.nrow(); i++)
					{
						double chisq;
						size_t df;
						vector<vector<int> > tt1row;
						for (size_t j = 0; j < tts.size(); j++)
						{
							tt1row.push_back(tts[j].getTransitionTable()[i]);
						}
						ChisqDirTest(tt1row, chisq, df, "");
						chisq_d += chisq;
						df_d += df;
					}
	}
	default:
		df_d = ttTot.getDf() - ttHom.getDf(); // degrees of freedom for chisq_d
		break;
	}

	if (marginal != CPFUNX2)
	{
		chisq_d = HeteroChisqTest(associationChisq, ttHom.getChisq(), df_d, p_d);
	}
	else
	{
		p_d = ChisqPvalue(chisq_d, (int)df_d);
	}

	switch (pValMode) {
	case 8:
	case 9:
		p_d = exact_heterogeneity_test(tts, pValMode, marginal);
		break;

	default:
		break;
	}

	ttHet.setpValue(p_d);
	ttHet.setChisq(chisq_d);
	ttHet.setDf(df_d);
	ttHet.setParents(ttHom.getParents());

	return ttHet;
}

TransitionTable fillPooledTable(int child, const vector< int > & parents,
                                const vector< int > & delays,
                                const vector<TrajectoryCollection> & trajCols)
{
    TransitionTable ttPooled(child, parents, delays, trajCols[0]);
    
    if(parents.size()>0) {
        ttPooled.fill(trajCols[0]);
    }
    
    for(size_t k=1; k<trajCols.size(); k++) {
        TransitionTable	tt_k(child, parents, delays, trajCols[k]);
        if(parents.size()>0) {
            tt_k.fill(trajCols[k]);
        }
        ttPooled  += tt_k;
    }
    
    return ttPooled;
}

TransitionTable pool(const vector<TransitionTable> & tts)
// Compute pooled transition table and transition table with unique parents
//   for each condition
// OUTPUT:
//   pooled transition table
// ASSUMPTIONS:
//   1. The same parent ID will have the same delay!
//      Needs to be fixed!!!!!!!!!!!
{
    TransitionTable pooled;
    
    vector<int> intersecParents = tts[0].getParents();
    for(size_t i=1; i<tts.size(); i++) {
        if(intersecParents.empty()) {
            break;
        }
        intersecParents = intersectionSet(intersecParents, tts[i].getParents());
    }
    
    for(size_t i=0; i<tts.size(); i++) {
        pooled += tts[i].subtable(intersecParents);
    }
    
    return pooled;
}

TransitionTable & computeHomChisq(TransitionTable & ttPooled, int pValueMode)
{
    if(ttPooled.getParents().size() == 0) {
        
        ttPooled.setChisq(0);
        ttPooled.setDf(0);
        ttPooled.setpValue(1);
        
    } else {
        ttPooled.applyChisquareTest(pValueMode);
    }
    
    return ttPooled;
}

/*
TransitionTable computeHomChisq(int child, const vector< int > & parents,
                                const vector< int > & delays,
                                const vector<TrajectoryCollection> & trajCols,
                                int pValueMode)
{
    TransitionTable ttPooled(child, parents, delays, trajCols[0]);
    
    if(parents.size()>0) {
        ttPooled.fill(trajCols[0]);
    }
    
    for(size_t k=1; k<trajCols.size(); k++) {
        TransitionTable	tt_k(child, parents, delays, trajCols[k]);
        if(parents.size()>0) {
            tt_k.fill(trajCols[k]);
        }
        ttPooled  += tt_k;
    }
    
    if(ttPooled.getParents().size() == 0) {
        
        ttPooled.setChisq(0);
        ttPooled.setDf(0);
        ttPooled.setpValue(1);
        
    } else {
        ttPooled.applyChisquareTest(pValueMode);
    }
    
    return ttPooled;
}
*/

double HeteroChisqTest(double sumChisq, double pooledChisq, size_t df, double & p_value)
{
	// double heteroChisq = abs(sumChisq - pooledChisq); MS 2/14/2013
    double heteroChisq = sumChisq - pooledChisq;
    
	// p_value = qchisq((int)df, heteroChisq);
	
    p_value = ChisqPvalue(abs(heteroChisq), (int)df);
    
	return heteroChisq;
}

double HeteroChisqTest(const vector< vector<int> > & table1_obs,
                       const vector< vector<int> > & table2_obs,
                       size_t & df, double & p_value)
{
	vector< vector<int> > tablePooledObs(table1_obs.size());
	
	for(int rowindex=0; rowindex < (int) table1_obs.size(); rowindex++)
	{
		tablePooledObs[rowindex].resize(table1_obs[0].size());
	}
    
	for(int r=0; r<(int)tablePooledObs.size(); r++) {
		for(int c=0; c<(int)table1_obs[0].size(); c++) {
			tablePooledObs[r][c] = table1_obs[r][c] + table2_obs[r][c];
		}
	}
    
	double chisq1, chisq2, chisqPooled;
    
	int K = (int) table1_obs[0].size();
    
	ChisquareTest(table1_obs, chisq1, df);
	ChisquareTest(table2_obs, chisq2, df);
	ChisquareTest(tablePooledObs, chisqPooled, df);
	
	// chisquare = abs(chisq1 + chisq2 - chisqPooled);
    
	df = (table1_obs.size() - 1) * (K-1) + (table2_obs.size() - 1) * (K-1)
    - (tablePooledObs.size() - 1) * (K-1);
    
	// double p_value = 1 - pchisq((int)df, chisquare);
	double heteroChisq = HeteroChisqTest(chisq1 + chisq2, chisqPooled, df, p_value);
    
	return heteroChisq;
}

void test_Hetero()
{
	// Test HeteroChisquareTest()
	cout << ">>>> Testing HeteroChisquareTest() ... " << endl;
    
	// int x_obs[][4] = {{1, 3, 2, 4}, {10, 3, 20, 4}};
	// double p_null[] = {0.2, 0.3, 0.4, 0.1};
    
	//vector< vector<int> > x1_obs ( 2 );
	// x1_obs[0].resize ( 2 ), x1_obs[1].resize ( 2 );
	//        vector< vector<int> > x1_obs (4);
    
    
	//      cout << "  elment capacity " << (x1_obs.at(i)).capacity() << endl;
	//x1_obs[0][1] = 1, x1_obs[0][2] = 3, x1_obs[0][3] = 2, x1_obs[0][4] = 4;
	//x1_obs[1][1] = 10, x1_obs[1][2] = 3, x1_obs[1][3] = 20, x1_obs[1][4] = 4;
	//vector< vector<int> > x1_obs (2,4);
	//x1_obs[0][0] = 1, x1_obs[0][1] = 3, x1_obs[0][2] = 2, x1_obs[0][3] = 4;
	//x1_obs[1][0] = 10, x1_obs[1][1] = 2, x1_obs[1][2] = 20, x1_obs[1][3] = 4;
	vector< vector<int> > x1_obs (2);
	for(int rowindex=0;rowindex<2;rowindex++)
	{
		x1_obs[rowindex].resize(2);
	}
    
	x1_obs[0][0] = 4, x1_obs[0][1] = 7;
	x1_obs[1][0] = 25, x1_obs[1][1] =13;
    
	//vector< vector<int> > x2_obs ( 2 );
	//x2_obs[0].resize ( 2 ), x2_obs[1].resize ( 2 );
    
	//vector< vector<int> > x2_obs ( 4 );
	//x2_obs[0].resize ( 4 ), x2_obs[1].resize ( 4 );
    
	//x2_obs[0][1] = 1, x2_obs[0][2] = 3, x2_obs[0][3] = 2, x2_obs[0][4] = 4;
	//x2_obs[1][1] = 10, x2_obs[1][2] = 3, x2_obs[1][3] = 20, x2_obs[1][4] = 4;
	//vector< vector<int> > x2_obs (2,4);
	//x2_obs[0][0] = 1, x2_obs[0][1] = 3, x2_obs[0][2] = 2, x2_obs[0][3] = 4;
	//x2_obs[1][0] = 10, x2_obs[1][1] = 3, x2_obs[1][2] = 20, x2_obs[1][3] = 4;
	vector< vector<int> > x2_obs (2);
	for(int rowindex=0;rowindex<2;rowindex++)
	{
		x2_obs[rowindex].resize(2);
	}
    
	x2_obs[0][0] = 7, x2_obs[0][1] = 10;
	x2_obs[1][0] = 27, x2_obs[1][1] =18;
	
	// double chisquare1=0;
	// double chisquare2=0;
    
	double chisquare, pval;
	size_t df;
    
	chisquare = HeteroChisqTest(x1_obs, x2_obs, df, pval);
    
	cout << "The chi square value is " << chisquare << endl;
	cout << "The p-value is " << pval << endl << endl;
    
}

double TwoSampleChisqTest(const vector<int> & freq1, const vector<int> & freq2,
                          size_t & df, double & p_value)
{
	double chisq=0;
    
	df = freq1.size() - 1;
    
	int n1=0, n2=0;
    
	assert(freq1.size() == freq2.size());
	
	if(freq1.size() != freq2.size()) {
		cerr << "Distributions to be compared has unequal number of bins!\n";
		exit(EXIT_FAILURE);
	}
    
	for(int i=0; i<(int) freq1.size(); i++) {
		n1 += freq1[i];
		n2 += freq2[i];
	}
    
	double sqrt_n1_over_n2 = sqrt(n1 / (double) n2);
    
	for (int j=0; j < (int) freq1.size(); j++) {
		if (freq1[j] == 0.0 && freq2[j] == 0.0) {
			--df;
            
		} else {
            
			double temp = freq1[j] / sqrt_n1_over_n2 - freq2[j] * sqrt_n1_over_n2;
            
			chisq += temp * temp / (freq1[j] + freq2[j]);
            
		}
        
	}
    
	// cout << "df=" << df << ", " << "chisq=" << chisq << endl;
    
	// p_value = qchisq((int)df, chisq);
    
	p_value = ChisqPvalue(chisq, (int)df);
    
	return chisq;
}

void test_TwoSampleChisqTest()
{
	cout << ">>>> Testing TwoSampleChisqTest() ..." << endl;
	
	vector<int> f1(3);
	vector<int> f2(3);
    
	//f1[0] = 3, f1[1] = 4, f1[2] = 1;
	//f2[0] = 3*2, f2[1] = 4*2, f2[2] = 1*2;
    
	f1[0] = 3, f1[1] = 4, f1[2] = 1;
	f2[0] = 3*2-5, f2[1] = 4*2-3, f2[2] = 1*2+10;
    
	cout << "Sample 1: ";
	for(int i=0; i < (int) f1.size(); i++) {
		cout << f1[i] << "\t";
	}
	cout << endl;
    
	cout << "Sample 2: ";
	for(int i=0; i < (int) f2.size(); i++) {
		cout << f2[i] << "\t";
	}
	cout << endl;
    
	double chisq, pval;
	size_t df;
    
	chisq = TwoSampleChisqTest(f2, f1, df, pval);
    
	cout << "The chi square value is " << chisq << endl;
	cout << "The degrees of freedom is " << df << endl;
	cout << "The p-value is " << pval << endl << endl;
    
}

double KSampleChisqTest(const vector< vector<int> > & freq,
                        size_t & df, double & p_value)
{
	size_t K = freq.size(); // number of samples K
	size_t nB = freq[0].size();  // number of bins
	
	vector<int> sumObs(nB,0);
	vector<int> n(K,0);
	int N = 0;
    
	for(size_t k=0; k<K; k++) {
        
		if(freq[k].size() != nB) {
			cerr << "Distributions to be compared has unequal number of bins!\n";
			exit(EXIT_FAILURE);
		}
        
		for(size_t b=0; b<nB; b++) {
			sumObs[b] += freq[k][b];
			n[k] += freq[k][b];
			N += freq[k][b];
		}
	}
    
	double chisq=0;
    
	df=(nB-1)*K; // MS 6/15/2013 Replace (K-1) by K; // needs to be adjusted for
    // empty row and empty columns in freq
    
	for(size_t b=0; b<nB; b++) {
        
		if(sumObs[b] == 0) {
			continue;
		}
		
		double expectedRatio = sumObs[b]/(double) N;
        
		for (size_t k=0; k<K; k++) {
            
			double expectedCount = expectedRatio * n[k];
            
			double d = freq[k][b] - expectedCount;
            
            if(expectedCount > 0) { // Condition added by MS on May 5, 2014
                chisq += d * d / expectedCount;
            }
            
		}
        
	}
    
	// cout << "df=" << df << ", " << "chisq=" << chisq << endl;
    
	// p_value = qchisq((int)df, chisq);
    
	p_value = ChisqPvalue(chisq, (int)df);
    
	return chisq;
}

void test_KSampleChisqTest()
{
	cout << ">>>> Testing KSampleChisqTest() ..." << endl;
	
	vector<int> f1(3);
	vector<int> f2(3);
    
	//f1[0] = 3, f1[1] = 4, f1[2] = 1;
	//f2[0] = 3*2, f2[1] = 4*2, f2[2] = 1*2;
    
	f1[0] = 3, f1[1] = 4, f1[2] = 1;
	f2[0] = 3*2-5, f2[1] = 4*2-3, f2[2] = 1*2+10;
    
	cout << "Sample 1: ";
	for(int i=0; i<(int) f1.size(); i++) {
		cout << f1[i] << "\t";
	}
	cout << endl;
    
	cout << "Sample 2: ";
	for(int i=0; i< (int) f2.size(); i++) {
		cout << f2[i] << "\t";
	}
	cout << endl;
    
	vector< vector<int> > f(2);
	f[0] = f1;
	f[1] = f2;
    
	double chisq, pval;
	size_t df;
    
	chisq = KSampleChisqTest(f, df, pval);
    
	cout << "The chi square value is " << chisq << endl;
	cout << "The degrees of freedom is " << df << endl;
	cout << "The p-value is " << pval << endl << endl;
    
}

TransitionTable applyHeteroChisqTest(const TransitionTable & tt1,
                                     const TransitionTable & tt2)
{
	TransitionTable tt = tt1 - tt2;
    
	double heteroChisq;
	double pValue;
	size_t df;
    
	heteroChisq = HeteroChisqTest(tt1.getTransitionTable(),
                                  tt2.getTransitionTable(),
                                  df, pValue);
    
	tt.setpValue(pValue);
	tt.setChisq(heteroChisq);
	tt.setDf((int) df);
    
	return tt;
}

void applyHeteroChisqTest(const vector<string> & files, int pValueMode,
                          const CMPX_MARGINAL & marginal)
{
    vector<TransitionTable> tts(files.size());
    
    for (size_t k=0; k<files.size(); k++) {
        tts[k].scan4by2(files[k]);
    }
    
    TransitionTable ttHom, ttHet, ttTot;
    
    ttHom = pool(tts);
    computeHomChisq(ttHom, pValueMode);
    
    // Compute invdividual chisqs using expected value from *pooled* marginals
    switch (marginal) {
        case POOLED:
        case POOLED_NO_ROWCOL:
            computeIndividualChisq(tts, pValueMode, ttHom);
            break;
            
        default:
            computeIndividualChisq(tts, pValueMode);
            break;
    }
    
    // Total strength of all interactions
    ttTot = computeTotalChisq(tts, pValueMode);
    
    // statistic for heterogeneity
    ttHet = computeHetChisq(tts, ttTot, ttHom, pValueMode, marginal);
    
    
    cout << "pd\tchisqd\tdf" << endl;
    cout << ttHet.getpValue() << "\t" << ttHet.getChisq()
    << "\t" << ttHet.getDf() << endl;
}

double applyKSampleRowChisqTest(const vector<TransitionTable> & tts,
                                size_t & df, double & pValue)
{
	double chisq = 0;
    
	vector< vector<int> > rowCounts(tts.size());
    
	for(size_t i=0; i < tts.size(); i++) {
		rowCounts[i] = tts[i].getRowSums();
	}
    
	chisq = KSampleChisqTest(rowCounts, df, pValue);
    
	return chisq;
}

double applyKSampleColChisqTest(const vector<TransitionTable> & tts,
                                size_t & df, double & pValue)
{
	double chisq = 0;
    
	vector< vector<int> > colCounts(tts.size());
    
	for(size_t i=0; i < tts.size(); i++) {
		colCounts[i] = tts[i].getColSums();
	}

	chisq = KSampleChisqTest(colCounts, df, pValue);
    
	return chisq;
}

void computeParentWorkingZoneChange(const vector<TransitionTable> & tts,
                                    double & chisq_z, size_t & df_z,
                                    double & p_z)
{
    // Compute work zone change statistics
    
    bool sameParents = true; // true if parents are identical
    
    for(size_t i=1; i<tts.size(); i++) {
        if( tts[i].getParents() != tts[0].getParents() ) {
            sameParents = false;
            break;
        }
    }
    
    if(sameParents) {
        
        chisq_z = applyKSampleRowChisqTest(tts, df_z, p_z);
        
    } else {
        
        // if parents are different always return
        chisq_z = 0.0;
        p_z = 0;
        df_z = 0;
        
        /* instead of
         chisq_z = 0.0;
         p_z = 1.0;
         df_z = 0;
         */
        
    }
}

void computeChildWorkingZoneChange(const vector<TransitionTable> & tts,
                                   double & chisq_z, size_t & df_z,
                                   double & p_z)
{
    chisq_z = applyKSampleColChisqTest(tts, df_z, p_z);
}

/*
 //add pvaluemode to adapt to fisher exact test
 //by Yang Zhang 2013.1.17
 TransitionTable applyHeteroChisqTest(vector<TransitionTable> & tts,
 TransitionTable & ttPooled,
 int pValueMode)
 {
 TransitionTable ttDiff;
 
 double sumChisq = 0;
 double pooledChisq;
 double heteroChisq;
 double pValue;
 size_t df;
 
 ttDiff = ttPooled;
 
 for(size_t i=0; i < tts.size(); i++) {
 
 // Apply chisquare test on each table
 tts[i].applyChisquareTest(pValueMode);
 sumChisq += tts[i].getChisq();
 
 // Accumulate ttPooled
 ttPooled += tts[i];
 }
 
 // Compute the deviation from ttPooled and save them in ttDiff
 for(size_t i=0; i < tts.size(); i++) {
 TransitionTable ttd = tts[i] - ttPooled;
 ttDiff += abs( ttd );
 }
 
 // Apply chisquare test on pooled table
 ttPooled.applyChisquareTest(pValueMode);
 pooledChisq = ttPooled.getChisq();
 
 // Compute degrees of freedom
 df = ( tts.size() - 1 ) * (tts[0].m_transitionTable.size() - 1)
 * (tts[0].m_transitionTable[0].size() - 1);
 
 // Apply the heterogeneity chisquare test
 heteroChisq = HeteroChisqTest(sumChisq, pooledChisq, df, pValue);
 
 ttDiff.setpValue(pValue);
 ttDiff.setChisq(heteroChisq);
 ttDiff.setDf((int) df);
 
 if(pValueMode == 8 || pValueMode == 9)
 {
 double pd = 1.0;
 double pc = 1.0;
 double pt = 1.0;
 exact_comparative_chisq_3(tts, pd, pc, pt, pValueMode);
 ttDiff.setpValue(pd);
 }
 
 return ttDiff;
 }
 
 TransitionTable applyHeteroChisqTest(vector<TransitionTable> & tts,
 TransitionTable & ttPooled)
 {
 TransitionTable ttDiff;
 
 double sumChisq = 0;
 double pooledChisq;
 double heteroChisq;
 double pValue;
 size_t df;
 
 ttDiff = ttPooled;
 
 for(size_t i=0; i < tts.size(); i++) {
 
 // Apply chisquare test on each table
 tts[i].applyChisquareTest(3);
 sumChisq += tts[i].getChisq();
 
 // Accumulate ttPooled
 ttPooled += tts[i];
 }
 
 // Compute the deviation from ttPooled and save them in ttDiff
 for(size_t i=0; i < tts.size(); i++) {
 TransitionTable ttd = tts[i] - ttPooled;
 ttDiff += abs( ttd );
 }
 
 // Apply chisquare test on pooled table
 ttPooled.applyChisquareTest(3);
 pooledChisq = ttPooled.getChisq();
 
 // Compute degrees of freedom
 df = ( tts.size() - 1 ) * (tts[0].m_transitionTable.size() - 1)
 * (tts[0].m_transitionTable[0].size() - 1);
 
 // Apply the heterogeneity chisquare test
 heteroChisq = HeteroChisqTest(sumChisq, pooledChisq, df, pValue);
 
 ttDiff.setpValue(pValue);
 ttDiff.setChisq(heteroChisq);
 ttDiff.setDf((int) df);
 
 return ttDiff;
 }
 */
