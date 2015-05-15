//
//  ComparativeChisq.h
//  gln
//
//  Created by Joe Song on 1/24/13.
//
//

#ifndef __gln__ComparativeChisq__
#define __gln__ComparativeChisq__

#include <iostream>
#include "TransitionTable.h"

// enum CMPX_MARGINAL {INDEP, POOLED, STACKED, POOLED_NO_ROWCOL};
enum CMPX_MARGINAL {POOLED, STACKED, POOLED_NO_ROWCOL,
    INDEP_SUB, INDEP_SUP, DIFF_COR_SUP, CPFUNX2};

//void exact_comparative_chisq_1(const vector<TransitionTable> & tts,
//                             double &pd, double &pc, double &pt, int pvaluemode);

void exact_comparative_chisq_2(const vector<TransitionTable> & tts,
                               double &pd, double &pc, double &pt, int pvaluemode);


void exact_comparative_chisq_3(const vector<TransitionTable> & tts,
                               double &pd, double &pc, double &pt,
                               int pvaluemode,
                               const CMPX_MARGINAL & marginal = INDEP_SUB);

double exact_heterogeneity_test(const vector< TransitionTable > & Cs,
                                const string & discrepancy_measure,
                                const CMPX_MARGINAL & marginal = INDEP_SUB);

double exact_heterogeneity_test(const vector<TransitionTable> & tts,
                                int pvaluemode,
                                const CMPX_MARGINAL & marginal = INDEP_SUB);

double exact_total_strength_test(const vector<TransitionTable> & tts,
                                 int pvaluemode);

vector<TransitionTable> fillIndividualTables
(int child, const vector< vector<int> > & parents,
 const vector< vector<int> > & delays,
 const vector<TrajectoryCollection> & trajCols);

const vector<TransitionTable> & computeIndividualChisq
(vector<TransitionTable> & tts, int pValMode,
 const TransitionTable & ttHom = TransitionTable());

TransitionTable computeTotalChisq(const vector<TransitionTable> & tts, int pValMode);

TransitionTable computeHetChisq(const vector<TransitionTable> & tts,
                                const TransitionTable & ttTot,
                                const TransitionTable & ttHom,
                                int pValMode, const CMPX_MARGINAL & marginal);

TransitionTable pool(const vector<TransitionTable> & tts);

TransitionTable fillPooledTable(int child, const vector< int > & parents,
                                const vector< int > & delays,
                                const vector<TrajectoryCollection> & trajCols);

TransitionTable & computeHomChisq(TransitionTable & ttPooled, int pValueMode);


double applyKSampleRowChisqTest(const vector<TransitionTable> & tts,
                                size_t & df, double & pValue);

double applyKSampleColChisqTest(const vector<TransitionTable> & tts,
                                size_t & df, double & pValue);

void computeParentWorkingZoneChange(const vector<TransitionTable> & tts,
                                    double & chisq_z, size_t & df_z,
                                    double & p_z);

void computeChildWorkingZoneChange(const vector<TransitionTable> & tts,
                                   double & chisq_z, size_t & df_z,
                                   double & p_z);

void applyHeteroChisqTest(const vector<string> & files, int pValueMode,
                          const CMPX_MARGINAL & marginal);

double HeteroChisqTest(double sumChisq, double pooledChisq,
                       size_t df, double & p_value);

double HeteroChisqTest(const vector< vector<int> > & table1_obs,
                       const vector< vector<int> > & table2_obs,
                       size_t & df, double & p_value);


#endif /* defined(__gln__ComparativeChisq__) */
