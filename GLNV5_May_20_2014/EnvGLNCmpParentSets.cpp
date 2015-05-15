//
//  EnvGLNCmpParentSets.cpp -- Compare the parent sets and determine if a given set is better
//                             than the current parent set in comparative analysis of
//                             interactions 
//  gln
//
//  Created by Joe Song on 10/1/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//

#include <iostream>
#include <cmath>

#include "EnvGLNCmp.h"
#include "StatDistributions.h"
#include "ComparativeChisq.h"

bool compare_chisq(double p_1, int df_1, double chisquare_1,
				   double p_2, int df_2, double chisquare_2);

bool compare_chisq_imp(double p1, int df1, double chisq1, 
                       double p2, int df2, double chisq2);

bool compare_F_imp_over_chisq(double p1, int df1, double chisq1, 
                              double p2, int df2, double chisq2);

bool compare_F_imp(int df1, double chisq1, int df2, double chisq2);
double compute_F_imp(int df1, double chisq1, int df2, double chisq2);

bool EnvGLNCmp::select(int child, 
                       const vector<TransitionTable> & tts,
                       const TransitionTable & ttTot, 
                       const TransitionTable & ttHom, 
                       const TransitionTable & ttHet)
{ // compare the current parent sets with the best parent sets and decide if 
    // the current one is better. If so, select the current one and update the 
    // best parent sets.
    
    bool selected = false;
    
    switch (m_gc.m_methodCmpxParCmp) {
            
        case BY_BEST_HOM:
            
            selected = select_BY_BEST_HOM(child, tts, ttTot, ttHom);
            
            break;
            
        default:
            
            // something is not correct with this one for the BY_BEST_HOM mode
            CmpIntxType type = labelIntxType(ttTot.getpValue(), ttHom.getpValue(), ttHet.getpValue(), m_gc.m_alpha);
            
            // Compare if the current parent sets have improved over the current best parent sets
            bool improved = compareParentSets(child, type, tts, ttTot, ttHom, ttHet);
            
            if(improved) {
                
                double chisq_parent_z;
                size_t df_parent_z;
                double p_parent_z;
                
                computeParentWorkingZoneChange(tts, chisq_parent_z, df_parent_z, p_parent_z);
                
                updateEnvironment(child, type, tts, ttTot, ttHom, ttHet, 
                                  p_parent_z, (int) df_parent_z, chisq_parent_z);
            }
            selected = improved;
    }
    
    return selected;
}

void EnvGLNCmp::collectStats(int child,
                             const vector<vector< int > > & parents,
                             const vector<vector< int > > & delays,
                             vector<TransitionTable> & tts,
                             TransitionTable & ttTot,
                             TransitionTable & ttHom,
                             TransitionTable & ttHet) const
{
    switch (m_gc.m_methodMarginal) {
        case INDEP_SUB:  // INDEP
            collectStatsIndep(child, parents, delays, tts, ttTot, ttHom, ttHet);
            break;
        case INDEP_SUP:
        case POOLED:
        case POOLED_NO_ROWCOL:
        case STACKED:
			//added by Yang Zhang 6.3.2014
		case CPFUNX2:
            collectStatsPooled(child, parents, delays, tts, ttTot, ttHom, ttHet);
            break;
        case DIFF_COR_SUP:
            collectDiffCorPooled(child, parents, delays, tts, ttTot, ttHom, ttHet);
            break;
        default:
            break;
    }
}

void EnvGLNCmp::collectStatsIndep(int child,
                             const vector<vector< int > > & parents, 
                             const vector<vector< int > > & delays, 
                             vector<TransitionTable> & tts,
                             TransitionTable & ttTot, 
                             TransitionTable & ttHom, 
                             TransitionTable & ttHet) const
{ // Collect statistics related to total strength, homegeneity, and heterogeneity
    // for a given parent sets across condition.
    
    switch (m_gc.m_methodCmpxParCmp) {
            
        case BY_BEST_HOM:
            
            size_t k;
            
            for (k=1; k<parents.size(); k++) {
                if( parents[k] != parents[k-1] )
                    break;
                else if( delays[k] != delays[k-1] )
                    break;
            }
            
            if (k == parents.size()) { // parents are the same across all conditions
                
                // Homoegeneity test
                
                // ttHom = computeHomChisq(child, parents[0], delays[0],
                //                        m_trajCols, m_gc.m_pValueMode);
                ttHom = fillPooledTable(child, parents[0], delays[0], m_trajCols);
                computeHomChisq(ttHom, m_gc.m_pValueMode);
            
            }
            
            break;
            
        default:
            
            // Homoegeneity test
            vector<int> commonParents, commonParentDelays;
            
            // Commnted out by MS Nov 21, 2011 findCommonParents(parents, delays, commonParents, commonParentDelays);
            
            findCommonParentsWithSameDelays(parents, delays, commonParents, commonParentDelays); // Added MS Nov 21, 2011
            
            //ttHom = computeHomChisq(child, commonParents, commonParentDelays,
            //                        m_trajCols, m_gc.m_pValueMode);
            ttHom = fillPooledTable(child, commonParents, commonParentDelays,
                                    m_trajCols);
            computeHomChisq(ttHom, m_gc.m_pValueMode);

            break;
    }

    // Fill in tts using the exact parent set for each condition
    // tts = computeIndividualChisq(child, parents, delays, m_trajCols, m_gc.m_pValueMode);

    tts = fillIndividualTables(child, parents, delays, m_trajCols);
    computeIndividualChisq(tts, m_gc.m_pValueMode);
    
    // Total strength of all interactions
    ttTot = computeTotalChisq(tts, m_gc.m_pValueMode);

    ttHet = computeHetChisq(tts, ttTot, ttHom, m_gc.m_pValueMode,
                            m_gc.m_methodMarginal); // statistic for heterogeneity
    
}

void computerDifferentialCorrelation(vector<TransitionTable> & tts,
                                     TransitionTable & ttTot,
                                     TransitionTable & ttHom,
                                     TransitionTable & ttHet)
{
    vector<double> corCoeffs(tts.size());
    double avgCor=0.0, diffCor=0.0, avgAbsCor=0.0;
    for (size_t k=0; k<tts.size(); k++) {
        corCoeffs[k] = tts[k].cor();
        avgCor += corCoeffs[k];
        avgAbsCor += abs(corCoeffs[k]);
        tts[k].setChisq(corCoeffs[k]);
        tts[k].setpValue(1-abs(corCoeffs[k]));
    }
    avgCor /= tts.size();
    avgAbsCor /= tts.size();
    
    for (size_t k=0; k<tts.size(); k++) {
        diffCor += abs(corCoeffs[k] - avgCor);
    }
    diffCor /= (2*tts.size()); // normalize to between 0 and 1
    
    // Save correlation related statistics into the chisq spots
    //   in the data structure
    ttTot.setChisq(avgAbsCor);
    ttTot.setpValue(1-avgAbsCor);
    ttHom.setChisq(abs(avgCor));
    ttHom.setpValue(1-abs(avgCor));
    ttHet.setChisq(diffCor);
    ttHet.setpValue(1-diffCor);
}

void EnvGLNCmp::collectDiffCorPooled(int child,
                                   const vector<vector< int > > & parents,
                                   const vector<vector< int > > & delays,
                                   vector<TransitionTable> & tts,
                                   TransitionTable & ttTot,
                                   TransitionTable & ttHom,
                                   TransitionTable & ttHet) const
{
    // Fill in tts using the exact parent set for each condition
    tts = fillIndividualTables(child, parents, delays, m_trajCols);
    
    // Compute invdividual chisqs using expected value from *individual* marginals
    computeIndividualChisq(tts, m_gc.m_pValueMode);
    
    // Find parent super set
    vector< int > parentSuperset, delaySuperset;
    vector<TransitionTable> ttsSuperset;
    
    mergeParentsWithSameDelays(parents, delays,
                               parentSuperset, delaySuperset);
    
    vector< vector< int > > parentSupersets(parents.size(), parentSuperset);
    vector< vector< int > > delaySupersets(delays.size(), delaySuperset);
    
    ttsSuperset = fillIndividualTables(child, parentSupersets, delaySupersets,
                                       m_trajCols);
    
    computerDifferentialCorrelation(ttsSuperset, ttTot, ttHom, ttHet);
}


void EnvGLNCmp::collectStatsPooled(int child,
                             const vector<vector< int > > & parents,
                             const vector<vector< int > > & delays,
                             vector<TransitionTable> & tts,
                             TransitionTable & ttTot,
                             TransitionTable & ttHom,
                             TransitionTable & ttHet) const
{
    // Fill in tts using the exact parent set for each condition
    tts = fillIndividualTables(child, parents, delays, m_trajCols);
    
    // Compute invdividual chisqs using expected value from *individual* marginals
    computeIndividualChisq(tts, m_gc.m_pValueMode);
    
    // Find parent super set
    vector< int > parentSuperset, delaySuperset;
    vector<TransitionTable> ttsSuperset;
    
    mergeParentsWithSameDelays(parents, delays,
                               parentSuperset, delaySuperset);
    
    vector< vector< int > > parentSupersets(parents.size(), parentSuperset);
    vector< vector< int > > delaySupersets(delays.size(), delaySuperset);
    
    // Fill in tts using the exact parent set for each condition
    // ttHom = computeHomChisq(child, parentSuperset, delaySuperset,
    //                         m_trajCols, m_gc.m_pValueMode);
    ttHom = fillPooledTable(child, parentSuperset, delaySuperset, m_trajCols);
    computeHomChisq(ttHom, m_gc.m_pValueMode);
    
    // ttsSuperset = computeIndividualChisq(child, parentSupersets, delaySupersets,
    //                                     m_trajCols, m_gc.m_pValueMode, ttHom);

    ttsSuperset = fillIndividualTables(child, parentSupersets, delaySupersets,
                                       m_trajCols);

    if(m_gc.m_methodMarginal != INDEP_SUP) {
        // Compute invdividual chisqs using expected value from *pooled* marginals
        computeIndividualChisq(ttsSuperset, m_gc.m_pValueMode, ttHom);
    } else {
        computeIndividualChisq(ttsSuperset, m_gc.m_pValueMode);
    }
    
    // Total strength of all interactions
    ttTot = computeTotalChisq(ttsSuperset, m_gc.m_pValueMode);
    
    // statistic for heterogeneity
    ttHet = computeHetChisq(ttsSuperset, ttTot, ttHom, m_gc.m_pValueMode,
                            m_gc.m_methodMarginal);
}

bool EnvGLNCmp::select_BY_BEST_HOM(int child, 
                                   const vector<TransitionTable> & tts,
                                   const TransitionTable & ttTot, 
                                   const TransitionTable & ttHom)
{
    bool selected = false;
    
    bool improvedPooled = compare_chisq(ttHom.getpValue(), ttHom.getDf(), 
                                        ttHom.getChisq(), 
                                        m_pooledTransTables[child].getpValue(),
                                        m_pooledTransTables[child].getDf(),
                                        m_pooledTransTables[child].getChisq());

    const TransitionTable * p_ttHom_best;
    const vector<TransitionTable> * p_tts_best = & tts;
    
    if(improvedPooled) {
        
        p_ttHom_best = & ttHom;
        m_pooledTransTables[child] = ttHom;
        
    } else {
        
        p_ttHom_best = & m_pooledTransTables[child];
        
    }
        
    double p_t = ttTot.getpValue(); // total strength of interactions across all conditions
    int df_t = ttTot.getDf();
    double chisq_t = ttTot.getChisq();
   
    bool improvedTotal = compare_chisq(p_t, df_t, chisq_t, m_pchisqs_t[child], 
                                       m_dfs_t[child], m_chisqs_t[child]);
        
    if(! improvedTotal) {
        chisq_t = m_chisqs_t[child];
        p_t = m_pchisqs_t[child];
        df_t = m_dfs_t[child];
        p_tts_best = & m_TransTables[child];
    }
    
    selected = improvedTotal || improvedPooled;
        
    if(selected) { // Heterogeneity test
        
        double F_a, p_a;
        
        if(p_ttHom_best -> getDf() > 0) {

            F_a = ( chisq_t / df_t ) / ( p_ttHom_best -> getChisq() / p_ttHom_best -> getDf() );
            
            p_a = FPvalue(F_a, df_t, p_ttHom_best -> getDf());
            
        } else {
            
            F_a = chisq_t;
            p_a = df_t;
            
        }
        
        TransitionTable ttHet; // statistic for heterogeneity

        ttHet.setpValue(p_a);
        ttHet.setChisq(F_a);
        ttHet.setDf(df_t);

        // Alternative way of testing heterogeneity:
        // int df_d = df_t - ttHom.getDf(); // degrees of freedom for chisq_d 
        // double p_d;	// p-value for chisq_d
        // double chisq_d = HeteroChisqTest(chisq_t, ttHom.getChisq(), df_d, p_d);
        // ttHet.setpValue(p_d);
        // ttHet.setChisq(chisq_d);
        // ttHet.setDf(df_d);
        
        double chisq_parent_z;
        size_t df_parent_z;
        double p_parent_z;
        
        computeParentWorkingZoneChange(* p_tts_best, chisq_parent_z, df_parent_z, p_parent_z);
        
        CmpIntxType type;
        
        // type = labelIntxType(ttTot.getpValue(), ttHom.getpValue(), ttHet.getpValue(), m_gc.m_alpha);
        
		type = labelIntxType(ttTot.getpValue(), p_ttHom_best -> getpValue(), ttHet.getpValue(), m_gc.m_alpha);

        updateEnvironment(child, type, * p_tts_best, ttTot, * p_ttHom_best, ttHet,  
                          p_parent_z, (int) df_parent_z, chisq_parent_z);
    }
    
    return selected;
}


bool EnvGLNCmp::compareParentSets(int child, int type, 
                                  const vector <TransitionTable> & tts,
                                  const TransitionTable & ttTot,
                                  const TransitionTable & ttHom,
                                  const TransitionTable & ttHet) const
{    
    //	return( ( p_t < m_pchisqs_t[child] )
    //			||	( ( p_t == m_pchisqs_t[child] ) && ( chi_t > m_chisqs_t[child] ) ) );
    
    bool improved = false;

    double chisq_t = ttTot.getChisq();
    int df_t = ttTot.getDf();
    double p_t = ttTot.getpValue();
    
	switch (m_gc.m_methodCmpxParCmp) {

        case BY_INTX_TYPE:
            improved = cmpParents_BY_INTX_TYPE(child, type, tts, ttTot, ttHom, ttHet);
            break;
            
		case BY_TOT:  // by chisq_t
            
			/* Should we use parent size or degrees of freedom?
             improved = compare(p_t, (int) pooled.getParents().size(), chi_t,
                                m_pchisqs_t[child], 
                                (int) m_pooledTransTables[child].getParents().size(),
                                m_chisqs_t[child]);
             */
			
			improved = compare_chisq(p_t, df_t, chisq_t, 
									 m_pchisqs_t[child], m_dfs_t[child], m_chisqs_t[child]);
			break;
            
		case BY_HET: // by chisq_d

            improved = compare_chisq(ttHet.getpValue(),
                                     ttHet.getDf(),
                                     ttHet.getChisq(),
                                     m_diffTransTables[child].getpValue(),
                                     m_diffTransTables[child].getDf(),
                                     m_diffTransTables[child].getChisq());

            break;

		case BY_HOM: // Somehow emphasizing conserved interactions

            improved = cmpParents_BY_HOM(child, tts, ttTot, ttHom, ttHet);
            
            break;

        case BY_EACH_COND: // by chisq_1, ..., chisq_k
            
            improved = cmpParents_BY_EACH_COND(child, tts, ttTot, ttHom, ttHet);

            break;
            
        case BY_T_IMP_PER_DF:
            
            improved = cmpParents_BY_T_IMP_PER_DF(child, tts, ttTot, ttHom, ttHet);
            break;
            
		case BY_T_C_IMP:
            
            improved = cmpParents_BY_T_C_IMP(child, tts, ttTot, ttHom, ttHet);
            
			break;
            
		default:
			break;
	}
	
    if(child == -1) { // for debugging purpose
        cout << child << ":";
        for(size_t k=0; k<tts.size(); k++) {
            cout << "\tCondition " << k << " chisq=" << tts[k].getChisq() << "--";
            cout << "\t";
            for(size_t i=0; i<tts[k].getParents().size(); i++) {
                cout << tts[k].getParents()[i] << ",";
            }
            cout << "\t\t";
        }
        cout << endl;
        for(size_t k=0; k<m_TransTables[child].size(); k++) {
            cout << "\tCondition " << k << " chisq=" << m_TransTables[child][k].getChisq() << "--";
            cout << "\t";
            for(size_t i=0; i<m_TransTables[child][k].getParents().size(); i++) {
                cout << m_TransTables[child][k].getParents()[i] << ",";
            }
            cout << "\t\t";
        }
        cout << endl << endl;
    }
    
	return improved;
}

bool EnvGLNCmp::cmpParents_BY_T_IMP_PER_DF(int child, 
                                           const vector <TransitionTable> & tts,
                                           const TransitionTable & ttTot,
                                           const TransitionTable & ttHom,
                                           const TransitionTable & ttHet) const
{
    bool improved;
    
    // encouraging conserved interactions
    
    //int df_c = pooled.getDf();
    
    int df_t_best = m_dfs_t[child];
    double chisq_t_best = m_chisqs_t[child];
    
    //int df_c_best = m_pooledTransTables[child].getDf();
    
    // improved = compare_F_imp(df_t - df_c, chisq_t, df_t_best - df_c_best, chisq_t_best);            
    improved = compare_F_imp(ttTot.getDf(), ttTot.getChisq(), df_t_best, chisq_t_best);            
    
    /*
     cout << "\tThe child " << chisq_t << "/" << df_t << " has ";
     if(improved) {
     cout << "improved over";
     } else {
     cout << "not improved over";
     }
     cout << " the current best " << chisq_t_best << "/" << df_t_best << endl;
     */
    
    return improved;
}
    
bool EnvGLNCmp::cmpParents_BY_T_C_IMP(int child, const vector <TransitionTable> & tts,
                                      const TransitionTable & ttTot,
                                      const TransitionTable & ttHom,
                                      const TransitionTable & ttHet) const
{
    bool improved;
    size_t k;
    
    double p_cur = ttTot.getpValue();
    int df_cur = ttTot.getDf();
    double chisq_cur = ttTot.getChisq();
    
    for (k=1; k<tts.size() ; k++) {
        // The comparison assumes parents are sorted in each tt table. Need to check.
        if (tts[k].getParents() != tts[k-1].getParents() ) {
            break;
        }
    }
    
    if(k == tts.size()) { // when parents are the same across conditions, use smaller of p_t and p_c
        double p_c_cur = ttHom.getpValue();
        int df_c_cur = ttHom.getDf();
        double chisq_c_cur = ttHom.getChisq();
        if( compare_chisq_imp(p_c_cur, df_c_cur, chisq_c_cur, p_cur, df_cur, chisq_cur) ) {
            p_cur = p_c_cur, df_cur = df_c_cur, chisq_cur = chisq_c_cur;
        }
    }
    
    for (k=1; k < m_TransTables[child].size(); k++) {
        // Parents are assumed to be sorted in each tt table.  Seems to be so.
        if (m_TransTables[child][k].getParents() != m_TransTables[child][k-1].getParents() ) {
            break;
        }
    }
    
    double p_best = m_pchisqs_t[child];
    int df_best = m_dfs_t[child];
    double chisq_best = m_chisqs_t[child];
    
    if( k == m_TransTables[child].size() ) { 
        // when parents are the same across conditions, use smaller of p_t and p_c
        double p_c_best = m_pooledTransTables[child].getpValue();
        int df_c_best = m_pooledTransTables[child].getDf();
        double chisq_c_best = m_pooledTransTables[child].getChisq();		
        
        if( compare_chisq_imp(p_c_best, df_c_best, chisq_c_best, p_best, df_best, chisq_best) ) {
            p_best = p_c_best, df_best = df_c_best, chisq_best = chisq_c_best;
        }
    }
    
    /*
     improved = compare_chisq_improvement(p_cur, df_cur, chisq_cur, p_best, df_best, chisq_best);
     */
    improved = compare_F_imp_over_chisq(p_cur, df_cur, chisq_cur, p_best, df_best, chisq_best);
    
    return improved;    
}

enum CMP_STAT {TOT, HET, HOM};

CMP_STAT ChooseChisq_BY_HOM(const vector <TransitionTable> & tts,
                            const TransitionTable & ttTot,
                            const TransitionTable & ttHom,
                            const TransitionTable & ttHet)
{
    CMP_STAT choice;
        
    size_t k;
    
    for (k=1; k<tts.size() ; k++) {
        // The comparison assumes parents are sorted in each tt table. Need to check.
        if (tts[k].getParents() != tts[k-1].getParents() ) {
            break;
        }
    }
    
    if(k == tts.size()) { // same parents across conditions, use better of p_t and p_c
                
        double p_c_cur = ttHom.getpValue();
        int df_c_cur = ttHom.getDf();
        double chisq_c_cur = ttHom.getChisq();
        
        if( compare_chisq(p_c_cur, df_c_cur, chisq_c_cur, 
                          ttTot.getpValue(), ttTot.getDf(), ttTot.getChisq()) ) {
            
            choice = HOM;
            
        } else {
            
            choice = TOT;
        }
        
    } else { // different parents
      
         choice = TOT;
        
    }
    return choice;
}


bool EnvGLNCmp::cmpParents_BY_HOM(int child, const vector <TransitionTable> & tts,
                                  const TransitionTable & ttTot,
                                  const TransitionTable & ttHom,
                                  const TransitionTable & ttHet) const
{
    bool improved;
    
    double p_cur=1;
    int df_cur=0;
    double chisq_cur=0;

    switch( ChooseChisq_BY_HOM(tts, ttTot, ttHom, ttHet) ) {
        case HOM:
            p_cur = ttHom.getpValue();
            df_cur = ttHom.getDf();
            chisq_cur = ttHom.getChisq();
            break;

        case HET:
            p_cur = ttHet.getpValue();
            df_cur = ttHet.getDf();
            chisq_cur = ttHet.getChisq();
            break;
            
        case TOT:
            p_cur = ttTot.getpValue();
            df_cur = ttTot.getDf();
            chisq_cur = ttTot.getChisq();
            break;
    }

    double p_best=1;
    int df_best=0;
    double chisq_best=0;
    
    TransitionTable ttTotChild;
    ttTotChild.setpValue(m_pchisqs_t[child]); 
    ttTotChild.setDf(m_dfs_t[child]); 
    ttTotChild.setChisq(m_chisqs_t[child]);
    
    switch( ChooseChisq_BY_HOM(m_TransTables[child],
                               ttTotChild,
                               m_pooledTransTables[child], 
                               m_diffTransTables[child] 
                               ) ) 
    {
        case HOM:
            p_best = m_pooledTransTables[child].getpValue();
            df_best = m_pooledTransTables[child].getDf();
            chisq_best = m_pooledTransTables[child].getChisq();
            break;
            
        case HET:
            p_best = m_diffTransTables[child].getpValue();
            df_best = m_diffTransTables[child].getDf();
            chisq_best = m_diffTransTables[child].getChisq();
            break;
            
        case TOT:
            p_best = m_pchisqs_t[child];
            df_best = m_dfs_t[child];
            chisq_best = m_chisqs_t[child];
            break;
    }
    
    improved = compare_chisq(p_cur, df_cur, chisq_cur, p_best, df_best, chisq_best);
    
    return improved;    
}

/*
bool EnvGLNCmp::cmpParents_BY_HOM2(int child, const TransitionTable & ttHom,
                                  const TransitionTable & ttHet, 
                                  const vector <TransitionTable> & tts, 
                                  double p_t, int df_t, double chisq_t) const
{
    bool improved;
    
    double p_cur;
    int df_cur;
    double chisq_cur;

    if( compare_chisq(ttHet.getpValue(), ttHet.getDf(), ttHet.getChisq(),
                      ttHom.getpValue(), ttHom.getDf(), ttHom.getChisq()
                      ) 
    {
        p_cur = ttHet.getpValue();
        df_cur = ttHet.getDf();
        chisq_cur = ttHet.getChisq();
        
    } else {
        p_cur = ttHom.getpValue();
        df_cur = ttHom.getDf();
        chisq_cur = ttHom.getChisq();
    }
    
    double p_best;
    int df_best;
    double chisq_best;
    
    switch( ChooseChisq_BY_HOM(m_pooledTransTables[child], 
                               m_diffTransTables[child], 
                               m_TransTables[child], 
                               m_pchisqs_t[child], 
                               m_dfs_t[child], 
                               m_chisqs_t[child]) ) 
    {
        case HOM:
            p_best = m_pooledTransTables[child].getpValue();
            df_best = m_pooledTransTables[child].getDf();
            chisq_best = m_pooledTransTables[child].getChisq();
            break;
            
        case HET:
            p_best = m_diffTransTables[child].getpValue();
            df_best = m_diffTransTables[child].getDf();
            chisq_best = m_diffTransTables[child].getChisq();
            break;
            
        case TOT:
            p_best = m_pchisqs_t[child];
            df_best = m_dfs_t[child];
            chisq_best = m_chisqs_t[child];
            break;
    }
    
    improved = compare_chisq(p_cur, df_cur, chisq_cur, p_best, df_best, chisq_best);
    
    return improved;    
}
*/

bool EnvGLNCmp::cmpParents_BY_EACH_COND(int child, 
                                        const vector <TransitionTable> & tts,
                                        const TransitionTable & ttTot,
                                        const TransitionTable & ttHom,
                                        const TransitionTable & ttHett) const
{
    bool improved = false;

    size_t k;
    
    /*
    double total_impr = 0;
 
    for(k=0; k<tts.size(); k++) {

        double improvement = compute_F_imp(tts[k].getDf(), 
                                           tts[k].getChisq(), 
                                           m_TransTables[child][k].getDf(),
                                           m_TransTables[child][k].getChisq());            
        if(improvement < 0) break;
        total_impr += improvement;
    }
    
    if(k == tts.size()) {
        
        improved
        = (total_impr > 0 || ttTot.getDf() < m_dfs_t[child]) ? true : false; 
        
    } else {
        
        improved = false;
        
    }
    */

    for(k=0; k<tts.size(); k++) {
        
        if (tts[k].getpValue() == m_TransTables[child][k].getpValue()
            && tts[k].getDf() == m_TransTables[child][k].getDf()
            && tts[k].getChisq() == m_TransTables[child][k].getChisq()) {
            continue;
        }
        
        improved = 
        compare_chisq(tts[k].getpValue(), 
                      tts[k].getDf(),
                      tts[k].getChisq(), 
                      m_TransTables[child][k].getpValue(), 
                      m_TransTables[child][k].getDf(), 
                      m_TransTables[child][k].getChisq() 
                      );
        
        if(! improved ) {
            break;
        }
        
    }
    
    return improved;
}


bool EnvGLNCmp::cmpParents_BY_INTX_TYPE(int child, int type, 
                                        const vector <TransitionTable> & tts,
                                        const TransitionTable & ttTot,
                                        const TransitionTable & ttHom,
                                        const TransitionTable & ttHet) const
// *** This function needs to be kept for the Drosophila study ***
{
    bool improved = false;
    
    switch ( m_types[child] ) {
            
        case NULL_INTX:
            
            improved = (type != NULL_INTX);
            
            break;
            
        case REL_DIFF:
            
            switch (type) {
                    
                case NULL_INTX:
                    improved = false;
                    break;
                    
                case REL_DIFF:
                    improved = // compare(m_diffTransTable, m_PooledTransTable[child]);
                    compare_chisq(ttHet.getpValue(), ttHet.getDf(), ttHet.getChisq(), 
                                  m_pooledTransTables[child].getpValue(), 
                                  m_pooledTransTables[child].getDf(),
                                  m_pooledTransTables[child].getChisq() 
                                  );
                    break;
                    
                case ABS_DIFF:
                case CONSERVED:
                    improved = true;
                    break;
                default:
                    break;
            }
            
            break;
            
        case ABS_DIFF:
            
            switch (type) {
                case NULL_INTX:
                case REL_DIFF:
                    improved = false;
                    break;
                    
                case ABS_DIFF:
                    improved = compare_chisq(ttTot.getpValue(), 
                                             (int) ttHom.getParents().size(),
                                             ttTot.getChisq(),
                                             m_pchisqs_t[child], 
                                             (int) m_pooledTransTables[child].getParents().size(),
                                             m_chisqs_t[child] 
                                             );
                    break;
                    
                case CONSERVED:
                    improved = compare_chisq(ttHom.getpValue(), 
                                             (int) ttHom.getParents().size(),
                                             ttHom.getChisq(),
                                             m_pchisqs_t[child], 
                                             (int) m_pooledTransTables[child].getParents().size(),
                                             m_chisqs_t[child]
                                             );
                    break;
                    
                default:
                    break;
            }
            
            break;
            
        case CONSERVED:
            
            switch (type) {
                case NULL_INTX:
                case REL_DIFF:
                    improved = false;
                    break;
                    
                case ABS_DIFF:
                    improved = compare_chisq(ttTot.getpValue(), 
                                             (int) ttHom.getParents().size(),
                                             ttTot.getChisq(),
                                             m_pooledTransTables[child].getpValue(), 
                                             (int) m_pooledTransTables[child].getParents().size(),
                                             m_pooledTransTables[child].getChisq()
                                             );
                    
                    break;
                case CONSERVED:
                    improved = // compare(m_pooledTransTable, m_PooledTransTable[child]);
                    compare_chisq(ttHom.getpValue(), 
                                  ttHom.getDf(), 
                                  ttHom.getChisq(),
                                  m_pooledTransTables[child].getpValue(), 
                                  m_pooledTransTables[child].getDf(),
                                  m_pooledTransTables[child].getChisq() 
                                  );
                    
                    break;
                default:
                    break;
            }
            
            break;
            
        default:
            break;
    }
    
    return improved;
}

/*
{
    bool improved = false;
    //so far there is no interaction detected, so any new discovery will update the environment
    if(m_types[child] == NULL_INTX && type != NULL_INTX) //if(m_types[child]==0&&type!=0)
    {
        improved = true;
        return improved;
    }
    if(m_types[child] == REL_DIFF) // if(m_types[child]==1)
    {
        //if(type==1&&ttHet.getpValue()<=m_pooledTransTables[child].getpValue())
        if(type == REL_DIFF && compare(diff, m_pooledTransTables[child])) // if(type==1&&compare(diff,m_pooledTransTables[child]))
        {
            improved = true;
            return improved;
        }
        if(type == ABS_DIFF || type == CONSERVED) // if(type==2||type==3)
        {
            improved = true;
            return improved;
        }
    }
    if(m_types[child] == ABS_DIFF) // if(m_types[child]==2)
    {
        if(type == ABS_DIFF  // if(type==2 
           && compare_chisq(p_t, (int) pooled.getParents().size(), chi_t,
                      m_pchisqs_t[child], (int) m_pooledTransTables[child].getParents().size(),m_chisqs_t[child]))
        {
            improved = true;
            return improved;
        }
        if(type== CONSERVED  // if(type== 3
           && compare_chisq(pooled.getpValue(), (int) pooled.getParents().size(),pooled.getChisq(),
                      m_pchisqs_t[child], (int) m_pooledTransTables[child].getParents().size(), 
                      m_chisqs_t[child]))
        {
            improved = true;
            return improved;
        }
    }
    //so far best parents show a homogeneity interaction
    if(m_types[child] == CONSERVED) // if(m_types[child]==3)
    {
        if(type==CONSERVED && compare(pooled,m_pooledTransTables[child])) // if(type==3 && compare(pooled,m_pooledTransTables[child]))
        {
            improved = true;
            return improved;
        }
        if(type==ABS_DIFF // if(type==2
           && compare_chisq(p_t, (int) pooled.getParents().size(), chi_t,
                      m_pooledTransTables[child].getpValue(), 
                      (int) m_pooledTransTables[child].getParents().size(),
                      m_pooledTransTables[child].getChisq()))
        {
            improved = true;
            return improved;
        }
    }
    return improved;
}
*/
