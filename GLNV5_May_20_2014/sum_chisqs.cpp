//
//  sum_chisqs.cpp
//  gln
//
//  Created by Joe Song on 3/10/13.
//
//
#include <vector>
using namespace std;

#include "StatDistributions.h"
#include "TrajectoryCollection.h"
#include "ChisqTests.h"

size_t cal_eff_dim(const vector< vector<int> > & X);

vector< vector<int> > stack(const vector<TrajectoryCollection> & trjCols,
                            const vector<int> & nodes)
{
    size_t N = 0; // total number of states on all trajectories
    for (size_t k=0; k<trjCols.size(); k++) {
        for (size_t i=0; i<trjCols[k].getNumTrajectories(); i++) {
            N += trjCols[k].getTrajectoryLength(i);
        }
    }
    
    size_t D = 0;
    
    if(nodes.empty()) {
        if ( N > 0) {
            D = trjCols[0].getTrajectory(0)[0].size();
        }
    } else {
        D = nodes.size();
    }
    
    vector< vector<int> > stacked(N, vector<int>(D));
    vector< vector<int> >::iterator itr = stacked.begin();
    
    for (size_t k=0; k<trjCols.size(); k++) {
        for (size_t i=0; i<trjCols[k].getNumTrajectories(); i++) {
            if(! nodes.empty() ) {
                for(size_t t=0; t<trjCols[k].getTrajectory(i).size(); t++) {
                    for(size_t m=0; m<nodes.size(); m++) {
                        (*itr)[m] = trjCols[k].getTrajectory(i)[t][nodes[m]-1];
                    }
                    itr ++;
                }
            } else {
                copy(trjCols[k].getTrajectory(i).cbegin(),
                     trjCols[k].getTrajectory(i).cend(),
                     itr);
                itr += trjCols[k].getTrajectory(i).size();
            }
        }
    }
    
    return stacked;
}

double scalar(const vector<TrajectoryCollection> & trjCols,
              const vector<int> & nodesInvolved)
{
    vector< vector<int> > stacked = stack(trjCols, nodesInvolved);
    
    size_t D = stacked[0].size();
    size_t D_eff = cal_eff_dim(stacked);
    
    double a = D / (double) D_eff;

    return a;
}

void scale_chisqs(Chisq & stat, double a)
{
    stat.m_pval = ChisqPvalue(stat.m_x / a, (int) (stat.m_df / a) );
}

/*
 vector<double>
 scale_chisqs(const vector< double > & chisqs,
 const vector< double > & dfs,
 const vector<TrajectoryCollection> & trjCols,
 const vector<int> & nodesInvolved)
 {
 double a = scalar(trjCols, nodesInvolved);
 vector<double> pvals(chisqs.size());
 
 for(size_t j=0; j<chisqs.size(); j++) {
 pvals[j] = ChisqPvalue(chisqs[j] / a, dfs[j] / a);
 }
 
 return pvals;
 }
 
 void scale_chisqs(vector<Chisq> & stats,
 const vector<TrajectoryCollection> & trjCols,
 const vector<int> & nodesInvolved)
 {
 double a = scalar(trjCols, nodesInvolved);
 
 for(size_t j=0; j<stats.size(); j++) {
 stats[j].m_pval = ChisqPvalue(stats[j].m_x / a, stats[j].m_df / a);
 }
 }
 */
