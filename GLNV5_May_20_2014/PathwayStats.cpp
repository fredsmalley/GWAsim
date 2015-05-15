//
//  PathwayStats.cpp
//  gln
//
//  Created by Joe Song on 3/11/13.
//
//
#include <fstream>
#include <cmath>

using namespace std;

#include "PathwayStats.h"
#include "Topology.h"
#include "TransitionTable.h"
#include "StatDistributions.h"

double scalar(const vector<TrajectoryCollection> & trjCols,
              const vector<int> & nodesInvolved);
void scale_chisqs(Chisq & stat, double a);

double gamma_pval(double x, double sum, double ss, size_t n)
{
    double mean = sum / n;
    double var = (ss - n * mean * mean) / ( n - 1 );
    
    if(sum==0||var==0)
    {
        return 1.0;
    }
    
    double scale = var / mean;
    double shape = mean * mean / var;
    
    double pval = GammaPvalue(x, shape, 1.0/scale);
    
    cout << "Gamma: x=" << x << ", mean=" << mean << ", " << "var=" << var
    << ", shape=" << shape << ", scale=" << scale << ", pval=" << pval << endl;
    
    return pval;
}

bool operator != (const MultiplePathwayStats & s1, const MultiplePathwayStats & s2)
{
    return !(s1 == s2);
}

bool operator == (const MultiplePathwayStats & s1, const MultiplePathwayStats & s2)
{
    if (s1.m_het != s2.m_het) {
        return false;
    } else if (s1.m_hom != s2.m_hom) {
        return false;
    } else if (s1.m_tot != s2.m_tot) {
        return false;
    } else if (s1.m_wzChildren != s2.m_wzChildren) {
        return false;
    } else if (s1.m_wzParents != s2.m_wzParents) {
        return false;
    } else {
        return true;
    }
}

const MultiplePathwayStats &
MultiplePathwayStats::approximate(const vector<MultiplePathwayStats> & statsBS)
{
    // MultiplePathwayStats stats = *this; // statsObs;
    
    cout << "Het-";
    double sum = 0.0, ss = 0.0;
    for (size_t i=0; i<statsBS.size(); i++) {
        sum += statsBS[i].m_het.m_x;
        ss += statsBS[i].m_het.m_x * statsBS[i].m_het.m_x;
    }
    m_het.m_pval = gamma_pval(m_het.m_x, sum, ss, statsBS.size());
    
    cout << "Hom-";
    sum = 0.0, ss = 0.0;
    for (size_t i=0; i<statsBS.size(); i++) {
        sum += statsBS[i].m_hom.m_x;
        ss += statsBS[i].m_hom.m_x * statsBS[i].m_hom.m_x;
    }
    m_hom.m_pval = gamma_pval(m_hom.m_x, sum, ss, statsBS.size());
    
    cout << "Tot-";
    sum = 0.0, ss = 0.0;
    for (size_t i=0; i<statsBS.size(); i++) {
        sum += statsBS[i].m_tot.m_x;
        ss += statsBS[i].m_tot.m_x * statsBS[i].m_tot.m_x;
    }
    m_tot.m_pval = gamma_pval(m_tot.m_x, sum, ss, statsBS.size());
    
    cout << "wzParents-";
    sum = 0.0, ss = 0.0;
    for (size_t i=0; i<statsBS.size(); i++) {
        sum += statsBS[i].m_wzParents.m_x;
        ss += statsBS[i].m_wzParents.m_x * statsBS[i].m_wzParents.m_x;
    }
    m_wzParents.m_pval = gamma_pval(m_wzParents.m_x, sum, ss, statsBS.size());
    
    cout << "wzChildren-";
    sum = 0.0, ss = 0.0;
    for (size_t i=0; i<statsBS.size(); i++) {
        sum += statsBS[i].m_wzChildren.m_x;
        ss += statsBS[i].m_wzChildren.m_x * statsBS[i].m_wzChildren.m_x;
    }
    m_wzChildren.m_pval = gamma_pval(m_wzChildren.m_x, sum, ss, statsBS.size());
    
    return *this;
}

void MultiplePathwayStats::compare
(const Topology & topo,
 const vector<TransitionTable> & tts,
 const vector<TransitionTable> & ttsPooled,
 const vector<double> & chisqs_t,
 const vector<int> & dfs_t,
 const vector<double> & chisqs_z,
 const vector<int> & dfs_z,
 const vector<double> & child_chisqs_z,
 const vector<int> & child_dfs_z)
{
    // Pathway interaction heterogeneity across all conditions
    m_het.m_df = 0;
    m_het.m_x = 0;
    for(size_t i=0; i < tts.size(); i++) {
        m_het.m_df += tts[i].getDf();
        m_het.m_x += abs(tts[i].getChisq());
    }
    m_het.m_pval = ChisqPvalue(m_het.m_x, m_het.m_df);

    // Pathway interaction homogeneity across all conditions
    m_hom.m_df = 0;
    m_hom.m_x = 0;
    for(size_t i=0; i < ttsPooled.size(); i++) {
        m_hom.m_df += ttsPooled[i].getDf();
        m_hom.m_x += ttsPooled[i].getChisq();
    }
    m_hom.m_pval = ChisqPvalue(m_hom.m_x, m_hom.m_df);
    
    // Pathway total interaction strength across all conditions
    m_tot.m_df = 0;
    m_tot.m_x = 0;
    for(size_t i=0; i < chisqs_t.size(); i++) {
        m_tot.m_df += dfs_t[i];
        m_tot.m_x += chisqs_t[i];
    }
    m_tot.m_pval = ChisqPvalue(m_tot.m_x, m_tot.m_df);
    
    // ??? Pathway parent working zone change across all conditions
    m_wzParents.m_df = 0;
    m_wzParents.m_x = 0;
    for(size_t i=0; i < chisqs_z.size(); i++) {
        m_wzParents.m_df += dfs_z[i];
        m_wzParents.m_x += chisqs_z[i];
    }
    m_wzParents.m_pval = ChisqPvalue(m_wzParents.m_x, m_wzParents.m_df);
    
    // Pathway child working zone change across all conditions
    m_wzChildren.m_df = 0;
    m_wzChildren.m_x = 0;
    vector<int> children = topo.getAllChildNodes();
    for(size_t i =0; i < children.size(); i ++ )
    {
        m_wzChildren.m_df += child_dfs_z[children[i]-1];
        m_wzChildren.m_x += child_chisqs_z[children[i]-1];
    }
    
    m_wzChildren.m_pval = ChisqPvalue(m_wzChildren.m_x, m_wzChildren.m_df);
    
}

void MultiplePathwayStats::scale(const vector<TrajectoryCollection> & trjCols,
           const vector<int> & nodeIDsOnPathway)
{
    double a = scalar(trjCols, nodeIDsOnPathway);
    cout << "Pathway chisq scale factor is " << a << endl;
    scale_chisqs(m_het, a);
    scale_chisqs(m_hom, a);
    scale_chisqs(m_tot, a);
    scale_chisqs(m_wzParents, a);
    scale_chisqs(m_wzChildren, a);
}

void MultiplePathwayStats::save(const string & file) const
{
    ofstream out(file, ios::out);
    out << "Chisq_d\t" << m_het.m_pval << "\t" << m_het.m_df << "\t" << m_het.m_x << endl;
    out << "Chisq_c\t" << m_hom.m_pval << "\t" << m_hom.m_df << "\t" << m_hom.m_x << endl;
    out << "Chisq_t\t" << m_tot.m_pval << "\t" << m_tot.m_df << "\t" << m_tot.m_x << endl;
    out << "Parents_Chisq_z\t" << m_wzParents.m_pval << "\t" << m_wzParents.m_df
        << "\t" << m_wzParents.m_x << endl;
    out << "Child_Chisq_z\t" << m_wzChildren.m_pval << "\t" << m_wzChildren.m_df
        << "\t" << m_wzChildren.m_x << endl;
    out.close();
}

void SinglePathwayStats::append(const string & file) const
{
    ofstream out(file, ios::app);
    //out.precision(2);
    out<<"GPEA pvalue:\t" << m_tot.m_pval <<endl;
    out<<"GPEA df:\t" << m_tot.m_df << endl;
    out<<"GPEA chisq:\t" << m_tot.m_x << endl;
    out.close();
}

void SinglePathwayStats::save(const string & file) const
{
    ofstream out(file);
    //out.precision(2);
    out<<"GPEA pvalue:\t" << m_tot.m_pval <<endl;
    out<<"GPEA df:\t" << m_tot.m_df << endl;
    out<<"GPEA chisq:\t" << m_tot.m_x << endl;
    out.close();
}
