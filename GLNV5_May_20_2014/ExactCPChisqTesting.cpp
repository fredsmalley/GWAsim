//
//  ExactCPChisqTesting.cpp
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//

#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

#include "ComparativeChisq.h"


void run_exact_cpchisq_test(const vector<TransitionTable> & tts, int pvalmode,
                            double pd_truth, double pc_truth, double pt_truth,
                            short version)
{
    /*
    void exact_comparative_chisq_1(const vector<TransitionTable> & tts,
                                   double &pd, double &pc, double &pt,
                                   int pvaluemode);
    */
    
    double pd=1.0, pc=1.0, pt=1.0;
    
    switch (version) {
        case 1:
            // exact_comparative_chisq_1(tts, pd, pc, pt, pvalmode);
            break;

        case 2:
            exact_comparative_chisq_2(tts, pd, pc, pt, pvalmode);
            break;

        case 3:
            exact_comparative_chisq_3(tts, pd, pc, pt, pvalmode);
            break;

        default:
            break;
    }
    
    std::ostringstream sd, sc, st;
    sd << pd_truth;
    sc << pc_truth;
    st << pt_truth;
    
    cout << "   pvalue mode = " << pvalmode << ":" << endl
    << "     pd = " << setprecision(3) << pd << " (truth "
    << (pd_truth==-1 ? "?" : sd.str()) << ")" << endl
    << "     pc = " << setprecision(3) << pc << " (truth "
    << (pc_truth==-1 ? "?" : sc.str()) << ")" << endl
    << "     pt = " << setprecision(3) << pt << " (truth "
    << (pt_truth==-1 ? "?" : st.str()) << ")"
    << endl << endl;
}


void test_exact_comparative(short version)
{
	cout << ">>>> Testing exact comparative analysis ... version " << version << endl;
    
    int mat1[3][3]= {
        {3, 0, 4},
        {1, 3, 0},
        {1, 3, 0}
    };
    
    int mat2[3][3]= {
        {0, 5, 0},
        {5, 1, 0},
        {0, 0, 4}
    };

    /*
    int mat3[3][3]= {
        {2, 1, 0},
        {3, 4, 1},
        {0, 1, 3}
    };
    */
    
    int mat4[9][3]= {
        {0, 0, 0},
        {3, 0, 0},
        {0, 0, 4},
        {0, 3, 0},
        {1, 0, 0},
        {0, 0, 0},
        {0, 2, 0},
        {1, 1, 0},
        {0, 0, 0}
    };
    
	vector<TransitionTable> tts(2);
	vector<vector<int> > parents(2, vector<int>(1, 2));
	vector<vector<int> > delays(2, vector<int>(1, 0));
    
    for( size_t k=0; k<2; k++) {
        tts[k]=TransitionTable(0, 3, parents[k], vector<int>(parents[k].size(),3), delays[k]);
        for(size_t i=0; i<tts[k].nrow(); i++)
            for(size_t j=0; j<tts[k].ncol(); j++)
                tts[k].set(i, j, mat1[i][j]);
    }
    
    if(1) {
        cout << "   ---- parents and transition tables are the same ....   " << endl << endl;
        
        run_exact_cpchisq_test(tts, 8, 1, 3.37e-5, -1, version);
        run_exact_cpchisq_test(tts, 9, 1, 2.25e-5, -1, version);
    }
    
    if(1) {
        cout << "   ---- parents are same but transition tables are different ....   " << endl << endl;
        
        for( size_t k=1; k<2; k++) {
            tts[k]=TransitionTable(0, 3, parents[k], vector<int>(parents[k].size(),3), delays[k]);
            for(size_t i=0; i<tts[k].nrow(); i++)
                for(size_t j=0; j<tts[k].ncol(); j++)
                    tts[k].set(i, j, mat2[i][j]);
        }
        
        run_exact_cpchisq_test(tts, 8, -1, 8.77e-2, -1, version);
        run_exact_cpchisq_test(tts, 9, -1, 9.84e-2, -1, version);
        
    }
    
	parents[1][0] = 3;
    
    if(1) {
        
        cout << "   ------ parents are different without intersection ....   " << endl << endl;
        
        tts[1]=TransitionTable(0, 3, parents[1], vector<int>(parents[1].size(),3), delays[1]);
        for(size_t i=0; i<tts[1].nrow(); i++)
            for(size_t j=0; j<tts[1].ncol(); j++)
                tts[1].set(i, j, mat1[i][j]);

        run_exact_cpchisq_test(tts, 8, -1, 1, -1, version);
        run_exact_cpchisq_test(tts, 9, -1, 1, -1, version);
    }
    
	parents[1].push_back(2);
	delays[1].push_back(0);
    
    if(1) {
        
        cout << "   ----- parents are different and have common parents ....   " << endl << endl;
                
        for( size_t k=1; k<2; k++) {
            tts[k]=TransitionTable(0, 3, parents[k], vector<int>(parents[k].size(),3), delays[k]);
            for(size_t i=0; i<tts[k].nrow(); i++)
                for(size_t j=0; j<tts[k].ncol(); j++)
                    tts[k].set(i, j, mat4[i][j]);
        }
        
        run_exact_cpchisq_test(tts, 8, -1, 8.77e-2, -1, version);
        run_exact_cpchisq_test(tts, 9, -1, -1 /*9.84e-2*/, -1, version);
    }
}


void test_fisher_comparative()
{
	cout << ">>>> Testing exact comparative analysis (OLD) ..." << endl;
    
    string tcf1 = "test_gln1.trj";
	string tcf2 = "test_gln2.trj";
    
	ofstream out1(tcf1.c_str(),ios::out);
	out1<<"TRAJECTORY_VER2"<<endl;
	out1<<"1\t3\t0\n";
	out1<<"3\t3\t3\n";
	out1<<"Node1\tNode2\tNode3\n";
	out1<<"15\n";
	out1<<"0\t0\t0\n";
	out1<<"0\t0\t1\n";
	out1<<"0\t0\t0\n";
	out1<<"0\t1\t1\n";
	out1<<"1\t1\t1\n";
	out1<<"1\t1\t2\n";
	out1<<"1\t1\t1\n";
	out1<<"2\t0\t2\n";
	out1<<"2\t0\t1\n";
	out1<<"2\t0\t2\n";
	out1<<"2\t0\t2\n";
	out1<<"0\t2\t1\n";
	out1<<"1\t2\t1\n";
	out1<<"1\t2\t0\n";
	out1<<"1\t2\t1\n";
	out1.close();
    
	ofstream out2(tcf2.c_str(),ios::out);
	out2<<"TRAJECTORY_VER2"<<endl;
	out2<<"1\t3\t0\n";
	out2<<"3\t3\t3\n";
	out2<<"Node1\tNode2\tNode3\n";
	out2<<"15\n";
	out2<<"0\t1\t0\n";
	out2<<"0\t1\t0\n";
	out2<<"0\t1\t0\n";
	out2<<"0\t1\t1\n";
	out2<<"1\t0\t1\n";
	out2<<"1\t0\t1\n";
	out2<<"1\t0\t1\n";
	out2<<"2\t2\t0\n";
	out2<<"2\t2\t0\n";
	out2<<"2\t2\t0\n";
	out2<<"2\t2\t0\n";
	out2<<"0\t1\t2\n";
	out2<<"1\t0\t2\n";
	out2<<"1\t0\t2\n";
	out2<<"1\t1\t2\n";
	out2.close();
    
    TrajectoryCollection trajCol1;
    trajCol1.scan(tcf1);
    
    TrajectoryCollection trajCol2;
    trajCol2.scan(tcf2);
    
	vector<TrajectoryCollection> m_trajCols(2);
	m_trajCols[0] = trajCol1;
	m_trajCols[1] = trajCol1;
    
	vector<TransitionTable> tts(2);
	vector<vector<int> > parents(2, vector<int>(1, 2));
	vector<vector<int> > delays(2, vector<int>(1, 0));
    
    tts[0] = TransitionTable(0, parents[0], delays[0], m_trajCols[0]);
	tts[0].fill(m_trajCols[0]);
    
    tts[1] = TransitionTable(0, parents[1], delays[1], m_trajCols[1]);
	tts[1].fill(m_trajCols[1]);
    
    if(1) {
        cout << "   ---- parents and transition tables are the same ....   " << endl << endl;
        
        run_exact_cpchisq_test(tts, 8, 1, 3.37e-5, 4.42e-2, 1);
        run_exact_cpchisq_test(tts, 9, 1, 2.25e-5, 5.40e-2, 1);
    }
    
    if(1) {
        cout << "   ---- parents are same but transition tables are different ....   " << endl << endl;
        
        m_trajCols[1] = trajCol2;
        tts[1] = TransitionTable(0, parents[1], delays[1], m_trajCols[1]);
        tts[1].fill(m_trajCols[1]);
        
        run_exact_cpchisq_test(tts, 8, 2.22e-1, 8.77e-2, 2.24e-2, 1);
        run_exact_cpchisq_test(tts, 9, 2.22e-1, 9.84e-2, 2.74e-2, 1);
        
    }
    
	parents[1][0] = 3;
    
    if(1) {
        
        cout << "   ------ parents are different without intersection ....   " << endl << endl;
        
        tts[1] = TransitionTable(0, parents[1], delays[1], m_trajCols[1]);
        tts[1].fill(m_trajCols[1]);
        
        run_exact_cpchisq_test(tts, 8, 5.40e-2, 1, 4.42e-2, 1);
        run_exact_cpchisq_test(tts, 9, 5.40e-2, 1, 5.40e-2, 1);
    }
    
	parents[1].push_back(2);
	delays[1].push_back(0);
    
    if(1) {
        
        cout << "   ----- parents are different and have common parents ....   " << endl << endl;
        
        tts[1] = TransitionTable(0, parents[1], delays[1], m_trajCols[1]);
        tts[1].fill(m_trajCols[1]);
        
        run_exact_cpchisq_test(tts, 8, 3.04e-1, 8.77e-2, 2.24e-2, 1);
        run_exact_cpchisq_test(tts, 9, 3.04e-1, 9.84e-2, 2.74e-2, 1);
        
    }
}
