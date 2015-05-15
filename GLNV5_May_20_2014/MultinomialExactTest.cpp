//
//  MultinomialExactTest.cpp
//  GLNV5_Apr_14_2014
//
//  Created by Hua Zhong on 4/17/14.
//  Copyright (c) 2014 ZH. All rights reserved.
//

#include "MultinomialExactTest.h"

double factorial (double n);

double factorial_vector(vector<int> v);

//double factorial_double (double p, int k){
//	if(k<0){
//		for(int i=1; i<=abs(k); i++) {
//			p /= i;
//		}
//	}else{
//		for(int i=1; i<=k; i++) {
//			p *= i;
//		}
//	}
//	return p;
//}

double marginalDistribution(const vector<TransitionTable> & Cs,  int index,  vector<vector<int>> & marginal, string mode){//index:0->row, 1->column
	double P = 1.0;
	if(mode == "uniform"){
		for (size_t k=0; k<Cs.size(); k++) {
            //condition: enumerating marginal
			P *= (double)1/(pow((int)marginal[k].size(), (int)Cs[k].getTotalSum()));

			vector<int> v = marginal[k];
			v.push_back(-Cs[k].getTotalSum());
			P /= factorial_vector(v);
			//
		}
	}else if(mode == "observed"){
		
		for (size_t k=0; k<Cs.size(); k++) {
            //condition: enumerating marginal
			//P *= (double)1/(pow((int)marginal[k].size(), (int)e.As[k].getTotalSum()));

			vector<double> observedMarginal;
			int total = Cs[k].getTotalSum();

			if(index == 0){
				vector<int> tmp = Cs[k].getRowSums();
				for (int j = 0; j < tmp.size(); j++)
				{
					observedMarginal.push_back(tmp[j]);
				}
				
			}else if(index == 1){
				vector<int> tmp = Cs[k].getColSums();
				for (int j = 0; j < tmp.size(); j++)
				{
					observedMarginal.push_back(tmp[j]);
				}
			}
			
			for (int i = 0; i < observedMarginal.size(); i++)
			{
				P *= pow((double)observedMarginal[i] / (double) total, (double) marginal[k][i]);
			}

			vector<int> v = marginal[k];
			v.push_back(-Cs[k].getTotalSum());
			P /= factorial_vector(v);
			//
		}
	}
	return P;
}

void oneWayEnumerator::initialize(const MultinomialExactTest & v, const vector<TransitionTable> & Cs){
    for(size_t k=0; k<Cs.size(); k++) {
        vector<int> parents(1);
        parents[0]=1;
        
        vector<int> parentBases(1);
        parentBases[0]=1;
    
        vector<string> nameParents(1);
        nameParents[0]="P1";
    
        vector<int> delays(1, 0);
    
        // child
        int child = 0, base = v.getNum()[k];
        string nameChild="Child";
        TransitionTable tmp;
        tmp.Transition::initialize(child, base, parents, parentBases, delays, nameChild, nameParents);
        As.push_back(tmp);
    }
    
    ARowsums.resize(Cs.size());
    AColsums.resize(Cs.size());
    
    for(size_t k=0; k<Cs.size(); k++) {
        ARowsums[k].resize(Cs[k].nrow());
        AColsums[k].resize(Cs[k].ncol());
    }
    
}

double
oneWayEnumerator::traverse(size_t k, size_t i, size_t j,
                        MultinomialExactTest & v,
                        const vector< TransitionTable > & Cs, const oneWayEnumerator & oe)
{
    double P;
    
    // Base case:
    if(k == As.size()) {
        // k being out of bound indciates (k,i,j) is beyond a leaf
        P = v.evaluate(*this, Cs, oe);
        return P;
    } else {
        P = 0;
    }
    
    // branch-and-bound: decide whether to skip the entire branch starting
    //   at and including Ak[i,j]
    BOUND_CHECK check = v.bound(k, i, j, *this, Cs);
    if( check == TO_SKIP_ENTIRE_BRANCH ) {
        return 0;
    } else if( check == TO_KEEP_ENTIRE_BRANCH ) {
        P = v.add(k, i, j, *this, Cs);
        return P;
    }
    
    size_t k_next=k, i_next=i, j_next=j;
    next(k_next, i_next, j_next);
    
    if (0) {
        cout << "k=" << k << ", i=" << i << ", j=" << j
        << "; next k=" << k_next << ", i=" << i_next << ", j=" << j_next
        << endl;
    }
    
    if(As[k].nrow() == 0 || As[k].ncol() == 0) {
        
        // As[k] is an empty matrix.  Move on to next one
        update(k_next, i_next, j_next, k, i, j, v, Cs);
        P += traverse(k_next, i_next, j_next, v, Cs, oe);
        
    } else {
        // lower and uppler limits for Ak[i,j]
        int lij, uij;
        limits(k, i, j, v, lij, uij);
        
        for( int a=lij; a <= uij; a++) { // will enter loop if and only if lij<=uij
            
            As[k].set(i, j, a); // As[k][i][j] = a;
            
            ARowsums[k][i] += a;
            AColsums[k][j] += a;
            
            //size_t k_next=k, i_next=i, j_next=j;
            //next(k_next, i_next, j_next);
            
            update(k_next, i_next, j_next, k, i, j, v, Cs);
            
            // cout << "k=" << k << ", i=" << i << ", j=" << j << ", a=" << a << endl;
            
            P += traverse(k_next, i_next, j_next, v, Cs, oe);
            
            ARowsums[k][i] -= a;
            AColsums[k][j] -= a;
            
            As[k].set(i, j, 0); // As[k][i][j] = 0; // clear up the entry to return
            
        }
    }
    
    return P;
}

void oneWayEnumerator::next(size_t & k, size_t & i, size_t & j)
{
    // k: contingency table index from 0 to K-1 (#matrices)
    // i: row index from 0 to I-1 (I: #rows)
    // j: column index from 0 to J-1 (J: #columns)
    
    if(As[k].nrow() == 0 || As[k].ncol() == 0) {
        k = k+1;
        i = 0;
        j = 0;
    } else if( j < As[k].ncol()-1 ) {  // # columns in Ak
        j = j+1;
    } else if ( i < As[k].nrow()-1 ) { // #rows in Ak
        j = 0;
        i = i+1;
    } else {
        k = k+1;
        i = 0;
        j = 0;
    }
}

void oneWayEnumerator::limits(size_t k, size_t i, size_t j, const MultinomialExactTest & v, int &lij, int &uij){//Only limit the row sum, because the table only have 1 row.
    
    // lower and upper bounds for Ak[i,j]
    int bri = v.m_requiredRowSums[k][i] - ARowsums[k][i];
    lij = 0;
    
    if(j == As[k].ncol() - 1) { // last column
        lij = max(lij, bri);
    }
    uij = bri;
}


void MultinomialExactTest::initialize(const vector<TransitionTable> & Cs){
    for (size_t i=0; i<Cs.size(); i++) {
        m_totalSampleSize.push_back(Cs[i].getTotalSum());
    }
}

void MultinomialExactTest::initializeRequirement(const oneWayEnumerator e, const vector<TransitionTable> & Cs){
    m_requiredRowSums.resize(Cs.size());
    m_requiredColSums.resize(Cs.size());
    
    if (e.m_index == 0) {
		m_num.clear();
        for (size_t k=0; k<Cs.size(); k++) {
            m_num.push_back(Cs[k].getTransitionTable().size());
            m_requiredRowSums[k].resize(1);
            m_requiredRowSums[k][0] = Cs[k].getTotalSum();
            m_requiredColSums[k] = vector<int> (m_num[k], Cs[k].getTotalSum());
        }
    }else if (e.m_index == 1){
		m_num.clear();
        for (size_t k=0; k<Cs.size(); k++) {
            m_num.push_back(Cs[k].getTransitionTable()[0].size());
            m_requiredRowSums[k].resize(1);
            m_requiredRowSums[k][0] = Cs[k].getTotalSum();
            m_requiredColSums[k] = vector<int> (m_num[k], Cs[k].getTotalSum());
        }
    }
}

double MultinomialExactTest::evaluate(const oneWayEnumerator & e, const vector<TransitionTable> & Cs, const oneWayEnumerator & oe){
    
    if (e.m_index==0) {
        oneWayEnumerator e2;
        e2.setIndex(1);
		MultinomialExactTest v2;
		v2.initialize(Cs);
        v2.initializeRequirement(e2, Cs);
        e2.initialize(v2, Cs);
        return e2.traverse(0, 0, 0, v2, Cs, e);
    }else if (e.m_index==1){
        EMTEnumerator EMTe;
        double P=1;
        vector<vector<int>> row;
        vector<vector<int>> col;
        for (size_t k=0; k<Cs.size(); k++) {
            vector<int> parents(1);
            parents[0]=1;
            
            vector<int> parentBases(1);
            parentBases[0]=oe.As[k].getTransitionTable()[0].size();
            
            vector<string> nameParents(1);
            nameParents[0]="P1";//, nameParents[1]="P2";
            
            vector<int> delays(1, 0);
            
            // child
            int child = 0, base = getNum()[k];
            string nameChild="Child";
            TransitionTable tmp;
            tmp.Transition::initialize(child, base, parents, parentBases, delays, nameChild, nameParents);
            EMTe.As.push_back(tmp);
            
            row.push_back(oe.As[k].getTransitionTable()[0]);
            col.push_back(e.As[k].getTransitionTable()[0]);
            
            //exact_functional_test(e.As[k], "chisq");
            
        }
        
        EMTFunctionalChisq EMTv;
        EMTv.setRowAndColSum(row, col, Cs);
        EMTv.initialize_customized_row_col_sum(Cs);
        EMTe.initialize(EMTv, Cs);
        P = EMTe.EMTEnumerator::traverse(0, 0, 0, EMTv, Cs);
//Hua commented, Apr 13 2015
//        if(EMTv.getExtremeness()==true){
//            P = 1-P;
//        }
        P *= marginalDistribution(Cs, 0, row, "observed");
		P *= marginalDistribution(Cs, 1, col, "observed");
        return P;
    }
    return 0;
}

//////
double multi_table_multinomial_exact_test(const vector<TransitionTable> & Cs,
                              const string & discrepancy_measure){
    oneWayEnumerator e;
    e.setIndex(0);
    if ( discrepancy_measure == "chisq" ) {
        MultinomialExactTest v;
        v.initialize(Cs);
        v.initializeRequirement(e, Cs);
        e.initialize(v, Cs);
        oneWayEnumerator e_no_use;
        double p = e.traverse(0, 0, 0, v, Cs, e_no_use);
        return p;
    } else {
        cerr << "ERROR: unknown discrepancy measure in multi_table_functional_test!"
        << endl;
        return 1;
    }
}

double multinomial_exact_test(const TransitionTable & C,
                             const string & discrepancy_measure)
{
    vector<TransitionTable> Cs(1, C);
    return multi_table_multinomial_exact_test(Cs, discrepancy_measure);
}

double multinomial_exact_test(const vector< vector<int> > & C,
                             const string & discrepancy_measure)
{
    double pval;
    
    if( discrepancy_measure == "chisq" ) {
        
        TransitionTable tt(0, C[0].size(), vector<int>(1,0),
                           vector<int>(1, C.size()), vector<int>(1,0));
        tt.setTransitionTable(C);
        pval = multinomial_exact_test(tt, discrepancy_measure);
    }
    
    return pval;
}

