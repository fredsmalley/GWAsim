//
//  EMTFisherTests.cpp
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//
#include <cmath>
#include <ctime>
#include <numeric>

using namespace std;

#include "ExactFunctionalTest.h"

double factorial_vector(vector<int>);//Hua Added, May 2 2014

bool sortAbs(double i, double j){
    //return (abs(i) < abs(j));
    return (i < j);
}
//----------------------------------------------------------------------------

void EMTFunctionalChisq::initialize_customized_row_col_sum(const vector<TransitionTable> & Cs){
		m_observedChisq.resize(Cs.size());
        m_nullChisq.resize(Cs.size());
        m_colSumChisq.resize(Cs.size());
        m_boundChisqs.resize(Cs.size());
        //        m_boundHeteroChisq = 0;
        //        m_row=0;
        
		vector<double> moreOrLessExtreme_tmp;//Hua added, May 17 2014
		moreOrLessExtreme_tmp.resize(Cs.size());//Hua added, May 17 2014
    
		m_totalObservedChisq = 0.0; // Line added by MS to prevent memory leak, Feb 23, 2015
        for (size_t k=0; k<Cs.size(); k++) {
            double chisq;
            size_t df;
            double p_tmp = ChisqDirTest(Cs[k].getTransitionTable(), chisq, df);
            m_observedChisq[k] = chisq;
            
			//Hua added, May 17 2014
            moreOrLessExtreme_tmp[k] = p_tmp;
			////
            
            m_totalObservedChisq += m_observedChisq[k];
        }
        
        //Hua added, Jun 5 2014, commented Apr 12, 2015
//        for (size_t k=0; k<Cs.size(); k++) {
//            double min=moreOrLessExtreme_tmp[0];
//            if(min < 1-moreOrLessExtreme_tmp[0])min = 1-moreOrLessExtreme_tmp[0];
//            int minIndex=0;
//            int index=0;
//            while(index < (int)moreOrLessExtreme_tmp.size()){
//                if(moreOrLessExtreme_tmp[index] < min || 1-moreOrLessExtreme_tmp[index] < min){
//                    minIndex = index;
//                }
//                index++;
//            }
//            if(moreOrLessExtreme_tmp[minIndex]<=0.5){
//                m_moreOrLessExtreme = false;
//            }else{
//                m_moreOrLessExtreme = true;
//            }
//        }
        ////
        
        ////Add by Hua Mar 20 2014
        for (size_t k=0; k<Cs.size(); k++) {
            size_t df=0;
            int nrow = (int) Cs[k].getTransitionTable().size();
            int K = (int) Cs[k].getTransitionTable()[0].size();
            
            vector<double> p_null(K, 1.0/K);
            vector<int> n_y(K);  // the histogram of child y
            
            for(int j = 0; j < K; j++){
                
                for(int i = 0; i < nrow; i++){
                    
                    n_y[j] += Cs[k].getTransitionTable()[i][j];
                    
                }//end for
                
            }//end for
            
            double chisq_y;
            ChisquareTest1DNoPValue(n_y, p_null, K, chisq_y, df);
            m_colSumChisq[k] = chisq_y;
        }
        /////
}

void EMTFunctionalChisq::initialize(const vector<TransitionTable> & Cs)
    {
        // compute required row sums
		m_requiredRowSums.resize(Cs.size());
		for(size_t k=0; k < Cs.size(); k++) {
			m_requiredRowSums[k] = Cs[k].getRowSums();
		}
    
		// compute required column sums
		m_requiredColSums.resize(Cs.size());
		for(size_t k=0; k < Cs.size(); k++) {
			m_requiredColSums[k] = Cs[k].getColSums();
		}

        initialize_customized_row_col_sum(Cs);
    }
    
void EMTFunctionalChisq::processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs)
    {
        double chisq;
        size_t df;
        ChisqDirTest(e.As[k].getTransitionTable(), chisq, df);
        m_nullChisq[k] =chisq;
    }
    
bool EMTFunctionalChisq::isMoreExtreme() const
    {
        double nullTotalChisq = 0, observedTotalChisq = 0;
        
        for(size_t k=0; k < m_nullChisq.size(); k++) {
            nullTotalChisq += m_nullChisq[k];
            observedTotalChisq += m_observedChisq[k];
        }
        
        //Commented by Hua, Apr 12 2015
//		if(m_moreOrLessExtreme ==false){//Hua added, May 17 2014
//			return nullTotalChisq >= observedTotalChisq;
//		}else{
//            return nullTotalChisq < observedTotalChisq;//No equal, modified by Hua, Jun 5 2014
//		}
        return nullTotalChisq >= observedTotalChisq;
    }

double EMTFunctionalChisq::evaluate(const EMTEnumerator & e, const vector<TransitionTable> & Cs)
    {
        // Multiple independent Fisher's exact tests
        double P;
        
		vector<vector<int> > row = e.ARowsums;
		vector<vector<int> > col = e.AColsums;
		
        if( isMoreExtreme() ) {
            P = 1;
            for(size_t k=0; k < Cs.size(); k++) {
				P *= FisherProb(e.As[k].getTransitionTable());
            }
        } else {
            P = 0;
        }
        
        return P;
    }
    
	vector<TransitionTable> EMTFunctionalChisq::generateTables(const vector<TransitionTable> & Cs) const
	{
		vector<TransitionTable> As(Cs);
    
		for (size_t k=0; k<Cs.size(); k++) {
			As[k].reset();
		}
    
		return As;
	}

//Added by Hua Apr 13 2015
double EMTFunctionalChisq::add(size_t k, size_t i, size_t j, const EMTEnumerator & e, const vector<TransitionTable> & Cs){
    
    vector<vector <int> > T;
    vector<vector <int> > e_table = e.As[k].getTransitionTable();
    if(i>0 && i<(int)e_table.size() && j==0){
        for (size_t t=0; t<i; t++) {
            T.push_back(e_table[t]);
        }
        
        vector<int> last_row(e_table[0].size());
        for (size_t c=0; c<e_table[0].size(); c++) {
            last_row[c] = m_requiredColSums[k][c] - e.AColsums[k][c];
        }
        T.push_back(last_row);
    }else{
        cerr<<"ExactFunctionalTest.cpp::add()::Not in correct corrdinates in table!"<<endl;
    }
    return FisherProb(T);
}

//Added by Hua Apr 13 2015
BOUND_CHECK EMTFunctionalChisq::bound(size_t k, size_t i, size_t j, const EMTEnumerator & e, const vector<TransitionTable> & Cs){
    BOUND_CHECK result = NOT_TO_SKIP; // "to keep entire branch"
    
    int rowNum = (int)e.As[k].getTransitionTable().size();
    int colNum = (int)e.As[k].getTransitionTable()[0].size();
    
    
    
    if(i>0 && i<rowNum && j==0){
        m_skip = false;
        double boundHeteroChisq_upper = 0.0;
        double boundHeteroChisq_lower = 0.0;
        
        for (size_t r=0; r<i; r++) {
            double chisq_tmp;
            size_t df=0;
            vector<double> p_null(colNum, 1.0/colNum);
            ChisquareTest1DNoPValue(e.As[k].getTransitionTable()[r], p_null, colNum, chisq_tmp, df);
            boundHeteroChisq_upper += chisq_tmp;
        }
        
        vector<int> colTotalSum_const = m_requiredColSums[k];//Cs[k].getColSums();
        
        //columnsum remaining
        for (size_t c=0; c<colTotalSum_const.size(); c++) {
            colTotalSum_const[c] = colTotalSum_const[c] - e.AColsums[k][c];
        }
        ////
        
        
        vector<int> rowTotalSum = m_requiredRowSums[k];//s[k].getRowSums();
        
        vector<int> tmpRow;
        
        //if(m_moreOrLessExtreme == false){//Hua added, May 17 2014
            //Upper bound
        if(result == NOT_TO_SKIP){
            bool skip_bound = false;
            vector<int> colTotalSum = colTotalSum_const;
            
            for (size_t r=i; r <rowNum; r++) {
                
                //|A11 - ave| < |A12 - ave|
                vector<double> colDistToAvg(colTotalSum.size());
                double everage = (double)rowTotalSum[r] / (double)colNum;//Could be double
                
                for (size_t c=0; c<colDistToAvg.size(); c++) {
                    colDistToAvg[c] = (double)colTotalSum[c] - everage;
                }
                
                sort(colDistToAvg.begin(), colDistToAvg.end(), sortAbs);
                for (size_t c=0; c<colDistToAvg.size(); c++) {
                    colDistToAvg[c] = colDistToAvg[c] + everage;
                }
                for (size_t c=0; c<colTotalSum.size(); c++) {
                    colTotalSum[c] = (int)colDistToAvg[c];
                }
                ////
                
                tmpRow.resize(colTotalSum.size());
                int balls = rowTotalSum[r];
                int index = colNum-1;//From last one to the first one
                while (balls > 0) {
                    if(balls !=0 && index == -1){
                        skip_bound = true;
                        cout<<"EFT::bound()::Failed calculating upper bound, use 0!"<<endl;
                        break;
                    }
                    if(balls <= colTotalSum[index]){
                        tmpRow[index] = balls;
                        balls = 0;
                    }else{
                        tmpRow[index] = colTotalSum[index];
                        balls -= tmpRow[index];
                        index--;
                    }
                    
                }
                
                if(skip_bound)break;
                
                double chisq_tmp;
                size_t df=0;
                vector<double> p_null(colNum, 1.0/colNum);
                ChisquareTest1DNoPValue(tmpRow, p_null, colNum, chisq_tmp, df);
                boundHeteroChisq_upper += chisq_tmp;
                tmpRow.clear();
            }
            
            if(skip_bound){
                result = NOT_TO_SKIP;
            }else{
                result = (boundHeteroChisq_upper - m_colSumChisq[k] +1e-07 < m_totalObservedChisq) ?
                TO_SKIP_ENTIRE_BRANCH : NOT_TO_SKIP;
            }
            
        }
        
        if(result == NOT_TO_SKIP){
            //Hua added, May 17 2014
            // compare lowest FunChisq and observed FunChisq
            //Lower bound
            bool skip_bound = false;
            vector<int> colTotalSum = colTotalSum_const;
            sort(colTotalSum.begin(), colTotalSum.end());
            
            for (size_t r=i; r <rowNum; r++) {
                tmpRow.resize(colTotalSum.size());
                int balls = rowTotalSum[r];
                int index = 0;
                double everage = (double)balls / (double)colNum;//Could be double
                
                
                
                while (balls > 0){
                    if(balls !=0 && index == (int) colTotalSum.size()){
                        skip_bound = true;
                        cout<<"EFT::bound()::Failed calculating lower bound, use 0!"<<endl;
                        break;
                    }
                    if(balls <= colTotalSum[index] && balls <= everage){
                        tmpRow[index] = balls;
                        balls = 0;
                        break;
                    }
                    if(colTotalSum[index] < everage){
                        tmpRow[index] = colTotalSum[index];
                        balls -= tmpRow[index];
                        if(balls==0)break;
                        index++;
                        everage = (double)balls / (double)(colNum - index);
                    }else{
                        
                        double chisq_tmp = 0;
                        
                        int N = rowTotalSum[r];
                        
                        if(N <= 0) {
                            skip_bound = true;
                            break;
                        }
                        
                        double x_exp = (double)N / (double)colNum;
                        for(int c=0; c<index; c++) {
                            
                            if(x_exp != 0) {
                                chisq_tmp += (tmpRow[c] - x_exp)*(tmpRow[c] - x_exp)/ x_exp;
                            } else if(tmpRow[c] != 0) {
                                cerr << "ERROR: expected is zero, but observed is not. Impossible!" << endl;
                                exit(EXIT_FAILURE);
                            }
                        }
                        
                        for(int k=index; k<colNum; k++) {
                            chisq_tmp += (everage - x_exp)*(everage - x_exp)/ x_exp;
                        }
                        
                        boundHeteroChisq_lower += chisq_tmp;
                        
                        balls = 0;
                        tmpRow.clear();
                        break;
                    }
                }
                
                if(skip_bound)break;
                
                
            }
            if(skip_bound){
                result = NOT_TO_SKIP;
            }else{
                result = (boundHeteroChisq_lower - m_colSumChisq[k] +1e-07 > m_totalObservedChisq) ?
                TO_KEEP_ENTIRE_BRANCH : NOT_TO_SKIP;
            }
            
        }
        if(result == TO_SKIP_ENTIRE_BRANCH) {
            // cout << "Suppose to skip";
            m_skip = true;
            // result = "not-to-skip";
        }
    }
    return result;
}
//////

    //Hua added, Apr 21 2014
	void EMTFunctionalChisq::setRowAndColSum(vector<vector<int> > &row, vector<vector<int> > &col, const vector<TransitionTable> & Cs){
        m_requiredRowSums.resize(row.size());
        m_requiredColSums.resize(col.size());
//m_nullFisherProb.resize(Cs.size());
        
        for (size_t k=0; k<row.size(); k++) {
            m_requiredRowSums[k].resize(row[k].size());
            m_requiredColSums[k].resize(col[k].size());
            for (size_t i=0; i<row[k].size(); i++) {
                m_requiredRowSums[k][i] = row[k][i];
            }
            for (size_t i=0; i<col[k].size(); i++) {
                m_requiredColSums[k][i] = col[k][i];
            }
        }
    }


////
//Add by Hua Apr 14 2014

double multi_table_functional_test(const vector<TransitionTable> & Cs,
                                   const string & discrepancy_measure)
{
    EMTEnumerator e;
    
    if ( discrepancy_measure == "chisq" ) {
        EMTFunctionalChisq v;
        double P = exact_multi_table_test(Cs, v, e);
//Hua commented, Apr 13 2015
//        if(v.getExtremeness()==true){
//            P = 1 - P;
//        }
        return P;
    } else {
        cerr << "ERROR: unknown discrepancy measure in multi_table_functional_test!"
        << endl;
        return 1;
    }
}

double exact_functional_test(const TransitionTable & C,
                             const string & discrepancy_measure)
{
    vector<TransitionTable> Cs(1, C);
    return multi_table_functional_test(Cs, discrepancy_measure);
}

double exact_functional_test(const vector< vector<int> > & C,
                             const string & discrepancy_measure)
{
    double pval=1.0;
    
    if( discrepancy_measure == "chisq" ) {
        
        TransitionTable tt(0, (int)C[0].size(), vector<int>(1,0),
                           vector<int>(1, (int)C.size()), vector<int>(1,0));
        tt.setTransitionTable(C);
        pval = exact_functional_test(tt, discrepancy_measure);
    }
    
    return pval;
}

////